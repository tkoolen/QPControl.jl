using RigidBodyDynamics: lower, upper, velocity_to_configuration_derivative_jacobian, velocity_to_configuration_derivative_jacobian!
using Compat: adjoint

struct MPCStage{C}
    Δt::Float64
    q::Vector{Variable}
    v::Vector{Variable}
    v̇::Vector{Variable}
    u::Vector{Variable}
    contacts::Vector{C}
    initialized::typeof(Ref(false))
end

struct MPCController{C, O <: MOI.AbstractOptimizer, M <: Mechanism, S <: MechanismState}
    mechanism::M
    state::S
    dynamicsresult::DynamicsResult{Float64, Float64}
    qpmodel::SimpleQP.Model{Float64, O}
    stages::Vector{MPCStage{C}}
    initialized::Base.RefValue{Bool}

    function MPCController{C}(mechanism::Mechanism, optimizer::O) where {C, O <: MOI.AbstractOptimizer}
        state = MechanismState(mechanism)
        dynamicsresult = DynamicsResult(mechanism)
        qpmodel = SimpleQP.Model(optimizer)
        stages = Vector{MPCStage{C}}()
        initialized = Ref(false)
        new{C, O, typeof(mechanism), typeof(state)}(mechanism, state, dynamicsresult, qpmodel, stages, initialized)
    end
end

horizon(c::MPCController) = length(c.stages)
stages(c::MPCController) = c.stages

# Needed to work around problems with infinite bounds in OSQP
const NEARLY_INFINITE = 1e9
make_finite(x::Real) = isfinite(x) ? x : NEARLY_INFINITE * sign(x)

function addstage!(controller::MPCController{C}, Δt::Real) where {C}
    @assert !controller.initialized[]
    model = controller.qpmodel
    mechanism = controller.mechanism
    state = controller.state

    q = [Variable(model) for _ in 1:num_positions(state)]
    v = [Variable(model) for _ in 1:num_velocities(state)]
    v̇ = [Variable(model) for _ in 1:num_velocities(state)]
    u = [Variable(model) for _ in 1:num_velocities(state)]

    u_limits = effort_bounds.(tree_joints(mechanism))
    u_lb = make_finite.(lower.(vcat(u_limits...)))
    u_ub = make_finite.(upper.(vcat(u_limits...)))
    @constraint model u >= u_lb
    @constraint model u <= u_ub

    contacts = Vector{C}()
    initialized = Ref(false)
    stage = MPCStage{C}(Δt, q, v, v̇, u, contacts, initialized)
    push!(controller.stages, stage)
    stage
end

function addstages!(c::MPCController, num_stages::Integer, Δt::Real)
    map(1:num_stages) do _
        addstage!(c, Δt)
    end
end

function position end
function wrench_world end


function addcontact!(stage::MPCStage{C}, contact::C) where C
    @assert !stage.initialized[]
    push!(stage.contacts, contact)
    contact
end

function addcontact!(controller::MPCController, stage::MPCStage{C}, position::Point3D, normal::FreeVector3D, μ::Float64) where C
    addcontact!(stage, C(position, normal, μ, controller.state, controller.qpmodel))
end

abstract type ContactModel end

struct BooleanContact <: ContactModel
    ρ_max::Float64
    separation_atol::Float64
    separation_max::Float64
end

BooleanContact(; ρ_max=10000, separation_atol=1e-6, separation_max=100) =
    BooleanContact(ρ_max, separation_atol, separation_max)

struct LCPContact <: ContactModel
    boolean_params::BooleanContact
    velocity_max::Float64
end

LCPContact(; ρ_max=10000, separation_atol=1e-6, separation_max=100, velocity_max=1000) =
    LCPContact(BooleanContact(ρ_max=ρ_max, separation_atol=separation_atol, separation_max=separation_max), velocity_max)


function addcontact!(controller::MPCController, stage::MPCStage{C}, position::Point3D, surface::HalfSpace3D, μ::Float64, contact_model::ContactModel=BooleanContact()) where C
    model = controller.qpmodel

    normal_in_body = let state = controller.state, surface = surface, position = position
        Parameter(model) do
            transform(state, surface.outward_normal, position.frame)
        end
    end

    contact_point = C(position, normal_in_body, μ, controller.state, controller.qpmodel)
    add_contact_indicators!(controller, stage, contact_point, surface, contact_model)
    addcontact!(stage, contact_point)
end

function add_contact_indicators!(controller::MPCController, stage::MPCStage, contact_point::ContactPoint, surface::HalfSpace3D, contact_model::BooleanContact)
    model = controller.qpmodel
    state = Parameter(identity, controller.state, model)
    mechanism = controller.state.mechanism
    position = contact_point.position
    body = body_fixed_frame_to_body(mechanism, position.frame)
    path_to_body = path(mechanism, root_body(mechanism), body)

    current_position_world = @expression transform(state, position, root_frame(mechanism))
    J_point = let state = state, position_world = current_position_world, path_to_body = path_to_body
        Parameter(point_jacobian(state(), path_to_body, position_world()), model) do J
            point_jacobian!(J, state(), path_to_body, position_world())
        end
    end
    position_world = @expression current_position_world.v + stage.Δt * (J_point.J * stage.v)
    surface_world = @expression(HalfSpace3D(
            transform(state, surface.point, root_frame(mechanism)),
            transform(state, surface.outward_normal, root_frame(mechanism))))

    separation = @expression(dot(surface_world.outward_normal.v, position_world) -
                             dot(surface_world.outward_normal.v, surface_world.point.v))

    indicator = Variable(model)
    @constraint(model, indicator ∈ {0, 1})

    separation_min = @expression(
        min((dot(surface_world.outward_normal.v, current_position_world.v) -
             dot(surface_world.outward_normal.v, surface_world.point.v)),
             -contact_model.separation_atol))

    @constraint(model, [separation] >= [separation_min])
    @constraint(model, [separation] <= [contact_model.separation_atol + contact_model.separation_max * (1 - indicator)])
    @constraint(model, contact_point.ρ <= fill(contact_model.ρ_max * indicator, length(contact_point.ρ)))

    indicator
end

function tangential_basis(num_basis_vectors::Val{N}) where N
    Δθ = 2 * π / N
    basisvectors = ntuple(num_basis_vectors) do i
        θ = (i - 1) * Δθ
        SVector(cos(θ), sin(θ), 0.0)
    end
    hcat(basisvectors...)
end

function add_contact_indicators!(controller::MPCController, stage::MPCStage, contact_point::ContactPoint{N}, surface::HalfSpace3D, contact_model::LCPContact) where N

    # separation ⫺ 0 ⟂ contact force ⫺ 0
    # Covers equations (7) and (10)
    add_contact_indicators!(controller, stage, contact_point, surface, contact_model.boolean_params)

    model = controller.qpmodel
    state = Parameter(identity, controller.state, model)
    mechanism = controller.state.mechanism
    position = contact_point.position
    body = body_fixed_frame_to_body(mechanism, position.frame)
    path_to_body = path(mechanism, root_body(mechanism), body)

    current_position_world = @expression transform(state, position, root_frame(mechanism))
    J_point = let state = state, position_world = current_position_world, path_to_body = path_to_body
        Parameter(point_jacobian(state(), path_to_body, position_world()), model) do J
            point_jacobian!(J, state(), path_to_body, position_world())
        end
    end
    velocity_world = @expression J_point.J * stage.v
    @assert contact_point.toroot().to == root_frame(mechanism)
    velocity_local = @expression rotation(inv(contact_point.toroot())) * velocity_world

    # Constraints to enforce sliding friction
    λ = Variable(model)
    @constraint(model, [λ] >= [0.0])
    @constraint(model, [λ] <= [contact_model.velocity_max])

    # The contact point defines variables ρ which multiply the basis vectors
    # but we need the variables β which multiply just the *tangential* basis.
    μ = contact_point.μ
    β = @expression((μ / sqrt(μ * μ + 1)) * contact_point.ρ)

    D = tangential_basis(Val(N))
    @constraint(model, fill(λ, N) + D' * velocity_local >= zeros(N))  # (8)

    λ_indicator = Variable(model)
    @constraint(model, λ_indicator ∈ {0, 1})
    # Equation (12)
    # λ_indicator == 1 => μc_n - sum(β) == 0
    # λ_indicator == 0 => λ == 0
    c_n = @expression(dot(contact_point.force_local, FreeVector3D(contact_point.normal_aligned_frame, 0.0, 0.0, 1.0)))
    @constraint(model, [μ * c_n - ones(N)' * β] <= [contact_model.boolean_params.ρ_max * (1 - λ_indicator)])
    @constraint(model, [λ] <= [λ_indicator * contact_model.velocity_max])

    β_indicator = [Variable(model) for _ in 1:N]
    for i in 1:N
        @constraint(model, β_indicator[i] ∈ {0, 1})
    end
    # Equation (11)
    # β_indicator[i] == 1 => λ_plus_D'v[i] == 0
    # β_indicator[i] == 0 => β[i] == 0
    @constraint(model, fill(λ, N) + D' * velocity_local <= contact_model.velocity_max * (ones(N) - β_indicator))
    @constraint(model, β <= contact_model.boolean_params.ρ_max * β_indicator)
end

function addjointlimit!(controller::MPCController, stage::MPCStage, joint::Joint, τ_max::Real)
    bounds = RigidBodyDynamics.position_bounds(joint)
    lb = RigidBodyDynamics.lower.(bounds)
    ub = RigidBodyDynamics.upper.(bounds)
    if any(isfinite, lb) || any(isfinite, ub)
        throw(ArgumentError("Cannot apply position bound enforcement torques to this joint type: $(joint.joint_type). You may set the position bounds on this joint to (-Inf, Inf) to avoid attempting to add such torques."))
    end
    τ = [Variable(controller.qpmodel) for _ in 1:num_velocities(joint)]
    if num_velocities(joint) > 0
        @constraint(controller.qpmodel, τ == zeros(num_velocities(joint)))
    end
    τ
end

function total_Δt(controller::MPCController, stage::MPCStage)
    # Compute the total time until the end of the given stage
    stage_index = findfirst(controller.stages, stage)
    if VERSION <= v"0.7-alpha"
        @assert stage_index != 0        # 0.6
    else
        @assert stage_index !== nothing # 0.7
    end
    sum(s -> s.Δt, controller.stages[1:stage_index])
end

function addjointlimit!(controller::MPCController, stage::MPCStage, joint::Joint{<:Any, <:OneDegreeOfFreedomFixedAxis}, τ_max::Real)
    @assert num_velocities(joint) == 1

    bounds = only(RigidBodyDynamics.position_bounds(joint))
    q_lb = RigidBodyDynamics.lower(bounds)
    q_ub = RigidBodyDynamics.upper(bounds)

    model = controller.qpmodel
    τ = Variable(model)
    if isfinite(q_lb) || isfinite(q_ub)
        (isfinite(q_lb) && isfinite(q_ub)) || throw(ArgumentError("Upper and lower bounds must be both finite or both infinite"))

        # indicator_lb is 1 iff this joint is at its lower limit
        indicator_lb = Variable(model)
        @constraint(model, indicator_lb ∈ {0, 1})
        # indicator_ub is 1 iff this joint is at its upper limit
        indicator_ub = Variable(model)
        @constraint(model, indicator_ub ∈ {0, 1})

        state = Parameter(identity, controller.state, model)
        q_current = @expression only(configuration(state, joint))
        q_next = stage.q[configuration_range(controller.state, joint)]
        # If the joint is *already* violating its limits (e.g. due to soft
        # enforcement in the simulator), then we adjust the lower and upper bounds
        # accordingly. Otherwise soft enforcement in the simulation can make the control
        # problem infeasible.
        q_lb_widened = @expression min(q_current, q_lb)
        q_ub_widened = @expression max(q_current, q_ub)
        q_range = @expression q_ub_widened - q_lb_widened

        @constraint(model, q_next <= [q_lb_widened + q_range * (1 - indicator_lb)])
        @constraint(model, q_next >= [q_ub_widened - q_range * (1 - indicator_ub)])
        @constraint(model, [indicator_lb + indicator_ub] <= 1)

        # If indicator_lb is zero, then τ <= 0
        @constraint(model, [τ] <= [τ_max * indicator_lb])
        # If indicator_ub is zero, then τ >= 0
        @constraint(model, [τ] >= [-τ_max * indicator_ub])

        # Most joints will be far from their upper or lower limit most of the time,
        # so we can simplify the optimization by checking if a joint could possibly
        # hit its limit during the given optimization and eliminate one or both
        # indicator variables if it cannot.
        v_lb = RigidBodyDynamics.lower(only(RigidBodyDynamics.velocity_bounds(joint)))
        v_ub = RigidBodyDynamics.upper(only(RigidBodyDynamics.velocity_bounds(joint)))

        Δt = total_Δt(controller, stage)

        q_after_max_v = @expression q_current + v_ub * Δt
        indicator_ub_restriction = @expression(ifelse(q_after_max_v < q_ub, 0.0, 1.0))
        @constraint(model, [indicator_ub] <= [indicator_ub_restriction])

        q_after_min_v = @expression q_current + v_lb * Δt
        indicator_lb_restriction = @expression(ifelse(q_after_min_v > q_lb, 0.0, 1.0))
        @constraint(model, [indicator_lb] <= [indicator_lb_restriction])
    else
        @constraint(model, [τ] == [0.0])
    end
    [τ]
end

function addjointlimits!(controller::MPCController, stage::MPCStage, force_max=10000.0)
    vcat([addjointlimit!(controller, stage, joint, force_max)
          for joint in tree_joints(controller.mechanism)]...)
end


"""
Holds a few useful dynamics parameters which are re-used for each stage in the MPC
optmization.
"""
struct DynamicsParams{TH, Tc, TJ}
    H::TH
    c::Tc
    J_qv::TJ
end

function generalized_torque(state, contact_point::ContactPoint, model::Model)
    mechanism = state.mechanism
    body = body_fixed_frame_to_body(mechanism, contact_point.position.frame)
    path_to_root = path(mechanism, body, root_body(mechanism))
    J_body_to_root = let state = state, path_to_root = path_to_root
        Parameter(geometric_jacobian(state, path_to_root), model) do J
            geometric_jacobian!(J, state, path_to_root)
        end
    end
    @expression adjoint(angular(J_body_to_root)) * angular(contact_point.wrench_world) + adjoint(linear(J_body_to_root)) * linear(contact_point.wrench_world)
end

function initialize!(controller::MPCController, stage::MPCStage,
                     params::DynamicsParams, q_prev, v_prev)
    model = controller.qpmodel
    state = controller.state
    @assert !stage.initialized[]
    Δt = stage.Δt
    H = params.H
    c = params.c
    J_qv = params.J_qv
    τ_ext = zeros(num_velocities(state))
    for contact in stage.contacts
        τ_ext = @expression τ_ext + generalized_torque(state, contact, model)
    end
    τ_joint_limit = addjointlimits!(controller, stage)
    @constraint(model, H * (stage.v - v_prev) == Δt * (stage.u + τ_joint_limit - c - τ_ext))
    q̇ = @expression(J_qv * stage.v)
    @constraint(model, q̇ == (1 / Δt) * (stage.q - q_prev))
    @constraint(model, stage.v̇ == (1 / Δt) * (stage.v - v_prev))
    stage.initialized[] = true
end


function initialize!(controller::MPCController)
    @assert !controller.initialized[]
    model = controller.qpmodel
    # terminal_state_cost = controller.terminal_state_cost

    H = let state = controller.state
        Parameter(mass_matrix(state), model) do H
            mass_matrix!(H, state)
        end
    end
    c = let state = controller.state, result = controller.dynamicsresult
        Parameter(result.dynamicsbias, model) do _
            dynamics_bias!(result, state)
        end
    end
    J_qv = let state = controller.state
        Parameter(velocity_to_configuration_derivative_jacobian(state), model) do J
            velocity_to_configuration_derivative_jacobian!(J, state)
        end
    end
    dynamics_params = DynamicsParams(H, c, J_qv)

    for i in eachindex(controller.stages)
        if i == 1
            q_prev = let state = controller.state
                Parameter(() -> configuration(state), model)
            end
            v_prev = let state = controller.state
                Parameter(() -> velocity(state), model)
            end
        else
            q_prev = controller.stages[i - 1].q
            v_prev = controller.stages[i - 1].v
        end
        stage = controller.stages[i]
        initialize!(controller, stage, dynamics_params, q_prev, v_prev)
    end
    controller.initialized[] = true
end

function (controller::MPCController)(τ::AbstractVector, t::Number, x::Union{<:Vector, <:MechanismState})
    if !controller.initialized[]
        initialize!(controller)
    end
    @assert controller.initialized[]
    # println("got:    ", Vector(configuration(x)), " ", Vector(velocity(x)))

    copyto!(controller.state, x)
    solve!(controller.qpmodel)
    τ .= SimpleQP.value.(controller.qpmodel, first(stages(controller)).u)
    if any(isnan, τ)
        @show SimpleQP.primalstatus(controller.qpmodel)
        @show t, Vector(x), τ
    end
    # println("expect: ", SimpleQP.value.(controller.qpmodel, controller.stages[1].q), " ", SimpleQP.value.(controller.qpmodel, controller.stages[1].q))
end
