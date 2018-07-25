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
    initialized::typeof(Ref(false))

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

function addcontact!(controller::MPCController, stage::MPCStage{C}, position::Point3D, surface::HalfSpace3D, μ::Float64, max_ρ=10000) where C
    model = controller.qpmodel

    normal_in_body = let state = controller.state, surface = surface, position = position
        Parameter(model) do
            transform(state, surface.outward_normal, position.frame)
        end
    end

    indicator = add_boolean_indicator!(controller, stage, position, surface)
    contact_point = C(position, normal_in_body, μ, controller.state, controller.qpmodel)
    @constraint(model, contact_point.ρ <= fill(max_ρ * indicator, length(contact_point.ρ)))

    addcontact!(stage, contact_point)
end

function add_boolean_indicator!(controller::MPCController, stage::MPCStage, position::Point3D, surface::HalfSpace3D, separation_atol=1e-2, separation_max=100)
    state = controller.state
    mechanism = state.mechanism
    body = body_fixed_frame_to_body(mechanism, position.frame)
    path_to_body = path(mechanism, root_body(mechanism), body)
    model = controller.qpmodel

    position_world = let state = state,  position = position
        Parameter(model) do
            transform(state, position, root_frame(state.mechanism))
        end
    end
    J_point = let state = state, position_world = position_world, path_to_body = path_to_body
        Parameter(point_jacobian(state, path_to_body, position_world()), model) do J
            point_jacobian!(J, state, path_to_body, position_world())
        end
    end
    linearized_position_world = @expression position_world.v + stage.Δt * (J_point.J * stage.v)

    surface_world = let surface = surface, state = state
        Parameter(model) do
            HalfSpace3D(
                transform(state, surface.point, root_frame(state.mechanism)),
                transform(state, surface.outward_normal, root_frame(state.mechanism)))
        end
    end
    separation = @expression(dot(surface_world.outward_normal.v, linearized_position_world) -
                             dot(surface_world.outward_normal.v, surface_world.point.v))

    indicator = Variable(model)
    @constraint(model, indicator ∈ {0, 1})
    @constraint(model, [separation] >= [-separation_atol])
    @constraint(model, [separation] <= [separation_atol + separation_max * (1 - indicator)])

    indicator
end

"""
Holds a few useful dynamics paramters which are re-used for each stage in the MPC
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

function initialize!(stage::MPCStage, model::Model, state::MechanismState,
                     params::DynamicsParams, q_prev, v_prev)
    @assert !stage.initialized[]
    Δt = stage.Δt
    H = params.H
    c = params.c
    J_qv = params.J_qv
    τ_ext = zeros(num_velocities(state))
    for contact in stage.contacts
        τ_ext = @expression τ_ext + generalized_torque(state, contact, model)
    end
    @show τ_ext
    @constraint(model, H * (stage.v - v_prev) == Δt * (stage.u - c - τ_ext))
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
        Parameter(result.dynamicsbias, model) do c
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
        initialize!(stage, controller.qpmodel, controller.state, dynamics_params, q_prev, v_prev)
    end
    controller.initialized[] = true
end

function (controller::MPCController)(τ::AbstractVector, t::Number, x::Union{<:Vector, <:MechanismState})
    if !controller.initialized[]
        initialize!(controller)
    end
    @assert controller.initialized[]

    copyto!(controller.state, x)
    solve!(controller.qpmodel)
    τ .= SimpleQP.value.(controller.qpmodel, first(stages(controller)).u)
    if any(isnan, τ)
        @show SimpleQP.primalstatus(controller.qpmodel)
        @show t, Vector(x), τ
    end
end
