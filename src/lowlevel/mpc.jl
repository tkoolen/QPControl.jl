using RigidBodyDynamics: lower, upper, velocity_to_configuration_derivative_jacobian, velocity_to_configuration_derivative_jacobian!
using Compat: adjoint

struct QuadraticCost{T}
    Q::Matrix{T}
    q::Vector{T}
    x0::Vector{T}
end

QuadraticCost{T}(n::Integer) where {T} = QuadraticCost(zeros(T, n, n), zeros(T, n), zeros(T, n))

all_effort_bounds(m::Mechanism) =
    collect(Base.Iterators.flatten(map(effort_bounds, joints(m))))

struct MPCStage{N, O<:MOI.AbstractOptimizer, S<:MechanismState}
    qpmodel::SimpleQP.Model{Float64, O}
    state::S
    result::DynamicsResult{Float64, Float64}
    Δt::Float64
    configuration::Vector{Variable}
    velocity::Vector{Variable}
    input::Vector{Variable}
    contacts::Dict{RigidBody{Float64},
                   Vector{Pair{ContactPoint{N},
                               HalfSpace3D{Float64}}}}
end


mutable struct MPCController{N, O<:MOI.AbstractOptimizer, M<:Mechanism, S<:MechanismState}
    mechanism::M
    qpmodel::SimpleQP.Model{Float64, O}
    stages::Vector{MPCStage{N, O, S}}
    running_state_cost::QuadraticCost{Float64}
    running_input_cost::QuadraticCost{Float64}
    terminal_state_cost::QuadraticCost{Float64}
    objective::SimpleQP.LazyExpression
    initialized::Bool

    function MPCController{N}(mechanism::Mechanism{Float64}, optimizer::O) where {N, O<:MOI.AbstractOptimizer}
        model = SimpleQP.Model(optimizer)
        M = typeof(mechanism)
        S = typeof(MechanismState(mechanism))
        stages = Vector{MPCStage{N, O, S}}()
        nq = num_positions(mechanism)
        nv = num_velocities(mechanism)
        running_state_cost = QuadraticCost{Float64}(nq + nv)
        running_input_cost = QuadraticCost{Float64}(nv)
        terminal_state_cost = QuadraticCost{Float64}(nq + nv)
        objective = SimpleQP.LazyExpression(identity, zero(QuadraticFunction{Float64}))

        new{N, O, M, S}(
            mechanism,
            model,
            stages,
            running_state_cost,
            running_input_cost,
            terminal_state_cost,
            objective,
            false
            )
    end
end

horizon(c::MPCController) = length(c.stages)
stages(c::MPCController) = c.stages

function addstage!(controller::MPCController{N, O, M, S}, Δt::Real = 0.001) where {N, O, M, S}
    @assert !controller.initialized
    model = controller.qpmodel
    mechanism = controller.mechanism
    state::S = MechanismState(mechanism)

    qnext = [Variable(model) for _ in 1:num_positions(state)]
    vnext = [Variable(model) for _ in 1:num_velocities(state)]
    u = [Variable(model) for _ in 1:num_velocities(state)]
    @constraint model u >= lower.(all_effort_bounds(mechanism))
    @constraint model u <= upper.(all_effort_bounds(mechanism))

    contacts = Dict{RigidBody{Float64}, Vector{Pair{ContactPoint{N}, HalfSpace3D{Float64}}}}()
    result = DynamicsResult(mechanism)
    stage = MPCStage{N, O, S}(
        controller.qpmodel, state, result, Δt,
        qnext, vnext, u, contacts)
    push!(controller.stages, stage)
    stage
end

function addstages!(c::MPCController, num_stages::Integer, Δt::Real = 0.001)
    for _ in 1:num_stages
        addstage!(c, Δt)
    end
end

function addcontact!(
        stage::MPCStage{N, O, S},
        body::RigidBody{Float64},
        position::Point3D,
        normal::FreeVector3D,
        μ::Float64,
        surface::HalfSpace3D,
        ) where {N, O, S}
    contact_point = ContactPoint{N}(position, normal,
                                    μ, stage.state,
                                    stage.qpmodel)
    push!(get!(Vector{Pair{ContactPoint{N}, HalfSpace3D{Float64}}}, stage.contacts, body), (contact_point => surface))
    # TODO: add objective term
    contact_point
end

function initialize_stage!(controller::MPCController, stage_index::Integer)
    stage = stages(controller)[stage_index]

    state = stage.state
    mechanism = state.mechanism
    result = stage.result
    model = stage.qpmodel
    Δt = stage.Δt
    qnext = stage.configuration
    vnext = stage.velocity
    u = stage.input

    H = let state = state
        Parameter(mass_matrix(state), model) do H
            mass_matrix!(H, state)
        end
    end
    c = let state = state, result = result
        Parameter(result.dynamicsbias, model) do c
            dynamics_bias!(result, state)
        end
    end
    J_qv = let state = state
        Parameter(velocity_to_configuration_derivative_jacobian(state), model) do J
            velocity_to_configuration_derivative_jacobian!(J, state)
        end
    end
    τ_ext = SimpleQP.LazyExpression(identity, zeros(AffineFunction{Float64}, num_velocities(state)))

    for (body, contact_points) in stage.contacts
        path_to_body = path(state.mechanism, root_body(state.mechanism), body)
        for (contact_point, surface) in contact_points
            point_in_world = let state = state, point = contact_point
                Parameter(model) do
                    (transform_to_root(state, point.position.frame) * point.position).v
                end
            end
            J_point = let state = state,
                          point = contact_point,
                          path_to_body = path_to_body,
                          J_point = point_jacobian(state, path_to_body, transform_to_root(state, point.position.frame) * point.position)
                Parameter(J_point.J, model) do _
                    p = transform_to_root(state, point.position.frame) * point.position
                    point_jacobian!(J_point, state, path_to_body, p)
                end
            end
            linearized_position_in_world = @expression point_in_world + Δt * (J_point * vnext)
            surface_origin = let surface = surface, mechanism = mechanism
                Parameter(model) do
                    @framecheck surface.point.frame root_frame(mechanism)
                    surface.point.v
                end
            end
            surface_normal = let surface = surface, mechanism = mechanism
                Parameter(model) do
                    @framecheck surface.outward_normal.frame root_frame(mechanism)
                    surface.outward_normal.v
                end
            end
            separation = @expression(dot(surface_normal,
                                         linearized_position_in_world) -
                                     dot(surface_normal, surface_origin))
            lb = -1e-2
            ub = let point = contact_point
                Parameter(model) do
                    isenabled(point) ? 1e-2 : 1e9  # TODO: Inf?
                end
            end
            @constraint(model, [separation] >= [lb])  # TODO: make scalar
            @constraint(model, [separation] <= [ub])

            J = let state = state, path_to_body = path_to_body
                Parameter(geometric_jacobian(state, path_to_body), model) do J
                    geometric_jacobian!(J, state, path_to_body)
                end
            end
            τ_ext = @expression τ_ext + adjoint(angular(J)) * angular(contact_point.wrench_world) + adjoint(linear(J)) * linear(contact_point.wrench_world)
        end
    end

    if stage_index == 1
        q_current = let state = state
            Parameter(model) do
                configuration(state)
            end
        end
        v_current = let state = state
            Parameter(model) do
                velocity(state)
            end
        end
    else
        q_current = controller.stages[stage_index - 1].configuration
        v_current = controller.stages[stage_index - 1].velocity
    end
    @constraint(model, H * (vnext - v_current) == Δt * (u - c - τ_ext))
    q̇ = @expression(J_qv * vnext)
    @constraint(model, qnext - q_current == Δt * q̇)

    objective = controller.objective
    objective = @expression(objective +
        objectiveterm(model,
                      controller.running_state_cost,
                      vcat(qnext, vnext)))
    objective = @expression(objective +
        objectiveterm(model,
                      controller.running_input_cost,
                      u))
    controller.objective = objective
end


function initialize!(controller::MPCController)
    @assert !controller.initialized
    model = controller.qpmodel
    terminal_state_cost = controller.terminal_state_cost

    for i in eachindex(controller.stages)
        initialize_stage!(controller, i)
    end

    objective = controller.objective
    xfinal = vcat(controller.stages[end].configuration,
                  controller.stages[end].velocity)
    objective = @expression(objective +
        objectiveterm(model, terminal_state_cost, xfinal))

    controller.objective = objective
    setobjective!(controller)
    controller.initialized = true
end

function setobjective!(controller::MPCController)
    @objective(controller.qpmodel, Minimize, controller.objective)
end

function objectiveterm(model, cost::QuadraticCost, x)
    x0 = let cost = cost
        Parameter(() -> cost.x0, model)
    end
    Q = let cost = cost
        Parameter(() -> cost.Q, model)
    end
    q = let cost = cost
        Parameter(() -> cost.q, model)
    end
    x̄ = [Variable(model) for _ in 1:length(cost.x0)]
    @constraint(model, x̄ == x - x0)
    @expression(x̄' * Q * x̄ + dot(q, x̄))
end

function (controller::MPCController)(τ::AbstractVector, t::Number, x::Union{<:Vector, <:MechanismState})
    if !controller.initialized
        initialize!(controller)
    end
    @assert controller.initialized

    # TODO: update controller contact normals
    for stage in stages(controller)
        copyto!(stage.state, x)
    end
    solve!(controller.qpmodel)
    τ .= SimpleQP.value.(controller.qpmodel, first(stages(controller)).input)
end

function addcontact!(
        controller::MPCController{N},
        body::RigidBody{Float64},
        position::Point3D,
        normal::FreeVector3D,
        μ::Float64,
        surface::HalfSpace3D,
        ) where {N}
    @assert !controller.initialized
    contact_point = ContactPoint{N}(position, normal,
                                    μ, controller.state,
                                    controller.qpmodel)
    push!(get!(Vector{Pair{ContactPoint{N}, HalfSpace3D{Float64}}}, controller.contacts, body), (contact_point => surface))
    # TODO: add objective term
    contact_point
end
