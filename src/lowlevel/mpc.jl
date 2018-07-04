using RigidBodyDynamics: lower, upper, velocity_to_configuration_derivative_jacobian, velocity_to_configuration_derivative_jacobian!

struct QuadraticCost{T}
    Q::Matrix{T}
    q::Vector{T}
    x0::Vector{T}
end

QuadraticCost{T}(n::Integer) where {T} = QuadraticCost(zeros(T, n, n), zeros(T, n), zeros(T, n))

all_effort_bounds(m::Mechanism) =
    collect(Base.Iterators.flatten(map(effort_bounds, joints(m))))

struct MPCController{N, O<:MOI.AbstractOptimizer, S<:MechanismState}
    state::S
    result::DynamicsResult{Float64, Float64}
    horizon::Int
    Δt::Float64
    contacts::Dict{RigidBody{Float64},
                   Vector{Pair{ContactPoint{N},
                               HalfSpace3D{Float64}}}}
    qpmodel::SimpleQP.Model{Float64, O}
    configurations::Vector{Vector{Variable}}
    velocities::Vector{Vector{Variable}}
    inputs::Vector{Vector{Variable}}
    running_state_cost::QuadraticCost{Float64}
    running_input_cost::QuadraticCost{Float64}
    terminal_state_cost::QuadraticCost{Float64}
    initialized::Bool

    function MPCController{N}(mechanism::Mechanism{Float64}, optimizer::O;
                              horizon::Integer = 1,
                              Δt = 0.001) where {N, O<:MOI.AbstractOptimizer}
        state = MechanismState(mechanism)
        result = DynamicsResult(mechanism)
        model = SimpleQP.Model(optimizer)

        configurations = Vector{Vector{Variable}}()
        velocities = Vector{Vector{Variable}}()
        inputs = Vector{Vector{Variable}}()
        running_state_cost = QuadraticCost{Float64}(num_positions(state) + num_velocities(state))
        running_input_cost = QuadraticCost{Float64}(num_velocities(state))
        terminal_state_cost = QuadraticCost{Float64}(num_positions(state) + num_velocities(state))



        new{N, O, typeof(state)}(
            state,
            result,
            horizon,
            Δt,
            Dict{RigidBody{Float64}, Vector{Pair{ContactPoint{N}, HalfSpace3D{Float64}}}}(),
            model,
            configurations,
            velocities,
            inputs,
            running_state_cost,
            running_input_cost,
            terminal_state_cost,
            false)
    end
end

function initialize!(controller::MPCController)
    horizon = controller.horizon
    Δt = controller.Δt
    model = controller.qpmodel
    state = controller.state
    mechanism = state.mechanism
    result = controller.result
    running_state_cost = controller.running_state_cost
    running_input_cost = controller.running_input_cost
    terminal_state_cost = controller.terminal_state_cost
    configurations = controller.configurations
    velocities = controller.velocities

    objective = SimpleQP.LazyExpression(identity, zero(QuadraticFunction{Float64}))

    for i in 1:horizon
        qnext = [Variable(model) for _ in 1:num_positions(state)]
        push!(controller.configurations, qnext)
        vnext = [Variable(model) for _ in 1:num_velocities(state)]
        push!(controller.velocities, vnext)
        u = [Variable(model) for _ in 1:num_velocities(state)]
        push!(controller.inputs, u)
        @constraint model u >= lower.(all_effort_bounds(mechanism))
        @constraint model u <= upper.(all_effort_bounds(mechanism))

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

        for (body, contact_points) in controller.contacts
            path_to_root = path(state.mechanism, body, root_body(state.mechanism))
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
                              J_point = point_jacobian(state, path_to_body, transform_to_root(state, point.positio.frame) * point.position)
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

                J = let state = state, path_to_root = path_to_root
                    Parameter(geometric_jacobian(state, path_to_root), model) do J
                        geometric_jacobian!(J, state, path_to_root)
                    end
                end
                τ_ext = @expression τ_ext + transpose(angular(J)) * angular(contact_point.wrench_world) + transpose(linear(J)) * linear(contact_point.wrench_world)
            end
        end

        if i == 1
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
            q_current = configurations[i - 1]
            v_current = velocities[i - 1]
        end
        @constraint(model, H * (vnext - v_current) == Δt * (u - c - τ_ext))
        q̇ = @expression(J_qv * vnext)
        @constraint(model, qnext - q_current == Δt * q̇)

        objective = @expression(objective +
            objectiveterm(model,
                          running_state_cost,
                          vcat(qnext, vnext)))
        objective = @expression(objective +
            objectiveterm(model,
                          running_input_cost,
                          u))
    end
    objective = @expression(objective +
        objectiveterm(model,
                      terminal_state_cost,
                      vcat(configurations[end],
                           velocities[end])))

    @objective(model, Minimize, objective)
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
    @expression(x̄' * Q * x̄ + q' * x̄)
end
