struct MomentumBasedController{O<:MOI.AbstractOptimizer, S<:MechanismState, N}
    # dynamics-related
    state::S
    result::DynamicsResult{Float64, Float64}
    centroidalframe::CartesianFrame3D
    externalwrenches::Dict{RigidBody{Float64}, Wrench{Float64}}
    contactsettings::Vector{Pair{ContactSettings{N}, Vector{Variable}}}
    qpmodel::SimpleQP.Model{O}
    objective::SimpleQP.LazyExpression # to incrementally build the objective function

    function MomentumBasedController{N}(mechanism::Mechanism{Float64}, optimizer::O) where {O<:MOI.AbstractOptimizer, N}
        state = MechanismState(mechanism)
        result = DynamicsResult(mechanism)
        centroidalframe = CartesianFrame3D("centroidal")
        rootframe = root_frame(mechanism)
        externalwrenches = Dict(b => zero(Wrench{Float64}, rootframe) for b in bodies(mechanism))
        contactsettings = Vector{Pair{ContactSettings{N}, Vector{Variable}}}()
        qpmodel = SimpleQP.Model(optimizer)
        v̇ = [Variable(model) for _ = 1 : num_velocities(state)]
        objective = SimpleQP.LazyExpression(identity, 0.0)
        # ρ = Vector{Float64}()
        new{O, typeof(state), N}(state, result, centroidalframe, externalwrenches, contactsettings, qpmodel, v̇, objective)
    end
end

centroidal_frame(controller::MomentumBasedController) = controller.centroidalframe

function (controller::MomentumBasedController)(τ::AbstractVector, t::Number, x::Union{<:Vector, <:MechanismState})
    state = controller.state
    result = controller.result
    qpmodel = controller.qpmodel
    externalwrenches = result.externalwrenches
    copyto!(state, x)
    solve!(qpmodel)
    # TODO: check status
    result.v̇ .= value.(qpmodel, controller.v̇)
    com = center_of_mass(state)
    centroidal_to_world = Transform3D(controller.centroidalframe, com.frame, com.v)
    map!(@closure(wrench -> transform(wrench, centroidal_to_world)), values(externalwrenches), values(externalwrenches))
    inverse_dynamics!(τ, result.jointwrenches, result.accelerations, state, result.v̇, externalwrenches)
    τ
end

function addtask!(controller::MomentumBasedController, task::AbstractMotionTask)
    model = controller.qpmodel
    state = controller.state
    v̇ = controller.v̇
    @constraint(model, task_error(task, model, state, v̇) == zeros(dimension(task)))
    nothing
end

function add_task_error_slack_variables!(controller::MomentumBasedController, task::AbstractMotionTask)
    model = controller.qpmodel
    state = controller.state
    v̇ = controller.v̇
    e = [Variable(model) for _ = 1 : dimension(task)]
    @constraint(model, e == task_error(task, model, state, v̇))
    e
end

function addtask!(controller::MomentumBasedController, task::AbstractMotionTask, weight::Number)
    e = add_task_error_slack_variables!(controller, task)
    controller.objective = @expression controller.objective + weight * (e ⋅ e)
    e
end

function addtask!(controller::MomentumBasedController, task::AbstractMotionTask, weight::AbstractMatrix)
    e = add_task_error_slack_variables!(controller, task)
    controller.objective = @expression controller.objective + e' * weight * e
    e
end

function initialize!(controller::MomentumBasedController)
    # TODO: add wrench balance constraint
end
