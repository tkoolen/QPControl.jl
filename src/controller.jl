struct MomentumBasedController{N, O<:MOI.AbstractOptimizer, S<:MechanismState}
    state::S
    result::DynamicsResult{Float64, Float64}
    centroidalframe::CartesianFrame3D
    externalwrenches::Dict{RigidBody{Float64}, Wrench{Float64}}
    contactsettings::Dict{ContactSettings{N}, SVector{N, Variable}}
    qpmodel::SimpleQP.Model{Float64, O}
    v̇::Vector{Variable}
    objective::SimpleQP.LazyExpression # buffer to incrementally build the objective function
    initialized::Base.RefValue{Bool}

    function MomentumBasedController{N}(mechanism::Mechanism{Float64}, optimizer::O) where {N, O<:MOI.AbstractOptimizer}
        state = MechanismState(mechanism)
        result = DynamicsResult(mechanism)
        centroidalframe = CartesianFrame3D("centroidal")
        rootframe = root_frame(mechanism)
        externalwrenches = Dict(b => zero(Wrench{Float64}, rootframe) for b in bodies(mechanism))
        contactsettings = Dict{ContactSettings{N}, SVector{N, Variable}}()
        qpmodel = SimpleQP.Model(optimizer)
        v̇ = [Variable(qpmodel) for _ = 1 : num_velocities(state)]
        objective = SimpleQP.LazyExpression(identity, zero(QuadraticFunction{Float64}))
        new{N, O, typeof(state)}(state, result, centroidalframe, externalwrenches, contactsettings, qpmodel, v̇, objective, Ref(false))
    end
end

centroidal_frame(controller::MomentumBasedController) = controller.centroidalframe

function (controller::MomentumBasedController)(τ::AbstractVector, t::Number, x::Union{<:Vector, <:MechanismState})
    if !controller.initialized[]
        initialize!(controller)
        controller.initialized[] = true
    end
    qpmodel = controller.qpmodel
    state = controller.state
    result = controller.result
    # externalwrenches = result.externalwrenches
    copyto!(state, x)
    solve!(qpmodel)
    # TODO: check status
    result.v̇ .= value.(qpmodel, controller.v̇)
    # com = center_of_mass(state)
    # centroidal_to_world = Transform3D(controller.centroidalframe, com.frame, com.v)
    # map!(@closure(wrench -> transform(wrench, centroidal_to_world)), values(externalwrenches), values(externalwrenches)) # TODO: desirable?
    inverse_dynamics!(τ, result.jointwrenches, result.accelerations, state, result.v̇)#, externalwrenches)
    τ
end

function addtask!(controller::MomentumBasedController, task::AbstractMotionTask)
    model = controller.qpmodel
    state = controller.state
    v̇ = controller.v̇
    @constraint(model, task_error(task, model, state, v̇) == zeros(dimension(task)))
    nothing
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

function add_task_error_slack_variables!(controller::MomentumBasedController, task::AbstractMotionTask)
    model = controller.qpmodel
    state = controller.state
    v̇ = controller.v̇
    e = [Variable(model) for _ = 1 : dimension(task)]
    @constraint(model, e == task_error(task, model, state, v̇))
    e
end

function addcontact!(controller::MomentumBasedController{N}, contactsettings::ContactSettings{N}) where N
    model = controller.qpmodel
    ρ = [Variable(model) for _ = 1 : N]
    @constraint(model, ρ >= zeros(N))
    push!(controller.contactsettings, contactsettings => ρ)
    nothing
end

@noinline function initialize!(controller::MomentumBasedController)
    @objective(controller.qpmodel, Minimize, controller.objective)
    # TODO: add wrench balance constraint
end
