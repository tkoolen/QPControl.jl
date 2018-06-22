mutable struct MomentumBasedController{N, O<:MOI.AbstractOptimizer, S<:MechanismState}
    state::S
    result::DynamicsResult{Float64, Float64}
    centroidalframe::CartesianFrame3D
    momentum_matrix::MomentumMatrix{Matrix{Float64}}
    contactdata::Dict{RigidBody{Float64}, Vector{ContactQPData{N}}}
    contactwrenches::Dict{BodyID, Wrench{Float64}}
    qpmodel::SimpleQP.Model{Float64, O}
    v̇::Vector{Variable}
    objective::SimpleQP.LazyExpression # buffer to incrementally build the objective function
    initialized::Bool

    function MomentumBasedController{N}(mechanism::Mechanism{Float64}, optimizer::O) where {N, O<:MOI.AbstractOptimizer}
        state = MechanismState(mechanism)
        result = DynamicsResult(mechanism)
        worldframe = root_frame(mechanism)
        centroidalframe = CartesianFrame3D("centroidal")
        nv = num_velocities(state)
        momentum_matrix = MomentumMatrix(worldframe, zeros(3, nv), zeros(3, nv))
        rootframe = root_frame(mechanism)
        contactdata = Dict{RigidBody{Float64}, Vector{ContactQPData{N}}}()
        contactwrenches = Dict{BodyID, Wrench{Float64}}()
        qpmodel = SimpleQP.Model(optimizer)
        v̇ = [Variable(qpmodel) for _ = 1 : nv]
        objective = SimpleQP.LazyExpression(identity, zero(QuadraticFunction{Float64}))
        new{N, O, typeof(state)}(
            state, result, centroidalframe, momentum_matrix, contactdata, contactwrenches, qpmodel, v̇, objective, false)
    end
end

centroidal_frame(controller::MomentumBasedController) = controller.centroidalframe

function (controller::MomentumBasedController)(τ::AbstractVector, t::Number, x::Union{<:Vector, <:MechanismState})
    if !controller.initialized
        initialize!(controller)
        controller.initialized = true
    end

    qpmodel = controller.qpmodel
    state = controller.state
    result = controller.result
    contactdata = controller.contactdata
    contactwrenches = controller.contactwrenches
    worldframe = root_frame(state.mechanism)

    copyto!(state, x)
    solve!(qpmodel)
    # @assert terminationstatus(qpmodel) == MOI.Success
    # @assert primalstatus(qpmodel) == MOI.FeasiblePoint
    @show terminationstatus(qpmodel) # FIXME
    @show primalstatus(qpmodel) # FIXME

    result.v̇ .= value.(qpmodel, controller.v̇)
    empty!(contactwrenches)
    for body in keys(controller.contactdata)
        contactwrench = zero(Wrench{Float64}, worldframe)
        for data in controller.contactdata[body]
            contactwrench += value(qpmodel, data.wrench_world)
        end
        contactwrenches[BodyID(body)] = contactwrench
    end

    inverse_dynamics!(τ, result.jointwrenches, result.accelerations, state, result.v̇, contactwrenches)
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

function regularize!(controller::MomentumBasedController, joint::Joint, weight)
    task = JointAccelerationTask(joint)
    setdesired!(task, zeros(num_velocities(joint)))
    addtask!(controller, task, weight)
    task
end

function addcontact!(controller::MomentumBasedController{N}, body::RigidBody{Float64}, point::ContactPoint) where N
    contactdata = ContactQPData{N}(point, controller.state, controller.qpmodel)
    push!(get!(Vector{ContactQPData{N}}, controller.contactdata, body), contactdata)
    objterm = objectiveterm(contactdata, controller.qpmodel)
    controller.objective = @expression controller.objective + objterm # TODO: currently kind of inefficient; would be better to have a single multi-arg addition
    contactdata
end

@noinline function initialize!(controller::MomentumBasedController)
    setobjective!(controller)
    mechanism = controller.state.mechanism
    if any(isfloating, out_joints(root_body(mechanism), mechanism))
        add_wrench_balance_constraint!(controller)
    end
end

function setobjective!(controller::MomentumBasedController)
    @objective(controller.qpmodel, Minimize, controller.objective)
end

function add_wrench_balance_constraint!(controller::MomentumBasedController{N}) where N
    # TODO: repeated computation of A if there's a MomentumRateTask
    qpmodel = controller.qpmodel
    state = controller.state
    mechanism = state.mechanism
    fg = mass(mechanism) * mechanism.gravitational_acceleration
    v̇ = controller.v̇
    A = Parameter(A -> momentum_matrix!(A, state), controller.momentum_matrix, qpmodel)
    Ȧv = Parameter{Wrench{Float64}}(() -> momentum_rate_bias(state), qpmodel)
    Wg = Parameter{Wrench{Float64}}(() -> Wrench(center_of_mass(state) × fg, fg), qpmodel)

    torque = @expression angular(Wg)
    force = @expression linear(Wg)
    for contact_data_vec in values(controller.contactdata)
        for contactdata in contact_data_vec
            torque = @expression torque + angular(contactdata.wrench_world)
            force = @expression force + linear(contactdata.wrench_world)
        end
    end

    @constraint(qpmodel, angular(A) * v̇ + angular(Ȧv) == torque)
    @constraint(qpmodel, linear(A) * v̇ + linear(Ȧv) == force)

    nothing
end
