mutable struct MomentumBasedController{N, O<:MOI.AbstractOptimizer, S<:MechanismState}
    state::S
    result::DynamicsResult{Float64, Float64}
    floatingjoint::Union{Nothing, Joint{Float64}}
    floating_joint_velocity_range::UnitRange{Int}
    centroidalframe::CartesianFrame3D
    momentum_matrix::MomentumMatrix{Matrix{Float64}}
    contacts::Dict{RigidBody{Float64}, Vector{ContactPoint{N}}}
    contactwrenches::Dict{BodyID, Wrench{Float64}}
    qpmodel::Parametron.Model{Float64, O}
    v̇::Vector{Variable}
    objective::Parametron.LazyExpression # buffer to incrementally build the objective function
    initialized::Bool

    function MomentumBasedController{N}(
            mechanism::Mechanism{Float64}, optimizer::O; floatingjoint=nothing) where {N, O<:MOI.AbstractOptimizer}
        state = MechanismState(mechanism)
        result = DynamicsResult(mechanism)
        floating_joint_velocity_range = floatingjoint == nothing ? (1 : 0) : velocity_range(state, floatingjoint)
        worldframe = root_frame(mechanism)
        centroidalframe = CartesianFrame3D("centroidal")
        nv = num_velocities(state)
        momentum_matrix = MomentumMatrix(worldframe, zeros(3, nv), zeros(3, nv))
        rootframe = root_frame(mechanism)
        contacts = Dict{RigidBody{Float64}, Vector{ContactPoint{N}}}()
        contactwrenches = Dict{BodyID, Wrench{Float64}}()
        qpmodel = Parametron.Model(optimizer)
        v̇ = [Variable(qpmodel) for _ = 1 : nv]
        objective = Parametron.LazyExpression(identity, zero(QuadraticFunction{Float64}))
        new{N, O, typeof(state)}(
            state, result, floatingjoint, floating_joint_velocity_range, centroidalframe,
            momentum_matrix, contacts, contactwrenches, qpmodel, v̇, objective, false)
    end
end

Base.show(io::IO, controller::MomentumBasedController{N, O, S}) where {N, O, S} =
    print(io, "MomentumBasedController{$N, $O, $S}(…)")

centroidal_frame(controller::MomentumBasedController) = controller.centroidalframe

function (controller::MomentumBasedController)(τ::AbstractVector, t::Number, x::Union{<:Vector, <:MechanismState})
    # initialize on first solve
    if !controller.initialized
        initialize!(controller)
        controller.initialized = true
    end

    # unpack
    qpmodel = controller.qpmodel
    state = controller.state
    result = controller.result
    contacts = controller.contacts
    contactwrenches = controller.contactwrenches
    worldframe = root_frame(state.mechanism)

    # set up and solve qp
    copyto!(state, x)
    solve!(qpmodel)
    checkstatus(qpmodel)

    # extract joint accelerations and contact wrenches
    @inbounds for i in eachindex(controller.v̇)
        result.v̇[i] = value(qpmodel, controller.v̇[i])
    end
    empty!(contactwrenches)
    for body in keys(controller.contacts)
        contactwrench = zero(Wrench{Float64}, worldframe)
        for data in controller.contacts[body]
            contactwrench += value(qpmodel, data.wrench_world)
        end
        contactwrenches[BodyID(body)] = contactwrench
    end

    # compute torques
    inverse_dynamics!(τ, result.jointwrenches, result.accelerations, state, result.v̇, contactwrenches)

    # make floating joint torques precisely zero (they are approximately zero already)
    zero_floating_joint_torques!(τ, state, controller.floating_joint_velocity_range)

    τ
end

function checkstatus(qpmodel::Parametron.Model)
    ok = terminationstatus(qpmodel) == MOI.Success && primalstatus(qpmodel) == MOI.FeasiblePoint
    ok = ok || (primalstatus(qpmodel) == MOI.FeasiblePoint && dualstatus(qpmodel) == MOI.FeasiblePoint)
    if !ok
        okish = terminationstatus(qpmodel) == MOI.AlmostSuccess && primalstatus(qpmodel) == MOI.UnknownResultStatus
        if !okish
            throw(QPSolveFailure(terminationstatus(qpmodel), primalstatus(qpmodel), dualstatus(qpmodel)))
        end
    end
end

function zero_floating_joint_torques!(τ::AbstractVector, state::MechanismState, velocity_range::UnitRange)
    for i in velocity_range
        τ[i] = 0
    end
end

function addtask!(controller::MomentumBasedController, task::AbstractMotionTask)
    model = controller.qpmodel
    err = task_error(task, model, controller.state, controller.v̇)
    taskdim = dimension(task)
    @constraint(model, err == zeros(taskdim))
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
    v̇joint = controller.v̇[velocity_range(controller.state, joint)]
    controller.objective = @expression controller.objective + weight * (v̇joint ⋅ v̇joint)
end

function addcontact!(
        controller::MomentumBasedController{N}, body::RigidBody{Float64},
        point::ContactPoint{N}) where N
    push!(get!(Vector{ContactPoint{N}}, controller.contacts, body), point)
    objterm = objectiveterm(point, controller.qpmodel)
    controller.objective = @expression controller.objective + objterm # TODO: currently kind of inefficient; would be better to have a single multi-arg addition
    point
end

function addcontact!(
        controller::MomentumBasedController{N}, body::RigidBody{Float64},
        position::Union{Point3D, Parameter{<:Point3D}},
        normal::Union{FreeVector3D, Parameter{<:FreeVector3D}},
        μ::Union{Float64, Parameter{Float64}}) where N
    addcontact!(controller, body, ContactPoint{N}(position, normal, μ, controller.state, controller.qpmodel))
end

@noinline function initialize!(controller::MomentumBasedController)
    setobjective!(controller)
    mechanism = controller.state.mechanism
    if controller.floatingjoint !== nothing
        add_wrench_balance_constraint!(controller, controller.floatingjoint)
    end
end

function setobjective!(controller::MomentumBasedController)
    @objective(controller.qpmodel, Minimize, controller.objective)
end

function add_wrench_balance_constraint!(controller::MomentumBasedController{N}, joint::Joint) where N
    qpmodel = controller.qpmodel
    state = controller.state
    mechanism = state.mechanism
    fg = mass(mechanism) * mechanism.gravitational_acceleration
    v̇ = controller.v̇
    A = Parameter(A -> momentum_matrix!(A, state), controller.momentum_matrix, qpmodel)
    Ȧv = Parameter{Wrench{Float64}}(() -> momentum_rate_bias(state), qpmodel)
    Wg = Parameter{Wrench{Float64}}(() -> Wrench(center_of_mass(state) × fg, fg), qpmodel)
    S = let state = state, joint = joint, floatingbody = successor(joint, mechanism)
        Parameter(qpmodel) do
            qjoint = configuration(state, joint)
            H = transform_to_root(state, floatingbody)
            transform(motion_subspace(joint, qjoint), H)
        end
    end


    torque = @expression angular(Wg)
    force = @expression linear(Wg)
    for contactvec in values(controller.contacts)
        for contact in contactvec
            torque = @expression torque + angular(contact.wrench_world)
            force = @expression force + linear(contact.wrench_world)
        end
    end

    nv = num_velocities(joint)
    @constraint(qpmodel, transpose(angular(S)) * (angular(A) * v̇ + angular(Ȧv) - torque) + transpose(linear(S)) * (linear(A) * v̇ + linear(Ȧv) - force) == zeros(nv))

    nothing
end
