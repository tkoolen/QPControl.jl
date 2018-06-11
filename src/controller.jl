# TODO:
# * max normal force
# * split up motion tasks into hard constraints and objective terms:
#   * constraints: just pass in a Vector{<:AbstractMotionTask}
#   * objective terms: Dict{<:AbstractMotionTask, SparseSymmetric64}
# * inner constructor: take TypeSortedCollection `motiontasks` (all tasks) +
#   Dict{<:AbstractMotionTask, SparseSymmetric64} for weights
# * convert motion tasks to TypeSortedCollection in outer constructor

struct MomentumBasedController{O<:MOI.AbstractOptimizer, S<:MechanismState, N, M}
    # dynamics-related
    state::S
    result::DynamicsResult{Float64, Float64}
    centroidalframe::CartesianFrame3D
    momentummatrix::MomentumMatrix{Matrix{Float64}}
    externalwrenches::Dict{RigidBody{Float64}, Wrench{Float64}}

    # contact info and motion tasks
    contactsettings::Vector{ContactSettings{N}}
    motiontasks::M

    # qp-related
    qpmodel::SimpleQP.Model{O}

    # TODO: remove:
    Ȧv_angular::Vector{Float64}
    Ȧv_linear::Vector{Float64}
    external_torque_sum::Vector{Float64}
    external_force_sum::Vector{Float64}

    function MomentumBasedController{O, N, M}(
            mechanism::Mechanism{Float64}, optimizer::O,
            contactsettings::Vector{ContactSettings{N}}, motiontasks::M,
            taskweights::Dict{<:AbstractMotionTask, SparseSymmetric64}) where {O<:MOI.AbstractOptimizer, N, M}
        nv = num_velocities(mechanism)

        # dynamics-related
        centroidalframe = CartesianFrame3D("centroidal")
        result = DynamicsResult(mechanism)
        totalmass = mass(mechanism)

        # contact info and motion tasks
        rootframe = root_frame(mechanism)
        externalwrenches = Dict(b => zero(Wrench{Float64}, rootframe) for b in bodies(mechanism))

        # qp-related
        qpmodel = SimpleQP.Model(optimizer)
        Ȧv_angular = zeros(3)
        Ȧv_linear = zeros(3)
        # ρ = Vector{Float64}()

        controller = new{O, N, M}(
            mechanism, centroidalframe, totalmass, result, momentummatrix, externalwrenches,
            contactsettings, motiontasks,
            qpmodel, Ȧv_angular, Ȧv_linear)
        set_up_qp!(controller, taskweights)
        controller
    end
end

function MomentumBasedController(
        mechanism::Mechanism,
        optimizer::O,
        contactsettings::Vector{ContactSettings{N}},
        motionconstraints::Vector{<:AbstractMotionTask},
        motionobjectiveterms::Dict{<:AbstractMotionTask, SparseSymmetric64}) where {O<:MOI.AbstractOptimizer, N}
    state = MechanismState(mechanism)
    motiontasks = TypeSortedCollection(vcat(motionconstraints, collect(keys(motionobjectiveterms))))
    MomentumBasedController{O, N, typeof(motiontasks)}(mechanism, optimizer, contactsettings, motiontasks, motionobjectiveterms)
end

centroidal_frame(controller::MomentumBasedController) = controller.centroidalframe

function set_up_qp!(controller::MomentumBasedController, taskweights::Dict{<:AbstractMotionTask, SparseSymmetric64})
    model = controller.qpmodel
    v̇ = [Variable(model) for _ = 1 : num_velocities(controller.mechanism)]

    # Newton-Euler
    # A = controller.momentummatrix


    # Motion tasks
    foreach(controller.motiontasks) do task
        taskdim = dimension(task)
        if haskey(taskweights, task)
            @constraint(model, task_error(task, v̇) == Constant(zeros(taskdim)))
        else
            e = [Variable(model) for _ = 1 : taskdim]
            @constraint(model, LinearTerm(eye(taskdim), e) == task_error(task, v̇)) # FIXME: no sparsity
            # FIXME: handle weight
        end
    end

    # Contacts


    nothing
end

