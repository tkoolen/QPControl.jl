abstract type AbstractMotionTask end

function set_task_vector!(dest::Vector, src::Union{<:SpatialAcceleration, <:Wrench}, task)
    angularrows = task.angularrows
    linearrows = task.linearrows
    copyto!(dest, 1, angular(src), first(angularrows), last(angularrows))
    copyto!(dest, length(angularrows) + 1, linear(src), first(linearrows), last(linearrows))
end

function set_task_matrix!(dest::Matrix, src::Union{<:GeometricJacobian, <:MomentumMatrix}, task)
    srcang = angular(src)
    srclin = linear(src)
    angularrows = task.angularrows
    linearrows = task.linearrows
    destindex = 1
    @inbounds for col = 1 : size(dest, 2)
        for srcrow in angularrows
            dest[destindex] = srcang[srcrow, col]
            destindex += 1
        end
        for srcrow in linearrows
            dest[destindex] = srclin[srcrow, col]
            destindex += 1
        end
    end
end

struct SpatialAccelerationTask <: AbstractMotionTask
    path::TreePath{RigidBody{Float64}, Joint{Float64}}
    jacobian::GeometricJacobian{Matrix{Float64}}
    angularrows::UnitRange{Int}
    linearrows::UnitRange{Int}
    coeffmatrix::Matrix{Float64}
    bias::Vector{Float64}
    desired::Vector{Float64}

    function SpatialAccelerationTask(
                mechanism::Mechanism,
                path::TreePath{RigidBody{Float64}, Joint{Float64}};
                frame::CartesianFrame3D = default_frame(target(path)),
                angularrows::UnitRange{Int} = 1 : 3,
                linearrows::UnitRange{Int} = 1 : 3)
        nv = num_velocities(mechanism)
        bodyframe = default_frame(target(path))
        baseframe = default_frame(source(path))
        jacobian = GeometricJacobian(bodyframe, baseframe, frame, zeros(3, nv), zeros(3, nv))
        taskdim = length(angularrows) + length(linearrows)
        coeffmatrix = zeros(taskdim, nv)
        bias = zeros(taskdim)
        desired = zeros(taskdim)
        new(path, jacobian, angularrows, linearrows, coeffmatrix, bias, desired)
    end
end

dimension(task::SpatialAccelerationTask) = size(task.coeffmatrix, 1)

function set_desired!(task::SpatialAccelerationTask, desired::SpatialAcceleration)
    @framecheck task.jacobian.body desired.body
    @framecheck task.jacobian.base desired.base
    @framecheck task.jacobian.frame desired.frame
    set_task_vector!(task.desired, desired, task)
    nothing
end

function set_desired!(task::SpatialAccelerationTask, desired::AbstractVector{Float64})
    copyto!(task.desired, desired)
    nothing
end

function update!(task::SpatialAccelerationTask, state::MechanismState)
    path = task.path
    jacobian = task.jacobian
    taskframe = jacobian.frame
    J̇v = transform(state, -bias_acceleration(state, source(path)) + bias_acceleration(state, target(path)), taskframe)
    set_task_vector!(task.bias, J̇v, task)
    world_to_task = inv(transform_to_root(state, taskframe))
    geometric_jacobian!(jacobian, state, path, world_to_task)
    set_task_matrix!(task.coeffmatrix, task.jacobian, task)
    nothing
end

function task_error(task::SpatialAccelerationTask, v̇)
    LinearTerm(task.coeffmatrix, v̇) + Constant(task.bias) - Constant(task.desired)
end


struct JointAccelerationTask{JT<:JointType{Float64}} <: AbstractMotionTask
    joint::Joint{Float64, JT}
    v̇range::UnitRange{Int}
    desired::Vector{Float64}

    function JointAccelerationTask(mechanism::Mechanism, joint::Joint{Float64, JT}) where {JT<:JointType{Float64}}
        v̇range = velocity_range(MechanismState(mechanism), joint) # FIXME
        new{JT}(joint, v̇range, zeros(num_velocities(joint)))
    end
end

dimension(task::JointAccelerationTask) = length(task.desired)

function set_desired!(task::JointAccelerationTask, desired)
    set_velocity!(task.desired, task.joint, desired)
    nothing
end

update!(task::JointAccelerationTask, state::MechanismState) = nothing

function task_error(task::JointAccelerationTask, v̇)
    nv = num_velocities(task.joint)
    LinearTerm(eye(nv), v̇[task.v̇range]) - Constant(task.desired)
end


struct CentroidalMomentumRateTask
    momentummatrix::MomentumMatrix{Matrix{Float64}}
    angularrows::UnitRange{Int}
    linearrows::UnitRange{Int}
    coeffmatrix::Matrix{Float64}
    bias::Vector{Float64}
    desired::Vector{Float64}

    function CentroidalMomentumRateTask(
            mechanism::Mechanism,
            centroidalframe::CartesianFrame3D = CartesianFrame3D(),
            angularrows::UnitRange{Int} = 1 : 3,
            linearrows::UnitRange{Int} = 1 : 3)
        nv = num_velocities(mechanism)
        momentummatrix = MomentumMatrix(centroidalframe, zeros(3, nv), zeros(3, nv))
        taskdim = length(angularrows) + length(linearrows)
        coeffmatrix = zeros(taskdim, nv)
        bias = zeros(taskdim)
        desired = zeros(taskdim)
        new(momentummatrix, angularrows, linearrows, coeffmatrix, bias, desired)
    end
end

dimension(task::CentroidalMomentumRateTask) = size(task.coeffmatrix, 2)

function set_desired!(task::CentroidalMomentumRateTask, desired::Wrench)
    @framecheck task.momentummatrix.frame desired.frame
    set_task_vector!(task.desired, desired, task)
end

function set_desired!(task::CentroidalMomentumRateTask, desired::AbstractVector{Float64})
    copyto!(task.desired, desired)
end

function update!(task::CentroidalMomentumRateTask, state::MechanismState)
    A = task.momentummatrix
    centroidalframe = A.frame
    com = center_of_mass(state)
    centroidal_to_world = Transform3D(centroidalframe, com.frame, com.v)
    world_to_centroidal = inv(centroidal_to_world)
    Ȧv = transform(momentum_rate_bias(state), world_to_centroidal)
    set_task_vector!(task.bias, Ȧv, task)
    momentum_matrix!(A, state, world_to_centroidal)
    set_task_matrix!(task.coeffmatrix, A, task)
    nothing
end

function task_error(task::CentroidalMomentumRateTask, v̇)
    LinearTerm(task.coeffmatrix, v̇) + Constant(task.bias) - Constant(task.desired)
end


# TODO:
# const SparseSymmetric64 = Symmetric{Float64,SparseMatrixCSC{Float64,Int}}

# struct MotionTaskSpecification
#     objectiveterms::Vector{Pair{AbstractMotionTask, SparseSymmetric64}}
#     constraints::Vector{AbstractMotionTask}
#     MotionTaskSpecification() = new([], [])
# end

# Base.push!(spec::MotionTaskSpecification, task::AbstractMotionTask, weight::SparseSymmetric64) = push!(spec.objectiveterms, task => weight)
# Base.push!(spec::MotionTaskSpecification, task::AbstractMotionTask) = push!(spec.constraints, task)

# struct WeightedMotionTask{T<:AbstractMotionTask} <: AbstractMotionTask
#     task::T
#     weight::SparseSymmetric64
# end

# dimension(weighted::WeightedMotionTask) = dimension(weighted.task)
# set_desired!(weighted::WeightedMotionTask, desired) = set_desired!(weighted.task, desired)
# update!(weighted::WeightedMotionTask, state::MechanismState) = update!(weighted.task, state)
# task_error(weighted::WeightedMotionTask, v̇) = task_error(weighted.task, v̇)

# disable!(weighted::WeightedMotionTask) = weighted.weight.data[:] = 0.0

