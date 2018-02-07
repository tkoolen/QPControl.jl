abstract type MotionTask{T} end
abstract type MutableMotionTask{T} <: MotionTask{T} end


# SpatialAccelerationTask
mutable struct SpatialAccelerationTask{T} <: MutableMotionTask{T}
    path::TreePath{RigidBody{T}, GenericJoint{T}}
    jacobian::GeometricJacobian{Matrix{T}}
    desired::SpatialAcceleration{T}
    angularselectionmatrix::Matrix{T}
    linearselectionmatrix::Matrix{T}
    weight::T

    function SpatialAccelerationTask(nv::Int, path::TreePath{RigidBody{T}, GenericJoint{T}}, frame::CartesianFrame3D,
            angularselectionmatrix::Matrix{T}, linearselectionmatrix::Matrix{T}) where {T}
        @assert size(angularselectionmatrix, 2) == 3
        @assert size(linearselectionmatrix, 2) == 3
        bodyframe = default_frame(target(path))
        baseframe = default_frame(source(path))
        jacobian = GeometricJacobian(bodyframe, baseframe, frame, Matrix{T}(3, nv), Matrix{T}(3, nv))
        desired = zero(SpatialAcceleration{T}, bodyframe, baseframe, frame)
        weight = zero(T)
        new{T}(path, jacobian, desired, angularselectionmatrix, linearselectionmatrix, weight)
    end
end

function set!(task::SpatialAccelerationTask, desired::SpatialAcceleration, weight::Number)
    @framecheck task.desired.body desired.body
    @framecheck task.desired.base desired.base
    @framecheck task.desired.frame desired.frame
    @boundscheck weight >= 0 || error("Weight must be nonnegative.")
    task.desired = desired
    task.weight = weight
end

RigidBodyDynamics.zero!(task::SpatialAccelerationTask, weight::Number) = set!(task, zero(task.desired), weight)

# JointAccelerationTask
mutable struct JointAccelerationTask{T} <: MutableMotionTask{T}
    joint::GenericJoint{T}
    desired::Vector{T}
    weight::T

    JointAccelerationTask(joint::GenericJoint{T}) where {T} = new{T}(joint, zeros(T, num_velocities(joint)), zero(T))
end

function set!(task::JointAccelerationTask, desired::Union{Number, AbstractVector}, weight::Number)
    @boundscheck weight >= 0 || error("Weight must be nonnegative.")
    task.desired .= desired
    task.weight = weight
end


# MomentumRateTask
struct MomentumRateTask{T} <: MotionTask{T}
    desired::Wrench{T}
    angularselectionmatrix::SMatrix{3, 3, T, 9}
    linearselectionmatrix::SMatrix{3, 3, T, 9}
    weight::T
end

MomentumRateTask(::Type{T}, frame::CartesianFrame3D) where {T} = MomentumRateTask(zero(Wrench{T}, frame), eye(SMatrix{3, 3}), eye(SMatrix{3, 3}), T(0))

disable!(task::MutableMotionTask) = task.weight = 0
isenabled(task::MotionTask) = task.weight > 0
isconstraint(task::MotionTask) = task.weight == Inf
