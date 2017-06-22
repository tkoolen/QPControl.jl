@compat abstract type MotionTask{T} end
@compat abstract type MutableMotionTask{T} <: MotionTask{T} end


# SpatialAccelerationTask
type SpatialAccelerationTask{T} <: MutableMotionTask{T}
    path::TreePath{RigidBody{T}, Joint{T}}
    jacobian::GeometricJacobian{Matrix{T}}
    desired::SpatialAcceleration{T}
    angularselectionmatrix::Matrix{T}
    linearselectionmatrix::Matrix{T}
    weight::T

    function (::Type{SpatialAccelerationTask{T}}){T}(path::TreePath{RigidBody{T}, Joint{T}}, frame::CartesianFrame3D, angularselectionmatrix::Matrix{T}, linearselectionmatrix::Matrix{T})
        @assert size(angularselectionmatrix, 2) == 3
        @assert size(linearselectionmatrix, 2) == 3
        nv = num_velocities(path)
        bodyframe = default_frame(target(path))
        baseframe = default_frame(source(path))
        jacobian = GeometricJacobian(bodyframe, baseframe, frame, Matrix{T}(3, nv), Matrix{T}(3, nv))
        desired = zero(SpatialAcceleration{T}, bodyframe, baseframe, frame)
        weight = zero(T)
        new{T}(path, jacobian, desired, angularselectionmatrix, linearselectionmatrix, weight)
    end
end

function SpatialAccelerationTask{T}(path::TreePath{RigidBody{T}, Joint{T}}, frame::CartesianFrame3D, angularselectionmatrix::Matrix{T}, linearselectionmatrix::Matrix{T})
    SpatialAccelerationTask{T}(path, frame, angularselectionmatrix, linearselectionmatrix)
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
type JointAccelerationTask{T} <: MutableMotionTask{T}
    joint::Joint{T}
    desired::Vector{T}
    weight::T

    (::Type{JointAccelerationTask{T}}){T}(joint::Joint{T}) = new{T}(joint, zeros(T, num_velocities(joint)), zero(T))
end

JointAccelerationTask{T}(joint::Joint{T}) = JointAccelerationTask{T}(joint)

function set!(task::JointAccelerationTask, desired::Union{Number, AbstractVector}, weight::Number)
    @boundscheck weight >= 0 || error("Weight must be nonnegative.")
    task.desired .= desired
    task.weight = weight
end


# MomentumRateTask
type MomentumRateTask{T} <: MotionTask{T}
    desired::Wrench{T}
    angularselectionmatrix::SMatrix{3, 3, T, 9}
    linearselectionmatrix::SMatrix{3, 3, T, 9}
    weight::T
end

MomentumRateTask{T}(::Type{T}, frame::CartesianFrame3D) = MomentumRateTask(zero(Wrench{T}, frame), eye(SMatrix{3, 3}), eye(SMatrix{3, 3}), T(0))

disable!(task::MutableMotionTask) = task.weight = 0
isenabled(task::MotionTask) = task.weight > 0
isconstraint(task::MotionTask) = task.weight == Inf
