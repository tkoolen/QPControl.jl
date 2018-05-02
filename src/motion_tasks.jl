abstract type MotionTask end
abstract type MutableMotionTask <: MotionTask end # TODO: get rid of this distinction; make all MotionTasks mutable


# SpatialAccelerationTask
mutable struct SpatialAccelerationTask <: MutableMotionTask
    path::TreePath{RigidBody{Float64}, Joint{Float64}} # TODO
    jacobian::GeometricJacobian{Matrix{Float64}}
    desired::SpatialAcceleration{Float64}
    angularselectionmatrix::Matrix{Float64}
    linearselectionmatrix::Matrix{Float64}
    weight::Float64

    function SpatialAccelerationTask(nv::Int, path::TreePath{RigidBody{Float64}, Joint{Float64}}, frame::CartesianFrame3D,
            angularselectionmatrix::Matrix{Float64}, linearselectionmatrix::Matrix{Float64})
        @assert size(angularselectionmatrix, 2) == 3
        @assert size(linearselectionmatrix, 2) == 3
        bodyframe = default_frame(target(path))
        baseframe = default_frame(source(path))
        jacobian = GeometricJacobian(bodyframe, baseframe, frame, zeros(3, nv), zeros(3, nv))
        desired = zero(SpatialAcceleration{Float64}, bodyframe, baseframe, frame)
        weight = 0.0
        new(path, jacobian, desired, angularselectionmatrix, linearselectionmatrix, weight)
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
mutable struct JointAccelerationTask <: MutableMotionTask
    joint::Joint{Float64}
    desired::Vector{Float64}
    weight::Float64

    JointAccelerationTask(joint::Joint{Float64}) = new(joint, zeros(num_velocities(joint)), 0.0)
end

function set!(task::JointAccelerationTask, desired::Union{Number, AbstractVector}, weight::Number)
    @boundscheck weight >= 0 || error("Weight must be nonnegative.")
    task.desired .= desired
    task.weight = weight
end


# MomentumRateTask
struct MomentumRateTask <: MotionTask
    desired::Wrench{Float64}
    angularselectionmatrix::SMatrix{3, 3, Float64, 9}
    linearselectionmatrix::SMatrix{3, 3, Float64, 9}
    weight::Float64
end

MomentumRateTask(frame::CartesianFrame3D) = MomentumRateTask(zero(Wrench{Float64}, frame), eye(SMatrix{3, 3}), eye(SMatrix{3, 3}), 0.0)

disable!(task::MutableMotionTask) = task.weight = 0
isenabled(task::MotionTask) = task.weight > 0
isconstraint(task::MotionTask) = task.weight == Inf
