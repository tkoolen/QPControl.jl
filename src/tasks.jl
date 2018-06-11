abstract type AbstractMotionTask end

struct SpatialAccelerationTask <: AbstractMotionTask
    path::TreePath{RigidBody{Float64}, Joint{Float64}}
    jacobian::GeometricJacobian{Matrix{Float64}}
    desired::Base.RefValue{SpatialAcceleration{Float64}}

    function SpatialAccelerationTask(
                mechanism::Mechanism, # TODO: would be nice to get rid of this; possible if compact Jacobians were available
                path::TreePath{RigidBody{Float64}, Joint{Float64}};
                frame::CartesianFrame3D = default_frame(target(path)))
        nv = num_velocities(mechanism)
        bodyframe = default_frame(target(path))
        baseframe = default_frame(source(path))
        jacobian = GeometricJacobian(bodyframe, baseframe, frame, zeros(3, nv), zeros(3, nv))
        desired = Ref{SpatialAcceleration{Float64}}(zero(SpatialAcceleration{Float64}, bodyframe, baseframe, frame))
        new(path, jacobian, desired)
    end
end

dimension(task::SpatialAccelerationTask) = 6

function set_desired!(task::SpatialAccelerationTask, desired::SpatialAcceleration)
    @framecheck task.desired.body desired.body
    @framecheck task.desired.base desired.base
    @framecheck task.desired.frame desired.frame
    task.desired[] = desired
    nothing
end

function task_error(task::SpatialAccelerationTask, model::SimpleQP.Model, state::MechanismState, v̇::AbstractVector{SimpleQP.Variable})
    J = let task = task, state = state
        Parameter(model) do
            world_to_desired = inv(transform_to_root(state, task.desired[].frame))
            geometric_jacobian!(task.jacobian, state, task.path, world_to_desired)
        end
    end
    J̇v = let task = task, state = state
        Parameter(model) do
            world_to_desired = inv(transform_to_root(state, task.desired[].frame))
            -bias_acceleration(state, source(task.path)) + bias_acceleration(state, target(task.path))
        end
    end
    desired = Parameter(@closure(() -> task.desired[]), model)
    @expression [
        angular(desired) - (angular(J) * v̇ + angular(J̇v));
        linear(desired) - (linear(J) * v̇ + linear(J̇v))]
end


struct JointAccelerationTask{JT<:JointType{Float64}} <: AbstractMotionTask
    joint::Joint{Float64, JT}
    desired::Vector{Float64}

    function JointAccelerationTask(joint::Joint{Float64, JT}) where {JT<:JointType{Float64}}
        new{JT}(joint, zeros(num_velocities(joint)))
    end
end

dimension(task::JointAccelerationTask) = length(task.desired)
set_desired!(task::JointAccelerationTask, desired) = set_velocity!(task.desired, task.joint, desired)

function task_error(task::JointAccelerationTask, model::SimpleQP.Model, state::MechanismState, v̇::AbstractVector{SimpleQP.Variable})
    desired = Parameter(@closure(() -> task.desired), model)
    v̇joint = v̇[velocity_range(state, task.joint)]
    @expression desired - v̇joint
end

function momentum_rate_task_params(task, model, state, v̇)
    world_to_centroidal = let task = task, state = state
        Parameter(model) do
            A = task.momentum_matrix
            centroidalframe = A.frame
            com = center_of_mass(state)
            centroidal_to_world = Transform3D(centroidalframe, com.frame, com.v)
            inv(centroidal_to_world)
        end
    end
    A = Parameter(@closure(() -> momentum_matrix!(task.momentum_matrix, state, world_to_centroidal())), model)
    Ȧv = Parameter(@closure(() -> transform(momentum_rate_bias(state), world_to_centroidal())), model)
    A, Ȧv
end

struct MomentumRateTask
    momentum_matrix::MomentumMatrix{Matrix{Float64}}
    desired::Base.RefValue{Wrench{Float64}}

    function MomentumRateTask(mechanism::Mechanism, centroidalframe::CartesianFrame3D)
        nv = num_velocities(mechanism)
        momentum_matrix = MomentumMatrix(centroidalframe, zeros(3, nv), zeros(3, nv))
        desired = Ref(zero(Wrench{Float64}, centroidalframe))
        new(momentum_matrix, desired)
    end
end

dimension(task::MomentumRateTask) = 6

function set_desired!(task::MomentumRateTask, desired::Wrench)
    @framecheck task.momentum_matrix.frame desired.frame
    task.desired[] = desired
end

function task_error(task::MomentumRateTask, model::SimpleQP.Model, state::MechanismState, v̇::AbstractVector{SimpleQP.Variable})
    A, Ȧv = momentum_rate_task_params(task, model, state, v̇)
    desired = Parameter(@closure(() -> task.desired[]), model)
    @expression [
        angular(desired) - (angular(A) * v̇ + angular(Ȧv));
        linear(desired) - (linear(A) * v̇ + linear(Ȧv))]
end

struct LinearMomentumRateTask
    momentum_matrix::MomentumMatrix{Matrix{Float64}}
    desired::Base.RefValue{FreeVector3D{SVector{3, Float64}}}

    function LinearMomentumRateTask(mechanism::Mechanism, centroidalframe::CartesianFrame3D = CartesianFrame3D())
        nv = num_velocities(mechanism)
        momentum_matrix = MomentumMatrix(centroidalframe, zeros(3, nv), zeros(3, nv))
        desired = Ref(zero(FreeVector3D{SVector{3, Float64}}, centroidalframe))
        new(momentum_matrix, desired)
    end
end

dimension(task::LinearMomentumRateTask) = 3

function set_desired!(task::LinearMomentumRateTask, desired::FreeVector3D)
    @framecheck task.momentum_matrix.frame desired.frame
    task.desired[] = desired
end

function task_error(task::LinearMomentumRateTask, model::SimpleQP.Model, state::MechanismState, v̇::AbstractVector{SimpleQP.Variable})
    A, Ȧv = momentum_rate_task_params(task, model, state, v̇)
    desired = Parameter(@closure(() -> task.desired[].v), model)
    @expression desired - (angular(A) * v̇ + angular(Ȧv))
end
