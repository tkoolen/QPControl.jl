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
    @framecheck task.desired[].body desired.body
    @framecheck task.desired[].base desired.base
    @framecheck task.desired[].frame desired.frame
    task.desired[] = desired
    nothing
end

function task_error(task::SpatialAccelerationTask, qpmodel, state::MechanismState, v̇::AbstractVector{SimpleQP.Variable})
    J = Parameter(task.jacobian, qpmodel) do jac
        world_to_desired = inv(transform_to_root(state, task.desired[].frame))
        geometric_jacobian!(jac, state, task.path, world_to_desired)
    end
    J̇v = Parameter{SpatialAcceleration{Float64}}(qpmodel) do
        bias = -bias_acceleration(state, source(task.path)) + bias_acceleration(state, target(task.path))
        transform(state, bias, task.desired[].frame)
    end
    desired = Parameter{SpatialAcceleration{Float64}}(() -> task.desired[], qpmodel)
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

function task_error(task::JointAccelerationTask, qpmodel, state::MechanismState, v̇::AbstractVector{SimpleQP.Variable})
    desired = Parameter(() -> task.desired, qpmodel)
    v̇joint = v̇[velocity_range(state, task.joint)]
    @expression desired - v̇joint
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

function momentum_rate_task_params(task, qpmodel, state, v̇)
    # TODO: repeated computation of world_to_centroidal, but running into inference issues if that computation
    # is extracted out into its own Parameter
    centroidalframe = task.momentum_matrix.frame
    A = Parameter(task.momentum_matrix, qpmodel) do A
        com = center_of_mass(state)
        centroidal_to_world = Transform3D(centroidalframe, com.frame, com.v)
        world_to_centroidal = inv(centroidal_to_world)
        momentum_matrix!(A, state, world_to_centroidal)
    end
    Ȧv = Parameter{Wrench{Float64}}(qpmodel) do
        com = center_of_mass(state)
        centroidal_to_world = Transform3D(centroidalframe, com.frame, com.v)
        world_to_centroidal = inv(centroidal_to_world)
        transform(momentum_rate_bias(state), world_to_centroidal)
    end
    A, Ȧv
end

dimension(task::MomentumRateTask) = 6

function set_desired!(task::MomentumRateTask, desired::Wrench)
    @framecheck task.momentum_matrix.frame desired.frame
    task.desired[] = desired
end

function task_error(task::MomentumRateTask, qpmodel, state::MechanismState, v̇::AbstractVector{SimpleQP.Variable})
    A, Ȧv = momentum_rate_task_params(task, qpmodel, state, v̇)
    desired = Parameter{Wrench{Float64}}(() -> task.desired[], qpmodel)
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

function task_error(task::LinearMomentumRateTask, qpmodel, state::MechanismState, v̇::AbstractVector{SimpleQP.Variable})
    A, Ȧv = momentum_rate_task_params(task, qpmodel, state, v̇)
    desired = Parameter(@closure(() -> task.desired[].v), qpmodel)
    @expression desired - (angular(A) * v̇ + angular(Ȧv))
end
