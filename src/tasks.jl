abstract type MotionTask end

const ContiguousRowView{T} = SubArray{T,2,Array{T,2},Tuple{UnitRange{Int64},Base.Slice{Base.OneTo{Int64}}},false}

struct SpatialAccelerationTask <: MotionTask
    path::TreePath{RigidBody{Float64}, Joint{Float64}}
    Jmatrix::Matrix{Float64}
    J::GeometricJacobian{ContiguousRowView{Float64}}
    bias::Vector{Float64}
    desired::Vector{Float64}

    function SpatialAccelerationTask(
                mechanism::Mechanism,
                path::TreePath{RigidBody{Float64}, Joint{Float64}};
                frame::CartesianFrame3D = default_frame(target(path)))
        nv = num_velocities(mechanism)
        Jmatrix = zeros(6, nv)
        bodyframe = default_frame(target(path))
        baseframe = default_frame(source(path))
        J = GeometricJacobian(bodyframe, baseframe, frame, view(Jmatrix, 1 : 3, :), view(Jmatrix, 4 : 6, :))
        bias = zeros(6)
        desired = zeros(6)
        new(path, Jmatrix, J, bias, desired)
    end
end

dimension(::SpatialAccelerationTask) = 6

function set_desired!(task::SpatialAccelerationTask, desired::SpatialAcceleration)
    @framecheck task.J.body desired.body
    @framecheck task.J.base desired.base
    @framecheck task.J.frame desired.frame
    copyto!(task.desired, 1, angular(desired), 1, 3)
    copyto!(task.desired, 4, linear(desired), 1, 3)
    nothing
end

function update!(task::SpatialAccelerationTask, state::MechanismState)
    path = task.path
    J = task.J
    taskframe = J.frame
    world_to_task = inv(transform_to_root(state, taskframe))
    geometric_jacobian!(J, state, path, world_to_task)
    J̇v = transform(state, -bias_acceleration(state, source(path)) + bias_acceleration(state, target(path)), taskframe)
    copyto!(task.bias, 1, angular(J̇v), 1, 3)
    copyto!(task.bias, 4, linear(J̇v), 1, 3)
    nothing
end

function task_error(task::SpatialAccelerationTask, v̇)
    LinearTerm(task.Jmatrix, v̇) + Constant(task.bias) - Constant(task.desired)
end


struct AngularAccelerationTask <: MotionTask
    path::TreePath{RigidBody{Float64}, Joint{Float64}}
    J::GeometricJacobian{Matrix{Float64}}
    bias::Vector{Float64}
    desired::Vector{Float64}

    function AngularAccelerationTask(
                mechanism::Mechanism,
                path::TreePath{RigidBody{Float64}, Joint{Float64}};
                frame::CartesianFrame3D = default_frame(target(path)))
        nv = num_velocities(mechanism)
        bodyframe = default_frame(target(path))
        baseframe = default_frame(source(path))
        J = GeometricJacobian(bodyframe, baseframe, frame, zeros(3, nv), zeros(3, nv))
        bias = zeros(3)
        desired = zeros(3)
        new(path, J, bias, desired)
    end
end

dimension(::AngularAccelerationTask) = 3

function set_desired!(task::AngularAccelerationTask, desired::FreeVector3D)
    @framecheck task.J.frame desired.frame
    copyto!(task.desired, desired.v)
    nothing
end

function update!(task::AngularAccelerationTask, state::MechanismState)
    path = task.path
    J = task.J
    taskframe = J.frame
    world_to_task = inv(transform_to_root(state, taskframe))
    geometric_jacobian!(J, state, path, world_to_task)
    J̇v = transform(state, -bias_acceleration(state, source(path)) + bias_acceleration(state, target(path)), taskframe)
    copyto!(task.bias, angular(J̇v))
    nothing
end

function task_error(task::AngularAccelerationTask, v̇)
    LinearTerm(angular(task.J), v̇) + Constant(task.bias) - Constant(task.desired)
end


struct JointAccelerationTask{JT<:JointType{Float64}} <: MotionTask
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


struct LinearMomentumRateTask
    momentummatrix::MomentumMatrix{Matrix{Float64}}
    desired::Base.RefValue{FreeVector3D{SVector{3, Float64}}}
    biased_rate::FreeVector3D{Vector{Float64}}

    function LinearMomentumRateTask(mechanism::Mechanism, centroidalframe::CartesianFrame3D)
        nv = num_velocities(mechanism)
        momentummatrix = MomentumMatrix(centroidalframe, zeros(3, nv), zeros(3, nv))
        desired = Ref(FreeVector3D(zeros(SVector{3, Float64})))
        biased_rate = FreeVector3D(zeros(3))
        new(momentummatrix, desired, biased_rate)
    end
end

dimension(::LinearMomentumRateTask) = 3

function set_desired!(task::LinearMomentumRateTask, desired::FreeVector3D)
    @framecheck task.desired[].frame desired.frame
    task.desired[] = desired
    nothing
end

function update!(task::LinearMomentumRateTask, state::MechanismState)
    A = task.momentummatrix
    centroidalframe = A.frame
    com = center_of_mass(state)
    centroidal_to_world = Transform3D(centroidalframe, com.frame, com.v)
    world_to_centroidal = inv(centroidal_to_world)
    momentum_matrix!(A, state, world_to_centroidal)
    Ȧv = transform(momentum_rate_bias(state), world_to_centroidal)
    @framecheck Ȧv.frame task.biased_rate.frame
    task.biased_rate.v .= task.desired[].v .- Ȧv.linear
    nothing
end

function task_error(task::LinearMomentumRateTask, v̇)
    LinearTerm(task.momentummatrix.linear, v̇) - Constant(task.biased_rate.v)
end


const SparseSymmetric64 = Symmetric{Float64,SparseMatrixCSC{Float64,Int}}

struct Weighted{T<:MotionTask} <: MotionTask
    task::T
    weight::SparseSymmetric64
end

dimension(weighted::Weighted) = dimension(weighted.task)
set_desired!(weighted::Weighted, desired) = set_desired!(weighted.task, desired)
update!(weighted::Weighted, state::MechanismState) = update!(weighted.task, state)
task_error(weighted::Weighted, v̇) = task_error(weighted.task, v̇)

disable!(weighted::Weighted) = weighted.weight.data[:] = 0.0

