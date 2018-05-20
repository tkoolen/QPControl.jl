abstract type MotionTask end

struct SpatialAccelerationTask <: MotionTask
    desired::Vector{Float64}
    path::TreePath{RigidBody{Float64}, Joint{Float64}}
    J::GeometricJacobian{Matrix{Float64}}
    Sangular::Matrix{Float64}
    Slinear::Matrix{Float64}
    SJ::Matrix{Float64}
    SJangular::
    SJ̇v::Vector{Float64}

    function SpatialAccelerationTask(
            mechanism::Mechanism, path::TreePath{RigidBody{Float64}, Joint{Float64}},
            frame::CartesianFrame3D, Sangular::Matrix{Float64} = eye(3), Slinear::Matrix{Float64} = eye(3))
        #nv = mapreduce(num_velocities, +, 0, path)
        nv = num_velocities(mechanism)
        bodyframe = default_frame(target(path))
        baseframe = default_frame(source(path))
        desired = Ref(zero(SpatialAcceleration{Float64}, bodyframe, baseframe, frame))
        J = GeometricJacobian(bodyframe, baseframe, frame, zeros(3, nv), zeros(3, nv))
        @assert size(Sangular, 2) == 3
        @assert size(Slinear, 2) == 3
        nrows = size(Sangular, 2) + size(Sangular, 2)
        SJ = zeros(nrows, nv)
        SJ̇v = zeros(nrows)
        new(desired, path, J, angularselection, linearselection, SJ, SJ̇v)
    end
end

function set_desired!(task::SpatialAccelerationTask, desired::SpatialAcceleration)
    @framecheck task.desired[].body desired.body
    @framecheck task.desired[].base desired.base
    @framecheck task.desired[].frame desired.frame
    task.desired[] = desired
    nothing
end

function update!(task::SpatialAccelerationTask, state::MechanismState)
    path = task.path
    desired = task.desired[]
    world_to_desired = inv(transform_to_root(state, desired.frame))
    geometric_jacobian!(task.J, state, path, world_to_desired)
    task.SJ[:] = 0

    A_mul_B!(task.SJ, task.Sangular, angular(J))

    J̇v = transform(state, -bias_acceleration(state, source(path)) + bias_acceleration(state, target(path)), desired.frame)
    @framecheck J.frame J̇v.frame
    @framecheck J.frame desired.frame


    p = desired - J̇v
    copyto!(task.pω, p.angular)
    copyto!(task.pν, p.angular)
    nothing
end

function errorterm(task::SpatialAccelerationTask, v̇)
    angular = LinearTerm(task.J.angular, v̇) - Constant(task.pω.v)
    linear = LinearTerm(task.J.linear, v̇) - Constant(task.pν.v)
    [angular, linear]
end


struct JointAccelerationTask{JT<:JointType{Float64}} <: MotionTask
    joint::Joint{Float64, JT}
    desired::Vector{Float64}

    JointAccelerationTask(joint::Joint{Float64, JT}) where {JT<:JointType{Float64}} = new{JT}(joint, zeros(num_velocities(joint)))
end

function set_desired!(task::JointAccelerationTask, desired)
    set_velocity!(task.desired, task.joint, desired)
    nothing
end

update!(task::JointAccelerationTask, state::MechanismState) = nothing

function errorterm(task::JointAccelerationTask, v̇)
    nv = num_velocities(task.joint)
    LinearTerm(eye(nv), v̇[task.joint]) - Constant(task.desired)
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

function errorterm(task::LinearMomentumRateTask, v̇)
    LinearTerm(task.momentummatrix.linear, v̇) - Constant(task.biased_rate.v)
end


const SparseSymmetric64 = Symmetric{Float64,SparseMatrixCSC{Float64,Int}}

struct Weighted{T<:MotionTask} <: MotionTask
    task::T
    weight::SparseSymmetric64
end
