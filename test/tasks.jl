@testset "JointAccelerationTask" begin
    mechanism = RBD.rand_tree_mechanism(Float64, [Revolute{Float64} for _ = 1 : 10]...)
    state = MechanismState(mechanism)
    rand!(state)
    nv = num_velocities(mechanism)
    v̇ = [SimpleQP.Variable(i) for i = 1 : nv]
    v̇vals = rand!(similar(velocity(state)))
    qpmodel = MockModel()
    tasks = Dict{Joint{Float64, Revolute{Float64}}, JointAccelerationTask{Revolute{Float64}}}()
    for joint in tree_joints(mechanism)
        task = JointAccelerationTask(joint)
        tasks[joint] = task
        setdesired!(task, rand())
    end
    for (joint, task) in tasks
        err = MBC.task_error(task, qpmodel, state, v̇)
        @test map(f -> f(Dict(zip(v̇, v̇vals))), err()) == v̇vals[joint] .- task.desired
        allocs = @allocated err()
        @test allocs == 0
    end
end

@testset "SpatialAccelerationTask" begin
    mechanism = RBD.rand_tree_mechanism(Float64, [Revolute{Float64} for _ = 1 : 10]...)
    state = MechanismState(mechanism)
    rand_configuration!(state)
    nv = num_velocities(mechanism)
    v̇ = [SimpleQP.Variable(i) for i = 1 : nv]
    qpmodel = MockModel()
    for testnum = 1 : 10
        base = rand(bodies(mechanism))
        body = rand(setdiff(bodies(mechanism), [base]))
        p = RBD.path(mechanism, base, body)
        task = SpatialAccelerationTask(mechanism, p)
        err = MBC.task_error(task, qpmodel, state, v̇)

        bodyframe = default_frame(body)
        baseframe = default_frame(base)
        desired = SpatialAcceleration(bodyframe, baseframe, bodyframe, SVector(1., 2., 3.), SVector(4., 5., 6.))
        MBC.setdesired!(task, desired)

        zero_velocity!(state)
        v̇0 = zeros(length(v̇))
        setdirty!(qpmodel)
        @test map(f -> f(Dict(zip(v̇, v̇0))), err()) == -Array(desired)

        MBC.setdesired!(task, zero(desired))
        setdirty!(qpmodel)
        @test map(f -> f(Dict(zip(v̇, v̇0))), err()) == zeros(6)
        rand_velocity!(state)
        biasaccel = transform(state, -RBD.bias_acceleration(state, base) + RBD.bias_acceleration(state, body), bodyframe)
        setdirty!(qpmodel)
        @test map(f -> f(Dict(zip(v̇, v̇0))), err()) == Array(biasaccel)

        zero_velocity!(state)
        setdirty!(qpmodel)
        v̇rand = rand(nv)
        expected = Array(transform(state, SpatialAcceleration(geometric_jacobian(state, p), v̇rand), bodyframe))
        @test map(f -> f(Dict(zip(v̇, v̇rand))), err()) ≈ expected atol = 1e-12

        allocs = @allocated err()
        @test allocs == 0
    end
end

@testset "MomentumRateTask" begin
    mechanism = RBD.rand_tree_mechanism(Float64, [Revolute{Float64} for _ = 1 : 10]...)
    state = MechanismState(mechanism)
    rand_configuration!(state)
    nv = num_velocities(mechanism)
    v̇ = [SimpleQP.Variable(i) for i = 1 : nv]
    qpmodel = MockModel()

    centroidalframe = CartesianFrame3D("centroidal")
    task = MomentumRateTask(mechanism, centroidalframe)
    err = MBC.task_error(task, qpmodel, state, v̇)

    angular, linear = SVector(1., 2., 3.), SVector(4., 5., 6.)
    desired = Wrench(centroidalframe, angular, linear)
    MBC.setdesired!(task, desired)
    zero_velocity!(state)
    setdirty!(qpmodel)
    v̇0 = zeros(length(v̇))
    @test map(f -> f(Dict(zip(v̇, v̇0))), err()) == -Array(desired)

    MBC.setdesired!(task, zero(desired))
    rand_velocity!(state)
    setdirty!(qpmodel)
    world_to_centroidal = Transform3D(root_frame(mechanism), centroidalframe, -center_of_mass(state).v)
    Ȧv = transform(momentum_rate_bias(state), world_to_centroidal)
    @test map(f -> f(Dict(zip(v̇, v̇0))), err()) == Array(Ȧv)

    zero_velocity!(state)
    setdirty!(qpmodel)
    v̇rand = rand(nv)
    expected = transform(Wrench(momentum_matrix(state), v̇rand), world_to_centroidal)
    @test map(f -> f(Dict(zip(v̇, v̇rand))), err()) ≈ Array(expected) atol = 1e-12

    allocs = @allocated err()
    @test allocs == 0
end
