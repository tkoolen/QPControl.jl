@testset "JointAccelerationTask" begin
    Random.seed!(1)
    mechanism = RBD.rand_tree_mechanism(Float64, [Revolute{Float64} for _ = 1 : 10]...)
    state = MechanismState(mechanism)
    rand!(state)
    nv = num_velocities(mechanism)
    v̇ = [Parametron.Variable(i) for i = 1 : nv]
    v̇vals = rand!(similar(velocity(state)))
    qpmodel = mock_model()
    tasks = Dict{Joint{Float64, Revolute{Float64}}, JointAccelerationTask{Revolute{Float64}}}()
    for joint in tree_joints(mechanism)
        task = JointAccelerationTask(joint)
        @test QPC.dimension(task) == num_velocities(joint)
        tasks[joint] = task
        setdesired!(task, rand())
    end
    for (joint, task) in tasks
        err = QPC.task_error(task, qpmodel, state, v̇)
        @test map(f -> f(Dict(zip(v̇, v̇vals))), err()) == v̇vals[joint] .- task.desired
        allocs = @allocated err()
        @test allocs == 0
    end
end

@testset "AngularAccelerationTask" begin
    Random.seed!(2)
    mechanism = RBD.rand_tree_mechanism(Float64, [Revolute{Float64} for _ = 1 : 10]...)
    state = MechanismState(mechanism)
    rand_configuration!(state)
    nv = num_velocities(mechanism)
    v̇ = [Parametron.Variable(i) for i = 1 : nv]
    qpmodel = mock_model()
    for testnum = 1 : 10
        base = rand(bodies(mechanism))
        body = rand(setdiff(bodies(mechanism), [base]))
        p = RBD.path(mechanism, base, body)
        task = AngularAccelerationTask(mechanism, p)
        @test QPC.dimension(task) == 3
        err = QPC.task_error(task, qpmodel, state, v̇)

        bodyframe = default_frame(body)
        baseframe = default_frame(base)
        desired = FreeVector3D(bodyframe, SVector(1., 2, 3))
        QPC.setdesired!(task, desired)

        zero_velocity!(state)
        v̇0 = zeros(length(v̇))
        setdirty!(qpmodel)
        @test map(f -> f(Dict(zip(v̇, v̇0))), err()) == -desired.v

        QPC.setdesired!(task, zero(desired))
        setdirty!(qpmodel)
        @test map(f -> f(Dict(zip(v̇, v̇0))), err()) == zeros(3)
        rand_velocity!(state)
        biasaccel = transform(state, -RBD.bias_acceleration(state, base) + RBD.bias_acceleration(state, body), bodyframe)
        setdirty!(qpmodel)
        @test map(f -> f(Dict(zip(v̇, v̇0))), err()) == angular(biasaccel)

        zero_velocity!(state)
        setdirty!(qpmodel)
        v̇rand = rand(nv)
        expected = angular(transform(state, SpatialAcceleration(geometric_jacobian(state, p), v̇rand), bodyframe))
        @test map(f -> f(Dict(zip(v̇, v̇rand))), err()) ≈ expected atol = 1e-12

        allocs = @allocated err()
        @test allocs == 0
    end
end

@testset "PointAccelerationTask" begin
    Random.seed!(3)
    mechanism = RBD.rand_tree_mechanism(Float64, [Revolute{Float64} for _ = 1 : 10]...)
    state = MechanismState(mechanism)
    rand_configuration!(state)
    nv = num_velocities(mechanism)
    v̇ = [Parametron.Variable(i) for i = 1 : nv]
    qpmodel = mock_model()
    for testnum = 1 : 10
        base = rand(bodies(mechanism))
        body = rand(setdiff(bodies(mechanism), [base]))
        path_to_body = RBD.path(mechanism, base, body)
        point = Point3D(default_frame(body), randn(SVector{3, Float64}))
        task = PointAccelerationTask(mechanism, path_to_body, point)
        @test QPC.dimension(task) == 3
        err = QPC.task_error(task, qpmodel, state, v̇)

        bodyframe = default_frame(body)
        baseframe = default_frame(base)
        desired = FreeVector3D(baseframe, SVector(1., 2, 3))
        QPC.setdesired!(task, desired)

        zero_velocity!(state)
        v̇0 = zeros(length(v̇))
        setdirty!(qpmodel)
        @test map(f -> f(Dict(zip(v̇, v̇0))), err()) == -desired.v

        QPC.setdesired!(task, zero(desired))
        setdirty!(qpmodel)
        @test map(f -> f(Dict(zip(v̇, v̇0))), err()) == zeros(3)

        rand_velocity!(state)
        T = transform(state, relative_twist(state, body, base), baseframe)
        ω = angular(T)
        ṗ = point_velocity(T, transform(state, point, baseframe))
        bias = transform(state,
            -RBD.bias_acceleration(state, base) + RBD.bias_acceleration(state, body),
            baseframe)
        expected = ω × ṗ.v + (angular(bias) × transform(state, point, baseframe).v + linear(bias))
        setdirty!(qpmodel)
        @test map(f -> f(Dict(zip(v̇, v̇0))), err()) == expected

        zero_velocity!(state)
        setdirty!(qpmodel)
        v̇rand = rand(nv)
        J_point = point_jacobian(state, path_to_body, transform(state, point, baseframe))
        expected = point_velocity(J_point, v̇rand).v
        setdirty!(qpmodel)
        @test map(f -> f(Dict(zip(v̇, v̇rand))), err()) ≈ expected atol=1e-12

        allocs = @allocated err()
        @test allocs == 0
    end
end


@testset "LinearAccelerationTask" begin
    Random.seed!(4)
    mechanism = RBD.rand_tree_mechanism(Float64, [Revolute{Float64} for _ = 1 : 10]...)
    state = MechanismState(mechanism)
    rand_configuration!(state)
    nv = num_velocities(mechanism)
    v̇ = [Parametron.Variable(i) for i = 1 : nv]
    qpmodel = mock_model()
    for testnum = 1 : 10
        base = rand(bodies(mechanism))
        body = rand(setdiff(bodies(mechanism), [base]))
        p = RBD.path(mechanism, base, body)
        task = LinearAccelerationTask(mechanism, p)
        @test QPC.dimension(task) == 3
        err = QPC.task_error(task, qpmodel, state, v̇)

        bodyframe = default_frame(body)
        baseframe = default_frame(base)
        desired = FreeVector3D(bodyframe, SVector(1., 2, 3))
        QPC.setdesired!(task, desired)

        zero_velocity!(state)
        v̇0 = zeros(length(v̇))
        setdirty!(qpmodel)
        @test map(f -> f(Dict(zip(v̇, v̇0))), err()) == -desired.v

        QPC.setdesired!(task, zero(desired))
        setdirty!(qpmodel)
        @test map(f -> f(Dict(zip(v̇, v̇0))), err()) == zeros(3)
        rand_velocity!(state)
        biasaccel = transform(state, -RBD.bias_acceleration(state, base) + RBD.bias_acceleration(state, body), bodyframe)
        setdirty!(qpmodel)
        @test map(f -> f(Dict(zip(v̇, v̇0))), err()) == linear(biasaccel)

        zero_velocity!(state)
        setdirty!(qpmodel)
        v̇rand = rand(nv)
        expected = linear(transform(state, SpatialAcceleration(geometric_jacobian(state, p), v̇rand), bodyframe))
        @test map(f -> f(Dict(zip(v̇, v̇rand))), err()) ≈ expected atol = 1e-12

        allocs = @allocated err()
        @test allocs == 0
    end
end

@testset "SpatialAccelerationTask" begin
    Random.seed!(5)
    mechanism = RBD.rand_tree_mechanism(Float64, [Revolute{Float64} for _ = 1 : 10]...)
    state = MechanismState(mechanism)
    rand_configuration!(state)
    nv = num_velocities(mechanism)
    v̇ = [Parametron.Variable(i) for i = 1 : nv]
    qpmodel = mock_model()
    for testnum = 1 : 10
        base = rand(bodies(mechanism))
        body = rand(setdiff(bodies(mechanism), [base]))
        p = RBD.path(mechanism, base, body)
        task = SpatialAccelerationTask(mechanism, p)
        @test QPC.dimension(task) == 6
        err = QPC.task_error(task, qpmodel, state, v̇)

        bodyframe = default_frame(body)
        baseframe = default_frame(base)
        desired = SpatialAcceleration(bodyframe, baseframe, bodyframe, SVector(1., 2., 3.), SVector(4., 5., 6.))
        QPC.setdesired!(task, desired)

        zero_velocity!(state)
        v̇0 = zeros(length(v̇))
        setdirty!(qpmodel)
        @test map(f -> f(Dict(zip(v̇, v̇0))), err()) == -SVector(desired)

        QPC.setdesired!(task, zero(desired))
        setdirty!(qpmodel)
        @test map(f -> f(Dict(zip(v̇, v̇0))), err()) == zeros(6)
        rand_velocity!(state)
        biasaccel = transform(state, -RBD.bias_acceleration(state, base) + RBD.bias_acceleration(state, body), bodyframe)
        setdirty!(qpmodel)
        @test map(f -> f(Dict(zip(v̇, v̇0))), err()) == SVector(biasaccel)

        zero_velocity!(state)
        setdirty!(qpmodel)
        v̇rand = rand(nv)
        expected = SVector(transform(state, SpatialAcceleration(geometric_jacobian(state, p), v̇rand), bodyframe))
        @test map(f -> f(Dict(zip(v̇, v̇rand))), err()) ≈ expected atol = 1e-12

        allocs = @allocated err()
        @test allocs == 0
    end
end

@testset "MomentumRateTask" begin
    Random.seed!(6)
    mechanism = RBD.rand_tree_mechanism(Float64, [Revolute{Float64} for _ = 1 : 10]...)
    state = MechanismState(mechanism)
    rand_configuration!(state)
    nv = num_velocities(mechanism)
    v̇ = [Parametron.Variable(i) for i = 1 : nv]
    qpmodel = mock_model()

    centroidalframe = CartesianFrame3D("centroidal")
    task = MomentumRateTask(mechanism, centroidalframe)
    @test QPC.dimension(task) == 6
    err = QPC.task_error(task, qpmodel, state, v̇)

    angular, linear = SVector(1., 2., 3.), SVector(4., 5., 6.)
    desired = Wrench(centroidalframe, angular, linear)
    QPC.setdesired!(task, desired)
    zero_velocity!(state)
    setdirty!(qpmodel)
    v̇0 = zeros(length(v̇))
    @test map(f -> f(Dict(zip(v̇, v̇0))), err()) == -SVector(desired)

    QPC.setdesired!(task, zero(desired))
    rand_velocity!(state)
    setdirty!(qpmodel)
    world_to_centroidal = Transform3D(root_frame(mechanism), centroidalframe, -center_of_mass(state).v)
    Ȧv = transform(momentum_rate_bias(state), world_to_centroidal)
    @test map(f -> f(Dict(zip(v̇, v̇0))), err()) == SVector(Ȧv)

    zero_velocity!(state)
    setdirty!(qpmodel)
    v̇rand = rand(nv)
    expected = transform(Wrench(momentum_matrix(state), v̇rand), world_to_centroidal)
    @test map(f -> f(Dict(zip(v̇, v̇rand))), err()) ≈ SVector(expected) atol = 1e-12

    allocs = @allocated err()
    @test allocs == 0
end
