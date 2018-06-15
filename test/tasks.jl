@testset "SpatialAccelerationTask" begin
    mechanism = RigidBodyDynamics.rand_tree_mechanism(Float64, [Revolute{Float64} for _ = 1 : 10]...)
    state = MechanismState(mechanism)
    rand_configuration!(state)
    nv = num_velocities(mechanism)
    v̇ = [SimpleQP.Variable(i) for i = 1 : nv]
    qpmodel = MockModel()
    for testnum = 1 : 10
        base = rand(bodies(mechanism))
        body = rand(setdiff(bodies(mechanism), [base]))
        p = RigidBodyDynamics.path(mechanism, base, body)
        task = SpatialAccelerationTask(mechanism, p)
        err = MBC.task_error(task, qpmodel, state, v̇)

        bodyframe = default_frame(body)
        baseframe = default_frame(base)
        desired = SpatialAcceleration(bodyframe, baseframe, bodyframe, SVector(1., 2., 3.), SVector(4., 5., 6.))
        MBC.set_desired!(task, desired)

        zero_velocity!(state)
        v̇0 = zeros(length(v̇))
        setdirty!(qpmodel)
        @test map(f -> f(Dict(zip(v̇, v̇0))), err()) == Array(desired)

        MBC.set_desired!(task, zero(desired))
        setdirty!(qpmodel)
        @test map(f -> f(Dict(zip(v̇, v̇0))), err()) == zeros(6)
        rand_velocity!(state)
        biasaccel = transform(state, -RBD.bias_acceleration(state, base) + RBD.bias_acceleration(state, body), bodyframe)
        setdirty!(qpmodel)
        @test map(f -> f(Dict(zip(v̇, v̇0))), err()) == -Array(biasaccel)

        zero_velocity!(state)
        setdirty!(qpmodel)
        v̇rand = rand(nv)
        expected = Array(transform(state, -SpatialAcceleration(geometric_jacobian(state, p), v̇rand), bodyframe))
        @test map(f -> f(Dict(zip(v̇, v̇rand))), err()) ≈ expected atol = 1e-12

        allocs = @allocated err()
        @test allocs == 0
    end
end

@testset "CentroidalMomentumRateTask" begin
    mechanism = RigidBodyDynamics.rand_tree_mechanism(Float64, [Revolute{Float64} for _ = 1 : 10]...)
    state = MechanismState(mechanism)
    rand_configuration!(state)
    nv = num_velocities(mechanism)
    v̇ = [SimpleQP.Variable(i) for i = 1 : nv]

    centroidalframe = CartesianFrame3D("centroidal")
    task = CentroidalMomentumRateTask(mechanism, centroidalframe)
    err = MBC.task_error(task, v̇)

    angular, linear = SVector(1., 2., 3.), SVector(4., 5., 6.)
    desired = Wrench(centroidalframe, angular, linear)
    MBC.set_desired!(task, desired)
    zero_velocity!(state)
    MBC.update!(task, state)
    v̇0 = zeros(length(v̇))
    @test err(Dict(zip(v̇, v̇0))) == -Array(desired)

    MBC.set_desired!(task, zero(desired))
    rand_velocity!(state)
    MBC.update!(task, state)
    world_to_centroidal = Transform3D(root_frame(mechanism), centroidalframe, -center_of_mass(state).v)
    Ȧv = transform(momentum_rate_bias(state), world_to_centroidal)
    @test err(Dict(zip(v̇, v̇0))) == Array(Ȧv)

    zero_velocity!(state)
    MBC.update!(task, state)
    v̇rand = rand(nv)
    @test err(Dict(zip(v̇, v̇rand))) ≈ Array(transform(Wrench(momentum_matrix(state), v̇rand), world_to_centroidal)) atol = 1e-12
end
