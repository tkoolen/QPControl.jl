@testset "parameterized contacts" begin
    # Construct a mechanism consisting of a single body which can
    # rotate about its origin
    world = RigidBody{Float64}("world")
    mechanism = Mechanism(world, gravity=SVector(0, 0, -9.81))
    frame = CartesianFrame3D("body")
    inertia = SpatialInertia(frame, SDiagonal(1., 1, 1), SVector(0., 0, 0), 10.0)
    body = RigidBody(inertia)
    joint = Joint("rx", Revolute(SVector(1., 0, 0)))
    attach!(mechanism, world, body, joint)

    contactmodel = SoftContactModel(hunt_crossley_hertz(k = 500e3), ViscoelasticCoulombModel(0.8, 20e3, 100.))
    add_contact_point!(body, Contact.ContactPoint(Point3D(default_frame(body), 0., 0, 0), contactmodel))

    # Make the contact normal a parameter, so that it updates its representation
    # in body frame to match a fixed orientation in world frame
    state = MechanismState(mechanism)
    model = SimpleQP.Model(defaultoptimizer())
    position = Point3D(default_frame(body), 0., 0, 0)
    μ = 1.0
    normal = let state = state, body = body
        Parameter(model) do
            transform(state, FreeVector3D(root_frame(state.mechanism), 0., 0, 1), default_frame(body))
        end
    end
    controller_contact = QPControl.ContactPoint{4}(position, normal, μ, state, model)
    controller_contact.maxnormalforce[] = 1e3

    # Constrain the contact force to be non-zero for testing
    @constraint(model, controller_contact.force_local.v == [0.0, 0.0, 1.0])

    @objective(model, Minimize, 0 * dot(controller_contact.force_local.v, controller_contact.force_local.v)) # Need to have an objective to avoid throwing "unsupported objective" errors

    solve!(model)

    # Since the body is currently aligned to the world, the local force and linear world wrench should match:
    @test value.(model, linear(controller_contact.wrench_world)) ≈ [0.0, 0.0, 1.0]
    @test value.(model, controller_contact.force_local.v) ≈ [0.0, 0.0, 1.0]

    # Now we rotate the robot 180 degrees, so the z axis of the body is opposite the z axis of the world
    set_configuration!(state, [π])
    # Re-solve the model, which should update the parameters
    solve!(model)
    # Now that the body z is opposite world z, the linear part of wrench_world should
    # point in the opposite direction from the force_local
    @test value.(model, linear(controller_contact.wrench_world)) ≈ [0.0, 0.0, -1.0]
    @test value.(model, controller_contact.force_local.v) ≈ [0.0, 0.0, 1.0]
end

@testset "fixed base joint space control, constrained = $constrained" for constrained in [true, false]
    srand(42)
    mechanism = rand_tree_mechanism(Float64, Prismatic{Float64}, Revolute{Float64}, Revolute{Float64})
    N = 4
    controller = MomentumBasedController{N}(mechanism, defaultoptimizer())
    tasks = Dict{Joint{Float64}, JointAccelerationTask}()
    for joint in tree_joints(mechanism)
        task = JointAccelerationTask(joint)
        tasks[joint] = task
        if constrained
            addtask!(controller, task)
        else
            weight = 1.0
            addtask!(controller, task, weight)
        end
        setdesired!(task, rand(num_velocities(joint)))
    end
    state = MechanismState(mechanism)
    rand!(state)
    τ = similar(velocity(state))
    controller(τ, 0.0, state)
    allocs = @allocated controller(τ, 0.0, state)
    @test allocs == 0

    result = DynamicsResult(mechanism)
    dynamics!(result, state, τ)
    for joint in tree_joints(mechanism)
        @test result.v̇[joint] ≈ tasks[joint].desired atol = 1e-10
    end
end

"""
Automatically load contact points from each body in the mechanism and add them
to the controller. If `parametric_contact_surface=true`, then the `normal` and
`μ` associated with each contact will be Parameters set to random, but fixed,
values.
"""
function set_up_valkyrie_contacts!(controller::MomentumBasedController; parametric_contact_surface=false)
    valmechanism = controller.state.mechanism
    for body in bodies(valmechanism)
        for point in RBD.contact_points(body)
            position = RBD.Contact.location(point)
            if parametric_contact_surface
                normal = let frame = position.frame, direction = normalize(randn(SVector{3}))
                    Parameter(controller.qpmodel) do
                        FreeVector3D(frame, direction)
                    end
                end
                μ = let μ = rand()
                    Parameter(() -> μ, controller.qpmodel)
                end
            else
                normal = FreeVector3D(position.frame, 0.0, 0.0, 1.0)
                μ = point.model.friction.μ
            end
            addcontact!(controller, body, position, normal, μ)
        end
    end
end

@testset "zero velocity free fall" begin
    val = Valkyrie()
    mechanism = val.mechanism
    floatingjoint = val.basejoint
    N = 4
    controller = MomentumBasedController{N}(mechanism, defaultoptimizer())
    set_up_valkyrie_contacts!(controller)
    state = MechanismState(mechanism)
    τ = similar(velocity(state))

    srand(5354)
    zero!(state)
    rand_configuration!(state)
    for joint in tree_joints(mechanism)
        joint == floatingjoint && continue
        regularize!(controller, joint, 1.0)
    end
    controller(τ, 0., state)
    result = DynamicsResult(mechanism)
    accels = result.accelerations
    RBD.spatial_accelerations!(accels, state, controller.result.v̇)
    for joint in tree_joints(mechanism)
        if joint == floatingjoint
            baseaccel = relative_acceleration(accels, val.pelvis, root_body(mechanism))
            baseaccel = transform(state, baseaccel, frame_after(floatingjoint))
            angularaccel = FreeVector3D(baseaccel.frame, baseaccel.angular)
            @test isapprox(angularaccel, FreeVector3D(frame_after(floatingjoint), zeros(SVector{3})), atol = 1e-4)
            linearaccel = FreeVector3D(baseaccel.frame, baseaccel.linear)
            linearaccel = transform_to_root(state, linearaccel.frame) * linearaccel
            @test isapprox(linearaccel, mechanism.gravitational_acceleration; atol = 1e-4)
        else
            v̇joint = controller.result.v̇[velocity_range(state, joint)]
            @test isapprox(v̇joint, zeros(num_velocities(joint)); atol = 1e-4)
        end
    end
    allocs = @allocated controller(τ, 0., state)
    @test allocs == 0
end

const MAX_NORMAL_FORCE_FIXME = 1e9

@testset "achievable momentum rate" begin
    val = Valkyrie()
    mechanism = val.mechanism
    floatingjoint = val.basejoint
    state = MechanismState(mechanism)
    τ = similar(velocity(state))

    N = 4
    controller = MomentumBasedController{N}(mechanism, GurobiOptimizer(OutputFlag=0))

    set_up_valkyrie_contacts!(controller; parametric_contact_surface=false)
    ḣtask = MomentumRateTask(mechanism, centroidal_frame(controller))
    addtask!(controller, ḣtask, 1000.0)

    for joint in tree_joints(mechanism)
        regularize!(controller, joint, 1e-6)
    end

    srand(1)
    for p in linspace(0., 1., 5)
        rand!(state)
        com = center_of_mass(state)
        centroidal_to_world = Transform3D(centroidal_frame(controller), com.frame, com.v)
        world_to_centroidal = inv(centroidal_to_world)

        # set random active contacts and random achievable wrench
        fg = world_to_centroidal * (mass(mechanism) * mechanism.gravitational_acceleration)
        ḣdes = Wrench(zero(fg), fg)
        for body in keys(controller.contacts)
            for contact in controller.contacts[body]
                active = rand() < p
                if active
                    normal = contact.normal
                    μ = contact.μ
                    contact.weight[] = 1e-6
                    contact.maxnormalforce[] = MAX_NORMAL_FORCE_FIXME
                    fnormal = 50. * rand()
                    μreduced = sqrt(2) / 2 * μ # due to polyhedral inner approximation; assumes 4 basis vectors or more
                    ftangential = μreduced * fnormal * rand() * cross(normal, FreeVector3D(normal.frame, normalize(randn(SVector{3}))))
                    f = fnormal * normal + ftangential
                    @show normal μ f
                    @assert isapprox(ftangential ⋅ normal, 0., atol = 1e-12)
                    @assert norm(f - (normal ⋅ f) * normal) ≤ μreduced * (normal ⋅ f)
                    wrench = Wrench(contact.position, f)
                    ḣdes += transform(wrench, world_to_centroidal * transform_to_root(state, wrench.frame))
                else
                    disable!(contact)
                end
            end
        end
        setdesired!(ḣtask, ḣdes)

        controller(τ, 0., state)

        # Ensure that desired momentum rate is achieved.
        ḣ = Wrench(momentum_matrix(state), controller.result.v̇) + momentum_rate_bias(state)
        ḣ = transform(ḣ, world_to_centroidal)
        @show ḣ ḣdes
        @test isapprox(ḣdes, ḣ; atol = 1e-3)

        allocs = @allocated controller(τ, 0., state)
        # @test allocs == 0
    end
end

@testset "spatial acceleration, constrained = $constrained" for constrained in [true, false]
    val = Valkyrie()
    mechanism = val.mechanism
    floatingjoint = val.basejoint
    state = MechanismState(mechanism)
    τ = similar(velocity(state))

    srand(533)
    rand!(state)
    N = 4
    controller = MomentumBasedController{N}(mechanism, defaultoptimizer())
    body = val.feet[left]
    base = val.palms[right]
    frame = default_frame(base)
    task = SpatialAccelerationTask(mechanism, path(mechanism, base, body), frame=frame)

    result = DynamicsResult(mechanism)
    accels = result.accelerations
    if constrained
        addtask!(controller, task)
        regularize!.(controller, tree_joints(mechanism), 1.0)
    else
        addtask!(controller, task, 1.0)
    end
    desiredaccel = rand(SpatialAcceleration{Float64}, default_frame(body), default_frame(base), frame)
    setdesired!(task, desiredaccel)
    controller(τ, 0., state)
    v̇ = controller.result.v̇
    spatial_accelerations!(accels, state, v̇)
    accel = relative_acceleration(accels, body, base)
    accel = transform(state, accel, frame)
    @test isapprox(accel, desiredaccel, atol = 1e-8)

    allocs = @allocated controller(τ, 0., state)
    @test allocs == 0
end
