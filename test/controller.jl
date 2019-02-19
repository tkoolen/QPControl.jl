@testset "SE3PDController" begin
    val = Valkyrie()
    mechanism = val.mechanism
    state = MechanismState(mechanism)
    rng = MersenneTwister(1)
    base = rand(rng, bodies(mechanism))
    body = rand(rng, bodies(mechanism))
    angular = let angle = π / 2, axis = SVector(1.0, 0.0, 0.0), y0 = one(Quat), yf = Quat(AngleAxis(angle, axis...))
        Interpolated(0.0, 1.0, y0, yf)
    end
    linear = Interpolated(0.0, 1.0, SVector(0.0, 1.0, 2.0), SVector(2.0, 3.0, 4.0))
    frame = default_frame(body)
    gains = SE3PDGains(FramePDGains(frame, PDGains(100.0, 20.0)), FramePDGains(frame, PDGains(1000.0, 200.0)))
    trajectory = SE3Trajectory(default_frame(body), default_frame(base), angular, linear)
    weight = Diagonal(vcat(zeros(3), fill(10.0, 3)))
    controller = SE3PDController(BodyID(base), BodyID(body), trajectory, weight, gains)
    controller(0.5, state)

    # FIXME: proper tests
end

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
    model = Parametron.Model(defaultoptimizer())
    position = Point3D(default_frame(body), 0., 0, 0)
    μ = 1.0
    normal = let state = state, body = body
        Parameter(model) do
            transform(state, FreeVector3D(root_frame(state.mechanism), 0., 0, 1), default_frame(body))
        end
    end
    controller_contact = QPControl.ContactPoint{4}(position, normal, μ, state, model)
    controller_contact.maxnormalforce[] = 1e3
    @test QPControl.isenabled(controller_contact)

    # Constrain the contact force to be non-zero for testing
    @constraint(model, controller_contact.force_local.v == [0.0, 0.0, 1.0])

    for θ in range(-π, stop=π, length=10)
        set_configuration!(state, [θ])
        solve!(model)
        # Sanity check our constraint
        @test value.(Ref(model), controller_contact.force_local.v) ≈ [0.0, 0.0, 1.0]

        # No matter how we rotate the robot, the contact-aligned frame will still
        # be aligned to the contact normal, which is fixed in world frame. So the
        # linear component of the contact wrench in world frame will always be
        # along [0, 0, 1].
        @test value.(Ref(model), linear(controller_contact.wrench_world)) ≈ [0.0, 0.0, 1.0]
    end
end

@testset "fixed base joint space control, constrained = $constrained" for constrained in [true, false]
    Random.seed!(42)
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
    Random.seed!(5354)
    val = Valkyrie()
    mechanism = val.mechanism
    floatingjoint = val.basejoint
    N = 4
    controller = MomentumBasedController{N}(mechanism, defaultoptimizer(), floatingjoint=floatingjoint)
    set_up_valkyrie_contacts!(controller)
    state = MechanismState(mechanism)
    τ = similar(velocity(state))

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
    Random.seed!(533454)
    val = Valkyrie()
    mechanism = val.mechanism
    floatingjoint = val.basejoint
    state = MechanismState(mechanism)
    τ = similar(velocity(state))

    N = 4
    controller = MomentumBasedController{N}(mechanism, defaultoptimizer(), floatingjoint=floatingjoint)

    set_up_valkyrie_contacts!(controller; parametric_contact_surface=true)
    ḣtask = MomentumRateTask(mechanism, centroidal_frame(controller))
    addtask!(controller, ḣtask) #, 1.0)

    for joint in tree_joints(mechanism)
        regularize!(controller, joint, 1e-6)
    end

    for p in range(0., stop=1., length=5)
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
                    normal = Parametron.evalarg(contact.normal)
                    μ = Parametron.evalarg(contact.μ)
                    contact.weight[] = 1e-6
                    contact.maxnormalforce[] = MAX_NORMAL_FORCE_FIXME
                    fnormal = 50. * rand()
                    μreduced = sqrt(2) / 2 * μ # due to polyhedral inner approximation; assumes 4 basis vectors or more
                    ftangential = μreduced * fnormal * rand() * cross(normal, FreeVector3D(normal.frame, normalize(randn(SVector{3}))))
                    f = fnormal * normal + ftangential
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
        @test isapprox(ḣdes, ḣ; atol = 1e-3)

        allocs = @allocated controller(τ, 0., state)
        @test allocs == 0
    end
end

@enum TaskMode begin
    CONSTRAINT
    OBJECTIVE_TERM_SCALAR_WEIGHT
    OBJECTIVE_TERM_SCALAR_WEIGHT_PARAMETER
    OBJECTIVE_TERM_MATRIX_WEIGHT
    OBJECTIVE_TERM_MATRIX_WEIGHT_PARAMETER
end

@testset "spatial acceleration, mode = $mode" for mode in instances(TaskMode)
    Random.seed!(533)
    val = Valkyrie()
    mechanism = val.mechanism
    floatingjoint = val.basejoint
    state = MechanismState(mechanism)
    τ = similar(velocity(state))

    rand!(state)
    N = 4
    controller = MomentumBasedController{N}(mechanism, defaultoptimizer(), floatingjoint=floatingjoint)
    body = val.feet[left]
    base = val.palms[right]
    frame = default_frame(base)
    task = SpatialAccelerationTask(mechanism, path(mechanism, base, body), frame=frame)

    result = DynamicsResult(mechanism)
    accels = result.accelerations
    if mode == CONSTRAINT
        addtask!(controller, task)
        regularize!.(Ref(controller), tree_joints(mechanism), 1.0)
    elseif mode == OBJECTIVE_TERM_SCALAR_WEIGHT
        addtask!(controller, task, 1.0)
    elseif mode == OBJECTIVE_TERM_SCALAR_WEIGHT_PARAMETER
        weight = Parameter(controller.qpmodel, val=1.0)
        addtask!(controller, task, weight)
    elseif mode == OBJECTIVE_TERM_MATRIX_WEIGHT
        addtask!(controller, task, Matrix(I, 6, 6))
    elseif mode == OBJECTIVE_TERM_MATRIX_WEIGHT_PARAMETER
        weight = Parameter(controller.qpmodel, val=Matrix(I, 6, 6))
        addtask!(controller, task, weight)
    else
        error()
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
