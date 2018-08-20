using Gurobi

function add_base_joint!(mechanism, label, direction)
    frame = CartesianFrame3D("base_$(label)_dummy")
    inertia = SpatialInertia(frame, SDiagonal(0., 0, 0), SVector(0., 0, 0), 0.0)
    body = RigidBody(inertia)
    joint = Joint("base_$(label)", Prismatic(direction))
    effort_bounds(joint) .= RigidBodyDynamics.Bounds(0, 0)
    parent = last(bodies(mechanism))
    attach!(mechanism, parent, body, joint)
end

function three_dof_point_mass()
    world = RigidBody{Float64}("world")
    mechanism = Mechanism(world, gravity=SVector(0, 0, -9.81))
    add_base_joint!(mechanism, "x", SVector(1., 0, 0))
    add_base_joint!(mechanism, "y", SVector(0., 1, 0))
    add_base_joint!(mechanism, "z", SVector(0., 0, 1))

    frame = CartesianFrame3D("core")
    inertia = SpatialInertia(frame, SDiagonal(0.01, 0.01, 0.01), SVector(0., 0, 0), 1.0)
    body = RigidBody(inertia)
    joint = Joint("core_to_base", Fixed{Float64}())
    parent = last(bodies(mechanism))
    attach!(mechanism, parent, body, joint)
    mechanism
end

function add_floor!(mechanism, orientation, μ, contacting_body)
    world = root_body(mechanism)
    floor = HalfSpace3D(Point3D(default_frame(world), 0., 0, 0), FreeVector3D(default_frame(world), orientation * SVector(0., 0, 1)))
    add_environment_primitive!(mechanism, floor)

    contactmodel = SoftContactModel(hunt_crossley_hertz(k = 500e3), ViscoelasticCoulombModel(μ, 20e3, 100.))
    add_contact_point!(contacting_body, Contact.ContactPoint(Point3D(default_frame(contacting_body), 0., 0, 0), contactmodel))
    floor
end

function mpc_lcp_feasibility_problem(mechanism, optimizer, floor, Δt)
    mpc = QPControl.MPCController{QPControl.ContactPoint{8}}(mechanism, optimizer)
    stage = QPControl.addstage!(mpc, 0.01)
    for body in bodies(mechanism)
        for point in RigidBodyDynamics.contact_points(body)
            position = location(point)
            μ = point.model.friction.μ
            contact = addcontact!(mpc, stage, position, floor, μ, QPControl.LCPContact())
            contact.maxnormalforce[] = 1e6
        end
    end
    mpc
end

function simulate_lcp!(state::MechanismState, mpc::MPCController, t_final=1.0)
    # Simulate the system by repeatedly running the MPC controller
    # and then updating the state based on the expected state
    # at the end of the first MPC stage.
    ts = Float64[]
    qs = Vector{Float64}[]
    vs = Vector{Float64}[]
    τ = similar(velocity(state))
    for t in 0 : (mpc.stages[1].Δt) : t_final
        push!(ts, t)
        push!(qs, configuration(state))
        push!(vs, velocity(state))
        mpc(τ, t, state)
        # Assume the dynamics behave exactly as modeled in the optimization
        set_configuration!(state, Parametron.value.(mpc.qpmodel, mpc.stages[1].q))
        set_velocity!(state, Parametron.value.(mpc.qpmodel, mpc.stages[1].v))
    end
    ts, qs, vs
end

@testset "LCP: sticking and sliding" begin
    # This test is meant to verify the LCP contact physics by
    # simulating a point with friction on a variety of sloped planes.
    # We will test both that for a sufficiently shallow plane, the point
    # does not slide, and that for a steeper plane it does slide.

    optimizer = GurobiOptimizer(OutputFlag=0, TimeLimit=1)

    @testset "μ = $μ" for μ in [0.25, 0.5, 0.75]
        @testset "stick = $stick" for stick in [false, true]
            @testset "axis = $axis" for axis in [RotZ(ϕ) * SVector(1., 0, 0) for ϕ in Compat.range(0, stop=2π-π/4, step=π/4)]
                θ = atan(μ)
                if stick
                    θ -= 0.05
                else
                    θ += 0.05
                end
                orientation = AngleAxis(θ, axis...)
                Δt = 0.01

                # Create a mechanism with a 3-DoF, unactuated, floating base.
                mechanism = three_dof_point_mass()
                floor = add_floor!(mechanism, orientation, μ, findbody(mechanism, "core"))

                # Create an MPC "controller", which we will use as a simple
                # one-step LCP simulator
                mpc = mpc_lcp_feasibility_problem(mechanism, optimizer, floor, Δt)

                state = MechanismState(mechanism)
                # Start off the ground a bit, to be sure we can dissipate
                # velocity correctly
                set_configuration!(state, findjoint(mechanism, "base_z"), 0.1)

                ts, qs, vs = simulate_lcp!(state, mpc)

                # The point should be on the floor
                @test floor.outward_normal.v' * qs[end] ≈ floor.outward_normal.v' * floor.point.v atol=1e-4

                if stick
                    # The point should have dropped straight down and stopped when it
                    # hit the floor, without sliding
                    @test qs[end] ≈ [0, 0, 0] atol=1e-3
                    @test vs[end] ≈ [0, 0, 0] atol=1e-3
                else
                    # The point should be sliding downhill
                    expected_direction = orientation * RotZ(-π/2) * axis
                    @test dot(normalize(vs[end]), expected_direction) ≈ 1

                    # The point should not have moved perpendicular to the slope
                    @test dot(qs[end], axis) ≈ 0 atol=1e-2
                end
            end
        end
    end
end


"""
A two-link mechanism consisting of a body and a foot, both moving along the
world Z axis.
"""
function hopper()
    world = RigidBody{Float64}("world")
    mechanism = Mechanism(world; gravity=SVector(0, 0, -9.81))

    frame = CartesianFrame3D("core")
    inertia = SpatialInertia(frame, SDiagonal(1., 1, 1), SVector(0., 0, 0), 10.)
    core = RigidBody(inertia)
    joint = Joint("floating_base", Prismatic(SVector(0., 0, 1)))
    attach!(mechanism, world, core, joint)
    position_bounds(joint) .= RigidBodyDynamics.Bounds(0, 2)
    velocity_bounds(joint) .= RigidBodyDynamics.Bounds(-100, 100)
    effort_bounds(joint) .= RigidBodyDynamics.Bounds(0, 0)

    frame = CartesianFrame3D("foot")
    inertia = SpatialInertia(frame, SDiagonal(0.05, 0.05, 0.05), SVector(0., 0, 0), 1.0)
    foot = RigidBody(inertia)
    joint = Joint("foot_extension", Prismatic(SVector(0., 0, 1)))
    attach!(mechanism, core, foot, joint)
    position_bounds(joint) .= RigidBodyDynamics.Bounds(-1., -0.5)
    velocity_bounds(joint) .= RigidBodyDynamics.Bounds(-100, 100)
    effort_bounds(joint) .= RigidBodyDynamics.Bounds(0, 0)
    mechanism
end

@testset "LCP with joint limits" begin
    optimizer = GurobiOptimizer(OutputFlag=0, TimeLimit=1)
    mechanism = hopper()
    floor = add_floor!(mechanism, RotX(0.0), 1.0, findbody(mechanism, "foot"))
    # Create an MPC "controller", which we will use as a simple
    # one-step LCP simulator
    Δt = 0.01
    mpc = mpc_lcp_feasibility_problem(mechanism, optimizer, floor, Δt)

    state = MechanismState(mechanism)

    # We're going to start the hopper slightly off the ground, with both
    # bodies moving upward. The expected behavior is:
    # 1. Both bodies move up together, with v1 decelerating at 9.81m/s2 and
    #    v2 constant at 0
    # 2. The upper link hits its position limit at q1 = 2.0 and immediately
    #    stops.
    # 3. The lower link continues at the same velocity in world frame, while
    #    the upper link starts to fall
    # 4. The joint from upper to lower link hits its position limit at q2 = -0.5,
    #    causing v2 to drop to 0 immediately. v1 increases slightly due to the impact
    # 5. Both bodies fall together, with q2 fixed at -0.5 and v2 fixed at 0
    # 6. The foot hits the ground when q1 = 0.5 and q2 = -0.5, at which point
    #    v1 also drops to 0 immediately
    # 7. Both bodies remain at rest at q1 = 0.5, q2 = -0.5

    set_configuration!(state, [1.5, -1.0])
    set_velocity!(state, [5.0, 0.0])

    ts, qs, vs = simulate_lcp!(state, mpc)
    @test qs[1] ≈ [1.5, -1.0]
    @test vs[1] ≈ [5.0, 0.0]
    @test ts[1] ≈ 0

    # 1. Both bodies move up together, with v1 decelerating at 9.81m/s2 and
    #    v2 constant at 0
    #
    # q1(t) = 1.5 + 5.0 t - 1/2 * 9.81 * t^2 = 2.0
    # So the first collision occurs at t = 0.11239194149021249
    for i in 1:12
        t = Δt * (i - 1)
        @test ts[i] ≈ t
        @test vs[i] ≈ [5.0 - 9.81 * t, 0.0]
        @test qs[i] ≈ [1.5 + 5.0 * t - 1/2 * 9.81 * t^2, -1.0] rtol=1e-2
    end

    # 2. The upper link hits its position limit at q1 = 2.0 and immediately
    #    stops.
    @test qs[13][1] ≈ 2.0
    @test vs[14][1] ≈ 0 atol=1e-4

    # 3. The lower link continues at the same velocity in world frame, while
    #    the upper link starts to fall
    for i in 15:25
        @test vs[i] ≈ [-9.81 * Δt * (i - 14), vs[14][2]]
    end

    # 4. The joint from upper to lower link hits its position limit at q2 = -0.5,
    #    causing v2 to drop to 0 immediately. v1 increases slightly due to the impact
    # 5. Both bodies fall together, with q2 fixed at -0.5 and v2 fixed at 0
    for i in 27:length(ts)
        @test qs[i][2] ≈ -0.5
        @test vs[i][2] ≈ 0 atol=1e-4
    end

    # 6. The foot hits the ground when q1 = 0.5 and q2 = -0.5, at which point
    #    v1 also drops to 0 immediately
    #
    # this takes about 0.45 seconds after the previous impact
    @test qs[27 + 44] ≈ [0.531, -0.5] atol=1e-3
    @test vs[27 + 44] ≈ [-5.253, 0.0] atol=1e-3
    @test qs[27 + 45] ≈ [0.5, -0.5] atol=1e-4

    # 7. Both bodies remain at rest at q1 = 0.5, q2 = -0.5
    for i in 73:length(ts)
        @test qs[i] ≈ [0.5, -0.5] atol=1e-4
        @test vs[i] ≈ [0, 0] atol=1e-4
    end
end
