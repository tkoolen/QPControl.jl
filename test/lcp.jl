using Gurobi

function add_base_joint!(mechanism, label, direction)
    frame = CartesianFrame3D("base_$(label)_dummy")
    inertia = SpatialInertia(frame, SDiagonal(0., 0, 0), SVector(0., 0, 0), 0.0)
    body = RigidBody(inertia)
    joint = Joint("base_$(label)", Prismatic(direction))
    effort_bounds(joint) .= RigidBodyDynamics.Bounds(0, 0)
    attach!(mechanism, last(bodies(mechanism)), body, joint)
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
    attach!(mechanism, last(bodies(mechanism)), body, joint)
    mechanism
end

function add_floor!(mechanism, orientation, μ)
    world = root_body(mechanism)
    floor = HalfSpace3D(Point3D(default_frame(world), 0., 0, 0), FreeVector3D(default_frame(world), orientation * SVector(0., 0, 1)))
    add_environment_primitive!(mechanism, floor)

    body = findbody(mechanism, "core")
    contactmodel = SoftContactModel(hunt_crossley_hertz(k = 500e3), ViscoelasticCoulombModel(μ, 20e3, 100.))
    add_contact_point!(body, Contact.ContactPoint(Point3D(default_frame(body), 0., 0, 0), contactmodel))
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

using GLPK

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
                floor = add_floor!(mechanism, orientation, μ)

                # Create an MPC "controller", which we will use as a simple
                # one-step LCP simulator
                mpc = mpc_lcp_feasibility_problem(mechanism, optimizer, floor, Δt)

                state = MechanismState(mechanism)
                # Start off the ground a bit, to be sure we can dissipate
                # velocity correctly
                set_configuration!(state, findjoint(mechanism, "base_z"), 0.1)

                # Simulate the system by repeatedly running the MPC controller
                # and then updating the state based on the expected state
                # at the end of the first MPC stage.
                ts = Float64[]
                qs = Vector{Float64}[]
                vs = Vector{Float64}[]
                τ = similar(velocity(state))
                t_final = 1.0
                for t in 0:Δt:t_final
                    push!(ts, t)
                    push!(qs, configuration(state))
                    push!(vs, velocity(state))
                    mpc(τ, t, state)
                    # Assume the dynamics behave exactly as modeled in the optimization
                    set_configuration!(state, SimpleQP.value.(mpc.qpmodel, mpc.stages[1].q))
                    set_velocity!(state, SimpleQP.value.(mpc.qpmodel, mpc.stages[1].v))
                end

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



