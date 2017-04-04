using MomentumBasedControl
using RigidBodyDynamics
using BipedControlUtil
using ValkyrieRobot
using StaticArrays
using Rotations
using Base.Test

val = Valkyrie()
mechanism = val.mechanism
controller = MomentumBasedController{Float64}(mechanism; num_basis_vectors_per_contact = 4)
state = MechanismState(Float64, mechanism)

@testset "zero velocity free fall" begin
    zero!(state)
    rand_configuration!(state)
    set_joint_accel_regularization!(controller, 1.)
    set_joint_accel_regularization!(controller, val.floatingjoint, 0.)
    for side in instances(Side)
        foot = val.feet[side]
        for point in contact_points(foot)
            set_contact_regularization!(controller, point, 1.)
        end
    end
    control(controller, 0., state)

    for joint in tree_joints(mechanism)
        if joint == val.floatingjoint
            baseaccel = relative_acceleration(state, val.pelvis, root_body(mechanism), controller.result.v̇)
            baseaccel = transform(state, baseaccel, frame_after(val.floatingjoint))
            angularaccel = FreeVector3D(baseaccel.frame, baseaccel.angular)
            @test isapprox(angularaccel, FreeVector3D(frame_after(val.floatingjoint), zeros(SVector{3})), atol = 1e-7)
            linearaccel = FreeVector3D(baseaccel.frame, baseaccel.linear)
            linearaccel = transform_to_root(state, linearaccel.frame) * linearaccel
            @test isapprox(linearaccel, mechanism.gravitationalAcceleration; atol = 1e-7)
        else
            v̇joint = controller.result.v̇[velocity_range(state, joint)]
            @test isapprox(v̇joint, zeros(num_velocities(joint)); atol = 1e-7)
        end
    end
end

@testset "achievable momentum rate" begin
    srand(1)
    @assert controller.num_basis_vectors_per_contact >= 4 # test relies on this fact

    for p in linspace(0., 1., 5)
        rand!(state)
        com = center_of_mass(state)
        centroidal_to_world = Transform3D(centroidal_frame(controller), com.frame, com.v)
        world_to_centroidal = inv(centroidal_to_world)
        reset!(controller)

        # set random active contacts and random achievable wrench
        fg = world_to_centroidal * (mass(mechanism) * mechanism.gravitationalAcceleration)
        ḣdes = Wrench(zero(fg), fg)
        for foot in values(val.feet)
            for point in contact_points(foot)
                active = rand() < p
                if active
                    r = RigidBodyDynamics.Contact.location(point)
                    μ = rand()
                    normal = FreeVector3D(r.frame, normalize(randn(SVector{3})))
                    enable_contact!(controller, point, μ, normal)
                    fnormal = 0.#50. * rand()
                    μreduced = sqrt(2) / 2 * μ # due to polyhedral inner approximation; assumes 4 basis vectors or more
                    ftangential = μreduced * fnormal * rand() * cross(normal, FreeVector3D(normal.frame, normalize(randn(SVector{3}))))
                    f = fnormal * normal + ftangential
                    @assert isapprox(ftangential ⋅ normal, 0., atol = 1e-12)
                    @assert norm(f - (normal ⋅ f) * normal) ≤ μreduced * (normal ⋅ f)
                    wrench = Wrench(r, f)
                    ḣdes += transform(wrench, world_to_centroidal * transform_to_root(state, wrench.frame))
                else
                    disable_contact!(controller, point)
                end
            end
        end
        set_desired_momentum_rate!(controller, ḣdes, eye(SMatrix{6, 6}))

        control(controller, 0., state)

        # Ensure that desired momentum rate is achieved.
        ḣ = Wrench(momentum_matrix(state), controller.result.v̇) + momentum_rate_bias(state)
        ḣ = transform(ḣ, world_to_centroidal)
        @test isapprox(ḣdes, ḣ; atol = 1e-3)
    end
end
