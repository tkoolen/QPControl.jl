using MomentumBasedControl
using RigidBodyDynamics
using BipedControlUtil
using ValkyrieRobot
using StaticArrays
using Base.Test

val = Valkyrie()
mechanism = val.mechanism
controller = MomentumBasedController{Float64}(mechanism)
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
