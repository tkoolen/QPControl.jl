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
    @test_broken allocs == 0

    result = DynamicsResult(mechanism)
    dynamics!(result, state, τ)
    for joint in tree_joints(mechanism)
        @test result.v̇[joint] ≈ tasks[joint].desired atol = 1e-10
    end
end

function set_up_valkyrie_contacts!(controller::MomentumBasedController)
    valmechanism = controller.state.mechanism
    for body in bodies(valmechanism)
        for point in RBD.contact_points(body)
            position = RBD.Contact.location(point)
            normal = FreeVector3D(position.frame, 0.0, 0.0, 1.0)
            μ = point.model.friction.μ
            addcontact!(controller, body, ContactPoint(position, normal, μ))
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
            @test isapprox(angularaccel, FreeVector3D(frame_after(floatingjoint), zeros(SVector{3})), atol = 1e-6)
            linearaccel = FreeVector3D(baseaccel.frame, baseaccel.linear)
            linearaccel = transform_to_root(state, linearaccel.frame) * linearaccel
            @test isapprox(linearaccel, mechanism.gravitational_acceleration; atol = 1e-6)
        else
            v̇joint = controller.result.v̇[velocity_range(state, joint)]
            @test isapprox(v̇joint, zeros(num_velocities(joint)); atol = 1e-6)
        end
    end
end

# @testset "achievable momentum rate" begin
#     val = Valkyrie()
#     mechanism = val.mechanism
#     floatingjoint = val.basejoint
#     controller = MomentumBasedController(mechanism, defaultoptimizer())
#     contacts = add_mechanism_contacts!(controller)
#     jointacceltasks = add_mechanism_joint_accel_tasks!(controller)
#     state = MechanismState(mechanism)
#     τ = similar(velocity(state))

#     srand(1)
#     for p in linspace(0., 1., 5)
#         rand!(state)
#         com = center_of_mass(state)
#         centroidal_to_world = Transform3D(centroidal_frame(controller), com.frame, com.v)
#         world_to_centroidal = inv(centroidal_to_world)
#         reset!(controller)

#         # set random active contacts and random achievable wrench
#         fg = world_to_centroidal * (mass(mechanism) * mechanism.gravitational_acceleration)
#         ḣdes = Wrench(zero(fg), fg)
#         for foot in values(val.feet)
#             for contactsettings in contacts[foot]
#                 active = rand() < p
#                 if active
#                     @assert num_basis_vectors(contactsettings) >= 4 # test relies on this fact
#                     r = contactsettings.point
#                     μ = rand()
#                     normal = FreeVector3D(r.frame, normalize(randn(SVector{3})))
#                     set!(contactsettings, 0., μ, normal)
#                     fnormal = 50. * rand()
#                     μreduced = sqrt(2) / 2 * μ # due to polyhedral inner approximation; assumes 4 basis vectors or more
#                     ftangential = μreduced * fnormal * rand() * cross(normal, FreeVector3D(normal.frame, normalize(randn(SVector{3}))))
#                     f = fnormal * normal + ftangential
#                     @assert isapprox(ftangential ⋅ normal, 0., atol = 1e-12)
#                     @assert norm(f - (normal ⋅ f) * normal) ≤ μreduced * (normal ⋅ f)
#                     wrench = Wrench(r, f)
#                     ḣdes += transform(wrench, world_to_centroidal * transform_to_root(state, wrench.frame))
#                 else
#                     disable!(contactsettings)
#                 end
#             end
#         end
#         add!(controller, CentroidalMomentumRateTask(ḣdes, eye(SMatrix{3, 3}), eye(SMatrix{3, 3}), 1.))
#         controller(τ, 0., state)

#         # Ensure that desired momentum rate is achieved.
#         ḣ = Wrench(momentum_matrix(state), controller.result.v̇) + momentum_rate_bias(state)
#         ḣ = transform(ḣ, world_to_centroidal)
#         @test isapprox(ḣdes, ḣ; atol = 1e-3)
#     end
# end

# @testset "spatial acceleration" begin
#     val = Valkyrie()
#     mechanism = val.mechanism
#     floatingjoint = val.basejoint
#     state = MechanismState(mechanism)
#     τ = similar(velocity(state))

#     srand(533)
#     rand!(state)
#     controller = MomentumBasedController(mechanism, defaultoptimizer())
#     body = val.feet[left]
#     base = val.palms[right]
#     frame = default_frame(base)
#     task = SpatialAccelerationTask(num_velocities(mechanism), path(mechanism, base, body), frame, eye(3), eye(3))
#     add!(controller, task)

#     result = DynamicsResult(mechanism)
#     accels = result.accelerations
#     for weight in [1., Inf]
#         reset!(controller)
#         if isinf(weight)
#             regularize_joint_accels!(controller, 1.)
#         end
#         desiredaccel = rand(SpatialAcceleration{Float64}, default_frame(body), default_frame(base), frame)
#         set!(task, desiredaccel, weight)
#         controller(τ, 0., state)
#         v̇ = controller.result.v̇
#         spatial_accelerations!(accels, state, v̇)
#         accel = relative_acceleration(accels, body, base)
#         accel = transform(state, accel, frame)
#         @test isapprox(accel, desiredaccel, atol = 1e-8)
#     end
# end
