using MomentumBasedControl
using MomentumBasedControl.PDControl
using RigidBodyDynamics
using RigidBodyDynamics.OdeIntegrators
using BipedControlUtil
using ValkyrieRobot
using StaticArrays
using Rotations
using Base.Test

@testset "pd" begin
    srand(512)
    @test pd(PDGains(1, 2), 3, 4) == -11

    let
        gains = PDGains(1, 2)
        e = [3; 4; 5]
        ė = [6; 7; 8]
        @test pd(gains, e, ė) == ((e, ė) -> pd(gains, e, ė)).(e, ė)
    end

    let
        gains = PDGains(5, 6)
        x = SVector(1, 2)
        ẋ = SVector(5, 6)
        @test pd(gains, x, ẋ) == pd(gains, x, zero(x), ẋ, zero(ẋ))
    end

    let
        gains = PDGains(2, 3)
        e = RotX(rand())
        ė = rand() * SVector(1, 0, 0)
        @test pd(gains, e, ė)[1] ≈ pd(gains, e.theta, ė[1])
    end

    let
        mechanism = rand_floating_tree_mechanism(Float64) # single floating body
        joint = first(tree_joints(mechanism))
        body = successor(joint, mechanism)
        base = root_body(mechanism)

        state = MechanismState(Float64, mechanism)
        rand!(state)

        Rdes = rand(RotMatrix{3})
        ωdes = zeros(SVector{3})
        gains = PDGains(100, 20)

        result = DynamicsResult(Float64, mechanism)
        dyn! = (vd, sd, t, state) -> begin
            H = transform_to_root(state, body)
            T = transform(twist_wrt_world(state, body), inv(H))
            R = rotation(H)
            ω = T.angular
            ωddes = pd(gains, R, Rdes, ω, ωdes)
            v̇des = Array([ωddes; zeros(ωddes)])
            τ = inverse_dynamics(state, v̇des)
            dynamics!(result, state, τ)
            copy!(vd, result.v̇)
            copy!(sd, result.ṡ)
        end

        finaltime = 3.
        Δt = 1e-3
        tableau = runge_kutta_4(Float64)
        storage = ExpandingStorage{Float64}(ceil(Int64, finaltime / Δt * 1.001))
        integrator = MuntheKaasIntegrator(dyn!, tableau, storage)
        integrate(integrator, state, finaltime, Δt)

        H = transform_to_root(state, body)
        T = transform(twist_wrt_world(state, body), inv(H))
        R = rotation(H)
        ω = T.angular

        @test isapprox(R * Rdes', eye(Rdes); atol = 1e-8)
        @test isapprox(ω, zero(ω); atol = 1e-8)
    end

    let
        mechanism = rand_floating_tree_mechanism(Float64) # single floating body
        joint = first(tree_joints(mechanism))
        body = successor(joint, mechanism)
        base = root_body(mechanism)

        baseframe = default_frame(base)
        desiredframe = CartesianFrame3D("desired")
        actualframe = frame_after(joint)

        state = MechanismState(Float64, mechanism)
        rand!(state)

        xdes = rand(Transform3D, desiredframe, baseframe)
        vdes = zero(Twist{Float64}, desiredframe, baseframe, actualframe)
        gains = DoubleGeodesicPDGains(baseframe, PDGains(100, 20), PDGains(100, 20)) # world-fixed gains

        result = DynamicsResult(Float64, mechanism)
        dyn! = (vd, sd, t, state) -> begin
            x = transform_to_root(state, body)
            invx = inv(x)
            v = transform(twist_wrt_world(state, body), invx)
            v̇des = pd(transform(gains, invx), x, xdes, v, vdes)
            τ = inverse_dynamics(state, Array(v̇des))
            dynamics!(result, state, τ)
            copy!(vd, result.v̇)
            copy!(sd, result.ṡ)
            nothing
        end

        finaltime = 3.
        Δt = 1e-3
        tableau = runge_kutta_4(Float64)
        storage = ExpandingStorage{Float64}(ceil(Int64, finaltime / Δt * 1.001))
        integrator = MuntheKaasIntegrator(dyn!, tableau, storage)
        integrate(integrator, state, finaltime, Δt)

        x = transform_to_root(state, body)
        v = transform(twist_wrt_world(state, body), inv(x))
        @test isapprox(x, xdes * eye(Transform3D, actualframe, desiredframe), atol = 1e-6)
        @test isapprox(v, vdes + zero(Twist{Float64}, actualframe, desiredframe, actualframe), atol = 1e-6)
    end
end

val = Valkyrie()
mechanism = val.mechanism
controller = MomentumBasedController{Float64}(mechanism)
contacts = add_mechanism_contacts!(controller)
jointacceltasks = add_mechanism_joint_accel_tasks!(controller)
state = MechanismState(Float64, mechanism)
controllerstate = MomentumBasedControllerState(state)

@testset "zero velocity free fall" begin
    srand(5354)
    zero!(state)
    rand_configuration!(state)
    regularize_joint_accels!(controller, 1.)
    disable!(jointacceltasks[val.floatingjoint])
    control(controller, 0., controllerstate)

    for joint in tree_joints(mechanism)
        if joint == val.floatingjoint
            baseaccel = relative_acceleration(state, val.pelvis, root_body(mechanism), controller.result.v̇)
            baseaccel = transform(state, baseaccel, frame_after(val.floatingjoint))
            angularaccel = FreeVector3D(baseaccel.frame, baseaccel.angular)
            @test isapprox(angularaccel, FreeVector3D(frame_after(val.floatingjoint), zeros(SVector{3})), atol = 1e-6)
            linearaccel = FreeVector3D(baseaccel.frame, baseaccel.linear)
            linearaccel = transform_to_root(state, linearaccel.frame) * linearaccel
            @test isapprox(linearaccel, mechanism.gravitationalAcceleration; atol = 1e-6)
        else
            v̇joint = controller.result.v̇[velocity_range(state, joint)]
            @test isapprox(v̇joint, zeros(num_velocities(joint)); atol = 1e-6)
        end
    end
end

@testset "achievable momentum rate" begin
    srand(1)
    for p in linspace(0., 1., 5)
        rand!(state)
        com = center_of_mass(state)
        centroidal_to_world = Transform3D(centroidal_frame(controller), com.frame, com.v)
        world_to_centroidal = inv(centroidal_to_world)
        reset!(controller)
        reset!(controllerstate)

        # set random active contacts and random achievable wrench
        fg = world_to_centroidal * (mass(mechanism) * mechanism.gravitationalAcceleration)
        ḣdes = Wrench(zero(fg), fg)
        for foot in values(val.feet)
            for contactsettings in contacts[foot]
                active = rand() < p
                if active
                    @assert num_basis_vectors(contactsettings) >= 4 # test relies on this fact
                    r = contactsettings.point
                    μ = rand()
                    normal = FreeVector3D(r.frame, normalize(randn(SVector{3})))
                    set!(contactsettings, 0., μ, normal)
                    fnormal = 50. * rand()
                    μreduced = sqrt(2) / 2 * μ # due to polyhedral inner approximation; assumes 4 basis vectors or more
                    ftangential = μreduced * fnormal * rand() * cross(normal, FreeVector3D(normal.frame, normalize(randn(SVector{3}))))
                    f = fnormal * normal + ftangential
                    @assert isapprox(ftangential ⋅ normal, 0., atol = 1e-12)
                    @assert norm(f - (normal ⋅ f) * normal) ≤ μreduced * (normal ⋅ f)
                    wrench = Wrench(r, f)
                    ḣdes += transform(wrench, world_to_centroidal * transform_to_root(state, wrench.frame))
                else
                    disable!(contactsettings)
                end
            end
        end
        add!(controller, MomentumRateTask(ḣdes, eye(SMatrix{3, 3}), eye(SMatrix{3, 3}), 1.))
        control(controller, 0., controllerstate)

        # Ensure that desired momentum rate is achieved.
        ḣ = Wrench(momentum_matrix(state), controller.result.v̇) + momentum_rate_bias(state)
        ḣ = transform(ḣ, world_to_centroidal)
        @test isapprox(ḣdes, ḣ; atol = 1e-3)
    end
end

@testset "spatial acceleration" begin
    srand(533)
    rand!(state)
    controller = MomentumBasedController{Float64}(mechanism)
    body = val.feet[left]
    base = val.palms[right]
    frame = default_frame(base)
    task = SpatialAccelerationTask(path(mechanism, base, body), frame, eye(3), eye(3))
    add!(controller, task)

    for weight in [1., Inf]
        reset!(controller)
        reset!(controllerstate)
        if isinf(weight)
            regularize_joint_accels!(controller, 1.)
        end
        desiredaccel = rand(SpatialAcceleration{Float64}, default_frame(body), default_frame(base), frame)
        set!(task, desiredaccel, weight)
        control(controller, 0., controllerstate)
        v̇ = controller.result.v̇
        accel = relative_acceleration(state, body, base, v̇)
        accel = transform(state, accel, frame)
        @test isapprox(accel, desiredaccel, atol = 1e-8)
    end
end

@testset "Δt" begin
    srand(2533)
    rand!(state)
    Δt = 1.
    controller = MomentumBasedController{Float64}(mechanism, Δt)
    contacts = add_mechanism_contacts!(controller)
    tasks = add_mechanism_joint_accel_tasks!(controller)

    randomize = () -> begin
        for task in values(tasks)
            set!(task, rand(), 1.)
        end
        for (body, contactsettings) in contacts
            for contactpointsetting in contactsettings
                frame = default_frame(body)
                set!(contactpointsetting, 1e-3, 100., rand(FreeVector3D, Float64, frame))
            end
        end
        rand!(state)
    end

    controllerstate = MomentumBasedControllerState(state)
    t0 = 2.
    randomize()
    τ1 = copy(control(controller, t0, controllerstate))
    for t in linspace(t0, t0 + Δt - 1e-3, 10)
        randomize()
        τ = control(controller, t, controllerstate)
        @test all(τ1 .== τ)
    end
    randomize()
    @test any(τ1 .!= control(controller, t0 + Δt + 1e-3, controllerstate))
end
