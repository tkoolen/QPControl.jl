type ContactSettings{T}
    weight::T
    μ::T
    active::Bool
end

type MomentumBasedController{T}
    mechanism::Mechanism{T}
    num_basis_vectors_per_contact::Int64
    centroidalframe::CartesianFrame3D
    τ::Vector{T}
    ρ::Vector{T}
    result::DynamicsResult{T, T}
    momentummatrix::MomentumMatrix{Matrix{T}}
    desiredjointaccels::Dict{Joint{T}, Maybe{Vector{T}}}
    desiredspatialaccels::Dict{Pair{Pair{RigidBody{T}, RigidBody{T}}, CartesianFrame3D}, Maybe{SpatialAcceleration{T}}}
    paths::Dict{Pair{RigidBody{T}, RigidBody{T}}, TreePath{RigidBody{T}, Joint{T}}}
    jacobians::Dict{Pair{Pair{RigidBody{T}, RigidBody{T}}, CartesianFrame3D}, GeometricJacobian{Matrix{T}}}
    v̇weights::Dict{Joint{T}, Vector{T}}
    desiredmomentumrate::Wrench{T}
    momentumrateweight::SMatrix{6, 6, T, 36}
    contactingbodies::Vector{RigidBody{T}}
    contactsettings::Dict{ContactPoint, ContactSettings{T}}
    wrenchmatrix::WrenchMatrix{Matrix{T}}
    externalwrenches::Vector{Wrench{T}}
    mass::T

    function MomentumBasedController(mechanism::Mechanism; num_basis_vectors_per_contact::Int64 = 4)
        centroidalframe = CartesianFrame3D("centroidal")
        nv = num_velocities(mechanism)
        τ = zeros(num_velocities(mechanism))
        result = DynamicsResult(Float64, mechanism)
        momentummatrix = MomentumMatrix(centroidalframe, Matrix{T}(3, nv), Matrix{T}(3, nv))
        desiredjointaccels = Dict()
        desiredspatialaccels = Dict()
        paths = Dict()
        jacobians = Dict()
        v̇weights = Dict(joint => zeros(num_velocities(joint)) for joint in tree_joints(mechanism))
        desiredmomentumrate = zero(Wrench{T}, centroidalframe)
        momentumrateweight = zero(SMatrix{6, 6, T})
        contactingbodies = RigidBody{T}[]
        contactsettings = Dict{ContactPoint, ContactSettings{T}}()
        nρ = 0
        for body in bodies(mechanism)
            if !isempty(contact_points(body))
                push!(contactingbodies, body)
            end
            for point in contact_points(body)
                contactsettings[point] = ContactSettings(0., 0.8, false)
                nρ += num_basis_vectors_per_contact
            end
        end
        ρ = zeros(T, nρ)
        wrenchmatrix = WrenchMatrix(centroidalframe, Matrix{T}(3, nρ), Matrix{T}(3, nρ))
        externalwrenches = Vector{Wrench{T}}(num_bodies(mechanism))
        new(mechanism, num_basis_vectors_per_contact, centroidalframe, τ, ρ, result, momentummatrix, desiredjointaccels, desiredspatialaccels, paths, jacobians, v̇weights,
            desiredmomentumrate, momentumrateweight, contactingbodies, contactsettings, wrenchmatrix, externalwrenches, mass(mechanism))
    end
end

Base.eltype{T}(::Type{MomentumBasedController{T}}) = T
Base.eltype{T}(controller::MomentumBasedController{T}) = eltype(typeof(controller))

centroidal_frame(controller::MomentumBasedController) = controller.centroidalframe

function clear_contacts!(controller::MomentumBasedController)
    for settings in values(controller.contactsettings)
        settings.active = false
    end
end

function clear_desireds!{T}(controller::MomentumBasedController{T})
    for accel in values(controller.desiredjointaccels)
        invalidate!(accel)
    end
    for accel in controller.desiredspatialaccels
        invalidate!(accel)
    end
    controller.momentumrateweight = zeros(SMatrix{6, 6, T})
end

reset!(controller::MomentumBasedController) = (clear_contacts!(controller); clear_desireds!(controller))

function set_desired_accel!(controller::MomentumBasedController, joint::Joint, accel)
    nv = num_velocities(joint)
    copydata!(get!(() -> Maybe(Vector{eltype(accel)}(nv)), controller.desiredjointaccels, joint), accel)
    nothing
end

function set_desired_accel!(controller::MomentumBasedController, base::RigidBody, body::RigidBody, accel::SpatialAcceleration)
    @boundscheck begin
        @assert RigidBodyDynamics.is_fixed_to_body(body, accel.body)
        @assert RigidBodyDynamics.is_fixed_to_body(base, accel.base)
    end
    frame = accel.frame
    pathkey = base => body
    accelkey = pathkey => frame
    T = eltype(controller)
    p = get!(() -> path(controller.mechanism, base, body), controller.paths, pathkey)
    if !haskey(controller.jacobians, accelkey)
        nv = num_velocities(p)
        controller.jacobians[accelkey] = GeometricJacobian(default_frame(body), default_frame(base), frame, Matrix{T}(3, nv), Matrix{T}(3, nv))
    end
    maybe_accel = get!(Maybe{SpatialAcceleration{T}}, controller.desiredspatialaccels, accelkey)
    setdata!(maybe_accel, accel)
    nothing
end

set_contact_regularization!(controller::MomentumBasedController, point::ContactPoint, weight) = (controller.contactsettings[point].weight = weight; nothing)
set_joint_accel_regularization!(controller::MomentumBasedController, joint::Joint, weights) = (controller.v̇weights[joint] .= weights; nothing)
function set_joint_accel_regularization!(controller::MomentumBasedController, weight)
    for weights in values(controller.v̇weights)
        weights .= weight
    end
end

set_contact_active!(controller::MomentumBasedController, point::ContactPoint, active::Bool) = controller.contactsettings[point].active = active
set_friction_coefficient!(controller::MomentumBasedController, point::ContactPoint, μ) = controller.contactsettings[point] = μ

function set_desired_momentum_rate!(controller::MomentumBasedController, rate::Wrench, weight::SMatrix{6, 6})
    @framecheck centroidal_frame(controller) rate.frame
    controller.desiredmomentumrate = rate
    controller.momentumrateweight = weight
end

function update_wrench_matrix!(Q::WrenchMatrix, contactingbodies, contactsettings, num_basis_vectors_per_contact::Int64, state::MechanismState, world_to_centroidal::Transform3D)
    @framecheck Q.frame world_to_centroidal.to
    rootframe = world_to_centroidal.from
    T = eltype(typeof(Q))
    col = 1
    Δθ = 2 * π / num_basis_vectors_per_contact
    for body in contactingbodies
        for point in contact_points(body)
            settings = contactsettings[point]
            if settings.active
                r = location(point)
                frame = r.frame
                r = world_to_centroidal * transform(state, r, rootframe)
                μ = settings.μ
                for i = 0 : num_basis_vectors_per_contact - 1
                    θ = i * Δθ
                    # TODO: would be better not to use cos and sin
                    β = FreeVector3D(location(point).frame, normalize(SVector(μ * cos(θ), μ * sin(θ), one(T)))) # FIXME: assumes normal is z axis of frame in which point is expressed
                    β = world_to_centroidal * transform(state, β, rootframe)
                    Qcol = Wrench(r, β)
                    @framecheck Q.frame Qcol.frame
                    view(Q.angular, :, col)[:] = Qcol.angular
                    view(Q.linear, :, col)[:] = Qcol.linear
                    col += 1
                end
            else
                for i = 0 : num_basis_vectors_per_contact - 1
                    view(Q.angular, :, col)[:] = 0
                    view(Q.linear, :, col)[:] = 0
                    col += 1
                end
            end
        end
    end
end

function back_out_external_wrenches!(externalwrenches::Vector, Q::WrenchMatrix, ρ::Vector, num_basis_vectors_per_contact::Int64, contactingbodies)
    T = eltype(ρ)
    for i in eachindex(externalwrenches)
        externalwrenches[i] = zero(Wrench{T}, Q.frame)
    end
    col = 1
    for body in contactingbodies
        for i = 1 : num_basis_vectors_per_contact * length(contact_points(body))
            # TODO: allocations
            angular = SVector{3}(view(Q.angular, :, col)) * ρ[col]
            linear = SVector{3}(view(Q.linear, :, col)) * ρ[col]
            externalwrenches[vertex_index(body)] += Wrench(Q.frame, angular, linear)
            col += 1
        end
    end
end

function control(controller::MomentumBasedController, t, state)
    # Create model
    solver = Gurobi.GurobiSolver()
    Gurobi.setparameters!(solver, Silent = true)
    model = Model(solver = solver)

    mechanism = controller.mechanism
    nv = num_velocities(state)
    nρ = num_cols(controller.wrenchmatrix)

    @variable(model, v̇[1 : nv])
    @variable(model, ρ[1 : nρ] >= 0)

    # Compute centroidal frame transform
    com = center_of_mass(state)
    centroidal_to_world = Transform3D(controller.centroidalframe, com.frame, com.v)
    world_to_centroidal = inv(centroidal_to_world)

    # Momentum stuff
    A = controller.momentummatrix
    momentum_matrix!(A, state, world_to_centroidal)
    Ȧv = transform(momentum_rate_bias(state), world_to_centroidal)
    ḣd = controller.desiredmomentumrate
    b = ḣd - Ȧv
    ḣerror = Wrench(A, v̇) - b
    ḣerrorvec = [ḣerror.angular; ḣerror.linear]

    # Gravitational wrench
    m = controller.mass
    fg = transform(m * mechanism.gravitationalAcceleration, world_to_centroidal)
    Wg = Wrench(zero(fg), fg)

    # Wrench matrix
    Q = controller.wrenchmatrix
    update_wrench_matrix!(Q, controller.contactingbodies, controller.contactsettings, controller.num_basis_vectors_per_contact, state, world_to_centroidal)

    # Desired joint accelerations
    for (joint, maybe_accel) in controller.desiredjointaccels
        if isvalid(maybe_accel)
            accel = get(maybe_accel)
            v̇range = velocity_range(state, joint)
            @constraint(model, view(v̇, v̇range) .== accel)
        end
    end

    # Desired spatial accelerations
    for (accelkey, maybe_accel) in controller.desiredspatialaccels
        if isvalid(maybe_accel)
            pair, frame = accelkey
            base, body = pair
            accel = get(maybe_accel)
            path = controller.paths[pair]
            tf = inv(transform_to_root(state, accel.frame))
            J = controller.jacobians[accelkey]
            J̇v = transform(state, -bias_acceleration(state, base) + bias_acceleration(state, body), tf.to)
            geometric_jacobian!(J, state, path, tf)
            @framecheck J.frame J̇v.frame
            @framecheck J.frame accel.frame
            # @constraint(model, jacobian.angular * ) # TODO: need velocities for joints on path
        end
    end

    # Newton; TODO: add external wrenches to compensate for
    @framecheck A.frame Ȧv.frame
    @framecheck A.frame Wg.frame
    @framecheck A.frame Q.frame
    @constraint(model, A.angular * v̇ + Ȧv.angular .== Wg.angular + Q.angular * ρ)
    @constraint(model, A.linear * v̇ + Ȧv.linear .== Wg.linear + Q.linear * ρ)

    # Objective
    momentumrateterm = @expression(model, dot(ḣerrorvec, controller.momentumrateweight * ḣerrorvec))
    v̇regularization = @expression(model, sum(sum(controller.v̇weights[joint] .* view(v̇, velocity_range(state, joint)).^2) for joint in tree_joints(mechanism)))
    ρregularization = zero(typeof(momentumrateterm))
    index = 1
    for body in controller.contactingbodies
        for point in contact_points(body)
            weight = controller.contactsettings[point].weight
            for i = 1 : controller.num_basis_vectors_per_contact
                ρregularization += weight * ρ[index]^2
                index += 1
            end
        end
    end
    @objective(model, Min, momentumrateterm + v̇regularization + ρregularization)

    # Solve
    status = solve(model)

    # Inverse dynamics
    controller.ρ .= getvalue(ρ)
    result = controller.result
    result.v̇ .= getvalue(v̇)
    externalwrenches = controller.externalwrenches
    τ = controller.τ
    back_out_external_wrenches!(externalwrenches, Q, controller.ρ, controller.num_basis_vectors_per_contact, controller.contactingbodies)
    map!(wrench -> transform(wrench, centroidal_to_world), externalwrenches)
    inverse_dynamics!(τ, result.jointwrenches, result.accelerations, state, result.v̇, externalwrenches)

    τ
end
