type MomentumBasedController{T}
    mechanism::Mechanism{T}
    num_basis_vectors_per_contact::Int64
    centroidalframe::CartesianFrame3D
    τ::Vector{T}
    result::DynamicsResult{T, T}
    momentummatrix::MomentumMatrix{Matrix{T}}
    desiredjointaccels::Dict{Joint{T}, Maybe{Vector{T}}}
    desiredspatialaccels::Dict{Pair{RigidBody{T}, RigidBody{T}}, Maybe{SpatialAcceleration{T}}}
    paths::Dict{Pair{RigidBody{T}, RigidBody{T}}, TreePath{RigidBody{T}, Joint{T}}}
    jacobians::Dict{Pair{RigidBody{T}, RigidBody{T}}, GeometricJacobian{Matrix{T}}}
    v̇weights::Vector{T}
    desiredmomentumrate::Wrench{T}
    momentumrateweight::SMatrix{6, 6, T, 36}
    contactweights::OrderedDict{ContactPoint, T}
    wrenchmatrix::WrenchMatrix{Matrix{T}}
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
        v̇weights = zeros(nv)
        desiredmomentumrate = zero(Wrench{T}, centroidalframe)
        momentumrateweight = zero(SMatrix{6, 6, T})
        contactweights = OrderedDict{ContactPoint, T}()
        nρ = 0
        for body in bodies(mechanism)
            for point in contact_points(body)
                contactweights[point] = 0
                nρ += num_basis_vectors_per_contact
            end
        end
        wrenchmatrix = WrenchMatrix(centroidalframe, Matrix{T}(3, nρ), Matrix{T}(3, nρ))
        new(mechanism, num_basis_vectors_per_contact, centroidalframe, τ, result, momentummatrix, desiredjointaccels, desiredspatialaccels, paths, jacobians, v̇weights,
            desiredmomentumrate, momentumrateweight, contactweights, wrenchmatrix, mass(mechanism))
    end
end

Base.eltype{T}(::Type{MomentumBasedController{T}}) = T
Base.eltype{T}(controller::MomentumBasedController{T}) = eltype(typeof(controller))

function clear_contacts!(controller::MomentumBasedController)
    for point in keys(controller.contactweights)
        controller.contactweights[point] = 0
    end
end

function clear_desireds!{T}(controller::MomentumBasedController{T})
    map!(invalidate!, values(controller.desiredjointaccels))
    map!(invalidate!, values(controller.desiredspatialaccels))
    controller.momentumrateweight = zeros(SMatrix{6, 6, T})
end

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
    pair = base => body
    T = eltype(controller)
    if !haskey(controller.paths, pair)
        p = path(controller.mechanism, base, body)
        controller.paths[pair] = p
        nv = num_velocities(p)
         # TODO: what if accel's frames are different but body and base are the same?
        controller.jacobians[pair] = GeometricJacobian(default_frame(body), default_frame(base), default_frame(body), Matrix{T}(3, nv), Matrix{T}(3, nv))
    end
    maybe_accel = get!(Maybe{SpatialAcceleration{T}}, controller.desiredspatialaccels, pair)
    setdata!(maybe_accel, accel)
    nothing
end

function set_desired_momentum_rate!(controller::MomentumBasedController, rate::Wrench, weight::AbstractMatrix)
    @framecheck controller.rate.frame rate.frame
    controller.rate = rate
    controller.momentumrateweight = SMatrix{6, 6}(weight)
end

set_contact_weight(controller::MomentumBasedController, point::ContactPoint, weight::Number) = (controller.contactweights[point] = weight)
set_joint_accel_weights(controller::MomentumBasedController, weight::Number) = (controller.v̇weights[:] = weight)

centroidal_frame(controller::MomentumBasedController) = controller.centroidalframe

function update_wrench_matrix!(controller::MomentumBasedController, world_to_centroidal::Transform3D)
    # TODO
end

function control(controller::MomentumBasedController, t, state)
    # Create model
    solver = Gurobi.GurobiSolver()
    Gurobi.setparameters!(solver, Silent = true)
    model = Model(solver = solver)

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
    g = RigidBodyDynamics.gravitational_spatial_acceleration(controller.mechanism)
    Wg = Wrench(g.frame, zero(g.linear), m * g.linear)
    Wg = transform(Wg, world_to_centroidal)

    # Wrench matrix
    update_wrench_matrix!(controller, world_to_centroidal)
    Q = controller.wrenchmatrix

    # Desired joint accelerations
    for (joint, maybe_accel) in controller.desiredjointaccels
        if isvalid(maybe_accel)
            accel = get(maybe_accel)
            v̇range = velocity_range(state, joint)
            @constraint(model, view(v̇, v̇range) .== accel)
        end
    end

    # Desired spatial accelerations
    for (pair, maybe_accel) in controller.desiredspatialaccels
        if isvalid(maybe_accel)
            base, body = pair
            accel = get(maybe_accel)
            path = controller.paths[pair]
            tf = inv(transform_to_root(state, accel.frame))
            J = controller.jacobians[pair]
            J̇v = transform(state, -bias_acceleration(state, base) + bias_acceleration(state, body), tf.to)
            geometric_jacobian!(J, state, path, tf)
            @framecheck J.frame J̇v.frame
            @framecheck J.frame accel.frame
            # @constraint(model, jacobian.angular * ) # TODO: need velocities for joints on path
        end
    end

    # Newton; TODO: external wrenches to compensate for
    @framecheck A.frame Ȧv.frame
    @framecheck A.frame Wg.frame
    @framecheck A.frame Q.frame
    @constraint(model, A.angular * v̇ + Ȧv.angular .== Wg.angular + Q.angular * ρ)
    @constraint(model, A.linear * v̇ + Ȧv.linear .== Wg.linear + Q.linear * ρ)

    # Objective
    v̇regularization = @expression(model, sum(controller.v̇weights[i] * v̇[i]^2 for i = 1 : nv))
    momentumrateterm = @expression(model, dot(ḣerrorvec, controller.momentumrateweight * ḣerrorvec))
    @objective(model, Min, momentumrateterm + v̇regularization) # TODO: regularize ρ

    # Solve
    status = solve(model)

    # Inverse dynamics
    result = controller.result
    result.v̇ = getvalue(v̇)
    # TODO: external wrenches!

    inverse_dynamics!(controller.τ, result.jointwrenches, result.accelerations, state, result.v̇)
    controller.τ
end
