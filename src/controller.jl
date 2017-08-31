# TODO: max normal force

type MomentumBasedControllerState{T<:Number, S<:MechanismState}
    next_control_time::T
    τ::Vector{T}
    mechstate::S
end

function MomentumBasedControllerState(mechstate::MechanismState)
    τ = fill(NaN, num_velocities(mechstate))
    MomentumBasedControllerState(-Inf, τ, mechstate)
end

reset!(state::MomentumBasedControllerState) = (state.next_control_time = -Inf)

type MomentumBasedController{T}
    mechanism::Mechanism{T}
    Δt::T
    centroidalframe::CartesianFrame3D
    mass::T

    # Cache variables and tasks. TODO: consider moving these to state
    ρ::Vector{T}
    result::DynamicsResult{T, T}
    momentummatrix::MomentumMatrix{Matrix{T}}
    wrenchmatrix::WrenchMatrix{Matrix{T}}
    externalwrenches::RigidBodyDynamics.BodyDict{T, Wrench{T}}
    contactsettings::Vector{ContactSettings{T}}
    spatialacceltasks::Vector{SpatialAccelerationTask{T}}
    jointacceltasks::Vector{JointAccelerationTask{T}}
    momentumratetask::MomentumRateTask{T}

    function MomentumBasedController(mechanism::Mechanism{T}, Δt = 5e-3) where {T}
        centroidalframe = CartesianFrame3D("centroidal")
        m = mass(mechanism)
        nv = num_velocities(mechanism)
        ρ = zeros(T, 0)
        result = DynamicsResult{T}(mechanism)
        momentummatrix = MomentumMatrix(centroidalframe, Matrix{T}(3, nv), Matrix{T}(3, nv))
        wrenchmatrix = WrenchMatrix(centroidalframe, Matrix{T}(3, 0), Matrix{T}(3, 0))
        rootframe = root_frame(mechanism)
        externalwrenches = RigidBodyDynamics.BodyDict{T, Wrench{T}}(b => zero(Wrench{T}, rootframe) for b in bodies(mechanism))
        contactsettings = Vector{ContactSettings{T}}()
        spatialacceltasks = Vector{SpatialAccelerationTask{T}}()
        jointacceltasks = Vector{JointAccelerationTask{T}}()
        momentumratetask = MomentumRateTask(T, centroidalframe)

        new{T}(mechanism, Δt, centroidalframe, m, ρ, result, momentummatrix, wrenchmatrix, externalwrenches,
            contactsettings, spatialacceltasks, jointacceltasks, momentumratetask)
    end
end

Base.eltype(::Type{MomentumBasedController{T}}) where {T} = T
Base.eltype(controller::MomentumBasedController{T}) where {T} = eltype(typeof(controller))
centroidal_frame(controller::MomentumBasedController) = controller.centroidalframe

function add_contact!(controller::MomentumBasedController, body::RigidBody, point::Point3D, num_basis_vectors::Int64)
    T = eltype(controller)
    nρ = length(controller.ρ)
    nρnew = nρ + num_basis_vectors
    ρrange = nρ + 1 : nρnew
    settings = ContactSettings(body, point, ρrange)
    push!(controller.contactsettings, settings)
    resize!(controller.ρ, nρnew)
    controller.wrenchmatrix = WrenchMatrix(controller.centroidalframe, Matrix{T}(3, nρnew), Matrix{T}(3, nρnew))
    settings
end

function add_contacts!(controller::MomentumBasedController, body::RigidBody, num_basis_vectors::Int64)
    T = eltype(controller)
    ret = Vector{ContactSettings{T}}()
    for point in contact_points(body)
        push!(ret, add_contact!(controller, body, RigidBodyDynamics.Contact.location(point), num_basis_vectors))
    end
    ret
end

function add_mechanism_contacts!(controller::MomentumBasedController, num_basis_vectors::Int64 = 4)
    T = eltype(controller)
    ret = Dict{RigidBody{T}, Vector{ContactSettings{T}}}()
    for body in bodies(controller.mechanism)
        ret[body] = add_contacts!(controller, body, num_basis_vectors)
    end
    ret
end

add!(controller::MomentumBasedController, task::SpatialAccelerationTask) = push!(controller.spatialacceltasks, task)
add!(controller::MomentumBasedController, task::JointAccelerationTask) = push!(controller.jointacceltasks, task)
add!(controller::MomentumBasedController, task::MomentumRateTask) = controller.momentumratetask = task

function add_mechanism_joint_accel_tasks!(controller::MomentumBasedController)
    T = eltype(controller)
    ret = Dict{GenericJoint{T}, JointAccelerationTask{T}}()
    for joint in tree_joints(controller.mechanism)
        ret[joint] = task = JointAccelerationTask(joint)
        add!(controller, task)
    end
    ret
end

clear_contacts!(controller::MomentumBasedController) = foreach(disable!, controller.contactsettings)

function clear_desireds!(controller::MomentumBasedController{T}) where {T}
    foreach(disable!, controller.spatialacceltasks)
    foreach(disable!, controller.jointacceltasks)
    controller.momentumratetask = MomentumRateTask(T, controller.centroidalframe)
end

reset!(controller::MomentumBasedController) = (clear_contacts!(controller); clear_desireds!(controller))

function regularize_joint_accels!(controller::MomentumBasedController, weight)
    for task in controller.jointacceltasks
        set!(task, 0, weight)
    end
end

function pd_center_of_mass!(controller::MomentumBasedController, gains::PDGains, state::MechanismState,
        desired_com::Point3D, desired_com_velocity::FreeVector3D, weight::Number)
    m = controller.mass
    c = center_of_mass(state)
    world_to_centroidal = Transform3D(c.frame, centroidal_frame(controller), -c.v)
    h = transform(momentum(state), world_to_centroidal)

    cdes = world_to_centroidal * transform(state, desired_com, world_to_centroidal.from)
    ċdes = world_to_centroidal * transform(state, desired_com_velocity, world_to_centroidal.from)

    c = Point3D(centroidal_frame(controller), zero(c.v))
    ċ = FreeVector3D(h.frame, h.linear / m)

    c̈des = pd(gains, c, cdes, ċ, ċdes)
    linear_momentum_rate_des = m * c̈des

    momentum_rate_des = Wrench(zero(linear_momentum_rate_des), linear_momentum_rate_des)
    add!(controller, MomentumRateTask(momentum_rate_des, zero(SMatrix{3, 3}), eye(SMatrix{3, 3}), weight))
end

function update_wrench_matrix!(controller::MomentumBasedController, state::MechanismState, world_to_centroidal::Transform3D)
    Q = controller.wrenchmatrix
    @framecheck Q.frame world_to_centroidal.to
    rootframe = world_to_centroidal.from
    T = eltype(typeof(Q))
    for settings in controller.contactsettings
        ρrange = settings.ρrange
        if isenabled(settings)
            # TODO: consider making the following a method of ContactSettings
            body_to_centroidal = world_to_centroidal * transform_to_root(state, settings.body)
            r = body_to_centroidal * settings.point
            n = body_to_centroidal * settings.normal
            μ = settings.μ
            rot = Rotations.rotation_between(SVector(0., 0., 1.), n.v)
            Δθ = 2 * π / length(ρrange)
            for i = 1 : length(ρrange)
                θ = (i - 1) * Δθ
                β = FreeVector3D(n.frame, rot * normalize(SVector(μ * cos(θ), μ * sin(θ), one(T))))
                Qcol = Wrench(r, β)
                @framecheck Q.frame Qcol.frame
                view(Q.angular, :, ρrange[i])[:] = Qcol.angular
                view(Q.linear, :, ρrange[i])[:] = Qcol.linear
            end
        else
            view(Q.angular, :, ρrange) .= 0
            view(Q.linear, :, ρrange) .= 0
        end
    end
end

function back_out_external_wrenches!(controller::MomentumBasedController)
    Q = controller.wrenchmatrix
    ρ = controller.ρ
    externalwrenches = controller.externalwrenches
    T = eltype(controller)
    for body in keys(externalwrenches)
        externalwrenches[body] = zero(Wrench{T}, Q.frame)
    end
    for settings in controller.contactsettings
        body = settings.body
        ρrange = settings.ρrange
        for col in ρrange
            angular = SVector{3}(view(Q.angular, :, col)) * ρ[col]
            linear = SVector{3}(view(Q.linear, :, col)) * ρ[col]
            externalwrenches[body] += Wrench(Q.frame, angular, linear)
        end
    end
end

function control(controller::MomentumBasedController, t, controllerstate::MomentumBasedControllerState)
    t < controllerstate.next_control_time - controller.Δt - eps(typeof(t)) && reset!(controllerstate)
    if t >= controllerstate.next_control_time
        controllerstate.next_control_time = t + controller.Δt
        state = controllerstate.mechstate

        T = eltype(controller)
        mechanism = controller.mechanism
        nv = num_velocities(state)
        nρ = length(controller.ρ)

        # Create model
        solver = Gurobi.GurobiSolver()
        Gurobi.setparameters!(solver, Silent = true)
        model = Model(solver = solver)
        @variable(model, v̇[1 : nv])
        @variable(model, ρ[1 : nρ] >= 0)

        # Compute centroidal frame transform
        com = center_of_mass(state)
        centroidal_to_world = Transform3D(controller.centroidalframe, com.frame, com.v)
        world_to_centroidal = inv(centroidal_to_world)

        # Momentum-related quantities
        A = controller.momentummatrix
        momentum_matrix!(A, state, world_to_centroidal)
        Ȧv = transform(momentum_rate_bias(state), world_to_centroidal)

        # Gravitational wrench
        m = controller.mass
        fg = world_to_centroidal * (m * mechanism.gravitational_acceleration)
        Wg = Wrench(zero(fg), fg)

        # Wrench matrix
        update_wrench_matrix!(controller, state, world_to_centroidal)
        Q = controller.wrenchmatrix

        # Newton; TODO: add external wrenches to compensate for
        @framecheck A.frame Ȧv.frame
        @framecheck A.frame Wg.frame
        @framecheck A.frame Q.frame
        @constraint(model, A.angular * v̇ + Ȧv.angular .== Wg.angular + Q.angular * ρ)
        @constraint(model, A.linear * v̇ + Ȧv.linear .== Wg.linear + Q.linear * ρ)

        # Initialize objective
        obj = zero(JuMP.GenericQuadExpr{T,JuMP.Variable})

        # TODO: warm-start-enabled version of the following
        # Handle desired spatial accelerations
        for task in controller.spatialacceltasks
            if isenabled(task)
                J = task.jacobian
                desired = task.desired
                path = task.path
                world_to_desired = inv(transform_to_root(state, desired.frame))
                geometric_jacobian!(J, state, path, world_to_desired)
                J̇v = transform(state, -bias_acceleration(state, source(path)) + bias_acceleration(state, target(path)), desired.frame)
                @framecheck J.frame J̇v.frame
                @framecheck J.frame desired.frame
                error = -(SpatialAcceleration(J, v̇) + J̇v) + task.desired
                angularerror = task.angularselectionmatrix * error.angular
                linearerror = task.linearselectionmatrix * error.linear
                if isconstraint(task)
                    @constraint(model, angularerror .== 0)
                    @constraint(model, linearerror .== 0)
                else
                    obj += task.weight * (dot(angularerror, angularerror) + dot(linearerror, linearerror))
                end
            end
        end

        # Handle desired joint accelerations
        for task in controller.jointacceltasks
            if isenabled(task)
                joint = task.joint
                vrange = velocity_range(state, joint)
                error = task.desired - view(v̇, vrange)
                if isconstraint(task)
                    @constraint(model, error .== 0)
                else
                    obj += task.weight * dot(error, error)
                end
            end
        end

        # Handle desired momentum rate.
        momentumratetask = controller.momentumratetask
        if isenabled(momentumratetask)
            error = momentumratetask.desired - (Wrench(A, v̇) + Ȧv)
            angularerror = Array(momentumratetask.angularselectionmatrix * error.angular)
            linearerror = Array(momentumratetask.linearselectionmatrix * error.linear)
            if isconstraint(momentumratetask)
                @constraint(model, angularerror .== 0)
                @constraint(model, linearerror .== 0)
            else
                obj += momentumratetask.weight * (dot(angularerror, angularerror) + dot(linearerror, linearerror))
            end
        end

        # Handle contact weights
        for settings in controller.contactsettings
            ρcontact = view(ρ, settings.ρrange)
            obj += settings.weight * dot(ρcontact, ρcontact)
        end

        @objective(model, Min, obj)

        # Solve
        status = solve(model)

        # Inverse dynamics
        controller.ρ .= getvalue(ρ)
        result = controller.result
        result.v̇ .= getvalue(v̇)
        externalwrenches = controller.externalwrenches
        τ = controllerstate.τ
        back_out_external_wrenches!(controller)
        map!(wrench -> transform(wrench, centroidal_to_world), values(externalwrenches), values(externalwrenches))
        inverse_dynamics!(τ, result.jointwrenches, result.accelerations, state, result.v̇, externalwrenches)
    end

    controllerstate.τ
end
