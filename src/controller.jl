# TODO: max normal force

mutable struct MomentumBasedController
    # dynamics-related
    mechanism::Mechanism{Float64} # TODO: possibly remove
    centroidalframe::CartesianFrame3D
    mass::Float64
    result::DynamicsResult{Float64, Float64}
    momentummatrix::MomentumMatrix{Matrix{Float64}}
    wrenchmatrix::WrenchMatrix{Matrix{Float64}}

    # contact info and motion tasks
    externalwrenches::BodyDict{Wrench{Float64}}
    contactsettings::Vector{ContactSettings}
    spatialacceltasks::Vector{SpatialAccelerationTask}
    jointacceltasks::Vector{JointAccelerationTask}
    momentumratetask::MomentumRateTask

    # qp-related
    ρ::Vector{Float64}

    function MomentumBasedController(mechanism::Mechanism{Float64})
        nv = num_velocities(mechanism)

        # dynamics-related
        centroidalframe = CartesianFrame3D("centroidal")
        totalmass = mass(mechanism)
        result = DynamicsResult(mechanism)
        momentummatrix = MomentumMatrix(centroidalframe, zeros(3, nv), zeros(3, nv))
        wrenchmatrix = WrenchMatrix(centroidalframe, zeros(3, 0), zeros(3, 0))

        # contact info and motion tasks
        rootframe = root_frame(mechanism)
        externalwrenches = BodyDict{Wrench{Float64}}(BodyID(b) => zero(Wrench{Float64}, rootframe) for b in bodies(mechanism))
        contactsettings = Vector{ContactSettings}()
        spatialacceltasks = Vector{SpatialAccelerationTask}()
        jointacceltasks = Vector{JointAccelerationTask}()
        momentumratetask = MomentumRateTask(centroidalframe)

        # qp-related
        ρ = Vector{Float64}()

        new(
            mechanism, centroidalframe, totalmass, result, momentummatrix, wrenchmatrix,
            externalwrenches, contactsettings, spatialacceltasks, jointacceltasks, momentumratetask,
            ρ)
    end
end

centroidal_frame(controller::MomentumBasedController) = controller.centroidalframe

function add_contact!(controller::MomentumBasedController, body::RigidBody, point::Point3D, num_basis_vectors::Int64)
    nρ = length(controller.ρ)
    nρnew = nρ + num_basis_vectors
    ρrange = nρ + 1 : nρnew
    settings = ContactSettings(body, point, ρrange)
    push!(controller.contactsettings, settings)
    resize!(controller.ρ, nρnew)
    controller.wrenchmatrix = WrenchMatrix(controller.centroidalframe, zeros(3, nρnew), zeros(3, nρnew))
    settings
end

function add_contacts!(controller::MomentumBasedController, body::RigidBody, num_basis_vectors::Int64)
    ret = Vector{ContactSettings}()
    for point in contact_points(body)
        push!(ret, add_contact!(controller, body, RigidBodyDynamics.Contact.location(point), num_basis_vectors))
    end
    ret
end

function add_mechanism_contacts!(controller::MomentumBasedController, num_basis_vectors::Int64 = 4)
    ret = Dict{RigidBody{Float64}, Vector{ContactSettings}}() # TODO: BodyIDDict
    for body in bodies(controller.mechanism)
        ret[body] = add_contacts!(controller, body, num_basis_vectors)
    end
    ret
end

add!(controller::MomentumBasedController, task::SpatialAccelerationTask) = push!(controller.spatialacceltasks, task)
add!(controller::MomentumBasedController, task::JointAccelerationTask) = push!(controller.jointacceltasks, task)
add!(controller::MomentumBasedController, task::MomentumRateTask) = controller.momentumratetask = task

function add_mechanism_joint_accel_tasks!(controller::MomentumBasedController)
    ret = Dict{Joint{Float64}, JointAccelerationTask}() # TODO: JointIDDict
    for joint in tree_joints(controller.mechanism)
        ret[joint] = task = JointAccelerationTask(joint)
        add!(controller, task)
    end
    ret
end

clear_contacts!(controller::MomentumBasedController) = foreach(disable!, controller.contactsettings)

function clear_desireds!(controller::MomentumBasedController)
    foreach(disable!, controller.spatialacceltasks)
    foreach(disable!, controller.jointacceltasks)
    controller.momentumratetask = MomentumRateTask(controller.centroidalframe)
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
                β = FreeVector3D(n.frame, rot * normalize(SVector(μ * cos(θ), μ * sin(θ), 1.0)))
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
    for body in keys(externalwrenches)
        externalwrenches[body] = zero(Wrench{Float64}, Q.frame)
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

function (controller::MomentumBasedController)(τ::AbstractVector, t::Number, state::MechanismState)
    mechanism = controller.mechanism
    nv = num_velocities(state)
    nρ = length(controller.ρ)

    # Create model
    solver = OSQPMathProgBaseInterface.OSQPSolver(eps_abs = 1e-8, eps_rel = 1e-16, max_iter = 10000)
    setparameters!(solver, Silent = true)
    model = Model(solver = solver)
    @variable(model, v̇[1 : nv])
    @variable(model, ρ[1 : nρ])
    ρcon = @constraint(model, ρ .>= 0)

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
    obj = zero(JuMP.GenericQuadExpr{Float64,JuMP.Variable})

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
    back_out_external_wrenches!(controller)
    map!(wrench -> transform(wrench, centroidal_to_world), values(externalwrenches), values(externalwrenches))
    inverse_dynamics!(τ, result.jointwrenches, result.accelerations, state, result.v̇, externalwrenches)

    τ
end
