# TODO: max normal force

function contactsettings end
function motiontasks end

struct MomentumBasedController{O<:MOI.AbstractOptimizer, N, M}
    # dynamics-related
    mechanism::Mechanism{Float64} # TODO: possibly remove
    centroidalframe::CartesianFrame3D
    totalmass::Float64
    result::DynamicsResult{Float64, Float64}
    momentummatrix::MomentumMatrix{Matrix{Float64}}
    wrenchmatrix::WrenchMatrix{Matrix{Float64}}

    # contact info and motion tasks
    externalwrenches::Dict{RigidBody{Float64}, Wrench{Float64}}
    contactsettings::Vector{ContactSettings{N}}
    motiontasks::M

    # qp-related
    qpmodel::SimpleQP.Model{O}

    function MomentumBasedController(
            mechanism::Mechanism{Float64}, optimizer::O, contactsettings::Vector{ContactSettings{N}}, motiontasks::M) where {O<:MOI.AbstractOptimizer, N, M}
        nv = num_velocities(mechanism)

        # dynamics-related
        centroidalframe = CartesianFrame3D("centroidal")
        totalmass = mass(mechanism)
        result = DynamicsResult(mechanism)
        momentummatrix = MomentumMatrix(centroidalframe, zeros(3, nv), zeros(3, nv))
        wrenchmatrix = WrenchMatrix(centroidalframe, zeros(3, 0), zeros(3, 0))

        # contact info and motion tasks
        rootframe = root_frame(mechanism)
        externalwrenches = Dict(b => zero(Wrench{Float64}, rootframe) for b in bodies(mechanism))

        # qp-related
        qpmodel = SimpleQP.Model(optimizer)
        ρ = Vector{Float64}()

        new{O, N, M}(
            mechanism, centroidalframe, totalmass, result, momentummatrix, wrenchmatrix,
            externalwrenches, contactsettings, motiontasks,
            qpmodel)
    end
end

function MomentumBasedController(mechanism::Mechanism, optimizer::MOI.AbstractOptimizer, highlevelcontroller)
    MomentumBasedController(mechanism, optimizer, contactsettings(highlevelcontroller), motiontasks(highlevelcontroller))
end

centroidal_frame(controller::MomentumBasedController) = controller.centroidalframe

