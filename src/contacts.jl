struct ContactPoint
    position::Point3D{SVector{3, Float64}}
    normal::FreeVector3D{SVector{3, Float64}}
    μ::Float64
    localtransform::Transform3D{Float64}

    function ContactPoint(
            position::Point3D{SVector{3, Float64}},
            normal::FreeVector3D{SVector{3, Float64}},
            μ::Float64)
        @framecheck position.frame normal.frame
        point_origin_normal_z = CartesianFrame3D()
        z = SVector(0., 0., 1.)
        rot = Rotations.rotation_between(z, normal.v)
        # `localtransform`: transform from frame with:
        # * origin at `position`
        # * `normal` as z-axis
        # to frame in which `position` and `normal` are expressed:
        localtransform = Transform3D(point_origin_normal_z, position.frame, rot, position.v)
        new(position, normal, μ, localtransform)
    end
end

function wrenchbasis(
        point::ContactPoint,
        num_basis_vectors::Val{N},
        body_to_desired::Transform3D = eye(Transform3D{Float64}, point.localtransform.to)) where N
    Δθ = 2 * π / N
    μ = point.μ
    basis_vectors = ntuple(num_basis_vectors) do i
        θ = (i - 1) * Δθ
        normalize(SVector(μ * cos(θ), μ * sin(θ), 1.0))
    end
    linear = hcat(basis_vectors...)
    angular = zero(linear)
    wrenchmatrix = WrenchMatrix(point.localtransform.from, angular, linear)
    transform(wrenchmatrix, body_to_desired * point.localtransform)
end

mutable struct ContactConfiguration{N}
    point::ContactPoint
    ρ::SVector{N, Variable} # basis vector multipliers
    weight::Float64
    maxnormalforce::Float64

    function ContactConfiguration(point::ContactPoint, ρ::SVector{N, Variable}) where N
        new{N}(point, ρ, 0.0, 0.0)
    end
end

disable!(config::ContactConfiguration) = config.maxnormalforce = 0
isenabled(config::ContactConfiguration) = config.maxnormalforce > 0
num_basis_vectors(::ContactConfiguration{N}) where {N} = N

function wrenchbasis(config::ContactConfiguration{N}, body_to_desired::Transform3D) where N
    wrenchbasis(config.point, Val(N), body_to_desired)
end

function wrenchbasis(config::ContactConfiguration, state::MechanismState, qpmodel::SimpleQP.Model)
    Parameter(qpmodel) do # TODO: inference?
        to_root = transform_to_root(state, config.point.localtransform.to) # TODO: make nicer
        wrenchbasis(config, to_root)
    end
end
