struct ContactInfo
    point::Point3D{SVector{3, Float64}}
    normal::FreeVector3D{SVector{3, Float64}}
    μ::Float64
    localtransform::Transform3D{Float64} # from frame with origin at point and normal as z-axis to frame in which point and normal are expressed

    function ContactInfo(point::Point3D{SVector{3, Float64}}, normal::FreeVector3D{SVector{3, Float64}}, μ::Float64)
        @framecheck point.frame normal.frame
        point_origin_normal_z = CartesianFrame3D()
        z = SVector(0., 0., 1.)
        rot = Rotations.rotation_between(z, normal.v)
        localtransform = Transform3D(point_origin_normal_z, point.frame, rot, point.v)
        new(point, normal, μ, localtransform)
    end
end

function wrenchbasis(info::ContactInfo, num_basis_vectors::Val{N},
        body_to_desired::Transform3D = eye(Transform3D{Float64}, info.localtransform.to)) where N
    Δθ = 2 * π / N
    μ = info.μ
    basis_vectors = ntuple(num_basis_vectors) do i
        θ = (i - 1) * Δθ
        normalize(SVector(μ * cos(θ), μ * sin(θ), 1.0))
    end
    linear = hcat(basis_vectors...)
    angular = zero(linear)
    wrenchmatrix = WrenchMatrix(info.localtransform.from, angular, linear)
    transform(wrenchmatrix, body_to_desired * info.localtransform)
end

mutable struct ContactSettings{N}
    contactinfo::ContactInfo
    weight::Float64
    maxnormalforce::Float64

    function ContactSettings{N}(contactinfo::ContactInfo) where N
        new{N}(contactinfo, 0.0, 0.0)
    end
end

disable!(settings::ContactSettings) = settings.maxnormalforce = 0
isenabled(settings::ContactSettings) = settings.maxnormalforce > 0
num_basis_vectors(::ContactSettings{N}) where {N} = N

function wrenchbasis(settings::ContactSettings{N}, body_to_desired::Transform3D) where N
    wrenchbasis(settings.contactinfo, Val(N), body_to_desired)
end
