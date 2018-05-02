mutable struct ContactSettings # TODO: store contact force basis
    body::RigidBody{Float64}
    point::Point3D{SVector{3, Float64}}
    weight::Float64
    μ::Float64
    normal::FreeVector3D{SVector{3, Float64}}
    maxnormalforce::Float64
    ρrange::UnitRange{Int64}

    function ContactSettings(body::RigidBody{Float64}, point::Point3D, ρrange::UnitRange{Int64})
        @framecheck point.frame default_frame(body)
        weight = 0.0
        μ = 0.0
        normal = FreeVector3D(point.frame, zeros(SVector{3, Float64}))
        maxnormalforce = 0.0
        new(body, point, weight, μ, normal, maxnormalforce, ρrange)
    end
end

num_basis_vectors(settings::ContactSettings) = length(settings.ρrange)
isenabled(settings::ContactSettings) = settings.maxnormalforce > 0
disable!(settings::ContactSettings) = settings.maxnormalforce = 0

function set!(settings::ContactSettings, weight::Float64, μ::Float64, normal::FreeVector3D, maxnormalforce::Float64 = Inf)
    settings.weight = weight
    settings.μ = μ
    settings.normal = normal
    settings.maxnormalforce = maxnormalforce
end
