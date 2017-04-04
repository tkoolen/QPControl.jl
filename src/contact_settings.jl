type ContactSettings{T}
    body::RigidBody{T}
    point::Point3D{SVector{3, T}}
    weight::T
    μ::T
    normal::FreeVector3D{SVector{3, T}}
    maxnormalforce::T
    ρrange::UnitRange{Int64}

    function ContactSettings(body::RigidBody{T}, point::Point3D, ρrange::UnitRange{Int64})
        @framecheck point.frame default_frame(body)
        weight = zero(T)
        μ = zero(T)
        normal = FreeVector3D(point.frame, zeros(SVector{3, T}))
        maxnormalforce = zero(T)
        new(body, point, weight, μ, normal, maxnormalforce, ρrange)
    end
end

ContactSettings{T}(body::RigidBody{T}, point::Point3D, ρrange::UnitRange{Int64}) = ContactSettings{T}(body, point, ρrange)
num_basis_vectors(settings::ContactSettings) = length(settings.ρrange)
isenabled(settings::ContactSettings) = settings.maxnormalforce > 0
disable!(settings::ContactSettings) = settings.maxnormalforce = 0

function set!{T}(settings::ContactSettings{T}, weight::T, μ::T, normal::FreeVector3D, maxnormalforce::T = Inf)
    settings.weight = weight
    settings.μ = μ
    settings.normal = normal
    settings.maxnormalforce = maxnormalforce
end
