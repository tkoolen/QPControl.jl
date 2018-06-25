"""
    z_up_transform(origin::Point3D, zaxis::FreeVector3D, from::CartesianFrame3D)

Return `Transform3D` from `from` to the frame in which both `origin` and `zaxis`
is expressed, such that z-axis of `from` is `zaxis`, and origin of `from` is
`origin`.
"""
function z_up_transform(origin::Point3D, zaxis::FreeVector3D, from::CartesianFrame3D)
    @framecheck origin.frame zaxis.frame
    to = origin.frame
    rotation = Rotations.rotation_between(SVector(0., 0., 1.), zaxis.v)
    translation = origin.v
    Transform3D(from, to, rotation, translation)
end

function forcebasis(μ::Float64, num_basis_vectors::Val{N}) where N
    Δθ = 2 * π / N
    basisvectors = ntuple(num_basis_vectors) do i
        θ = (i - 1) * Δθ
        normalize(SVector(μ * cos(θ), μ * sin(θ), 1.0))
    end
    hcat(basisvectors...)
end

nvars(model::SimpleQP.Model, ::Val{N}) where {N} = SVector(ntuple(_ -> Variable(model), Val(N)))

mutable struct ContactPoint{N}
    normal_aligned_frame::CartesianFrame3D
    ρ::SVector{N, Variable} # basis vector multipliers
    force_local::FreeVector3D{SVector{3, Variable}} # contact force expressed in contact point's normal-aligned frame
    wrench_world::Wrench{Variable}
    position::Point3D{SVector{3, Float64}}
    normal::FreeVector3D{SVector{3, Float64}}
    μ::Float64
    weight::Float64
    maxnormalforce::Float64

    function ContactPoint{N}(
            position::Point3D, normal::FreeVector3D, μ::Float64,
            state::MechanismState, model::SimpleQP.Model) where N
        # frames
        normal_aligned_frame = CartesianFrame3D() # assumed not to change
        worldframe = root_frame(state.mechanism)

        # variables
        ρ = nvars(model, Val(N))
        force_local = FreeVector3D(normal_aligned_frame, nvars(model, Val(3)))
        wrench_world = Wrench(worldframe, nvars(model, Val(3)), nvars(model, Val(3)))

        ret = new{N}(normal_aligned_frame, ρ, force_local, wrench_world, position, normal, μ, 0.0, 0.0)

        # constraints
        basis = Parameter(() -> forcebasis(ret.μ, Val(N)), model)
        maxρ = Parameter(x -> x .= ret.maxnormalforce / (N * sqrt(ret.μ^2 + 1)), zeros(N), model)
        toroot = Parameter{Transform3D{Float64}}(model) do
            localtransform = z_up_transform(ret.position, ret.normal, ret.normal_aligned_frame)
            transform_to_root(state, localtransform.to) * localtransform
        end
        @constraint(model, force_local.v == basis * ρ)
        @constraint(model, ρ >= zeros(N))
        @constraint(model, ρ <= maxρ)
        @constraint(model, linear(wrench_world) == rotation(toroot) * force_local.v)
        @constraint(model, angular(wrench_world) == translation(toroot) × linear(wrench_world))

        ret
    end
end

disable!(point::ContactPoint) = point.maxnormalforce = 0
isenabled(point::ContactPoint) = point.maxnormalforce > 0

function objectiveterm(point::ContactPoint, model::SimpleQP.Model)
    weight = Parameter{Float64}(() -> point.weight, model)
    f = point.force_local
    @expression weight * (f ⋅ f)
end

# TODO: type piracy:
function SimpleQP.value(model::SimpleQP.Model, wrench::Wrench{Variable})
    # TODO: Ref probably allocates on 0.6
    Wrench(wrench.frame, value.(Ref(model), angular(wrench)), value.(Ref(model), linear(wrench)))
end
