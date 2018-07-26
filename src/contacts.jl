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

struct ContactPoint{N}
    normal_aligned_frame::CartesianFrame3D
    ρ::SVector{N, Variable} # basis vector multipliers
    force_local::FreeVector3D{SVector{3, Variable}} # contact force expressed in contact point's normal-aligned frame
    wrench_world::Wrench{Variable}
    position::Any # not used in any methods, just here for debugging
    normal::Any # not used in any methods, just here for debugging
    μ::Any  # not used in any methods, just here for debugging
    weight::typeof(Ref(1.0)) # In 0.7 and up, this is just Ref{Float64}
    maxnormalforce::typeof(Ref(1.0)) # In 0.7 and up, this is just Ref{Float64}

    function ContactPoint{N}(
            position::Union{Point3D, Parameter{<:Point3D}}, normal::Union{FreeVector3D, Parameter{<:FreeVector3D}}, μ::Union{Float64, Parameter{Float64}},
            state::MechanismState, model::SimpleQP.Model) where N
        # frames
        normal_aligned_frame = CartesianFrame3D()
        worldframe = root_frame(state.mechanism)

        # variables
        ρ = nvars(model, Val(N))
        force_local = FreeVector3D(normal_aligned_frame, nvars(model, Val(3)))
        wrench_world = Wrench(worldframe, nvars(model, Val(3)), nvars(model, Val(3)))

        ret = new{N}(normal_aligned_frame, ρ, force_local, wrench_world, position, normal, μ, Ref(0.0), Ref(0.0))

        # constraints
        basis = @expression(forcebasis(μ, Val(N)))
        maxρ = let ret = ret
            Parameter(x -> x .= ret.maxnormalforce[] / (N * sqrt(μ^2 + 1)), zeros(N), model)
        end
        toroot = @expression(transform_to_root(state, position.frame) * z_up_transform(position, normal, normal_aligned_frame))
        # toroot = let state = state, ret = ret
        hat = RBD.Spatial.hat
        @constraint(model, force_local.v == basis * ρ)
        @constraint(model, ρ >= zeros(N))
        @constraint(model, ρ <= maxρ)
        @constraint(model, linear(wrench_world) == rotation(toroot) * force_local.v)
        @constraint(model, angular(wrench_world) == hat(translation(toroot)) * linear(wrench_world)) # TODO: ×

        ret
    end
end

disable!(point::ContactPoint) = point.maxnormalforce[] = 0
isenabled(point::ContactPoint) = point.maxnormalforce[] > 0

function objectiveterm(point::ContactPoint, model::SimpleQP.Model)
    weight = Parameter{Float64}(() -> point.weight[], model)
    f = point.force_local
    @expression weight * (f ⋅ f)
end

# TODO: type piracy:
function SimpleQP.value(m::SimpleQP.Model, wrench::Wrench{Variable})
    # Wrench(wrench.frame, value.(Ref(m), angular(wrench)), value.(Ref(m), linear(wrench))) # allocates on 0.6
    τ = angular(wrench)
    f = linear(wrench)
    @inbounds τval = SVector(value(m, τ[1]), value(m, τ[2]), value(m, τ[3]))
    @inbounds fval = SVector(value(m, f[1]), value(m, f[2]), value(m, f[3]))
    return Wrench(wrench.frame, τval, fval)
end
