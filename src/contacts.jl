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

struct ContactPoint
    position::Point3D{SVector{3, Float64}}
    normal::FreeVector3D{SVector{3, Float64}}
    μ::Float64
    localtransform::Transform3D{Float64}

    function ContactPoint(
            position::Point3D{SVector{3, Float64}},
            normal::FreeVector3D{SVector{3, Float64}},
            μ::Float64)
        normal_aligned_frame = CartesianFrame3D()
        localtransform = z_up_transform(position, normal, normal_aligned_frame)
        new(position, normal, μ, localtransform)
    end
end

nvars(model::SimpleQP.Model, ::Val{N}) where {N} = SVector(ntuple(_ -> Variable(model), Val(N)))

mutable struct ContactQPData{N}
    ρ::SVector{N, Variable} # basis vector multipliers
    force_local::FreeVector3D{SVector{3, Variable}} # contact force expressed in contact point's normal-aligned frame
    wrench_world::Wrench{Variable}
    point::ContactPoint
    weight::Float64
    maxnormalforce::Float64

    function ContactQPData{N}(point::ContactPoint, state::MechanismState, model::SimpleQP.Model) where N
        # frames
        normal_aligned_frame = point.localtransform.from # assumed not to change
        worldframe = root_frame(state.mechanism)

        # variables
        ρ = nvars(model, Val(N))
        force_local = FreeVector3D(normal_aligned_frame, nvars(model, Val(3)))
        wrench_world = Wrench(worldframe, nvars(model, Val(3)), nvars(model, Val(3)))

        ret = new{N}(ρ, force_local, wrench_world, point, 0.0, 0.0)

        # constraints
        basis = Parameter(() -> forcebasis(ret.point.μ, Val(N)), model)
        maxnormalforce = Parameter(x -> x[1] = ret.maxnormalforce, zeros(1), model)
        toroot = Parameter{Transform3D{Float64}}(model) do
            localtransform = point.localtransform
            transform_to_root(state, localtransform.to) * localtransform
        end
        normalforcevec = [force_local.v[3]] # TODO: would be nicer with a scalar constraint for max normal force
        hat = RigidBodyDynamics.Spatial.hat

        @constraint(model, force_local.v == basis * ρ)
        @constraint(model, ρ >= zeros(N))
        # @constraint(model, normalforcevec <= maxnormalforce) # FIXME
        @constraint(model, linear(wrench_world) == rotation(toroot) * force_local.v)
        @constraint(model, angular(wrench_world) == hat(translation(toroot)) * linear(wrench_world))

        # TODO: add objective term here?

        ret
    end
end

disable!(data::ContactQPData) = data.maxnormalforce = 0
isenabled(data::ContactQPData) = data.maxnormalforce > 0

function objectiveterm(data::ContactQPData, model::SimpleQP.Model)
    weight = Parameter{Float64}(() -> data.weight, model)
    f = data.force_local
    @expression weight * (f ⋅ f)
end

# TODO: type piracy:
function SimpleQP.value(model::SimpleQP.Model, wrench::Wrench{Variable})
    # TODO: Ref probably allocates on 0.6
    Wrench(wrench.frame, value.(Ref(model), angular(wrench)), value.(Ref(model), linear(wrench)))
end
