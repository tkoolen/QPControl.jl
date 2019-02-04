struct Piecewise{T, K, TV<:AbstractVector{T}, KV<:AbstractVector{K}}
    subtrajectories::TV
    knots::KV
    clamp::Bool

    function Piecewise(subtrajectories::TV, knots::KV; clamp::Bool=true) where {T, K, TV<:AbstractVector{T}, KV<:AbstractVector{K}}
        @assert issorted(knots)
        @assert length(knots) == length(subtrajectories) + 1
        new{T, K, TV, KV}(subtrajectories, knots, clamp)
    end
end

function (traj::Piecewise)(x, args...)
    knots = traj.knots
    subtrajectories = traj.subtrajectories
    x0, xf = first(knots), last(knots)
    x′ = if traj.clamp
        clamp(x, x0, xf)
    else
        (x < x0 || x > xf) && throw_trajectory_domain_error(x, x0, xf)
        x
    end
    index = min(searchsortedlast(knots, x′), length(subtrajectories))
    @inbounds xknot = knots[index]
    @inbounds subtraj = traj.subtrajectories[index]
    return subtraj(x′ - xknot, args...)
end
