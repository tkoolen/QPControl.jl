struct Piecewise{T, B, VT<:AbstractVector{T}, VB<:AbstractVector{B}}
    subfunctions::VT
    breaks::VB
    clamp::Bool

    function Piecewise(subfunctions::VT, breaks::VB; clamp::Bool=true) where {T, B, VT<:AbstractVector{T}, VB<:AbstractVector{B}}
        @assert issorted(breaks)
        @assert length(breaks) == length(subfunctions) + 1
        new{T, B, VT, VB}(subfunctions, breaks, clamp)
    end
end

function (traj::Piecewise)(x, args...)
    breaks = traj.breaks
    subfunctions = traj.subfunctions
    x0, xf = first(breaks), last(breaks)
    x′ = if traj.clamp
        clamp(x, x0, xf)
    else
        (x < x0 || x > xf) && throw_trajectory_domain_error(x, x0, xf)
        x
    end
    index = min(searchsortedlast(breaks, x′), length(subfunctions))
    @inbounds xknot = breaks[index]
    @inbounds subtraj = traj.subfunctions[index]
    return subtraj(x′ - xknot, args...)
end

subfunctions(traj::Piecewise) = traj.subfunctions
breaks(traj::Piecewise) = traj.breaks
