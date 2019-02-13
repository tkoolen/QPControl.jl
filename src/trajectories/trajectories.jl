module Trajectories

export
    fit_cubic,
    fit_quintic,
    Constant,
    Interpolated,
    Piecewise,
    BezierCurve

using StaticUnivariatePolynomials
using StaticArrays
using Rotations

using Base: tail
using StaticUnivariatePolynomials: constant, derivative

const SUP = StaticUnivariatePolynomials

@noinline function throw_trajectory_domain_error(x, x0, xf)
    throw(DomainError(x, "Trajectory evaluated outside of range [$x0, $xf]"))
end

tangent_type(::Type{T}) where {T<:Number} = T

include("fit_polynomial.jl")
include("constant.jl")
include("interpolated.jl")
include("piecewise.jl")
include("bezier.jl")

end
