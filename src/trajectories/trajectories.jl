module Trajectories

export
    fit_cubic,
    fit_quintic,
    Constant,
    Interpolated,
    Piecewise,
    BezierCurve,
    PointTrajectory,
    FreeVectorTrajectory,
    SE3Trajectory,
    PolyDerivEvaluator

using StaticUnivariatePolynomials
using LinearAlgebra
using StaticArrays
using Rotations
using RigidBodyDynamics

using Base: tail
using StaticUnivariatePolynomials: constant, derivative, exponential_integral, _map

const SUP = StaticUnivariatePolynomials

const BezierCurve = BernsteinPolynomial

@noinline function throw_trajectory_domain_error(x, x0, xf)
    throw(DomainError(x, "Trajectory evaluated outside of range [$x0, $xf]"))
end

zero_tangent(x::Number) = zero(x)
zero_tangent(x::AbstractVector) = zero(x)
zero_tangent(x::Rotation{3, T}) where {T} = zero(SVector{3, T})

include("fit_polynomial.jl")
include("constant.jl")
include("poly_deriv_evaluator.jl")
include("interpolated.jl")
include("piecewise.jl")
include("point_vector.jl")
include("se3.jl")

end
