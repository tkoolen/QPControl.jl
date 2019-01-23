module Trajectories

export
    fit_cubic,
    fit_quintic,
    InterpolationTrajectory

using StaticUnivariatePolynomials
using StaticArrays
using Rotations

const SUP = StaticUnivariatePolynomials

include("fit_polynomial.jl")
include("interpolation.jl")

end
