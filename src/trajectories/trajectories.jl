module Trajectories

export
    InterpolationTrajectory

using StaticUnivariatePolynomials
using Rotations

const SUP = StaticUnivariatePolynomials

include("interpolation.jl")

end
