__precompile__()

module QPControl

# Contact-related
export
    ContactPoint,
    disable!

# Task-related
export
    AbstractMotionTask,
    SpatialAccelerationTask,
    AngularAccelerationTask,
    LinearAccelerationTask,
    PointAccelerationTask,
    JointAccelerationTask,
    MomentumRateTask,
    LinearMomentumRateTask,
    setdesired!

# Low level
export
    MomentumBasedController,
    addtask!,
    addcontact!,
    regularize!,
    centroidal_frame

# High level
export
    StandingController

# Trajectories
export
    InterpolationTrajectory

using LinearAlgebra
using RigidBodyDynamics
using RigidBodyDynamics.Graphs
using RigidBodyDynamics.Contact
using RigidBodyDynamics.PDControl
using Parametron
using StaticArrays
using Rotations
using StaticUnivariatePolynomials

import MathOptInterface

const MOI = MathOptInterface
const RBD = RigidBodyDynamics
const SUP = StaticUnivariatePolynomials

include("contacts.jl")
include("tasks.jl")
include("exceptions.jl")
include("lowlevel/momentum.jl")
include("highlevel/standing.jl")
include("trajectories.jl")

end # module
