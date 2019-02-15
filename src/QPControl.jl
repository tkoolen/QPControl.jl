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

# Low-level
export
    MomentumBasedController,
    SE3PDController,
    addtask!,
    addcontact!,
    regularize!,
    centroidal_frame

# High-level
export
    StandingController

using LinearAlgebra
using RigidBodyDynamics
using RigidBodyDynamics.Graphs
using RigidBodyDynamics.Contact
using RigidBodyDynamics.PDControl
using Parametron
using StaticArrays
using Rotations

import MathOptInterface

const MOI = MathOptInterface
const RBD = RigidBodyDynamics

include("contacts.jl")
include("tasks.jl")
include("exceptions.jl")
include(joinpath("lowlevel", "momentum.jl"))
include(joinpath("lowlevel", "se3pdcontroller.jl"))
include(joinpath("highlevel", "standing.jl"))
include(joinpath("trajectories", "trajectories.jl"))

end # module
