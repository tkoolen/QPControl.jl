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
    addcontact!
# Low level: Momentum
export
    MomentumBasedController,
    addtask!,
    regularize!,
    centroidal_frame
# Low level: MPC
export
    MPCController,
    addstage!,
    addstages!

# High level
export
    StandingController

using Compat
using Compat.LinearAlgebra
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
include("lowlevel/momentum.jl")
include("lowlevel/mpc.jl")
include("highlevel/standing.jl")

end # module
