__precompile__()

module MomentumBasedControl

# Contact-related
export
    ContactPoint,
    disable!

# Task-related
export
    AbstractMotionTask,
    SpatialAccelerationTask,
    JointAccelerationTask,
    MomentumRateTask,
    LinearMomentumRateTask,
    setdesired!

# Controller
export
    MomentumBasedController,
    addtask!,
    addcontact!,
    regularize!,
    centroidal_frame

using Compat
using RigidBodyDynamics
using RigidBodyDynamics.Graphs
using RigidBodyDynamics.Contact
using RigidBodyDynamics.PDControl
using SimpleQP
using StaticArrays
using Rotations
using FastClosures # TODO: necessary?

import MathOptInterface

const MOI = MathOptInterface

include("contacts.jl")
include("tasks.jl")
include("controller.jl")
# include("util.jl")

end # module
