__precompile__()

module MomentumBasedControl

# Contact types
export
    ContactInfo,
    ContactSettings

# Task types
export
    AbstractMotionTask,
    SpatialAccelerationTask,
    JointAccelerationTask,
    CentroidalMomentumRateTask
    # WeightedMotionTask,

# Control types
export
    MomentumBasedController

# Utility functions
export
    contactsettings,
    motiontasks,
    disable!

using Compat
using RigidBodyDynamics
using RigidBodyDynamics.Graphs
using RigidBodyDynamics.Contact
using RigidBodyDynamics.PDControl
using SimpleQP
using StaticArrays
using Rotations
import MathOptInterface

const MOI = MathOptInterface

include("contacts.jl")
include("tasks.jl")
include("controller.jl")
# include("util.jl")

end # module
