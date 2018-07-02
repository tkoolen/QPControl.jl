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
    AngularAccelerationTask,
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
using Nullables

import MathOptInterface

const MOI = MathOptInterface
const RBD = RigidBodyDynamics

include("contacts.jl")
include("tasks.jl")
include("exceptions.jl")
include("controller.jl")

end # module
