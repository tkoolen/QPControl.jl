__precompile__()

module MomentumBasedControl

export
    ContactInfo,
    ContactSettings,
    MotionTask,
    SpatialAccelerationTask,
    AngularAccelerationTask,
    JointAccelerationTask,
    LinearMomentumRateTask,
    Weighted

export
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
# include("controller.jl")
# include("util.jl")

end # module
