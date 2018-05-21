module MomentumBasedControlTests

using MomentumBasedControl
using RigidBodyDynamics
using RigidBodyDynamics.OdeIntegrators
using ValkyrieRobot.BipedControlUtil
using ValkyrieRobot
using StaticArrays
using Rotations
using Base.Test
using OSQP.MathOptInterfaceOSQP
import SimpleQP

const MBC = MomentumBasedControl
const RBD = RigidBodyDynamics

include("contacts.jl")
include("tasks.jl")
# include("controller.jl")
# include("notebooks.jl")

end # module
