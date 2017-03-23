module MomentumBasedControl

using RigidBodyDynamics
using RigidBodyDynamics.Graphs
using RigidBodyDynamics.Contact
using JuMP
using StaticArrays
import DataStructures: OrderedDict
import Gurobi

export MomentumBasedController,
    control,
    clear_contacts!,
    clear_desireds!,
    set_desired_accel!,
    set_desired_momentum_rate!,
    set_contact_weight,
    set_joint_accel_weights


include("util.jl")
include("controller.jl")
include("pd.jl")

end # module
