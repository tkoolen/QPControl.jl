__precompile__()

module MomentumBasedControl

using RigidBodyDynamics
using RigidBodyDynamics.Graphs
using RigidBodyDynamics.Contact
using JuMP
using StaticArrays
using Rotations
import DataStructures: OrderedDict
import Gurobi

export MomentumBasedController,
    PDGains,
    control,
    centroidal_transform,
    centroidal_frame,
    clear_contacts!,
    clear_desireds!,
    enable_contact!,
    disable_contact!,
    reset!,
    set_desired_accel!,
    set_desired_momentum_rate!,
    set_contact_weight,
    set_joint_accel_weights,
    set_contact_regularization!,
    set_joint_accel_regularization!


include("util.jl")
include("controller.jl")
include("pd.jl")

end # module
