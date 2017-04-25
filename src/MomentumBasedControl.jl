__precompile__()

module MomentumBasedControl

using RigidBodyDynamics
using RigidBodyDynamics.Graphs
using RigidBodyDynamics.Contact
using JuMP
using StaticArrays
using Rotations
import Gurobi

import RigidBodyDynamics: set! # TODO

export MomentumBasedController,
    MomentumBasedControllerState,
    ContactSettings,
    SpatialAccelerationTask,
    JointAccelerationTask,
    MomentumRateTask,
    PDGains,
    centroidal_transform,
    centroidal_frame,
    add_contact!,
    add_contacts!,
    add_mechanism_contacts!,
    add!,
    add_mechanism_joint_accel_tasks!,
    clear_contacts!,
    clear_desireds!,
    reset!,
    regularize_joint_accels!,
    control,
    set!,
    disable!,
    num_basis_vectors

include("contact_settings.jl")
include("motion_tasks.jl")
include("controller.jl")
include("pd.jl")

end # module
