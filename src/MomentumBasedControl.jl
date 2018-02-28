__precompile__()

module MomentumBasedControl

using RigidBodyDynamics
using RigidBodyDynamics.Graphs
using RigidBodyDynamics.Contact
using RigidBodyDynamics.PDControl
using JuMP
using StaticArrays
using Rotations
using Compat
using ForwardDiff
import MathProgBase: setparameters!
import OSQP: OSQPMathProgBaseInterface

import RigidBodyDynamics: set! # TODO

export MomentumBasedController,
    ContactSettings,
    SpatialAccelerationTask,
    JointAccelerationTask,
    MomentumRateTask,
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
    pd_center_of_mass!,
    set!,
    disable!,
    num_basis_vectors,
    val_deriv_deriv2

include("contact_settings.jl")
include("motion_tasks.jl")
include("controller.jl")
include("util.jl")

end # module
