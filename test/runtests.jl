module QPControlTests

using Test
using LinearAlgebra
using Random
using QPControl
using RigidBodyDynamics
using RigidBodyDynamics.Contact
using ValkyrieRobot.BipedControlUtil
using ValkyrieRobot
using StaticArrays
using Rotations
using MathOptInterface
using OSQP
using OSQP.MathOptInterfaceOSQP: OSQPSettings
using Parametron

import Parametron: mock_model, setdirty!

const QPC = QPControl
const RBD = RigidBodyDynamics
const MOI = MathOptInterface

macro test_noalloc(expr)
    quote
        let
            $expr
            allocs = @allocated $expr
            @test allocs == 0
        end
    end |> esc
end

function defaultoptimizer()
    optimizer = OSQP.Optimizer()
    MOI.set(optimizer, OSQPSettings.Verbose(), false)
    MOI.set(optimizer, OSQPSettings.EpsAbs(), 1e-8)
    MOI.set(optimizer, OSQPSettings.EpsRel(), 1e-16)
    MOI.set(optimizer, OSQPSettings.MaxIter(), 20000)
    MOI.set(optimizer, OSQPSettings.AdaptiveRhoInterval(), 25) # required for deterministic behavior
    optimizer
end

include("tasks.jl")
include("notebooks.jl")
include("controller.jl")

end # module
