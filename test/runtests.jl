module MomentumBasedControlTests

using Compat
using Compat.Test
using Compat.Random
using Compat.LinearAlgebra
using QPControl
using RigidBodyDynamics
using ValkyrieRobot.BipedControlUtil
using ValkyrieRobot
using StaticArrays
using Rotations
using MathOptInterface
using OSQP.MathOptInterfaceOSQP
using SimpleQP

import SimpleQP: MockModel, setdirty!

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
    optimizer = OSQPOptimizer()
    MOI.set!(optimizer, OSQPSettings.Verbose(), false)
    MOI.set!(optimizer, OSQPSettings.EpsAbs(), 1e-8)
    MOI.set!(optimizer, OSQPSettings.EpsRel(), 1e-16)
    MOI.set!(optimizer, OSQPSettings.MaxIter(), 10000)
    MOI.set!(optimizer, OSQPSettings.AdaptiveRhoInterval(), 25) # required for deterministic behavior
    optimizer
end

include("tasks.jl")
include("controller.jl")
include("notebooks.jl")

end # module
