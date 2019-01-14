module TrajectoriesTest

using QPControl.Trajectories
using Test
using LinearAlgebra
using Random
using Rotations
using RigidBodyDynamics
using StaticArrays

@testset "InterpolationTrajectory" begin
    let y0 = 1, yf = 2
        traj = InterpolationTrajectory(0.0, 1.0, y0, yf)

        y = traj(0.5)
        @test y === 1.5

        y, yd, ydd = traj(0.5, Val(2))
        @test y === 1.5
        @test yd === 1.0
        @test ydd === 0.0
    end

    let f = CartesianFrame3D(), y0 = Point3D(f, 1, 2, 3), yf = Point3D(f, 0, 2, 4)
        traj = InterpolationTrajectory(0.0, 1.0, y0, yf)
        y, yd, ydd = traj(0.5, Val(2))
        @test y ≈ Point3D(f, 0.5, 2.0, 3.5) atol=1e-15
        @test yd ≈ FreeVector3D(f, -1.0, 0.0, 1.0) atol=1e-15
        @test ydd ≈ FreeVector3D(f, 0.0, 0.0, 0.0) atol=1e-15
    end

    let angle = π / 2, axis = SVector(1.0, 0.0, 0.0), y0 = one(Quat), yf = Quat(AngleAxis(angle, axis...))
        traj = InterpolationTrajectory(0.0, 1.0, y0, yf)

        y = traj(0.0)
        @test y ≈ y0 atol=1e-15

        y = traj(1.0)
        @test y ≈ yf atol=1e-15

        y, yd, ydd = traj(0.5, Val(2))
        @test y ≈ Quat(AngleAxis(0.5 * angle, axis...)) atol=1e-15
        @test normalize(yd) ≈ axis atol=1e-15
        @test norm(yd) ≈ abs(angle) atol=1e-15
        @test ydd ≈ SVector(0.0, 0.0, 0.0) atol=1e-15
    end
end

end
