module TrajectoriesTest

using QPControl.Trajectories
using Test
using LinearAlgebra
using Random
using Rotations
using RigidBodyDynamics
using StaticArrays
using StaticUnivariatePolynomials
using ForwardDiff

const SUP = StaticUnivariatePolynomials

@testset "fit_quintic" begin
    rng = MersenneTwister()
    for i = 1 : 10
        x0 = rand(rng)
        xf = rand(rng)
        y0 = rand(rng)
        yd0 = rand(rng)
        ydd0 = rand(rng)
        yf = rand(rng)
        ydf = rand(rng)
        yddf = rand(rng)

        p = fit_quintic(x0=x0, xf=xf, y0=y0, yd0=yd0, ydd0=ydd0, yf=yf, ydf=ydf, yddf=yddf)
        pd = SUP.derivative(p)
        pdd = SUP.derivative(pd)

        atol = 1e-6
        @test p(x0) ≈ y0 atol=atol
        @test pd(x0) ≈ yd0 atol=atol
        @test pdd(x0) ≈ ydd0 atol=atol
        @test p(xf) ≈ yf atol=atol
        @test pd(xf) ≈ ydf atol=atol
        @test pdd(xf) ≈ yddf atol=atol
    end
end

@testset "InterpolationTrajectory, identity interpolator function" begin
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

@testset "InterpolationTrajectory, polynomial interpolator function" begin
    interpolator = fit_quintic(x0=0.0, xf=1.0, y0=0.0, yd0=0.0, ydd0=0.0, yf=1.0, ydf=0.0, yddf=0.0)
    x0 = -1.0
    xf = 2.0
    y0 = 2.0
    yf = 3.0
    traj = InterpolationTrajectory(x0, xf, y0, yf, interpolator, min_num_derivs=Val(2))

    y, yd, ydd = traj(x0, Val(2))
    @test y ≈ y0 atol=1e-10
    @test yd ≈ 0.0 atol=1e-10
    @test ydd ≈ 0.0 atol=1e-10

    y, yd, ydd = traj(xf, Val(2))
    @test y ≈ yf atol=1e-10
    @test yd ≈ 0.0 atol=1e-10
    @test ydd ≈ 0.0 atol=1e-10

    for x in range(x0; stop=xf, length=10)
        _, yd, ydd = traj(x, Val(2))
        yd_forward = ForwardDiff.derivative(traj, x)
        ydd_forward = ForwardDiff.derivative(x -> ForwardDiff.derivative(traj, x), x)
        @test yd ≈ yd_forward atol=1e-4
        @test ydd ≈ ydd_forward atol=1e-4
    end
end

end
