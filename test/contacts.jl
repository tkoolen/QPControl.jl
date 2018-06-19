@testset "ContactPoint" begin
    frame = CartesianFrame3D()
    position = Point3D(frame, 1., 2., 3.)
    normal = FreeVector3D(frame, normalize(SVector(1., 3., 2.)))
    μ = 0.5
    point = ContactPoint(position, normal, μ)

    for N = 1 : 2 : 10
        basis = MBC.wrenchbasis(point, Val(N))
        for i = 1 : N
            ρ = zeros(N)
            ρ[i] = 1
            wrench = Wrench(basis, ρ)
            @test wrench.frame == frame
            force = FreeVector3D(frame, linear(wrench))
            torque = FreeVector3D(frame, angular(wrench))
            @test torque ≈ position × force atol = 1e-12
            normalforce = force ⋅ normal
            @test normalforce > 0
            @test norm(force - normalforce * normal) <= μ * normalforce + 1e-12
        end
    end
    let point = point
        @test_noalloc MBC.wrenchbasis(point, Val(4))
    end
end

@testset "ContactConfiguration" begin
    frame = CartesianFrame3D()
    position = Point3D(frame, 1., 2., 3.)
    normal = FreeVector3D(frame, normalize(SVector(1., 3., 2.)))
    μ = 0.5
    point = ContactPoint(position, normal, μ)
    ρ = SVector(ntuple(i -> Variable(i), Val(4)))
    settings = ContactConfiguration(point, ρ)
    @test !MBC.isenabled(settings)
    settings.maxnormalforce = 1.0
    @test MBC.isenabled(settings)
    basis1 = MBC.wrenchbasis(settings, eye(Transform3D{Float64}, frame))
    basis2 = MBC.wrenchbasis(point, Val(4))
    @test angular(basis1) == angular(basis2)
    @test linear(basis1) == linear(basis2)
end

@testset "Wrench expression allocations" begin
    model = MockModel()

    frame = CartesianFrame3D()
    position = Point3D(frame, 1., 2., 3.)
    normal = FreeVector3D(frame, normalize(SVector(1., 3., 2.)))
    μ = 0.5
    point = ContactPoint(position, normal, μ)
    ρ = SVector(ntuple(i -> Variable(i), Val(4)))
    settings = ContactConfiguration(point, ρ)

    tf = eye(Transform3D{Float64}, frame)
    wrenchbasis = Parameter(@closure(() -> MBC.wrenchbasis(settings, tf)), model)

    torque = @expression angular(wrenchbasis) * ρ
    force = @expression linear(wrenchbasis) * ρ

    @test_noalloc (setdirty!(model); torque())
    @test_noalloc (setdirty!(model); force())
end
