# test_transforms.jl
using Transforms
using Test

@testset "LinearMap" begin
    m = LinearMap(0.0, 1000.0)

    @test forward(m, 0.0) ≈ -1.0
    @test forward(m, 1000.0) ≈ 1.0
    @test forward(m, 500.0) ≈ 0.0

    @test inverse(m, -1.0) ≈ 0.0
    @test inverse(m, 1.0) ≈ 1000.0
    @test inverse(m, 0.0) ≈ 500.0

    @test dzdξ(m, 0.0) ≈ 500.0
    @test d2zdξ2(m, 0.0) ≈ 0.0
end

@testset "TanhMap" begin
    m = TanhMap(0.0, 2000.0, 2.5)

    # Boundary conditions
    @test forward(m, 0.0) ≈ -1.0 atol=1e-10
    @test forward(m, 2000.0) ≈ 1.0 atol=1e-10

    # Round-trip
    z_test = 500.0
    ξ_test = forward(m, z_test)
    z_recovered = inverse(m, ξ_test)
    @test z_recovered ≈ z_test atol=1e-10

    # Metric positivity
    ξ_range = range(-0.9999, 0.9999, length=50)
    J1_vals = dzdξ.(Ref(m), ξ_range)
    @test all(J1_vals .> 0)  # Monotonic
end

@testset "LogarithmicMap" begin
    # Valid: zmin > 0
    m = LogarithmicMap(1.0, 1000.0)

    # Boundary conditions
    @test forward(m, 1.0) ≈ -1.0 atol=1e-10
    @test forward(m, 1000.0) ≈ 1.0 atol=1e-10

    # Invalid: zmin ≤ 0
    @test_throws DomainError LogarithmicMap(0.0, 1000.0)
    @test_throws DomainError LogarithmicMap(-10.0, 1000.0)
end

@testset "profile_transform" begin
    m = TanhMap(0.0, 1000.0, 2.0)
    z_profile = [10.0, 100.0, 500.0, 1000.0]
    ξ_profile = profile_transform(m, z_profile)

    @test ξ_profile[1] < ξ_profile[2] < ξ_profile[3] < ξ_profile[4]  # Monotonic
    @test ξ_profile[end] ≈ 1.0 atol=1e-10  # Upper boundary
end

@testset "JacobianStack" begin
    m1 = TanhMap(0.0, 1000.0, 2.0)
    m2 = HyperbolicMap(-1.0, 1.0, 1.0)
    stack = JacobianStack([m1, m2], ["TanhMap", "HyperbolicMap"])

    z = 500.0
    ξ, J1, J2 = evaluate_jacobian_stack(stack, z)

    @test -1.0 < ξ < 1.0  # Within computational domain
    @test J1 > 0  # Monotonic
end