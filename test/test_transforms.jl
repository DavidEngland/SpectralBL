# test/test_transforms.jl
include("../src/transforms.jl")

using .Transforms
using Test

@testset "LinearMap" begin
    m = LinearMap(0.0, 1000.0)
    @test forward(m, 0.0) ≈ -1.0
    @test forward(m, 1000.0) ≈ 1.0
    @test forward(m, 500.0) ≈ 0.0
    @test inverse(m, -1.0) ≈ 0.0
    @test inverse(m, 1.0) ≈ 1000.0
    @test dzdξ(m, 0.0) ≈ 500.0
    @test d2zdξ2(m, 0.0) ≈ 0.0
end

@testset "HyperbolicMap" begin
    m = HyperbolicMap(0.0, 1000.0, 3.0)
    @test forward(m, 0.0) ≈ -1.0 atol=1e-10
    @test forward(m, 1000.0) ≈ 1.0 atol=1e-10
    @test inverse(m, -1.0) ≈ 0.0 atol=1e-10
    @test dzdξ(m, 0.0) > 0.0
    @test_throws DomainError HyperbolicMap(-10.0, 1000.0, 3.0)
end

@testset "LogarithmicMap" begin
    m = LogarithmicMap(1.0, 1000.0)
    @test forward(m, 1.0) ≈ -1.0 atol=1e-10
    @test forward(m, 1000.0) ≈ 1.0 atol=1e-10
    @test_throws DomainError LogarithmicMap(0.0, 1000.0)
end

@testset "TanhMap" begin
    m = TanhMap(0.0, 2000.0, 2.5)
    @test forward(m, 0.0) ≈ -1.0 atol=1e-10
    @test forward(m, 2000.0) ≈ 1.0 atol=1e-10
    @test inverse(m, forward(m, 500.0)) ≈ 500.0 atol=1e-10
end

@testset "CfdWallMap" begin
    m = CfdWallMap(500.0, 3.5)
    @test forward(m, 0.0) ≈ -1.0 atol=1e-10
    @test forward(m, 500.0) ≈ 1.0 atol=1e-10
    @test dzdξ(m, -1.0) > 0.0
end

@testset "profile_transform" begin
    m = TanhMap(0.0, 1000.0, 2.0)
    z_profile = [10.0, 100.0, 500.0, 1000.0]
    ξ_profile = profile_transform(m, z_profile)
    @test ξ_profile[1] < ξ_profile[2] < ξ_profile[3] < ξ_profile[4]
end

@testset "JacobianStack" begin
    m1 = TanhMap(0.0, 1000.0, 2.0)
    m2 = HyperbolicMap(0.0, 1000.0, 1.0)
    stack = JacobianStack([m1, m2], ["TanhMap", "HyperbolicMap"])
    ξ, J1, J2 = evaluate_jacobian_stack(stack, z = 500.0)
    @test -1.0 <= ξ <= 1.0
    @test J1 > 0.0
end

@testset "to_latex dispatch" begin
    lm = LinearMap(0.0, 1000.0)
    lm_doc = to_latex(lm)
    @test occursin("Uniform Linear Coordinate Transformation", lm_doc)
    @test occursin("\\mathcal{F}(z)", lm_doc)
    @test !occursin("\$\$", lm_doc)

    hm = HyperbolicMap(0.0, 1000.0, 3.0)
    hm_doc = to_latex(hm)
    @test occursin("Clustered Hyperbolic Stretched Grid Transformation", hm_doc)
    @test occursin("\\alpha = 3.0", hm_doc)

    tm = TanhMap(0.5, 50.0, 2.3)
    tm_doc = to_latex(tm)
    @test occursin("Symmetric Hyperbolic Tangent Transformation", tm_doc)
    @test occursin("2.3", tm_doc)

    custom = CustomMap(x -> x, x -> x, _ -> 1.0, _ -> 0.0)
    custom_doc = to_latex(custom)
    @test occursin("Arbitrary Numerical Transformation Layer", custom_doc)
end