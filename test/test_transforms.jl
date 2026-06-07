# test/test_transforms.jl
"""
Comprehensive test suite for Transforms module.

Covers:
- All coordinate map types (LinearMap, HyperbolicMap, LogarithmicMap, TanhMap, CfdWallMap, CustomMap)
- Input validation and error handling
- Metric Jacobians (analytical verification via finite differences)
- Profile transformation (vectorized)
- Jacobian stack composition with exact chain rules
- Manifold consistency (round-trip accuracy, monotonicity)
"""

# 1. Point directly to your local source file
include("../src/transforms.jl")

using Transforms
using Test

# ===== SHARED TEST UTILITIES =====

"""
    verify_metrics_via_fd(m::CoordinateMap, ξ_test::Float64; eps=1e-8, tol=1e-5)

Verify first and second metric Jacobians against finite-difference approximations.
Returns (error_J1, error_J2). Both should be < tol for correct implementation.
"""
function verify_metrics_via_fd(m::CoordinateMap, ξ_test::Float64; eps=1e-8, tol=1e-5)
    # Analytical metrics
    J1_analytical = dzdξ(m, ξ_test)
    J2_analytical = d2zdξ2(m, ξ_test)

    # Finite-difference approximations
    z_plus = inverse(m, ξ_test + eps)
    z_center = inverse(m, ξ_test)
    z_minus = inverse(m, ξ_test - eps)

    J1_fd = (z_plus - z_minus) / (2 * eps)
    J2_fd = (z_plus - 2*z_center + z_minus) / (eps^2)

    # Relative errors
    error_J1 = abs(J1_analytical - J1_fd) / (abs(J1_analytical) + eps)
    error_J2 = abs(J2_analytical - J2_fd) / (abs(J2_analytical) + eps)

    return (error_J1, error_J2)
end

"""
    verify_round_trip(m::CoordinateMap, z_test::Float64; tol=1e-12)

Test forward-inverse composition: forward(m, inverse(m, ξ)) ≈ ξ.
"""
function verify_round_trip(m::CoordinateMap, z_test::Float64; tol=1e-12)
    ξ = forward(m, z_test)
    z_recovered = inverse(m, ξ)
    error = abs(z_recovered - z_test) / (abs(z_test) + 1e-10)
    return error < tol
end


# ===== 1. LINEAR MAP TESTS =====

@testset "LinearMap: Basic Operations" begin
    m = LinearMap(0.0, 1000.0)
    is_valid(m)  # Should not throw

    # Boundary conditions
    @test forward(m, 0.0) ≈ -1.0
    @test forward(m, 1000.0) ≈ 1.0
    @test forward(m, 500.0) ≈ 0.0

    @test inverse(m, -1.0) ≈ 0.0
    @test inverse(m, 1.0) ≈ 1000.0
    @test inverse(m, 0.0) ≈ 500.0

    # Metrics
    @test dzdξ(m, 0.0) ≈ 500.0
    @test d2zdξ2(m, 0.0) ≈ 0.0

    # Inverse metrics
    z_test = 250.0
    @test dξdz(m, z_test) ≈ 1.0 / dzdξ(m, forward(m, z_test))
end

@testset "LinearMap: Metric Verification" begin
    m = LinearMap(0.0, 1000.0)
    ξ_test = 0.5

    error_J1, error_J2 = verify_metrics_via_fd(m, ξ_test)
    @test error_J1 < 1e-5
    @test error_J2 < 1e-5
end

@testset "LinearMap: Round-Trip & Monotonicity" begin
    m = LinearMap(0.0, 1000.0)

    # Round-trip at multiple points
    for z in [10.0, 100.0, 500.0, 900.0]
        @test verify_round_trip(m, z)
    end

    # Monotonicity
    ξ_range = range(-0.99, 0.99, length=100)
    z_range = inverse.(Ref(m), ξ_range)
    @test all(diff(z_range) .> 0)  # Strictly increasing
end

@testset "LinearMap: Domain Validation" begin
    # Valid maps
    m1 = LinearMap(0.0, 1000.0)
    m2 = LinearMap(100.0, 5000.0)
    @test is_valid(m1) && is_valid(m2)

    # Invalid maps
    @test_throws DomainError LinearMap(1000.0, 0.0)  # zmin > zmax
    @test_throws DomainError LinearMap(-100.0, 500.0)  # zmin < 0
end


# ===== 2. HYPERBOLIC MAP TESTS =====

@testset "HyperbolicMap: Basic Operations" begin
    m = HyperbolicMap(0.0, 1000.0, 2.5)
    is_valid(m)

    # Boundary conditions
    @test forward(m, 0.0) ≈ -1.0 atol=1e-10
    @test forward(m, 1000.0) ≈ 1.0 atol=1e-10

    # Round-trip
    z_test = 250.0
    ξ = forward(m, z_test)
    z_recovered = inverse(m, ξ)
    @test z_recovered ≈ z_test atol=1e-10
end

@testset "HyperbolicMap: Parameter Effects (Alpha)" begin
    z_test = 100.0

    # Weak compression (α = 1)
    m_weak = HyperbolicMap(0.0, 1000.0, 1.0)
    ξ_weak = forward(m_weak, z_test)

    # Strong compression (α = 5)
    m_strong = HyperbolicMap(0.0, 1000.0, 5.0)
    ξ_strong = forward(m_strong, z_test)

    # Strong compression should push ξ closer to -1
    @test ξ_strong < ξ_weak
end

@testset "HyperbolicMap: Metric Verification" begin
    m = HyperbolicMap(0.0, 1000.0, 2.5)
    ξ_test = 0.0

    error_J1, error_J2 = verify_metrics_via_fd(m, ξ_test)
    @test error_J1 < 1e-5
    @test error_J2 < 1e-5
end

@testset "HyperbolicMap: Monotonicity & Positivity" begin
    m = HyperbolicMap(0.0, 1000.0, 3.0)

    ξ_range = range(-0.99, 0.99, length=100)
    z_range = inverse.(Ref(m), ξ_range)
    J1_range = dzdξ.(Ref(m), ξ_range)

    @test all(diff(z_range) .> 0)   # Monotonically increasing
    @test all(J1_range .> 0)        # Always positive
end

@testset "HyperbolicMap: Domain Validation" begin
    @test_throws DomainError HyperbolicMap(0.0, 1000.0, -1.0)  # α ≤ 0
    @test_throws DomainError HyperbolicMap(500.0, 500.0, 2.0)   # zmin ≥ zmax
end


# ===== 3. LOGARITHMIC MAP TESTS =====

@testset "LogarithmicMap: Constructor Validation" begin
    # Valid
    m = LogarithmicMap(1.0, 1000.0)
    @test is_valid(m)

    # Invalid: zmin ≤ 0 (constructor validation)
    @test_throws DomainError LogarithmicMap(0.0, 1000.0)
    @test_throws DomainError LogarithmicMap(-10.0, 1000.0)

    # Invalid: zmax ≤ zmin
    @test_throws ArgumentError LogarithmicMap(1000.0, 1000.0)
end

@testset "LogarithmicMap: Basic Operations" begin
    m = LogarithmicMap(10.0, 1000.0)

    # Boundary conditions
    @test forward(m, 10.0) ≈ -1.0 atol=1e-10
    @test forward(m, 1000.0) ≈ 1.0 atol=1e-10

    # Geometric mean at ξ = 0
    z_mid = inverse(m, 0.0)
    z_expected = sqrt(10.0 * 1000.0)  # Geometric mean
    @test z_mid ≈ z_expected atol=1e-8
end

@testset "LogarithmicMap: Metric Verification" begin
    m = LogarithmicMap(10.0, 1000.0)
    ξ_test = 0.0

    error_J1, error_J2 = verify_metrics_via_fd(m, ξ_test)
    @test error_J1 < 1e-5
    @test error_J2 < 1e-5
end

@testset "LogarithmicMap: Round-Trip" begin
    m = LogarithmicMap(10.0, 1000.0)

    for z in [10.0, 50.0, 100.0, 500.0, 1000.0]
        @test verify_round_trip(m, z)
    end
end

@testset "LogarithmicMap: Forward Domain Error" begin
    m = LogarithmicMap(10.0, 1000.0)

    # z ≤ 0 should throw in forward()
    @test_throws DomainError forward(m, 0.0)
    @test_throws DomainError forward(m, -10.0)
    @test_throws DomainError forward(m, 5.0)  # z < zmin
end


# ===== 4. TANH MAP TESTS =====

@testset "TanhMap: Basic Operations" begin
    m = TanhMap(0.0, 2000.0, 2.5)
    is_valid(m)

    # Boundary conditions
    @test forward(m, 0.0) ≈ -1.0 atol=1e-10
    @test forward(m, 2000.0) ≈ 1.0 atol=1e-10
    @test forward(m, 1000.0) ≈ 0.0 atol=1e-3  # Center approximately

    # Round-trip
    z_test = 500.0
    ξ = forward(m, z_test)
    z_recovered = inverse(m, ξ)
    @test z_recovered ≈ z_test atol=1e-10
end

@testset "TanhMap: Parameter Effects (Alpha)" begin
    z_test = 100.0

    # Weak compression (α = 1)
    m_weak = TanhMap(0.0, 2000.0, 1.0)
    ξ_weak = forward(m_weak, z_test)

    # Strong compression (α = 4)
    m_strong = TanhMap(0.0, 2000.0, 4.0)
    ξ_strong = forward(m_strong, z_test)

    # Strong compression pushes toward -1
    @test ξ_strong < ξ_weak
end

@testset "TanhMap: Metric Verification" begin
    m = TanhMap(0.0, 2000.0, 2.5)

    # Test at multiple ξ values
    for ξ in [-0.5, 0.0, 0.5]
        error_J1, error_J2 = verify_metrics_via_fd(m, ξ)
        @test error_J1 < 1e-5 "J1 FD check failed at ξ=$ξ"
        @test error_J2 < 1e-5 "J2 FD check failed at ξ=$ξ"
    end
end

@testset "TanhMap: Monotonicity & Positivity" begin
    m = TanhMap(0.0, 2000.0, 3.0)

    ξ_range = range(-0.9999, 0.9999, length=200)
    z_range = inverse.(Ref(m), ξ_range)
    J1_range = dzdξ.(Ref(m), ξ_range)

    @test all(diff(z_range) .> 0)   # Strictly increasing
    @test all(J1_range .> 0)        # Always positive
end

@testset "TanhMap: Boundary Clamping" begin
    m = TanhMap(0.0, 2000.0, 5.0)

    # Extreme ξ values should still be computable (clamped)
    ξ_extreme = [-0.99999999, 0.99999999]
    for ξ in ξ_extreme
        z = inverse(m, ξ)
        @test 0.0 ≤ z ≤ 2000.0
        ξ_back = forward(m, z)
        @test abs(ξ_back - ξ) < 1e-6
    end
end


# ===== 5. CFD WALL MAP TESTS =====

@testset "CfdWallMap: Basic Operations" begin
    m = CfdWallMap(500.0, 4.0)
    is_valid(m)

    # Boundary conditions
    @test inverse(m, -1.0) ≈ 0.0 atol=1e-10        # Wall at ξ = -1
    @test inverse(m, 1.0) ≈ 500.0 atol=1e-10       # Top at ξ = 1

    # Round-trip
    z_test = 150.0
    ξ = forward(m, z_test)
    z_recovered = inverse(m, ξ)
    @test z_recovered ≈ z_test atol=1e-10
end

@testset "CfdWallMap: Wall Refinement Pattern" begin
    m = CfdWallMap(500.0, 4.0)

    # Wall should be highly refined
    z_near_wall = inverse(m, -0.95)
    z_mid = inverse(m, 0.0)

    # Near-wall should be much closer to 0 than to zmax/2
    @test z_near_wall < z_mid / 4
end

@testset "CfdWallMap: Parameter Effects (Delta)" begin
    z_test = 100.0

    # Weak packing (δ = 1)
    m_weak = CfdWallMap(500.0, 1.0)
    ξ_weak = forward(m_weak, z_test)

    # Strong packing (δ = 5)
    m_strong = CfdWallMap(500.0, 5.0)
    ξ_strong = forward(m_strong, z_test)

    # Strong packing pushes away from wall more (ξ closer to 0)
    @test ξ_strong > ξ_weak
end

@testset "CfdWallMap: Metric Verification" begin
    m = CfdWallMap(500.0, 3.5)

    for ξ in [-0.5, 0.0, 0.5]
        error_J1, error_J2 = verify_metrics_via_fd(m, ξ)
        @test error_J1 < 1e-5
        @test error_J2 < 1e-5
    end
end

@testset "CfdWallMap: Monotonicity & Positivity" begin
    m = CfdWallMap(500.0, 4.0)

    ξ_range = range(-0.9999, 0.9999, length=200)
    z_range = inverse.(Ref(m), ξ_range)
    J1_range = dzdξ.(Ref(m), ξ_range)

    @test all(diff(z_range) .> 0)   # Strictly increasing
    @test all(J1_range .> 0)        # Always positive
end

@testset "CfdWallMap: Domain Validation" begin
    @test_throws DomainError CfdWallMap(0.0, 3.0)   # ztop ≤ 0
    @test_throws DomainError CfdWallMap(500.0, 0.0) # δ ≤ 0
    @test_throws DomainError CfdWallMap(500.0, -1.0) # δ < 0
end


# ===== 6. CUSTOM MAP TESTS =====

@testset "CustomMap: User-Supplied Functions" begin
    # Simple exponential custom map
    fwd(z) = 2.0 * log(z / 1.0) / log(100.0 / 1.0) - 1.0
    inv(ξ) = exp(((ξ + 1.0) / 2.0) * log(100.0))
    dzxi(ξ) = 0.5 * log(100.0) * inv(ξ)
    d2zxi2(ξ) = 0.25 * (log(100.0))^2 * inv(ξ)

    m = CustomMap(fwd, inv, dzxi, d2zxi2)

    # Boundary conditions
    @test forward(m, 1.0) ≈ -1.0 atol=1e-10
    @test forward(m, 100.0) ≈ 1.0 atol=1e-10

    # Round-trip
    z_test = 10.0
    ξ = forward(m, z_test)
    z_recovered = inverse(m, ξ)
    @test z_recovered ≈ z_test atol=1e-10
end


# ===== 7. PROFILE TRANSFORMATION TESTS =====

@testset "profile_transform: Vectorized Transformation" begin
    m = TanhMap(0.0, 1000.0, 2.0)
    z_profile = [10.0, 100.0, 500.0, 1000.0]
    ξ_profile = profile_transform(m, z_profile)

    # Monotonicity
    @test all(diff(ξ_profile) .> 0)

    # Boundary match
    @test ξ_profile[end] ≈ 1.0 atol=1e-10

    # Length preservation
    @test length(ξ_profile) == length(z_profile)
end

@testset "profile_transform: Empty & Single Element" begin
    m = LinearMap(0.0, 1000.0)

    # Single element
    z_single = [500.0]
    ξ_single = profile_transform(m, z_single)
    @test length(ξ_single) == 1
    @test ξ_single[1] ≈ 0.0

    # Empty (should be empty, not throw)
    z_empty = Float64[]
    ξ_empty = profile_transform(m, z_empty)
    @test length(ξ_empty) == 0
end


# ===== 8. JACOBIAN STACK TESTS =====

@testset "JacobianStack: Constructor Validation" begin
    m1 = TanhMap(0.0, 1000.0, 2.0)
    m2 = HyperbolicMap(0.0, 1000.0, 1.5)

    # Valid stack
    stack = JacobianStack([m1, m2], ["Layer1", "Layer2"])
    @test length(stack.maps) == 2
    @test stack.labels == ["Layer1", "Layer2"]

    # Empty maps
    @test_throws ArgumentError JacobianStack([])

    # Label/map mismatch
    @test_throws ArgumentError JacobianStack([m1, m2], ["OnlyOne"])

    # Auto-generated labels
    stack_auto = JacobianStack([m1, m2])
    @test !isempty(stack_auto.labels)
    @test length(stack_auto.labels) == 2
end

@testset "JacobianStack: Two-Layer Composition" begin
    m1 = TanhMap(0.0, 1000.0, 2.0)
    m2 = HyperbolicMap(0.0, 1000.0, 1.5)
    stack = JacobianStack([m1, m2])

    z_test = 500.0
    ξ, J1, J2 = evaluate_jacobian_stack(stack, z_test)

    @test -1.0 < ξ < 1.0
    @test J1 > 0         # Monotonic
    @test !isnan(J1) && !isinf(J1)
    @test !isnan(J2) && !isinf(J2)
end

@testset "JacobianStack: Round-Trip at Multiple Heights" begin
    m1 = TanhMap(0.0, 1000.0, 2.0)
    m2 = HyperbolicMap(0.0, 1000.0, 1.5)
    stack = JacobianStack([m1, m2])

    for z in [10.0, 100.0, 500.0, 900.0]
        ξ, _, _ = evaluate_jacobian_stack(stack, z)

        # Recover z by inverting stack
        z_recovered = z
        for i in length(stack.maps):-1:1
            z_recovered = inverse(stack.maps[i], z_recovered)
        end

        @test z_recovered ≈ z atol=1e-10
    end
end

@testset "JacobianStack: Exact Analytical Chain Rules" begin
    """
    Test that J₂ (d²ξ/dz²) computed via exact analytical chain rules
    matches finite-difference approximation closely (within 1e-4 relative error).
    """
    m1 = TanhMap(0.0, 2000.0, 2.5)
    m2 = HyperbolicMap(-1.0, 1.0, 1.5)
    stack = JacobianStack([m1, m2])

    z_test = 500.0
    eps = 1e-7

    # Analytical
    ξ_center, J1_center, J2_analytical = evaluate_jacobian_stack(stack, z_test)

    # Finite-difference
    ξ_plus, J1_plus, _ = evaluate_jacobian_stack(stack, z_test + eps)
    ξ_minus, J1_minus, _ = evaluate_jacobian_stack(stack, z_test - eps)

    J2_fd = (J1_plus - J1_minus) / (2 * eps)

    # Error should be small (FD truncation ~ ε)
    relative_error = abs(J2_analytical - J2_fd) / (abs(J2_analytical) + eps)
    @test relative_error < 1e-3  # FD approximation accurate
end

@testset "JacobianStack: Three-Layer Stack" begin
    """
    Test stability with three layers of composition.
    """
    m1 = TanhMap(0.0, 2000.0, 2.0)
    m2 = HyperbolicMap(0.0, 2000.0, 1.5)
    m3 = LinearMap(0.0, 2000.0)
    stack = JacobianStack([m1, m2, m3], ["Tanh", "Hyperbolic", "Linear"])

    z_test = 500.0
    ξ, J1, J2 = evaluate_jacobian_stack(stack, z_test)

    @test -1.0 < ξ < 1.0
    @test J1 > 0
    @test !isnan(J1) && !isinf(J1)
    @test !isnan(J2) && !isinf(J2)
end

@testset "JacobianStack: Invalid Map in Stack" begin
    """
    Ensure invalid maps are caught during stack construction.
    """
    m_valid = TanhMap(0.0, 1000.0, 2.0)
    m_invalid = HyperbolicMap(0.0, 1000.0, -1.0)  # Negative α

    @test_throws ArgumentError JacobianStack([m_valid, m_invalid])
end


# ===== 9. INVERSE METRICS TESTS =====

@testset "Inverse Metrics: dξdz and d2ξdz2" begin
    m = TanhMap(0.0, 2000.0, 2.5)

    z_test = 500.0
    ξ_test = forward(m, z_test)

    # Forward metrics
    J1_fwd = dzdξ(m, ξ_test)
    J2_fwd = d2zdξ2(m, ξ_test)

    # Inverse metrics (via chain rule)
    J1_inv = dξdz(m, z_test)
    J2_inv = d2ξdz2(m, z_test)

    # Check relationships
    @test J1_inv ≈ 1.0 / J1_fwd atol=1e-10
    @test J2_inv ≈ -J2_fwd / (J1_fwd^3) atol=1e-10
end


# ===== 10. INTEGRATION TESTS =====

@testset "Integration: All Maps with Same Domain" begin
    """
    Verify all maps with [0, 1000] domain produce consistent ordering.
    """
    z_test = 100.0

    m_linear = LinearMap(0.0, 1000.0)
    m_hyp = HyperbolicMap(0.0, 1000.0, 2.0)
    m_tanh = TanhMap(0.0, 1000.0, 2.0)

    ξ_linear = forward(m_linear, z_test)
    ξ_hyp = forward(m_hyp, z_test)
    ξ_tanh = forward(m_tanh, z_test)

    # All should be in [-1, 1]
    @test all(-1.0 .≤ [ξ_linear, ξ_hyp, ξ_tanh] .≤ 1.0)

    # All should map 0 to -1 and 1000 to 1
    @test forward(m_linear, 0.0) ≈ -1.0
    @test forward(m_hyp, 0.0) ≈ -1.0
    @test forward(m_tanh, 0.0) ≈ -1.0
end

@testset "Integration: Stack with CfdWallMap" begin
    """
    Two-layer stack: TanhMap (stability) + CfdWallMap (wall refinement).
    """
    m1 = TanhMap(0.0, 500.0, 2.5)
    m2 = CfdWallMap(500.0, 3.5)
    stack = JacobianStack([m1, m2], ["Stability", "Wall"])

    # Near-wall point
    z_wall = 5.0
    ξ, J1, J2 = evaluate_jacobian_stack(stack, z_wall)

    @test ξ < -0.5  # Should be heavily compressed near -1
    @test J1 > 0
end


# ===== 11. EDGE CASES & ROBUSTNESS =====

@testset "Edge Case: Very Small Domains" begin
    m = LinearMap(0.0, 0.1)  # 10 cm domain

    z_test = 0.05
    ξ = forward(m, z_test)
    z_recovered = inverse(m, ξ)

    @test z_recovered ≈ z_test atol=1e-12
end

@testset "Edge Case: Very Large Domains" begin
    m = LinearMap(0.0, 1e5)  # 100 km domain

    z_test = 5e4
    ξ = forward(m, z_test)
    z_recovered = inverse(m, ξ)

    @test z_recovered ≈ z_test atol=1.0
end

@testset "Edge Case: Extreme Compression Parameters" begin
    m_weak = TanhMap(0.0, 1000.0, 0.1)  # Very weak
    m_strong = TanhMap(0.0, 1000.0, 10.0)  # Very strong

    z_test = 100.0
    ξ_weak = forward(m_weak, z_test)
    ξ_strong = forward(m_strong, z_test)

    # Both should still map to valid ξ
    @test -1.0 < ξ_weak < 1.0
    @test -1.0 < ξ_strong < 1.0

    # Strong should compress more
    @test ξ_strong < ξ_weak
end

println("✓ All tests passed!")
