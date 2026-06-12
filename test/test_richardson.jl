#test/test_richardson.jl
include("../src/Richardson.jl")

using .Richardson
using Test

@testset "Richardson operator assembly" begin
    N = 32
    z_min, z_max, alpha = 1.5, 60.0, 2.5

    z, J, D_z, xi, D_xi = build_physical_gradient_operator(N, z_min, z_max, alpha)

    @test length(z) == N + 1
    @test length(J) == N + 1
    @test size(D_z) == (N + 1, N + 1)
    @test size(D_xi) == (N + 1, N + 1)
    @test all(J .> 0.0)

    # Chebyshev endpoints
    @test isapprox(xi[1], 1.0; atol=1e-12)
    @test isapprox(xi[end], -1.0; atol=1e-12)
end

@testset "Richardson profile evaluation" begin
    N = 32
    z_min, z_max, alpha = 1.5, 60.0, 2.5
    z, J, D_z, _, _ = build_physical_gradient_operator(N, z_min, z_max, alpha)

    c_theta = zeros(Float64, N + 1)
    c_u = zeros(Float64, N + 1)

    # Smooth low-mode structure
    c_theta[1] = 295.0
    c_theta[2] = 2.5
    c_theta[3] = -0.7

    c_u[1] = 4.0
    c_u[2] = 3.1
    c_u[3] = 0.5

    out = calculate_pseudospectral_ri_g(c_theta, c_u, D_z; g=9.81, theta_ref=293.15)

    @test length(out.dtheta_dz) == N + 1
    @test length(out.du_dz) == N + 1
    @test length(out.Ri_g) == N + 1
    @test all(isfinite, out.Ri_g)

    full = build_and_evaluate_ri_g(N, z_min, z_max, alpha, c_theta, c_u)
    @test length(full.z) == N + 1
    @test all(isfinite, full.Ri_g)

    # Sanity: strongest grid compression near boundaries for tanh map.
    @test minimum(J) < maximum(J)
    @test z[1] > z[end]
end
