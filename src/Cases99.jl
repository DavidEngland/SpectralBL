# src/Cases99.jl
module UnifiedManifold

using LinearAlgebra
include("transforms.jl")
using .Transforms

export UnifiedManifoldWorkspace, physical_to_computational
export CoordinateMap, LinearMap, HyperbolicMap, LogarithmicMap, TanhMap, CfdWallMap, CustomMap
export to_latex

"""
    UnifiedManifoldWorkspace(N, z_0m, z_top, alpha_stretch; ...)

Constructs a metric-consistent Riemannian geometry using Chebyshev polynomials T_n(ξ)
as the spectral basis. The physical-to-computational mapping via a coordinate-map compactification
(alpha_stretch parameter) ensures dense nodal concentration near z_0m where CASES-99
inversions are sharpest.
"""
struct UnifiedManifoldWorkspace{T<:AbstractFloat}
    N::Int
    K_q::Int
    xi_target::Vector{T}
    xi_q::Vector{T}
    J_q::Vector{T}
    psi_M::Vector{T}
    psi_W::Vector{T}
    psi_T::Vector{T}
    psi_M_z::Vector{T}
    psi_W_z::Vector{T}
    psi_T_z::Vector{T}
    Manifold_Mass::Matrix{T}
    Mass_Factored::Cholesky{T, Matrix{T}}
    z_atm::Vector{T}
    Dz_atm::Matrix{T}
    map::Transforms.CoordinateMap
    alpha_stretch::T
    sigma::T

    # Clean Inner Constructor Block
    function UnifiedManifoldWorkspace(N::Int, z_0m::T, z_top::T, alpha_stretch::T;
                                     n_m=3, n_w=12, delta=1.2, K_q::Int=72,
                                     map::Transforms.CoordinateMap=TanhMap(z_0m, z_top, alpha_stretch),
                                     z_surf::T=T(10.0), z_top_anchor::T=T(40.0), delta_z::T=T(5.0)) where {T<:AbstractFloat}

        xi_q = [cos(pi * (2k - 1) / (2K_q)) for k in 1:K_q]

        # Canonical scale retained for compatibility with legacy diagnostics.
        sigma = (z_top - z_0m) / 2.0

        # Use the selected map object as the source of truth for metric Jacobians.
        J_q = T[dzdξ(map, x) for x in xi_q]

        xi_target = [cos(pi * i / N) for i in 0:N]
        z_atm = T[inverse(map, x) for x in xi_target]
        inv_J_target = T[dξdz(map, z) for z in z_atm]

        # Chebyshev Differentiation Matrix Assembly
        D1 = zeros(T, N+1, N+1)
        for i in 1:(N+1), j in 1:(N+1)
            if i == j
                if i == 1;       D1[i,j] = (2.0 * N^2 + 1.0) / 6.0
                elseif i == N+1; D1[i,j] = -(2.0 * N^2 + 1.0) / 6.0
                else;            D1[i,j] = -xi_target[i] / (2.0 * (1.0 - xi_target[i]^2))
                end
            else
                c_i = (i == 1 || i == N+1) ? 2.0 : 1.0
                c_j = (j == 1 || j == N+1) ? 2.0 : 1.0
                D1[i,j] = (c_i / c_j) * ((-1)^(i+j)) / (xi_target[i] - xi_target[j])
            end
        end
        Dz_atm = zeros(T, N+1, N+1)
        for i in 1:(N+1); Dz_atm[i, :] = inv_J_target[i] .* D1[i, :]; end

        # Consistent Mass Matrix under Physical Metric Weight Cancellation
        Manifold_Mass = zeros(T, N+1, N+1)
        for m in 0:N, n in 0:N
            val = 0.0
            for k in 1:K_q
                val += cos(m * acos(xi_q[k])) * cos(n * acos(xi_q[k])) * J_q[k] * sin(acos(xi_q[k]))
            end
            Manifold_Mass[m+1, n+1] = val * (pi / K_q)
        end
        Mass_Factored = cholesky(Hermitian(Manifold_Mass))

        # Smooth Sub-meso Partitioning Windows
        psi_M = zeros(T, N + 1); psi_W = zeros(T, N + 1); psi_T = zeros(T, N + 1)
        for i in 1:(N+1)
            n = i - 1
            psi_M[i] = 0.5 * (1.0 - tanh((n - n_m) / delta))
            psi_W[i] = 0.5 * (1.0 + tanh((n - n_m) / delta)) * 0.5 * (1.0 - tanh((n - n_w) / delta))
            psi_T[i] = 1.0 - psi_M[i] - psi_W[i]
        end

        psi_M_z, psi_W_z, psi_T_z = compute_coordinate_windows(z_atm, z_surf, z_top_anchor, delta_z)

        new{T}(N, K_q, xi_target, xi_q, J_q, psi_M, psi_W, psi_T, psi_M_z, psi_W_z, psi_T_z, Manifold_Mass, Mass_Factored, z_atm, Dz_atm, map, alpha_stretch, sigma)
    end
end

function compute_coordinate_windows(z_nodes::Vector{T}, z_surf::T, z_top_anchor::T, delta_z::T) where {T<:AbstractFloat}
    z_lo, z_hi = minimum(z_nodes), maximum(z_nodes)
    z_s = clamp(z_surf, z_lo, z_hi)
    z_t = clamp(z_top_anchor, z_lo, z_hi)
    if z_s > z_t
        z_s, z_t = z_t, z_s
    end
    δz = max(abs(delta_z), sqrt(eps(T)))

    psi_M_z = zeros(T, length(z_nodes))
    psi_W_z = zeros(T, length(z_nodes))
    psi_T_z = zeros(T, length(z_nodes))
    for i in eachindex(z_nodes)
        z = z_nodes[i]
        psi_T_z[i] = 0.5 * (1.0 - tanh((z - z_s) / δz))
        psi_M_z[i] = 0.5 * (1.0 + tanh((z - z_t) / δz))
        psi_W_z[i] = 1.0 - psi_M_z[i] - psi_T_z[i]
    end

    @assert maximum(abs.(psi_M_z .+ psi_W_z .+ psi_T_z .- 1.0)) < 1e-8 "Coordinate windows violate partition of unity"
    @assert minimum(vcat(psi_M_z, psi_W_z, psi_T_z)) > -1e-8 "Coordinate windows contain negative values"

    i_surface = argmin(z_nodes)
    i_top = argmax(z_nodes)
    @assert psi_T_z[i_surface] >= psi_M_z[i_surface] && psi_T_z[i_surface] >= psi_W_z[i_surface] "Expected turbulence dominance near surface not satisfied"
    @assert psi_M_z[i_top] >= psi_T_z[i_top] && psi_M_z[i_top] >= psi_W_z[i_top] "Expected synoptic dominance aloft not satisfied"

    i_peak_w = argmax(psi_W_z)
    @assert i_peak_w != i_surface && i_peak_w != i_top "Wave window peak should be interior, not at boundaries"

    return psi_M_z, psi_W_z, psi_T_z
end

# --- Robust Outer Method Definition ---
"""
    physical_to_computational(ws, z_phys)

Maps physical heights z directly to computational coordinates ξ ∈ [-1, 1]
by analytically inverting the selected coordinate map profile.
"""
function physical_to_computational(ws::UnifiedManifoldWorkspace{T}, z_phys::Vector{T}) where {T<:AbstractFloat}
    xi = T[forward(ws.map, z) for z in z_phys]
    return clamp.(xi, T(-1.0), T(1.0))
end

end # module