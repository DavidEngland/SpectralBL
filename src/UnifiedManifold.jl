# src/UnifiedManifold.jl
module UnifiedManifold

using LinearAlgebra

export UnifiedManifoldWorkspace, physical_to_computational

"""
    UnifiedManifoldWorkspace(N, z_min, z_max, alpha_stretch; ...)

Constructs a metric-consistent Riemannian geometry matching the TanhMap formulation.
Grid nodes are clustered symmetrically near boundaries to capture intense nocturnal SBL
inversions and low-level shear layers without numerical aliasing.
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
    Manifold_Mass::Matrix{T}
    Mass_Factored::Cholesky{T, Matrix{T}}
    z_atm::Vector{T}
    Dz_atm::Matrix{T}
    alpha_stretch::T
    L_domain::T
    z_center::T

    function UnifiedManifoldWorkspace(N::Int, z_min::T, z_max::T, alpha_stretch::T;
                                     n_m=3, n_w=12, delta=1.2, K_q::Int=72, invert_windows::Bool=false) where {T<:AbstractFloat}

        L_domain = z_max - z_min
        z_center = (z_max + z_min) / 2.0

        # Quadrature Nodes (Computational Space)
        xi_q = [cos(pi * (2k - 1) / (2K_q)) for k in 1:K_q]

        # --- RECONCILED TanhMap METRIC JACOBIAN: dz/dξ ---
        # J = (L / 2) * (α / tanh(α)) * sech²(α * ξ)
        J_q = [(L_domain / 2.0) * (alpha_stretch / tanh(alpha_stretch)) * (sech(alpha_stretch * x))^2 for x in xi_q]

        # Target Collocation Nodes
        xi_target = [cos(pi * i / N) for i in 0:N]

        # --- RECONCILED TanhMap INVERSE PROFILE: ξ -> z ---
        z_atm = [z_center + (L_domain / 2.0) * tanh(alpha_stretch * x) / tanh(alpha_stretch) for x in xi_target]

        # Inverse Jacobian Metric: dξ/dz = 1 / J
        inv_J_target = [1.0 / ((L_domain / 2.0) * (alpha_stretch / tanh(alpha_stretch)) * (sech(alpha_stretch * x))^2) for x in xi_target]

        # Chebyshev Differentiation Matrix Assembly (D1)
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

        # Transform Computational Derivatives into True Physical Metric Space
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

        # Smooth sub-meso partitioning windows in mode space.
        # Optional inversion swaps low/high windows so turbulence can occupy
        # lower-order manifold support when requested by pipeline configuration.
        psi_M = zeros(T, N + 1); psi_W = zeros(T, N + 1); psi_T = zeros(T, N + 1)
        for i in 1:(N+1)
            n = i - 1
            psi_low  = 0.5 * (1.0 - tanh((n - n_m) / delta))
            psi_band = 0.5 * (1.0 + tanh((n - n_m) / delta)) * 0.5 * (1.0 - tanh((n - n_w) / delta))
            psi_high = 1.0 - psi_low - psi_band

            if invert_windows
                psi_M[i] = psi_high
                psi_W[i] = psi_band
                psi_T[i] = psi_low
            else
                psi_M[i] = psi_low
                psi_W[i] = psi_band
                psi_T[i] = psi_high
            end
        end

        new{T}(N, K_q, xi_target, xi_q, J_q, psi_M, psi_W, psi_T, Manifold_Mass, Mass_Factored, z_atm, Dz_atm, alpha_stretch, L_domain, z_center)
    end
end

# --- Reconciled Forward Transformation ---
"""
    physical_to_computational(ws, z_phys)

Maps physical heights z directly to computational coordinates ξ ∈ [-1, 1]
by analytically inverting the TanhMap profile.
"""
function physical_to_computational(ws::UnifiedManifoldWorkspace{T}, z_phys::Vector{T}) where {T<:AbstractFloat}
    xi = zeros(T, length(z_phys))
    for i in eachindex(z_phys)
        # Numerical argument clamp to avoid DomainErrors at outer limits
        arg = clamp(((z_phys[i] - ws.z_center) * tanh(ws.alpha_stretch)) / (ws.L_domain / 2.0), -0.9999999999999, 0.9999999999999)
        xi[i] = atanh(arg) / ws.alpha_stretch
    end
    return clamp.(xi, -1.0, 1.0)
end

end # module