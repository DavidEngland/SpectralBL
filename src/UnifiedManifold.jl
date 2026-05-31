module UnifiedManifold

using LinearAlgebra

export UnifiedManifoldWorkspace

struct UnifiedManifoldWorkspace{T<:AbstractFloat}
    N::Int
    K_q::Int
    xi_target::Vector{T}
    xi_q::Vector{T}
    J_q::Vector{T}
    psi_M::Vector{T}; psi_W::Vector{T}; psi_T::Vector{T}
    Manifold_Mass::Matrix{T}
    Mass_Factored::Cholesky{T, Matrix{T}}
    z_atm::Vector{T}
    Dz_atm::Matrix{T}

    function UnifiedManifoldWorkspace(N::Int, z_0m::T, z_top::T, alpha_stretch::T; 
                                     n_m=3, n_w=12, delta=1.2, K_q::Int=72) where {T<:AbstractFloat}
        
        xi_q = [cos(pi * (2k - 1) / (2K_q)) for k in 1:K_q]
        sigma = (z_top - z_0m) * alpha_stretch / 2.0
        J_q = [sigma * (2.0 + alpha_stretch) / (1.0 - x + alpha_stretch)^2 for x in xi_q]
        
        xi_target = [cos(pi * i / N) for i in 0:N]
        z_atm = [z_0m + sigma * (1.0 + x) / (1.0 - x + alpha_stretch) for x in xi_target]
        inv_J_target = [1.0 / (sigma * (2.0 + alpha_stretch) / (1.0 - x + alpha_stretch)^2) for x in xi_target]
        
        D1 = zeros(T, N+1, N+1)
        for i in 1:(N+1), j in 1:(N+1)
            if i == j
                if i == 1;     D1[i,j] = (2.0 * N^2 + 1.0) / 6.0
                elseif i == N+1; D1[i,j] = -(2.0 * N^2 + 1.0) / 6.0
                else;          D1[i,j] = -xi_target[i] / (2.0 * (1.0 - xi_target[i]^2))
                end
            else
                c_i = (i == 1 || i == N+1) ? 2.0 : 1.0
                c_j = (j == 1 || j == N+1) ? 2.0 : 1.0
                D1[i,j] = (c_i / c_j) * ((-1)^(i+j)) / (xi_target[i] - xi_target[j])
            end
        end
        Dz_atm = zeros(T, N+1, N+1)
        for i in 1:(N+1); Dz_atm[i, :] = inv_J_target[i] .* D1[i, :]; end

        Manifold_Mass = zeros(T, N+1, N+1)
        for m in 0:N, n in 0:N
            val = 0.0
            for k in 1:K_q
                val += cos(m * acos(xi_q[k])) * cos(n * acos(xi_q[k])) * J_q[k] * sin(acos(xi_q[k]))
            end
            Manifold_Mass[m+1, n+1] = val * (pi / K_q)
        end
        Mass_Factored = cholesky(Hermitian(Manifold_Mass))

        psi_M = zeros(T, N + 1); psi_W = zeros(T, N + 1); psi_T = zeros(T, N + 1)
        for i in 1:(N+1)
            n = i - 1
            psi_M[i] = 0.5 * (1.0 - tanh((n - n_m) / delta))
            psi_W[i] = 0.5 * (1.0 + tanh((n - n_m) / delta)) * 0.5 * (1.0 - tanh((n - n_w) / delta))
            psi_T[i] = 1.0 - psi_M[i] - psi_W[i]
        end

        new{T}(N, K_q, xi_target, xi_q, J_q, psi_M, psi_W, psi_T, Manifold_Mass, Mass_Factored, z_atm, Dz_atm)
    end
end

end
