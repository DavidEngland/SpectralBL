module SpectralDiagnostics

using LinearAlgebra

export HighFidelityRecord, process_timestamp_metrics

struct HighFidelityRecord{T<:AbstractFloat}
    time_idx::Int
    Ri_f::T
    R_W::T
    F_W::T
    chi_N::T
    D_eff::T
    Status::String
end

function process_timestamp_metrics(time_idx::Int, c_theta::Vector{T}, c_u::Vector{T}, ws, status_str::String; g = 9.81, theta_ref = 265.0) where {T<:AbstractFloat}
    N = ws.N
    
    c_θ_W = c_theta .* ws.psi_W; c_θ_T = c_theta .* ws.psi_T
    c_u_W = c_u .* ws.psi_W;     c_u_T = c_u .* ws.psi_T
    
    E_tot = dot(c_theta, ws.Manifold_Mass * c_theta) + dot(c_u, ws.Manifold_Mass * c_u)
    E_W   = dot(c_θ_W, ws.Manifold_Mass * c_θ_W) + dot(c_u_W, ws.Manifold_Mass * c_u_W)
    E_T   = dot(c_θ_T, ws.Manifold_Mass * c_θ_T) + dot(c_u_T, ws.Manifold_Mass * c_u_T)
    
    R_W = E_T > 1e-9 ? E_W / E_T : 0.0
    F_W = E_tot > 1e-9 ? E_W / E_tot : 0.0
    
    num_chi = 0.0; den_chi = 0.0
    for n in 1:N
        scale_fact = (n / N)^2
        grad_energy = scale_fact * c_theta[n+1]^2
        den_chi += grad_energy
        if n >= 2; num_chi += scale_fact * ((n - 1) / N)^2 * c_theta[n+1]^2; end
    end
    chi_N = den_chi > 1e-9 ? num_chi / den_chi : 0.0

    entropy = 0.0
    sum_c = sum(c_theta.^2)
    for i in 1:(N+1)
        p_n = c_theta[i]^2 / (sum_c + 1e-12)
        if p_n > 1e-9; entropy -= p_n * log(p_n); end
    end
    D_eff = exp(entropy)

    c_θ_loc = c_theta .- c_θ_W
    c_u_loc = c_u .- c_u_W
    
    theta_loc_profile = zeros(T, N + 1)
    u_loc_profile     = zeros(T, N + 1)
    for i in 1:(N+1)
        xi_i = ws.xi_target[i]
        for n in 0:N
            theta_loc_profile[i] += c_θ_loc[n+1] * cos(n * acos(xi_i))
            u_loc_profile[i]     += c_u_loc[n+1] * cos(n * acos(xi_i))
        end
    end
    
    dtheta_dz_profile = ws.Dz_atm * theta_loc_profile
    du_dz_profile     = ws.Dz_atm * u_loc_profile
    
    shear_sq = du_dz_profile[end]^2
    Ri_f = shear_sq > 1e-6 ? (g / theta_ref) * dtheta_dz_profile[end] / shear_sq : 0.5

    return HighFidelityRecord{T}(time_idx, Ri_f, R_W, F_W, chi_N, D_eff, status_str)
end

end
