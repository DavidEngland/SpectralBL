module SpectralDiagnostics

using LinearAlgebra

export HighFidelityRecord, process_timestamp_metrics, calculate_adaptive_wave_fraction

struct HighFidelityRecord{T<:AbstractFloat}
    time_idx::Int
    Ri_f::T
    R_W::T
    F_W::T
    chi_N::T
    D_eff::T
    E_total::T
    E_wave::T
    E_turb::T
    peak_mode::Int
    wave_window_min::Int
    wave_window_max::Int
    peak_in_wave_window::Bool
    Status::String
end

function process_timestamp_metrics(time_idx::Int, c_theta_raw::Vector{T}, c_u_raw::Vector{T}, ws, status_str::String; g = 9.81, theta_ref = 265.0, wave_threshold = 0.1) where {T<:AbstractFloat}
    N = ws.N
    D_manifold = N + 1 # Dynamic manifold length target (33)
    
    # --- PROJECTION ENHANCEMENT: Zero-pad observational slices to full workspace dimensions ---
    c_theta = zeros(T, D_manifold)
    c_u     = zeros(T, D_manifold)

    # Safely copy available coefficients (up to length 7) into full 33-dimensional space
    len_in = min(length(c_theta_raw), D_manifold)
    c_theta[1:len_in] .= c_theta_raw[1:len_in]
    c_u[1:len_in]     .= c_u_raw[1:len_in]
    # ----------------------------------------------------------------------------------------

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

    # Wave-window coverage QA: verify dominant modal energy falls in active psi_W support.
    modal_energy = c_theta.^2 .+ c_u.^2
    peak_mode = argmax(modal_energy) - 1
    active_modes = findall(x -> x >= wave_threshold, ws.psi_W)
    wave_window_min = isempty(active_modes) ? -1 : minimum(active_modes) - 1
    wave_window_max = isempty(active_modes) ? -1 : maximum(active_modes) - 1
    peak_in_wave_window = !isempty(active_modes) && (peak_mode >= wave_window_min) && (peak_mode <= wave_window_max)
    status_out = peak_in_wave_window ? status_str : string(status_str, " | PeakOutsidePsiW")

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

    return HighFidelityRecord{T}(
        time_idx,
        Ri_f,
        R_W,
        F_W,
        chi_N,
        D_eff,
        E_tot,
        E_W,
        E_T,
        peak_mode,
        wave_window_min,
        wave_window_max,
        peak_in_wave_window,
        status_out,
    )
end

function calculate_adaptive_wave_fraction(ws, c_u::AbstractVector, d_eff; alpha_floor::Real = 1.5, wave_threshold::Real = 0.1)
    n_modes = min(length(c_u), length(ws.psi_W))
    if n_modes == 0
        return (0.0, 0, -1, false, "NoModes")
    end

    modal_energy = abs2.(c_u[1:n_modes])
    e_total = sum(modal_energy)
    if e_total <= eps(Float64)
        return (0.0, 0, -1, false, "ZeroEnergy")
    end

    peak_mode = argmax(modal_energy) - 1

    psi_w = ws.psi_W[1:n_modes]
    active_modes = findall(x -> x >= wave_threshold, psi_w)
    n_min_eff = isempty(active_modes) ? -1 : minimum(active_modes) - 1
    n_max_eff = isempty(active_modes) ? -1 : maximum(active_modes) - 1
    in_window = !isempty(active_modes) && (peak_mode >= n_min_eff) && (peak_mode <= n_max_eff)

    weighted_wave_energy = sum(modal_energy .* psi_w)
    f_w_adaptive = clamp(weighted_wave_energy / (e_total + eps(Float64)), 0.0, 1.0)

    # Keep a minimal adaptive floor for highly compressed modal states.
    if d_eff < alpha_floor
        f_w_adaptive = max(f_w_adaptive, 0.30)
    end

    return (f_w_adaptive, peak_mode, n_min_eff, in_window, "AdaptiveWaveFractionOK")
end

end