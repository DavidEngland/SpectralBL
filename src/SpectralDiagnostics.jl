# src/SpectralDiagnostics.jl
module SpectralDiagnostics

using LinearAlgebra
using Statistics

export HighFidelityRecord, process_timestamp_metrics, calculate_adaptive_wave_fraction

struct HighFidelityRecord{T<:AbstractFloat}
    time_idx::Int
    Ri_g::T   # Gradient Richardson number (averaged near-surface nodes)
    Ri_b::T   # Bulk Richardson number (profile-integrated finite difference)
    R_W::T
    F_W::T
    chi_N::T
    D_eff::T
    E_total::T
    E_meso::T
    E_wave::T
    E_turb::T
    E_interaction::T
    peak_mode::Int
    wave_window_min::Int
    wave_window_max::Int
    peak_in_wave_window::Bool
    Status::String
end

function process_timestamp_metrics(time_idx::Int, c_theta_raw::Vector{T}, c_u_raw::Vector{T}, ws, status_str::String;
    g = 9.81, theta_ref::Union{Nothing,Float64} = nothing, wave_threshold = 0.1,
    entropy_prob_floor = 1e-9, active_mode_floor = 1e-8) where {T<:AbstractFloat}
    N = ws.N
    D_manifold = N + 1 # Dynamic manifold length target (33)

    # Manifold_Mass includes the physical coordinate Jacobian J_q by construction
    # (UnifiedManifoldWorkspace: M[m,n] += T_m * T_n * J_q[k] * sin(acos(xi_q[k])) * π/K_q).
    # Energy inner products are therefore physically consistent without additional correction.

    # Verify spectral partition of unity: psi_M + psi_W + psi_T = 1 for all modes.
    @assert maximum(abs.(ws.psi_M .+ ws.psi_W .+ ws.psi_T .- 1)) < 1e-10 "Spectral windows violate partition of unity"

    # Embed available coefficients into the full manifold dimension.
    # If upstream inversion returns fewer modes than D_manifold, remaining modes are
    # physically unresolved for this timestamp and represented as zeros.
    c_theta = zeros(T, D_manifold)
    c_u     = zeros(T, D_manifold)

    len_theta = min(length(c_theta_raw), D_manifold)
    len_u = min(length(c_u_raw), D_manifold)
    c_theta[1:len_theta] .= c_theta_raw[1:len_theta]
    c_u[1:len_u]         .= c_u_raw[1:len_u]

    c_θ_M = c_theta .* ws.psi_M; c_θ_W = c_theta .* ws.psi_W; c_θ_T = c_theta .* ws.psi_T
    c_u_M = c_u .* ws.psi_M;     c_u_W = c_u .* ws.psi_W;     c_u_T = c_u .* ws.psi_T

    # Window energies under a non-diagonal mass metric are not strictly additive;
    # interaction terms are captured explicitly in E_int.
    E_tot = dot(c_theta, ws.Manifold_Mass * c_theta) + dot(c_u, ws.Manifold_Mass * c_u)
    E_M   = dot(c_θ_M, ws.Manifold_Mass * c_θ_M) + dot(c_u_M, ws.Manifold_Mass * c_u_M)
    E_W   = dot(c_θ_W, ws.Manifold_Mass * c_θ_W) + dot(c_u_W, ws.Manifold_Mass * c_u_W)
    E_T   = dot(c_θ_T, ws.Manifold_Mass * c_θ_T) + dot(c_u_T, ws.Manifold_Mass * c_u_T)
    E_int = E_tot - (E_M + E_W + E_T)
    
    R_W = E_T > 1e-9 ? E_W / E_T : 0.0
    F_W = E_tot > 1e-9 ? E_W / E_tot : 0.0
    
    # χ_N: normalized fourth-moment spectral roughness index.
    # χ_N = (Σ n⁴ cₙ²) / (N² · Σ n² cₙ²)  ∈ [0, 1]
    # High values indicate gradient energy concentrated at fine scales (sharp inversions).
    num_chi = 0.0; den_chi = 0.0
    for n in 1:N
        cn2 = c_theta[n+1]^2
        num_chi += Float64(n)^4 * cn2
        den_chi += Float64(n)^2 * cn2
    end
    chi_N = den_chi > 1e-9 ? num_chi / (Float64(N)^2 * den_chi) : 0.0

    # Shannon entropy over actively resolved modal support only.
    entropy = 0.0
    sum_c = sum(abs2, c_theta)
    active_modes = findall(i -> abs2(c_theta[i]) > active_mode_floor, eachindex(c_theta))
    if isempty(active_modes) || sum_c <= eps(T)
        entropy = 0.0
    else
        for i in active_modes
            p_n = abs2(c_theta[i]) / (sum_c + 1e-12)
            if p_n > entropy_prob_floor
                entropy -= p_n * log(p_n)
            end
        end
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
    if len_theta < D_manifold || len_u < D_manifold
        status_out = string(status_out, " | TruncatedInput(theta=", len_theta, ",u=", len_u, ")")
    end

    c_θ_loc = c_theta .- c_θ_W
    c_u_loc = c_u .- c_u_W

    # Build Chebyshev evaluation matrix T_eval[i, n+1] = Tₙ(ξᵢ) = cos(n·acos(ξᵢ)).
    # Single BLAS matrix-vector multiply replaces the O(N²) scalar loop.
    T_eval = zeros(T, N+1, N+1)
    for (i, xi_i) in enumerate(ws.xi_target)
        xi_c = clamp(xi_i, -one(T), one(T))
        acos_xi = acos(xi_c)
        for n in 0:N
            T_eval[i, n+1] = cos(T(n) * acos_xi)
        end
    end
    theta_loc_profile = T_eval * c_θ_loc
    u_loc_profile     = T_eval * c_u_loc

    theta_ref_eff = isnothing(theta_ref) ? max(mean(theta_loc_profile), 1.0) : theta_ref
    
    dtheta_dz_profile = ws.Dz_atm * theta_loc_profile
    du_dz_profile     = ws.Dz_atm * u_loc_profile

    # Gradient Ri: average over last 3 near-surface nodes to suppress Chebyshev endpoint noise.
    # xi_target[end] = -1 corresponds to the lower physical boundary (z_min).
    avg_idx = max(1, N-1):(N+1)
    dtdz_avg  = mean(dtheta_dz_profile[avg_idx])
    shear_avg = mean(abs2, du_dz_profile[avg_idx]) + 1e-8
    Ri_g = (g / theta_ref_eff) * dtdz_avg / shear_avg

    # Bulk Ri: finite difference between profile-averaged top and bottom layer means.
    # xi_target[1] = +1 → z_top; xi_target[end] = -1 → z_bottom.
    n_avg = max(1, min(3, (N+1) ÷ 4))
    z_top_idx = 1:n_avg
    z_bot_idx = (N+2-n_avg):(N+1)
    dtheta_bulk = mean(theta_loc_profile[z_top_idx]) - mean(theta_loc_profile[z_bot_idx])
    du_bulk     = mean(u_loc_profile[z_top_idx])     - mean(u_loc_profile[z_bot_idx])
    dz_bulk     = mean(ws.z_atm[z_top_idx])          - mean(ws.z_atm[z_bot_idx])
    Ri_b = abs(du_bulk) > 1e-3 ? (g / theta_ref_eff) * dtheta_bulk * dz_bulk / du_bulk^2 : 0.5

    return HighFidelityRecord{T}(
        time_idx,
        Ri_g,
        Ri_b,
        R_W,
        F_W,
        chi_N,
        D_eff,
        E_tot,
        E_M,
        E_W,
        E_T,
        E_int,
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
        return (0.0, 0, -1, false, "NoModes", 1.0)
    end

    modal_energy = abs2.(c_u[1:n_modes])
    e_total = sum(modal_energy)
    if e_total <= eps(Float64)
        return (0.0, 0, -1, false, "ZeroEnergy", 1.0)
    end

    peak_mode = argmax(modal_energy) - 1

    psi_w = ws.psi_W[1:n_modes]
    active_modes = findall(x -> x >= wave_threshold, psi_w)
    n_min_eff = isempty(active_modes) ? -1 : minimum(active_modes) - 1
    n_max_eff = isempty(active_modes) ? -1 : maximum(active_modes) - 1
    in_window = !isempty(active_modes) && (peak_mode >= n_min_eff) && (peak_mode <= n_max_eff)

    weighted_wave_energy = sum(modal_energy .* psi_w)
    f_w_adaptive = clamp(weighted_wave_energy / (e_total + eps(Float64)), 0.0, 1.0)

    # compression_factor: quantifies modal-state collapse without injecting artificial wave signal.
    # Approaches 1 when D_eff ≈ 1 (single dominant mode), decays to 0 for fully distributed states.
    # The classifier can combine f_w_adaptive and compression_factor without either being modified.
    compression_factor = exp(-(d_eff - 1.0))

    return (f_w_adaptive, peak_mode, n_min_eff, in_window, "AdaptiveWaveFractionOK", compression_factor)
end

end