module CasesIngestion

using NCDatasets, LinearAlgebra, Statistics, UnifiedManifold

export ingest_and_project_slice!, project_with_svd_truncation

function highest_active_mode(c::AbstractVector{<:AbstractFloat}; rel_threshold::Float64=1e-3)
    max_amp = maximum(abs.(c))
    if max_amp <= eps(Float64)
        return 0
    end
    active = findall(x -> abs(x) >= rel_threshold * max_amp, c)
    return isempty(active) ? 0 : (maximum(active) - 1)
end

"""
    ingest_and_project_slice!(nc_path, t_idx, ws)

Parses a 7-point vertical tower snapshot from the NetCDF data stream, validates
sensor bounds, and maps the observed fields into a stable manifold projection layer.
"""
function ingest_and_project_slice!(nc_path::String, t_idx::Int, ws)
    Dataset(nc_path, "r") do ds
        # Validated instrumentation heights from EOL tower layout
        z_obs = [1.5, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 55.0]
        M_obs = length(z_obs)
        
        θ_slice = zeros(Float64, M_obs)
        u_slice = zeros(Float64, M_obs)
        
        try
            for i in 1:M_obs
                h = z_obs[i]
                h_str = h == 1.5 ? "1_5m" : string(Int(h)) * "m"

                t_key = haskey(ds, "tc_" * h_str) ? "tc_" * h_str : (haskey(ds, "T_" * h_str) ? "T_" * h_str : error("Missing T"))
                u_key = haskey(ds, "u_" * h_str) ? "u_" * h_str : (haskey(ds, "U_" * h_str) ? "U_" * h_str : error("Missing U"))

                t_var = ds[t_key]
                u_var = ds[u_key]

                if ndims(t_var) == 1
                    raw_θ = t_var[t_idx]
                    raw_u = u_var[t_idx]
                    if ismissing(raw_θ) || ismissing(raw_u)
                        error("Dropout at $h_str")
                    end
                    θ_slice[i] = Float64(raw_θ)
                    u_slice[i] = Float64(raw_u)
                else
                    θ_col = t_var[:, t_idx]
                    u_col = u_var[:, t_idx]

                    valid_θ = filter(x -> !ismissing(x) && x != -1037.0 && !isnan(Float64(x)), θ_col)
                    valid_u = filter(x -> !ismissing(x) && x != -1037.0 && !isnan(Float64(x)), u_col)

                    if isempty(valid_θ) || isempty(valid_u)
                        error("Dropout at $h_str")
                    end

                    θ_slice[i] = mean(Float64.(valid_θ))
                    u_slice[i] = mean(Float64.(valid_u))
                end
            end
        catch e
            return nothing, nothing, "Sensor Dropout Detected"
        end

        if any(isnan, θ_slice) || any(isnan, u_slice)
            return nothing, nothing, "NaN Anomaly"
        end
        
        if any(x -> (x < -50.0 || x > 380.0), θ_slice) || any(x -> (abs(x) > 100.0), u_slice)
            return nothing, nothing, "Data Out of Bounds"
        end

        # Synthesize design matrix H (Size: 7 x 33)
        N_poly = ws.N
        H = zeros(Float64, M_obs, N_poly + 1)
        xi_obs = physical_to_computational(ws, z_obs)

        for i in 1:M_obs
            for n in 0:N_poly
                H[i, n+1] = cos(n * acos(xi_obs[i]))
            end
        end
        
        # Deconstruct modes utilizing a fixed precision barrier
        c_theta, rank_θ, κ_θ = project_with_svd_truncation(H, θ_slice)
        c_u,     rank_u, κ_u = project_with_svd_truncation(H, u_slice)

        n_active = highest_active_mode(c_theta)
        status = "Rank=$(rank_θ), Cond=$(round(κ_θ, sigdigits=2)), ActiveMode=$(n_active)"
        return c_theta, c_u, status
    end
end

"""
    project_with_svd_truncation(H, A_obs)

Performs a stabilized inversion using column-normalized least squares with
ridge and high-mode smoothness regularization to reduce artificial low-order
mode locking when observations are sparse.
"""
function project_with_svd_truncation(H::Matrix{T}, A_obs::Vector{T}) where {T<:AbstractFloat}
    H_scaled = copy(H)

    # 1. Coordinate normalization across high-order columns
    col_norms = [norm(H_scaled[:, j]) for j in 1:size(H_scaled, 2)]
    for j in 1:size(H_scaled, 2)
        if col_norms[j] > 1e-8
            H_scaled[:, j] ./= col_norms[j]
        end
    end

    # 2. Economy SVD Execution
    U, S, Vt = svd(H_scaled)

    # Dynamically select active modes from both singular decay and projected data energy.
    # This avoids a static rank mask when the forcing profile is weak in higher singular vectors.
    tol_s = max(T(1e-4) * S[1], T(eps(Float64)))
    beta = U' * A_obs
    beta_max = maximum(abs.(beta))
    tol_beta = max(T(1e-3) * beta_max, T(eps(Float64)))
    active_idx = findall(i -> (S[i] > tol_s) && (abs(beta[i]) > tol_beta), eachindex(S))
    if isempty(active_idx)
        active_idx = [1]
    end
    rank_eff = length(active_idx)

    # 3. Solve stabilized normal equations with smoothness prior in mode space.
    n_modes = size(H_scaled, 2)
    mode_idx = collect(0:(n_modes - 1))
    smooth_diag = T.(mode_idx .^ 4)
    smooth_diag[1] = zero(T) # Do not penalize mean mode.

    rank_scale = T(length(S) / rank_eff)
    λ_ridge = T(1e-5 * S[1]^2) * rank_scale
    λ_smooth = T(1e-8 * S[1]^2) * rank_scale

    A = (H_scaled' * H_scaled) .+
        λ_ridge .* Matrix{T}(I, n_modes, n_modes) .+
        λ_smooth .* Diagonal(smooth_diag)
    b = H_scaled' * A_obs
    c_scaled = A \ b

    # 4. Physical rescaling back to the original column space.
    c = copy(c_scaled)
    for j in 1:length(c)
        if col_norms[j] > 1e-8
            c[j] /= col_norms[j]
        end
    end

    kappa_eff = S[1] / max(S[active_idx[end]], T(eps(Float64)))
    return c, rank_eff, kappa_eff
end

end # module