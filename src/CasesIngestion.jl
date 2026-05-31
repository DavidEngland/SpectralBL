module CasesIngestion

using NCDatasets, LinearAlgebra, UnifiedManifold

export ingest_and_project_slice!, project_with_svd_truncation

"""
    ingest_and_project_slice!(nc_path, t_idx, ws; tol_frac=1e-3)

Ingests an instantaneous vertical profile from the NCAR ISFS CASES-99 netCDF schema,
resolves the station-flattened height variables, checks for missing entries,
and projects data onto the spectral coefficients via a rank-truncated SVD pseudo-inverse.
"""
function ingest_and_project_slice!(nc_path::String, t_idx::Int, ws; tol_frac=1e-3)
    Dataset(nc_path, "r") do ds
        # Exact main tower thermocouple and sonic instrument height arrays
        z_obs = [2.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0]
        M_obs = length(z_obs)
        
        θ_slice = zeros(Float64, M_obs)
        u_slice = zeros(Float64, M_obs)
        
        # Robustly parse station-flattened NCAR keys
        try
            for i in 1:M_obs
                h_str = string(Int(z_obs[i]))

                # Check for standard names or fallback to central-site naming variants
                t_key = haskey(ds, "tc_" * h_str * "m") ? "tc_" * h_str * "m" : "T_" * h_str * "m"
                u_key = haskey(ds, "u_" * h_str * "m") ? "u_" * h_str * "m" : "u_cs_" * h_str * "m"

                θ_slice[i] = ds[t_key][t_idx]
                u_slice[i] = ds[u_key][t_idx]
            end
        catch e
            return nothing, nothing, "Missing Variable/Sensor Drop"
        end

        # Defensive screening for missing data flags or unhandled NaNs
        if any(isnan, θ_slice) || any(isnan, u_slice) || any(x -> x < -500.0 || x > 5000.0, θ_slice)
            return nothing, nothing, "Data Dropout Detected"
        end
        
        # Build global structural H-matrix from the unified coordinate mapping engine
        N_poly = ws.N
        H = zeros(Float64, M_obs, N_poly + 1)
        xi_obs = physical_to_computational(ws, z_obs)

        for i in 1:M_obs
            for n in 0:N_poly
                H[i, n+1] = cos(n * acos(xi_obs[i]))
            end
        end
        
        # Unbiased singular-space projection
        c_theta, rank_θ, κ_θ = project_with_svd_truncation(H, θ_slice, tol_frac=tol_frac)
        c_u,     rank_u, κ_u = project_with_svd_truncation(H, u_slice, tol_frac=tol_frac)
        
        # Return a rich status string for runtime log aggregation
        status = "Rank=$(rank_θ), Cond=$(round(κ_θ, sigdigits=2))"

        return c_theta, c_u, status
    end
end

"""
    project_with_svd_truncation(H, A_obs; tol_frac=1e-3)

Solves the underdetermined or rank-deficient least-squares projection problem via
Singular Value Decomposition. Truncates all spectral modes whose singular values
fall below the structural data tolerance threshold to completely eliminate ghost modes.
"""
function project_with_svd_truncation(H::Matrix{T}, A_obs::Vector{T}; tol_frac=1e-3) where {T<:AbstractFloat}
    U, S, Vt = svd(H)
    tol = tol_frac * S[1]
    rank_eff = count(>(tol), S)

    # Invert only the resolved singular subspace
    S_inv = [S[i] > tol ? 1.0 / S[i] : 0.0 for i in 1:length(S)]

    # Exact projection: c = V * S_inv * U' * A_obs
    c = zeros(T, size(H, 2))
    for i in 1:rank_eff
        scalar_proj = dot(U[:, i], A_obs) * S_inv[i]
        c .+= scalar_proj .* Vt[i, :]
    end

    kappa_eff = S[1] / S[rank_eff]
    return c, rank_eff, kappa_eff
end

end # module