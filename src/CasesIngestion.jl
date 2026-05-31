module CasesIngestion

using NCDatasets, LinearAlgebra, UnifiedManifold

export ingest_and_project_slice!, project_with_svd_truncation

"""
    ingest_and_project_slice!(nc_path, t_idx, ws; tol_frac=1e-3)

Ingests an instantaneous vertical profile from CASES-99 netCDF, handles orientation,
builds the global H-matrix using UnifiedManifold's coordinate mapper, and projects the data
onto the spectral coefficients using rank-decoupled SVD truncation.
"""
function ingest_and_project_slice!(nc_path::String, t_idx::Int, ws; tol_frac=1e-3)
    Dataset(nc_path, "r") do ds
        z_obs = Array(ds["height"])[:]
        θ_raw = ds["theta"][:]
        u_raw = ds["u"][:]
        
        # 1. Orientation-Agnostic Array Slicing
        θ_slice = size(θ_raw, 1) == length(z_obs) ? θ_raw[:, t_idx] : θ_raw[t_idx, :]
        u_slice = size(u_raw, 1) == length(z_obs) ? u_raw[:, t_idx] : u_raw[t_idx, :]
        
        # 2. Quality Control & Dropout Isolation
        if any(isnan, θ_slice) || any(isnan, u_slice) || any(x -> x < -500.0 || x > 5000.0, θ_slice)
            return nothing, nothing, "Data Dropout Detected"
        end
        
        # 3. Build Global H-Matrix Using Unified Coordinate Map
        M_obs = length(z_obs)
        N_poly = ws.N
        H = zeros(Float64, M_obs, N_poly + 1)
        
        # Call the single-source-of-truth coordinate mapper from UnifiedManifold
        xi_obs = physical_to_computational(ws, z_obs)

        for i in 1:M_obs
            for n in 0:N_poly
                H[i, n+1] = cos(n * acos(xi_obs[i]))
            end
        end
        
        # 4. Rank-Decoupled Spectral Projection via Truncated SVD
        c_theta, rank_θ, κ_θ = project_with_svd_truncation(H, θ_slice, tol_frac=tol_frac)
        c_u,     rank_u, κ_u = project_with_svd_truncation(H, u_slice, tol_frac=tol_frac)
        
        # Format a dynamic performance status key for monitoring via DataFrame output
        status = "Rank=$(rank_θ), Cond=$(round(κ_θ, sigdigits=3))"
        
        return c_theta, c_u, status
    end
end

"""
    project_with_svd_truncation(H, A_obs; tol_frac=1e-3)

Solves the rank-deficient least-squares projection problem via Singular Value Decomposition.
Truncates modes whose singular values fall below the structural data tolerance threshold.
"""
function project_with_svd_truncation(H::Matrix{T}, A_obs::Vector{T}; tol_frac=1e-3) where {T<:AbstractFloat}
    U, S, Vt = svd(H)
    tol = tol_frac * S[1]
    rank_eff = count(>(tol), S)

    # Solve strictly over the resolved singular sub-space
    S_inv = [S[i] > tol ? 1.0 / S[i] : 0.0 for i in 1:length(S)]

    # Unbiased reconstruction: c = V * S_inv * U' * A_obs
    c = zeros(T, size(H, 2))
    for i in 1:rank_eff
        scalar_proj = dot(U[:, i], A_obs) * S_inv[i]
        c .+= scalar_proj .* Vt[i, :]
    end

    kappa_eff = S[1] / S[rank_eff]
    return c, rank_eff, kappa_eff
end

end # module