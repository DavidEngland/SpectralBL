module CasesIngestion

using NCDatasets, LinearAlgebra, UnifiedManifold

export ingest_and_project_slice!, project_with_svd_truncation

"""
    ingest_and_project_slice!(nc_path, t_idx, ws)

Parses a 7-point vertical tower snapshot from the NetCDF data stream, validates
sensor bounds, and maps the observed fields into a stable manifold projection layer.
"""
function ingest_and_project_slice!(nc_path::String, t_idx::Int, ws)
    Dataset(nc_path, "r") do ds
        # Validated instrumentation heights from EOL tower layout
        z_obs = [1.5, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0]
        M_obs = length(z_obs)
        
        θ_slice = zeros(Float64, M_obs)
        u_slice = zeros(Float64, M_obs)
        
        try
            for i in 1:M_obs
                h = z_obs[i]
                h_str = h == 1.5 ? "1_5m" : string(Int(h)) * "m"

                t_key = haskey(ds, "tc_" * h_str) ? "tc_" * h_str : (haskey(ds, "T_" * h_str) ? "T_" * h_str : error("Missing T"))
                u_key = haskey(ds, "u_" * h_str) ? "u_" * h_str : (haskey(ds, "U_" * h_str) ? "U_" * h_str : error("Missing U"))

                raw_θ = ds[t_key][t_idx]
                raw_u = ds[u_key][t_idx]

                if ismissing(raw_θ) || ismissing(raw_u)
                    error("Dropout at $h_str")
                end

                θ_slice[i] = Float64(raw_θ)
                u_slice[i] = Float64(raw_u)
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

        status = "Rank=$(rank_θ), Cond=$(round(κ_θ, sigdigits=2))"
        return c_theta, c_u, status
    end
end

"""
    project_with_svd_truncation(H, A_obs)

Decomposes and scales the linear design system via economy Singular Value Decomposition.
Implements a hardcoded floor to utilize all available physical degrees of freedom.
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

    # Hard-coded floor forces the projection to unlock up to the 7-level structural limit
    machine_floor = 1e-12 * S[1]
    rank_eff = count(>(machine_floor), S)

    # 3. Projection to singular space
    beta = U' * A_obs

    # 4. Inverse scaling singular coefficients
    scaled_beta = zeros(T, length(S))
    for i in 1:rank_eff
        scaled_beta[i] = beta[i] / S[i]
    end

    # 5. Mapping back to full polynomial array
    c = zeros(T, size(Vt, 2))
    for i in 1:rank_eff
        v_row = vec(Vt[i, :])
        c .+= scaled_beta[i] .* v_row
    end

    # 6. Physical rescaling
    for j in 1:length(c)
        if col_norms[j] > 1e-8
            c[j] /= col_norms[j]
        end
    end

    kappa_eff = S[1] / S[rank_eff]
    return c, rank_eff, kappa_eff
end

end # module