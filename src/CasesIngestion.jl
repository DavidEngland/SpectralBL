module CasesIngestion

using NCDatasets, LinearAlgebra

export ingest_and_project_slice!

function ingest_and_project_slice!(nc_path::String, t_idx::Int, ws; lambda_tikhonov = 1e-6)
    Dataset(nc_path, "r") do ds
        z_obs = Array(ds["height"])[:]
        θ_raw = ds["theta"][:]
        u_raw = ds["u"][:]
        
        θ_slice = size(θ_raw, 1) == length(z_obs) ? θ_raw[:, t_idx] : θ_raw[t_idx, :]
        u_slice = size(u_raw, 1) == length(z_obs) ? u_raw[:, t_idx] : u_raw[t_idx, :]
        
        if any(isnan, θ_slice) || any(isnan, u_slice) || any(x -> x < -500.0 || x > 5000.0, θ_slice)
            return nothing, nothing, "Data Dropout Detected"
        end
        
        M_obs = length(z_obs)
        N_poly = ws.N
        H = zeros(Float64, M_obs, N_poly + 1)
        
        sigma = (ws.z_atm[end] - ws.z_atm[1]) * 0.15 / 2.0
        for i in 1:M_obs
            xi_val = (z_obs[i] - ws.z_atm[1]) * (1.0 + 0.15) / (sigma + z_obs[i] - ws.z_atm[1]) - 1.0
            xi_val = clamp(xi_val, -1.0, 1.0)
            for n in 0:N_poly
                H[i, n+1] = cos(n * acos(xi_val))
            end
        end
        
        A_mat = H' * H
        κ = cond(A_mat)
        status = "Stable"
        if κ > 1e11
            A_mat += lambda_tikhonov * I
            status = "Regularized"
        end
        
        c_theta = A_mat \ (H' * θ_slice)
        c_u     = A_mat \ (H' * u_slice)
        
        return c_theta, c_u, status
    end
end

end
