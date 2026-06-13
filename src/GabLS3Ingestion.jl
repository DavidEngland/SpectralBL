# src/GabLS3Ingestion.jl
module GabLS3Ingestion

using NCDatasets
using Statistics
using Dates
using UnifiedManifold
using CasesIngestion: project_with_svd_truncation

export GABLS3VariableMap,
       default_variable_map,
       infer_hour_window,
       ingest_and_project_gabls3_slice!

struct GABLS3VariableMap
    u::String
    v::String
    theta::String
    z_mean::String
    heat_flux::String
    mom_u_flux::String
    mom_v_flux::String
    z_flux::String
    time::String
end

function default_variable_map()
    return GABLS3VariableMap(
        "u",
        "v",
        "th",
        "zf",
        "wt",
        "uw",
        "vw",
        "zh",
        "time",
    )
end

function infer_hour_window(time_hours::Real)
    # GABLS3 simulation starts at 12:00 UTC on 2006-07-01 (epoch hour 12).
    # "Hour 8" and "Hour 9" refer to elapsed hours 7-8 and 8-9 post-simulation-start,
    # corresponding to epoch hours 19-20 and 20-21 respectively.
    sim_hour = time_hours - 12.0
    if 7.0 <= sim_hour < 8.0
        return "HOUR8"
    elseif 8.0 <= sim_hour < 9.0
        return "HOUR9"
    end
    return "OTHER"
end

function _finite_pair_filter(z::Vector{Float64}, x::Vector{Float64})
    @assert length(z) == length(x)
    keep = map(i -> isfinite(z[i]) && isfinite(x[i]), eachindex(z))
    return z[keep], x[keep]
end

function _collapse_duplicate_heights(z::Vector{Float64}, x::Vector{Float64})
    @assert length(z) == length(x)
    groups = Dict{Float64, Vector{Float64}}()
    for i in eachindex(z)
        push!(get!(groups, z[i], Float64[]), x[i])
    end

    z_u = sort(collect(keys(groups)))
    x_u = [mean(groups[zi]) for zi in z_u]
    return z_u, x_u
end

function _extract_profile(ds, value_var::String, z_var::String, t_idx::Int)
    v = ds[value_var]
    z = ds[z_var]

    vals = Float64.(coalesce.(v[:, t_idx], NaN))
    zvals = if ndims(z) == 2
        Float64.(coalesce.(z[:, t_idx], NaN))
    else
        Float64.(coalesce.(z[:], NaN))
    end

    zf, xf = _finite_pair_filter(vec(zvals), vec(vals))
    zc, xc = _collapse_duplicate_heights(zf, xf)
    return zc, xc
end

function _project_profile(ws, z_obs::Vector{Float64}, profile::Vector{Float64})
    N_poly = ws.N
    M_obs = length(z_obs)

    xi_obs = physical_to_computational(ws, z_obs)
    H = zeros(Float64, M_obs, N_poly + 1)

    # Chebyshev 3-term recurrence: T_0=1, T_1=x, T_n = 2x*T_{n-1} - T_{n-2}.
    # Replaces M*N acos/cos transcendental calls with branchless FP arithmetic.
    for i in 1:M_obs
        xi_i = clamp(xi_obs[i], -1.0, 1.0)
        H[i, 1] = 1.0
        if N_poly >= 1
            H[i, 2] = xi_i
        end
        for n in 2:N_poly
            H[i, n + 1] = 2.0 * xi_i * H[i, n] - H[i, n - 1]
        end
    end

    c, rank_eff, cond_eff = project_with_svd_truncation(H, profile)
    return c, rank_eff, cond_eff
end

function _interp_linear(z_src::Vector{Float64}, x_src::Vector{Float64}, z_dst::Vector{Float64})
    @assert length(z_src) == length(x_src)
    if isempty(z_src)
        return fill(NaN, length(z_dst))
    end
    if length(z_src) == 1
        return fill(x_src[1], length(z_dst))
    end

    x_dst = similar(z_dst)
    for (i, zq) in enumerate(z_dst)
        if zq <= z_src[1]
            x_dst[i] = x_src[1]
            continue
        elseif zq >= z_src[end]
            x_dst[i] = x_src[end]
            continue
        end

        j = searchsortedlast(z_src, zq)
        j = clamp(j, 1, length(z_src) - 1)
        z1, z2 = z_src[j], z_src[j + 1]
        x1, x2 = x_src[j], x_src[j + 1]
        w = (zq - z1) / (z2 - z1 + eps(Float64))
        x_dst[i] = (1.0 - w) * x1 + w * x2
    end

    return x_dst
end

"""
    ingest_and_project_gabls3_slice!(nc_path, t_idx, ws; varmap=default_variable_map())

Read one GABLS3 time slice, project Set A mean profiles into spectral coefficients,
and return Set C flux targets with window tags.
"""
function ingest_and_project_gabls3_slice!(
    nc_path::String,
    t_idx::Int,
    ws;
    varmap::GABLS3VariableMap=default_variable_map(),
)
    Dataset(nc_path, "r") do ds
        required = [
            varmap.u,
            varmap.v,
            varmap.theta,
            varmap.z_mean,
            varmap.heat_flux,
            varmap.mom_u_flux,
            varmap.mom_v_flux,
            varmap.z_flux,
            varmap.time,
        ]
        for vname in required
            haskey(ds, vname) || error("Missing required variable: " * vname)
        end

        z_mean, theta_profile = _extract_profile(ds, varmap.theta, varmap.z_mean, t_idx)
        z_u, u_profile = _extract_profile(ds, varmap.u, varmap.z_mean, t_idx)
        z_v, v_profile = _extract_profile(ds, varmap.v, varmap.z_mean, t_idx)

        if length(z_mean) < 3 || length(z_u) < 3 || length(z_v) < 3
            return nothing
        end

        if z_mean != z_u
            u_profile = _interp_linear(z_u, u_profile, z_mean)
        end
        if z_mean != z_v
            v_profile = _interp_linear(z_v, v_profile, z_mean)
        end

        # Restrict profiles to the workspace physical domain [z_min, z_max].
        # Observations above z_max would be clamped to xi=1 in the Chebyshev
        # map, stacking at the boundary and producing artificial high-order
        # gradients. Filtering here eliminates that before projection.
        z_ws_min = minimum(ws.z_atm)
        z_ws_max = maximum(ws.z_atm)
        domain_mask = (z_mean .>= z_ws_min) .& (z_mean .<= z_ws_max)
        if count(domain_mask) < 3
            return nothing
        end
        z_mean   = z_mean[domain_mask]
        theta_profile = theta_profile[domain_mask]
        u_profile     = u_profile[domain_mask]
        v_profile     = v_profile[domain_mask]

        c_theta, rank_theta, cond_theta = _project_profile(ws, z_mean, theta_profile)
        c_u, rank_u, cond_u = _project_profile(ws, z_mean, u_profile)
        c_v, rank_v, cond_v = _project_profile(ws, z_mean, v_profile)

        z_flux, heat_flux = _extract_profile(ds, varmap.heat_flux, varmap.z_flux, t_idx)
        z_uw, mom_u_flux = _extract_profile(ds, varmap.mom_u_flux, varmap.z_flux, t_idx)
        z_vw, mom_v_flux = _extract_profile(ds, varmap.mom_v_flux, varmap.z_flux, t_idx)

        if z_flux != z_uw
            mom_u_flux = _interp_linear(z_uw, mom_u_flux, z_flux)
        end
        if z_flux != z_vw
            mom_v_flux = _interp_linear(z_vw, mom_v_flux, z_flux)
        end

        time_val = ds[varmap.time][t_idx]
        # Anchor to the GABLS3 simulation midnight epoch (2006-07-01 00:00:00 UTC).
        # Using the file's first timestamp as a reference would silently shift
        # HOUR8/HOUR9 tags if the file does not start exactly at t=0.
        gabls3_epoch = DateTime(2006, 7, 1, 0, 0, 0)
        t_hours = if time_val isa DateTime
            Dates.value(time_val - gabls3_epoch) / 3_600_000.0
        else
            Float64(time_val)
        end
        window_tag = infer_hour_window(t_hours)

        status = "RankTheta=$(rank_theta),RankU=$(rank_u),RankV=$(rank_v),CondTheta=$(round(cond_theta, sigdigits=3)),CondU=$(round(cond_u, sigdigits=3)),CondV=$(round(cond_v, sigdigits=3))"

        return (
            c_theta=c_theta,
            c_u=c_u,
            c_v=c_v,
            status=status,
            time_hours=t_hours,
            window_tag=window_tag,
            z_mean=z_mean,
            theta_profile=theta_profile,
            u_profile=u_profile,
            v_profile=v_profile,
            z_flux=z_flux,
            heat_flux=heat_flux,
            mom_u_flux=mom_u_flux,
            mom_v_flux=mom_v_flux,
        )
    end
end

end # module
