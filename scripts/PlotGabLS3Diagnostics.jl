# scripts/PlotGabLS3Diagnostics.jl

using CairoMakie
using CSV
using DataFrames
using NCDatasets
using LinearAlgebra
using Statistics

const DEFAULT_NC = joinpath("data", "gabs3", "gabls3_scm_cabauw_obs_v33.nc")
const DEFAULT_TRAJECTORY = joinpath("data", "trajectory_gabls3.csv")
const DEFAULT_PREDICTIONS = joinpath("reports", "gabs3", "metrics", "predictions.csv")
const DEFAULT_PLOT_DIR = joinpath("reports", "gabs3", "plots")

function _as_zt(a::AbstractArray, n_t::Int)
    if ndims(a) == 1
        v = [x isa Missing ? NaN : Float64(x) for x in a]
        return reshape(v, :, 1)
    end

    if ndims(a) != 2
        error("Expected a 1D/2D array, got ndims=$(ndims(a)).")
    end

    m = map(a) do x
        x isa Missing ? NaN : Float64(x)
    end
    if size(m, 2) == n_t
        return m
    elseif size(m, 1) == n_t
        return permutedims(m)
    elseif size(m, 2) == 1
        return repeat(m, 1, n_t)
    elseif size(m, 2) > n_t
        return m[:, 1:n_t]
    elseif size(m, 1) > n_t
        return permutedims(m[1:n_t, :])
    end

    error("Cannot orient matrix with shape $(size(m)) to (z, time=$n_t).")
end

function _safe_finite(v)
    out = Float64[]
    for x in v
        y = try
            Float64(x)
        catch
            NaN
        end
        isfinite(y) && push!(out, y)
    end
    return out
end

function _window_spans(df::DataFrame)
    spans = Dict{String, Tuple{Float64, Float64}}()
    for tag in ["HOUR8", "HOUR9"]
        g = df[df.WindowTag .== tag, :]
        nrow(g) == 0 && continue
        tvals = _safe_finite(g.SourceTimeHours)
        isempty(tvals) && continue
        spans[tag] = (minimum(tvals), maximum(tvals))
    end
    return spans
end

function _plot_time_height!(ax, t::Vector{Float64}, z::Matrix{Float64}, u::Matrix{Float64})
    nz, nt = size(u)
    size(z) == (nz, nt) || error("z and u must share shape (z, t).")

    t2d = repeat(reshape(t, 1, :), nz, 1)

    if isdefined(CairoMakie, :pcolormesh!)
        return CairoMakie.pcolormesh!(ax, t2d, z, u, colormap=:viridis)
    end

    z_mean = [begin
        vals = z[i, isfinite.(z[i, :])]
        isempty(vals) ? NaN : mean(vals)
    end for i in axes(z, 1)]
    row_mask = isfinite.(z_mean)

    zf = z_mean[row_mask]
    uf = u[row_mask, :]
    return contourf!(ax, t, zf, transpose(uf), levels=24, colormap=:viridis)
end

function _chebyshev_coefficients(z::Vector{Float64}, u::Vector{Float64}; nmax::Int=32)
    length(z) == length(u) || error("z and u length mismatch.")

    # NetCDF slices can contain NaNs/missings at some levels; filter to valid pairs.
    mask = isfinite.(z) .& isfinite.(u)
    zf = z[mask]
    uf = u[mask]
    length(zf) >= 4 || error("Need at least 4 finite points to estimate spectrum.")

    # Sort by height and collapse duplicate levels to keep the least-squares system stable.
    p = sortperm(zf)
    z_sorted = zf[p]
    u_sorted = uf[p]

    z_unique = Float64[]
    u_unique = Float64[]
    i = 1
    while i <= length(z_sorted)
        j = i
        while j < length(z_sorted) && isapprox(z_sorted[j + 1], z_sorted[i]; atol=1e-9, rtol=0.0)
            j += 1
        end
        push!(z_unique, z_sorted[i])
        push!(u_unique, mean(u_sorted[i:j]))
        i = j + 1
    end

    n = length(z_unique)
    n >= 4 || error("Need at least 4 unique finite levels to estimate spectrum.")

    zmin, zmax = extrema(z_unique)
    if zmax <= zmin
        error("Degenerate vertical coordinate for spectrum reconstruction.")
    end

    x = @. 2.0 * (z_unique - zmin) / (zmax - zmin) - 1.0
    ncoef = min(nmax + 1, n)

    T = zeros(Float64, n, ncoef)
    T[:, 1] .= 1.0
    if ncoef >= 2
        T[:, 2] .= x
    end
    for k in 3:ncoef
        @inbounds T[:, k] .= 2.0 .* x .* T[:, k - 1] .- T[:, k - 2]
    end

    c = T \ u_unique
    return c
end

function plot_gabls3_regimes(nc_path::String, trajectory_path::String, save_path::String)
    df_traj = CSV.read(trajectory_path, DataFrame)
    n_t = nrow(df_traj)
    n_t > 0 || error("Trajectory is empty: $trajectory_path")

    time_hours = Float64.(df_traj.SourceTimeHours)
    z_mesh, u_wind = Dataset(nc_path, "r") do ds
        z_raw = Array(ds["zf"])
        u_raw = Array(ds["u"])
        z = _as_zt(z_raw, n_t)
        u = _as_zt(u_raw, n_t)
        nt = min(size(z, 2), size(u, 2), n_t)
        return z[:, 1:nt], u[:, 1:nt]
    end

    t = time_hours[1:size(u_wind, 2)]
    set_theme!(theme_latexfonts())
    fig = Figure(size=(1000, 750), font="CMU Serif")

    z_vals = _safe_finite(vec(z_mesh))
    z_hi = isempty(z_vals) ? 800.0 : min(maximum(z_vals), 1000.0)

    ax1 = Axis(
        fig[1, 1],
        ylabel="Height z (m)",
        title="GABLS3 Jet Evolution and Effective-Dimension Tracking",
        limits=(nothing, nothing, 0, z_hi),
    )

    mesh_plot = _plot_time_height!(ax1, t, z_mesh, u_wind)
    Colorbar(fig[1, 2], mesh_plot, label="U wind (m s^-1)")

    d_vals = _safe_finite(df_traj.D_eff)
    d_lo, d_hi = if isempty(d_vals)
        (0.0, 1.0)
    else
        (minimum(d_vals), maximum(d_vals))
    end
    pad = max(1e-4, 0.05 * max(abs(d_lo), abs(d_hi), 1.0))

    ax2 = Axis(
        fig[2, 1],
        xlabel="Time (hours since 2006-07-01 00:00:00)",
        ylabel="Effective Dimension D_eff",
        limits=(extrema(t), (d_lo - pad, d_hi + pad)),
    )
    lines!(ax2, time_hours, Float64.(df_traj.D_eff), color=:black, linewidth=2.5)

    spans = _window_spans(df_traj)
    if haskey(spans, "HOUR8")
        s, e = spans["HOUR8"]
        vspan!(ax2, s, e, color=(:crimson, 0.15))
        text!(ax2, (s + e) / 2, d_hi + 0.6pad, text="HOUR8", align=(:center, :top), color=:crimson)
    end
    if haskey(spans, "HOUR9")
        s, e = spans["HOUR9"]
        vspan!(ax2, s, e, color=(:darkorange, 0.15))
        text!(ax2, (s + e) / 2, d_hi + 0.6pad, text="HOUR9", align=(:center, :top), color=:darkorange)
    end

    hidexdecorations!(ax1, grid=false, ticks=false)
    rowsize!(fig.layout, 1, Relative(0.6))

    mkpath(dirname(save_path))
    save(save_path, fig, px_per_unit=2)
    println("Regime contour plot saved to: $save_path")
end

function plot_profile_spectral_snapshot(nc_path::String, trajectory_path::String, t_idx::Int, save_path::String)
    df = CSV.read(trajectory_path, DataFrame)
    row = df[df.TimeIdx .== t_idx, :]
    nrow(row) == 0 && error("Time index $t_idx not found in trajectory dataset.")
    r = first(eachrow(row))

    z_prof, u_prof = Dataset(nc_path, "r") do ds
        n_t = nrow(df)
        z = _as_zt(Array(ds["zf"]), n_t)
        u = _as_zt(Array(ds["u"]), n_t)
        idx = clamp(t_idx, 1, min(size(z, 2), size(u, 2)))
        return z[:, idx], u[:, idx]
    end

    coeffs = _chebyshev_coefficients(z_prof, u_prof)
    cabs = abs.(coeffs) .+ eps(Float64)
    nmode = collect(0:length(coeffs)-1)
    spec_mask = isfinite.(cabs) .& (cabs .> 0.0)
    nmode_plot = nmode[spec_mask]
    spectrum_plot = log10.(cabs[spec_mask])

    fig = Figure(size=(1000, 450), font="CMU Serif")

    z_vals = _safe_finite(z_prof)
    z_hi = isempty(z_vals) ? 800.0 : min(maximum(z_vals), 1000.0)

    ax1 = Axis(
        fig[1, 1],
        xlabel="Zonal wind U (m s^-1)",
        ylabel="Height z (m)",
        title="Profile State: TimeIdx=$(t_idx) ($(r.WindowTag))",
        limits=(nothing, (0, z_hi)),
    )
    prof_mask = isfinite.(z_prof) .& isfinite.(u_prof)
    lines!(ax1, u_prof[prof_mask], z_prof[prof_mask], color=:navy, linewidth=2)
    scatter!(ax1, u_prof[prof_mask], z_prof[prof_mask], color=:navy, markersize=5)

    ax2 = Axis(
        fig[1, 2],
        xlabel="Chebyshev mode n",
        ylabel="log10(|c_n|)",
        title="Chebyshev Spectrum (reconstructed from U(z))",
    )
    if !isempty(nmode_plot)
        lines!(ax2, nmode_plot, spectrum_plot, color=:crimson, linewidth=2)
        scatter!(ax2, nmode_plot, spectrum_plot, color=:crimson, marker=:rect, markersize=7)
    else
        text!(ax2, 0.5, 0.5, text="No finite spectral coefficients at this snapshot", align=(:center, :center), space=:relative, color=:gray40)
    end

    mkpath(dirname(save_path))
    save(save_path, fig, px_per_unit=2)
    println("Profile/spectral snapshot saved to: $save_path")
end

function plot_model_comparison(preds_path::String, save_path::String)
    df = CSV.read(preds_path, DataFrame)

    hasproperty(df, :HeatFluxTruth) || error("Missing HeatFluxTruth in predictions CSV.")
    hasproperty(df, :HeatFluxTruth_Pred) || error("Missing HeatFluxTruth_Pred in predictions CSV.")
    hasproperty(df, :HeatFluxTruth_RegimePred) || error("Missing HeatFluxTruth_RegimePred in predictions CSV.")

    fig = Figure(size=(1000, 500), font="CMU Serif")
    colors = Dict("OTHER" => :dimgray, "HOUR8" => :crimson, "HOUR9" => :darkorange)

    obs = Float64.(df.HeatFluxTruth)
    pg = Float64.(df.HeatFluxTruth_Pred)
    pr = Float64.(df.HeatFluxTruth_RegimePred)
    finite_all = isfinite.(obs) .& (isfinite.(pg) .| isfinite.(pr))
    vals = _safe_finite(vcat(obs[finite_all], pg[isfinite.(pg)], pr[isfinite.(pr)]))
    xlo, xhi = isempty(vals) ? (-0.1, 0.1) : extrema(vals)

    ax1 = Axis(fig[1, 1], xlabel="Observed HeatFluxTruth", ylabel="Predicted", title="Global Weighted Linear Baseline", limits=((xlo, xhi), (xlo, xhi)))
    ax2 = Axis(fig[1, 2], xlabel="Observed HeatFluxTruth", title="Regime-Conditioned Sub-Model", limits=((xlo, xhi), (xlo, xhi)))

    lines!(ax1, [xlo, xhi], [xlo, xhi], color=:black, linestyle=:dash, alpha=0.6)
    lines!(ax2, [xlo, xhi], [xlo, xhi], color=:black, linestyle=:dash, alpha=0.6)

    for tag in ["OTHER", "HOUR8", "HOUR9"]
        sub = df[df.WindowTag .== tag, :]
        nrow(sub) == 0 && continue

        truth = Float64.(sub.HeatFluxTruth)
        pred_g = Float64.(sub.HeatFluxTruth_Pred)
        pred_r = Float64.(sub.HeatFluxTruth_RegimePred)

        mg = isfinite.(truth) .& isfinite.(pred_g)
        mr = isfinite.(truth) .& isfinite.(pred_r)

        any(mg) && scatter!(ax1, truth[mg], pred_g[mg], color=colors[tag], label=tag, markersize=8, alpha=0.75)
        any(mr) && scatter!(ax2, truth[mr], pred_r[mr], color=colors[tag], markersize=8, alpha=0.75)
    end

    axislegend(ax1, position=:rb)
    hideydecorations!(ax2, grid=false, ticks=false)

    mkpath(dirname(save_path))
    save(save_path, fig, px_per_unit=2)
    println("Model comparison plot saved to: $save_path")
end

function run_gabls3_plot_suite(
    ;
    nc_path::String=DEFAULT_NC,
    trajectory_path::String=DEFAULT_TRAJECTORY,
    predictions_path::String=DEFAULT_PREDICTIONS,
    plot_dir::String=DEFAULT_PLOT_DIR,
    snapshot_timeidx::Union{Nothing, Int}=nothing,
)
    isfile(nc_path) || error("Missing NetCDF file: $nc_path")
    isfile(trajectory_path) || error("Missing trajectory CSV: $trajectory_path")
    isfile(predictions_path) || error("Missing predictions CSV: $predictions_path")

    mkpath(plot_dir)

    df = CSV.read(trajectory_path, DataFrame)
    idx = if snapshot_timeidx === nothing
        cand = df[df.WindowTag .== "HOUR8", :]
        nrow(cand) > 0 ? Int(cand.TimeIdx[1]) : Int(df.TimeIdx[clamp(Int(round(nrow(df) / 2)), 1, nrow(df))])
    else
        snapshot_timeidx
    end

    regimes_path = joinpath(plot_dir, "llj_manifold_regimes.png")
    profile_path = joinpath(plot_dir, "profile_spectral_snapshot.png")
    model_path = joinpath(plot_dir, "model_flux_comparison.png")

    plot_gabls3_regimes(nc_path, trajectory_path, regimes_path)
    plot_profile_spectral_snapshot(nc_path, trajectory_path, idx, profile_path)
    plot_model_comparison(predictions_path, model_path)

    println("GABLS3 diagnostic plot suite complete.")
    println("  - $regimes_path")
    println("  - $profile_path")
    println("  - $model_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_gabls3_plot_suite()
end
