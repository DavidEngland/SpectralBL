# scripts/generate_coefficient_trajectory_figure.jl
using CSV
using DataFrames
using NCDatasets
using Plots
using LinearAlgebra
using Statistics
using Printf
using LaTeXStrings

include("../src/Cases99.jl")
using .UnifiedManifold: UnifiedManifoldWorkspace, physical_to_computational

const Z_OBS = [1.5, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 55.0]

function safe_float(x)
    if ismissing(x)
        return nothing
    elseif x isa Number
        v = Float64(x)
        return isfinite(v) ? v : nothing
    end
    sx = strip(string(x))
    if isempty(sx)
        return nothing
    end
    try
        v = parse(Float64, sx)
        return isfinite(v) ? v : nothing
    catch
        return nothing
    end
end

function infer_regime(d_eff::Float64, f_w::Float64)
    if d_eff >= 10.0 && f_w < 0.25
        return 1
    elseif d_eff <= 8.0 && f_w >= 0.30
        return 2
    else
        return 3
    end
end

function project_with_svd_truncation(H::Matrix{Float64}, a_obs::Vector{Float64})
    H_scaled = copy(H)

    col_norms = [norm(H_scaled[:, j]) for j in 1:size(H_scaled, 2)]
    for j in 1:size(H_scaled, 2)
        if col_norms[j] > 1e-8
            H_scaled[:, j] ./= col_norms[j]
        end
    end

    U, S, _ = svd(H_scaled)
    tol_s = max(1e-4 * S[1], eps(Float64))
    beta = U' * a_obs
    beta_max = maximum(abs.(beta))
    tol_beta = max(1e-3 * beta_max, eps(Float64))
    active_idx = findall(i -> (S[i] > tol_s) && (abs(beta[i]) > tol_beta), eachindex(S))
    if isempty(active_idx)
        active_idx = [1]
    end
    rank_eff = length(active_idx)

    n_modes = size(H_scaled, 2)
    mode_idx = collect(0:(n_modes - 1))
    smooth_diag = Float64.(mode_idx .^ 4)
    smooth_diag[1] = 0.0

    rank_scale = length(S) / rank_eff
    lambda_ridge = 1e-5 * S[1]^2 * rank_scale
    lambda_smooth = 1e-8 * S[1]^2 * rank_scale

    A = (H_scaled' * H_scaled) .+
        lambda_ridge .* Matrix{Float64}(I, n_modes, n_modes) .+
        lambda_smooth .* Diagonal(smooth_diag)
    b = H_scaled' * a_obs
    c_scaled = A \ b

    c = copy(c_scaled)
    for j in 1:length(c)
        if col_norms[j] > 1e-8
            c[j] /= col_norms[j]
        end
    end

    kappa_eff = S[1] / max(S[active_idx[end]], eps(Float64))
    return c, rank_eff, kappa_eff
end

function ingest_and_project_slice_local(nc_path::String, t_idx::Int, ws)
    z_obs = Float64[1.5, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 55.0]
    m_obs = length(z_obs)
    theta_slice = zeros(Float64, m_obs)
    u_slice = zeros(Float64, m_obs)

    Dataset(nc_path, "r") do ds
        for i in eachindex(z_obs)
            h = z_obs[i]
            h_str = h == 1.5 ? "1_5m" : string(Int(h)) * "m"

            t_key = haskey(ds, "tc_" * h_str) ? "tc_" * h_str : (haskey(ds, "T_" * h_str) ? "T_" * h_str : nothing)
            u_key = haskey(ds, "u_" * h_str) ? "u_" * h_str : (haskey(ds, "U_" * h_str) ? "U_" * h_str : nothing)
            if isnothing(t_key) || isnothing(u_key)
                return nothing, nothing, "Missing observation channels"
            end

            t_val = pull_scalar(ds[t_key], t_idx)
            u_val = pull_scalar(ds[u_key], t_idx)
            if isnothing(t_val) || isnothing(u_val)
                return nothing, nothing, "Sensor Dropout Detected"
            end

            theta_slice[i] = t_val
            u_slice[i] = u_val
        end
    end

    if any(isnan, theta_slice) || any(isnan, u_slice)
        return nothing, nothing, "NaN Anomaly"
    end

    n_poly = ws.N
    H = zeros(Float64, m_obs, n_poly + 1)
    xi_obs = physical_to_computational(ws, z_obs)
    for i in 1:m_obs
        for n in 0:n_poly
            H[i, n + 1] = cos(n * acos(xi_obs[i]))
        end
    end

    c_theta, rank_theta, kappa_theta = project_with_svd_truncation(H, theta_slice)
    c_u, _, _ = project_with_svd_truncation(H, u_slice)
    status = "Rank=$(rank_theta), Cond=$(round(kappa_theta, sigdigits=2))"
    return c_theta, c_u, status
end

function pull_scalar(var, t_idx::Int)
    if ndims(var) == 1
        if t_idx < 1 || t_idx > length(var)
            return nothing
        end
        raw = var[t_idx]
        if ismissing(raw)
            return nothing
        end
        val = Float64(raw)
        return isfinite(val) ? val : nothing
    end

    n_t = size(var, ndims(var))
    if t_idx < 1 || t_idx > n_t
        return nothing
    end

    slice = var[:, t_idx]
    vals = Float64[]
    for x in slice
        if ismissing(x)
            continue
        end
        xv = Float64(x)
        if xv == -1037.0 || !isfinite(xv)
            continue
        end
        push!(vals, xv)
    end
    return isempty(vals) ? nothing : mean(vals)
end

function fetch_observed_state(nc_path::String, t_idx::Int)
    if !isfile(nc_path)
        return nothing
    end

    Dataset(nc_path, "r") do ds
        theta = zeros(Float64, length(Z_OBS))
        u = zeros(Float64, length(Z_OBS))

        for i in eachindex(Z_OBS)
            h = Z_OBS[i]
            h_str = h == 1.5 ? "1_5m" : string(Int(h)) * "m"

            t_key = haskey(ds, "tc_" * h_str) ? "tc_" * h_str : (haskey(ds, "T_" * h_str) ? "T_" * h_str : nothing)
            u_key = haskey(ds, "u_" * h_str) ? "u_" * h_str : (haskey(ds, "U_" * h_str) ? "U_" * h_str : nothing)
            if isnothing(t_key) || isnothing(u_key)
                return nothing
            end

            t_val = pull_scalar(ds[t_key], t_idx)
            u_val = pull_scalar(ds[u_key], t_idx)
            if isnothing(t_val) || isnothing(u_val)
                return nothing
            end

            theta[i] = t_val
            u[i] = u_val
        end

        return theta, u
    end
end

function choose_rows(df::DataFrame; n_each_regime::Int=45)
    sort!(df, [:FileDate, :TimeIdx])
    selected_rows = DataFrame[]

    for r in 1:3
        sub = filter(row -> row.Regime == r, df)
        if nrow(sub) == 0
            continue
        end

        k = min(n_each_regime, nrow(sub))
        idx = unique(round.(Int, range(1, nrow(sub), length=k)))
        push!(selected_rows, sub[idx, :])
    end

    if isempty(selected_rows)
        return DataFrame()
    end

    out = reduce((a, b) -> vcat(a, b; cols=:union), selected_rows)
    unique!(out, [:FileDate, :TimeIdx])
    sort!(out, [:FileDate, :TimeIdx])
    return out
end

function safe_volume_proxy(X::Matrix{Float64})
    if size(X, 1) < 4
        return NaN
    end
    C = cov(X)
    C .+= 1e-10 * Matrix{Float64}(I, 3, 3)
    return sqrt(abs(det(C)))
end

function run()
    trajectory_path = joinpath("data", "diagnostic_trajectory.csv")
    if !isfile(trajectory_path)
        println("! Skipping attractor figure: missing ", trajectory_path)
        return
    end

    raw = CSV.read(trajectory_path, DataFrame)
    if nrow(raw) == 0
        println("! Skipping attractor figure: trajectory table is empty.")
        return
    end

    required = [:FileDate, :TimeIdx, :D_eff, :F_W]
    missing_cols = [c for c in required if !hasproperty(raw, c)]
    if !isempty(missing_cols)
        println("! Skipping attractor figure: missing columns ", join(string.(missing_cols), ", "))
        return
    end

    cleaned = DataFrame(FileDate=Int[], TimeIdx=Int[], D_eff=Float64[], F_W=Float64[], Regime=Int[])
    for row in eachrow(raw)
        d_eff = safe_float(row.D_eff)
        f_w = safe_float(row.F_W)
        file_date_f = safe_float(row.FileDate)
        time_idx_f = safe_float(row.TimeIdx)
        if isnothing(d_eff) || isnothing(f_w) || isnothing(file_date_f) || isnothing(time_idx_f)
            continue
        end

        file_date = Int(round(file_date_f))
        time_idx = Int(round(time_idx_f))
        reg = infer_regime(d_eff, f_w)
        if hasproperty(raw, :Regime)
            reg_f = safe_float(row.Regime)
            if !isnothing(reg_f)
                reg = Int(round(reg_f))
            end
        end
        push!(cleaned, (file_date, time_idx, d_eff, f_w, clamp(reg, 1, 3)))
    end

    if nrow(cleaned) == 0
        println("! Skipping attractor figure: no finite trajectory rows after cleaning.")
        return
    end

    selected = choose_rows(cleaned; n_each_regime=45)
    if nrow(selected) == 0
        println("! Skipping attractor figure: no sampled rows available.")
        return
    end

    ws = UnifiedManifoldWorkspace(32, 1.5, 55.0, 2.50; n_m=3, n_w=12, delta=1.2)

    p_phys = Float64[]
    p_coeff = Float64[]
    d_eff_vals = Float64[]
    f_w_vals = Float64[]
    reg_vals = Int[]
    seq_vals = Int[]

    seq = 0
    for row in eachrow(selected)
        nc_path = joinpath("data", "ncar_eol_dee0099881", "cases.$(row.FileDate).nc")
        obs = fetch_observed_state(nc_path, row.TimeIdx)
        if isnothing(obs)
            continue
        end

        theta, u = obs
        c_theta, c_u, status = ingest_and_project_slice_local(nc_path, row.TimeIdx, ws)
        if isnothing(c_theta) || isnothing(c_u) || occursin("Dropout", status)
            continue
        end

        # Physical space coordinates from direct tower observations.
        push!(p_phys, theta[1], theta[5], u[3])
        # Coefficient-space coordinates emphasizing resolved low-order structure.
        push!(p_coeff, c_theta[2], c_theta[3], c_u[2])

        push!(d_eff_vals, row.D_eff)
        push!(f_w_vals, row.F_W)
        push!(reg_vals, row.Regime)
        seq += 1
        push!(seq_vals, seq)
    end

    n_valid = length(d_eff_vals)
    if n_valid < 15
        println("! Skipping attractor figure: only ", n_valid, " valid trajectory points after data checks.")
        return
    end

    X_phys = Matrix(reshape(p_phys, 3, :)')
    X_coeff = Matrix(reshape(p_coeff, 3, :)')

    phys_volume = safe_volume_proxy(X_phys)
    coeff_volume = safe_volume_proxy(X_coeff)
    collapse_ratio = isfinite(phys_volume) && phys_volume > 0 ? coeff_volume / phys_volume : NaN

    regime_colors = Dict(1 => :blue, 2 => :red, 3 => :green)
    point_colors = [regime_colors[r] for r in reg_vals]

    p1 = scatter3d(
        X_phys[:, 1], X_phys[:, 2], X_phys[:, 3];
        markercolor=point_colors,
        markerstrokewidth=0,
        markersize=3.8,
        alpha=0.82,
        xlabel=L"\theta(1.5\,\mathrm{m})\,[\mathrm{K}]",
        ylabel=L"\theta(30\,\mathrm{m})\,[\mathrm{K}]",
        zlabel=L"u(10\,\mathrm{m})\,[\mathrm{m}\,\mathrm{s}^{-1}]",
        title="Physical Tower-State Trajectory",
        legend=false,
        left_margin=10Plots.mm,
        bottom_margin=8Plots.mm,
    )
    plot!(p1, X_phys[:, 1], X_phys[:, 2], X_phys[:, 3]; color=:black, linewidth=1.0, alpha=0.40, label=false)

    p2 = scatter3d(
        X_coeff[:, 1], X_coeff[:, 2], X_coeff[:, 3];
        markercolor=point_colors,
        markerstrokewidth=0,
        markersize=3.8,
        alpha=0.82,
        xlabel=L"c_{\theta,1}",
        ylabel=L"c_{\theta,2}",
        zlabel=L"c_{u,1}",
        title="Chebyshev Coefficient-Space Trajectory",
        legend=false,
        left_margin=10Plots.mm,
        bottom_margin=8Plots.mm,
    )
    plot!(p2, X_coeff[:, 1], X_coeff[:, 2], X_coeff[:, 3]; color=:black, linewidth=1.0, alpha=0.40, label=false)

    p3 = scatter(
        seq_vals, d_eff_vals;
        markercolor=point_colors,
        markerstrokewidth=0,
        markersize=3.6,
        alpha=0.86,
        xlabel="Sampled Trajectory Index",
        ylabel=L"D_{\mathrm{eff}}",
        title="Temporal Compression Trace",
        legend=false,
        left_margin=10Plots.mm,
        bottom_margin=8Plots.mm,
    )
    plot!(p3, seq_vals, d_eff_vals; color=:black, linewidth=1.0, alpha=0.35, label=false)

    p4 = scatter(
        d_eff_vals, f_w_vals;
        markercolor=point_colors,
        markerstrokewidth=0,
        markersize=3.6,
        alpha=0.86,
        xlabel=L"D_{\mathrm{eff}}",
        ylabel=L"F_W",
        title="State-Space Regime Envelope",
        legend=false,
        left_margin=10Plots.mm,
        bottom_margin=8Plots.mm,
    )

    ratio_text = isfinite(collapse_ratio) ? @sprintf("Volume ratio = %.3e", collapse_ratio) : "Volume ratio unavailable"
    fig = plot(p1, p2, p3, p4; layout=(2, 2), size=(1500, 1100), plot_title="Physical-to-Spectral Attractor Collapse ($ratio_text)")

    reports_dir = joinpath("reports", "ncar_eol_dee0099881")
    draft_fig_dir = joinpath("drafts", "figures")
    mkpath(reports_dir)
    mkpath(draft_fig_dir)

    out_name = "attractor_collapse_physical_vs_spectral"
    savefig(fig, joinpath(reports_dir, out_name * ".pdf"))
    savefig(fig, joinpath(draft_fig_dir, out_name * ".pdf"))

    println("✓ Attractor collapse figure saved: ", joinpath(reports_dir, out_name * ".pdf"))
    println("✓ Attractor collapse figure saved: ", joinpath(draft_fig_dir, out_name * ".pdf"))
    println("✓ Figure points used: ", n_valid, ", collapse ratio: ", ratio_text)
end

run()
