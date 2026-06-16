# scripts/GabLS3PredictiveBaseline.jl

using CSV
using DataFrames
using Statistics
using LinearAlgebra
using Printf

const DEFAULT_INPUT = joinpath("data", "trajectory_gabls3.csv")
const DEFAULT_MD = joinpath("reports", "gabs3", "gabls3_flux_skill.md")
const DEFAULT_METRICS_DIR = joinpath("reports", "gabs3", "metrics")

function _safe_mean(v::AbstractVector{<:Real})
    vals = filter(isfinite, Float64.(v))
    return isempty(vals) ? NaN : mean(vals)
end

function _to_float(x)
    if x isa Missing
        return NaN
    elseif x isa AbstractFloat || x isa Integer
        return Float64(x)
    elseif x isa AbstractString
        v = tryparse(Float64, x)
        return v === nothing ? NaN : v
    end
    return try
        Float64(x)
    catch
        NaN
    end
end

function _prepare_dataframe(df::DataFrame)
    keep_cols = [
        :WindowTag,
        :D_eff,
        :F_W,
        :chi_N,
        :Ri_g,
        :Ri_b,
        :UShearMean,
        :VShearMean,
        :ShearMagMean,
        :Deff_x_ShearMag,
        :HeatFluxTruth,
        :MomUFluxTruth,
        :MomVFluxTruth,
    ]

    # Lagged inertial features are optional — present only after the pipeline
    # has been re-run with the updated RunGabLS3Pipeline.jl.
    optional_cols = [:DeltaU_Dt, :DeltaV_Dt, :InertialMagDt]
    have_lag = all(c -> hasproperty(df, c), optional_cols)

    all_cols = have_lag ? vcat(keep_cols, optional_cols) : keep_cols
    for c in keep_cols
        hasproperty(df, c) || error("Missing required column: $(c)")
    end

    d = select(df, all_cols)

    numeric_cols = filter(c -> c != :WindowTag, all_cols)
    for c in numeric_cols
        d[!, c] = _to_float.(d[!, c])
    end

    mask = trues(nrow(d))
    for c in keep_cols
        c == :WindowTag && continue
        mask .&= isfinite.(d[!, c])
    end
    # Lag columns use NaN for the first row; mask them out separately so their
    # absence doesn't drop every row when lagged features are present.
    if have_lag
        lag_mask = trues(nrow(d))
        for c in optional_cols
            lag_mask .&= isfinite.(d[!, c])
        end
        d[!, :HasLag] = lag_mask
    end

    d = d[mask, :]
    nrow(d) > 0 || error("No finite rows available after filtering.")

    counts = combine(groupby(d, :WindowTag), nrow => :n_samples)
    count_map = Dict(row.WindowTag => row.n_samples for row in eachrow(counts))
    d[!, :SampleWeight] = [1.0 / count_map[tag] for tag in d.WindowTag]

    return d, counts, have_lag
end

# Column indices for the restricted model used when N < min_samples_full.
# Intercept (1), D_eff (2), ShearMagMean (10) — three degrees of freedom
# that span regime magnitude and spectral dimensionality with minimal parameters.
const REDUCED_COLS = [1, 2, 10]

function _build_design_matrix(df::DataFrame; have_lag::Bool=false)
    n = nrow(df)
    # Column 10 is always ShearMagMean; lag features start at 11 when present.
    ncols = have_lag ? 13 : 10
    X = zeros(Float64, n, ncols)
    X[:, 1] .= 1.0
    X[:, 2] .= df.D_eff
    X[:, 3] .= df.F_W
    X[:, 4] .= df.chi_N
    X[:, 5] .= df.Ri_g
    X[:, 6] .= df.Ri_b
    X[:, 7] .= df.UShearMean
    X[:, 8] .= df.VShearMean
    X[:, 9] .= df.Deff_x_ShearMag
    X[:, 10] .= df.ShearMagMean

    feature_names = [
        "Intercept", "D_eff", "F_W", "chi_N", "Ri_g", "Ri_b",
        "UShearMean", "VShearMean", "Deff_x_ShearMag", "ShearMagMean",
    ]

    if have_lag
        X[:, 11] .= _to_float.(df.DeltaU_Dt)
        X[:, 12] .= _to_float.(df.DeltaV_Dt)
        X[:, 13] .= _to_float.(df.InertialMagDt)
        push!(feature_names, "DeltaU_Dt", "DeltaV_Dt", "InertialMagDt")
    end

    return X, feature_names
end

function _standardise_features(X::Matrix{Float64})
    mu = vec(mean(X, dims=1))
    sigma = vec(std(X, dims=1))
    # Intercept column (index 1, all-ones) must not be centred or scaled.
    # Zero-variance columns (constant features) also left untouched.
    mu[1] = 0.0
    sigma[sigma .== 0.0] .= 1.0
    sigma[1] = 1.0
    Xs = (X .- mu') ./ sigma'
    return Xs, mu, sigma
end

function _fit_weighted_linear(X::Matrix{Float64}, y::Vector{Float64}, w::Vector{Float64}; ridge::Float64=1e-8)
    Xs, mu, sigma = _standardise_features(X)
    W = Diagonal(w)
    p = size(Xs, 2)
    # Normal equations are symmetric by construction; Cholesky exploits that and
    # guarantees a positive-definite solve without the extra overhead of a generic
    # factorisation path.
    A = Hermitian(Xs' * W * Xs + ridge * Matrix{Float64}(I, p, p))
    b = Xs' * W * y
    β_scaled = cholesky(A) \ b
    yhat = Xs * β_scaled
    # Recover coefficients on the original (unstandardised) feature scale so
    # reported values are directly interpretable in physical units.
    β = β_scaled ./ sigma
    β[1] -= dot(mu ./ sigma, β_scaled)
    return β, yhat
end

function _metrics(y::Vector{Float64}, yhat::Vector{Float64}, w::Vector{Float64})
    resid = y .- yhat
    wmse = sum(w .* resid.^2) / (sum(w) + eps(Float64))
    wrmse = sqrt(wmse)
    wmae = sum(w .* abs.(resid)) / (sum(w) + eps(Float64))
    wbias = sum(w .* resid) / (sum(w) + eps(Float64))

    ybar = sum(w .* y) / (sum(w) + eps(Float64))
    sst = sum(w .* (y .- ybar).^2)
    sse = sum(w .* resid.^2)
    wr2 = sst > eps(Float64) ? (1.0 - sse / sst) : NaN

    return (wrmse=wrmse, wmae=wmae, wr2=wr2, wbias=wbias)
end

function _metrics_by_tag(df::DataFrame, pred_col::Symbol, target_col::Symbol)
    out = DataFrame(WindowTag=String[], wrmse=Float64[], wmae=Float64[], wr2=Float64[], wbias=Float64[], n_samples=Int[])
    for g in groupby(df, :WindowTag)
        y = Float64.(g[!, target_col])
        yhat = Float64.(g[!, pred_col])
        w = ones(Float64, length(y))
        m = _metrics(y, yhat, w)
        push!(out, (String(g.WindowTag[1]), m.wrmse, m.wmae, m.wr2, m.wbias, nrow(g)))
    end
    return out
end

"""
Regime-conditioned sub-model: fit a separate weighted linear model per WindowTag.

When the effective sample count for a window is below `min_samples_full` (default 15),
the full feature set is over-determined (N ≤ P) and the ridge regulariser would
simply interpolate through those points. In that case, we fall back to a 3-column
restricted model [Intercept, D_eff, ShearMagMean] that is always well-posed and
still captures the primary magnitude/dimensionality regime contrast.
"""
function _fit_regime_conditioned(
    df::DataFrame,
    X::Matrix{Float64},
    y::Vector{Float64};
    ridge::Float64=1e-8,
    min_samples_full::Int=15,
)
    n = nrow(df)
    yhat = fill(NaN, n)
    betas = Dict{String, Vector{Float64}}()
    reduced_tags = Set{String}()   # tags that used the restricted feature set

    tags = unique(df.WindowTag)
    for tag in tags
        idx = findall(==(tag), df.WindowTag)
        length(idx) < 2 && continue

        Xi = X[idx, :]
        yi = y[idx]
        wi = Float64.(df.SampleWeight[idx])

        # For small windows, commit to the restricted column set *before* the
        # finite-row check so that lag NaN values in out-of-scope columns do
        # not falsely exclude rows that are perfectly valid for the restricted fit.
        use_reduced = length(idx) < min_samples_full
        Xi_check = use_reduced ? Xi[:, REDUCED_COLS] : Xi

        finite_rows = findall(i -> all(isfinite, Xi_check[i, :]) && isfinite(yi[i]), 1:length(idx))
        length(finite_rows) < 2 && continue

        Xi_f = Xi_check[finite_rows, :]
        yi_f = yi[finite_rows]
        wi_f = wi[finite_rows]

        if use_reduced
            push!(reduced_tags, String(tag))
        end

        βi, yhati = _fit_weighted_linear(Xi_f, yi_f, wi_f; ridge=ridge)
        betas[String(tag)] = βi
        yhat[idx[finite_rows]] .= yhati
    end

    return yhat, betas, reduced_tags
end

function _write_markdown_report(md_path::String, counts::DataFrame, target_results::Dict{String, Any}, have_lag::Bool)
    mkpath(dirname(md_path))
    open(md_path, "w") do io
        println(io, "# GABLS3 Flux Predictive Baseline")
        println(io)
        lag_note = have_lag ? " Lagged inertial-oscillation features (ΔU/Δt, ΔV/Δt, InertialMagDt) included." : " Lagged features not yet available — re-run `make gabls3-run` to add them."
        println(io, "This report fits (1) a global weighted linear baseline and (2) a regime-conditioned sub-model (one fit per WindowTag).$lag_note Sample weights are inversely proportional to window sample count.")
        println(io)
        println(io, "## Window Counts")
        println(io)
        println(io, "| WindowTag | n_samples |")
        println(io, "|---|---:|")
        for row in eachrow(counts)
            println(io, "| $(row.WindowTag) | $(row.n_samples) |")
        end

        for (target, res) in target_results
            println(io)
            println(io, "## Target: $(target)")
            println(io)

            # --- Global baseline ---
            println(io, "### Global Weighted Linear Baseline")
            println(io)
            println(io, "- wRMSE: ", @sprintf("%.6f", res.global_metrics.wrmse))
            println(io, "- wMAE: ", @sprintf("%.6f", res.global_metrics.wmae))
            println(io, "- wR2: ", @sprintf("%.6f", res.global_metrics.wr2))
            println(io, "- wBias: ", @sprintf("%.6f", res.global_metrics.wbias))
            println(io)
            println(io, "#### By-Window Metrics (Global Model)")
            println(io)
            println(io, "| WindowTag | wRMSE | wMAE | wR2 | wBias | n_samples |")
            println(io, "|---|---:|---:|---:|---:|---:|")
            for row in eachrow(res.by_tag)
                println(io, "| $(row.WindowTag) | ", @sprintf("%.6f", row.wrmse), " | ", @sprintf("%.6f", row.wmae), " | ", @sprintf("%.6f", row.wr2), " | ", @sprintf("%.6f", row.wbias), " | $(row.n_samples) |")
            end
            println(io)
            println(io, "#### Coefficients (Global Model)")
            println(io)
            println(io, "| Feature | Coefficient |")
            println(io, "|---|---:|")
            for (fname, beta) in zip(res.feature_names, res.beta)
                println(io, "| $(fname) | ", @sprintf("%.8f", beta), " |")
            end

            # --- Regime-conditioned sub-model ---
            println(io)
            println(io, "### Regime-Conditioned Sub-Model")
            println(io)
            println(io, "#### By-Window Metrics (Regime Model)")
            println(io)
            println(io, "| WindowTag | wRMSE | wMAE | wR2 | wBias | n_samples |")
            println(io, "|---|---:|---:|---:|---:|---:|")
            for row in eachrow(res.regime_by_tag)
                println(io, "| $(row.WindowTag) | ", @sprintf("%.6f", row.wrmse), " | ", @sprintf("%.6f", row.wmae), " | ", @sprintf("%.6f", row.wr2), " | ", @sprintf("%.6f", row.wbias), " | $(row.n_samples) |")
            end
            println(io)
            println(io, "#### Per-Tag Coefficients (Regime Model)")
            println(io)
            reduced_feature_names = [res.feature_names[i] for i in REDUCED_COLS]
            for (tag, betas) in sort(collect(res.regime_betas), by=first)
                is_reduced = tag in res.reduced_tags
                suffix = is_reduced ? " *(restricted: N < $(res.min_samples_full), features: $(join(reduced_feature_names, ", ")))*" : ""
                println(io, "**WindowTag = $(tag)**$(suffix)")
                println(io)
                println(io, "| Feature | Coefficient |")
                println(io, "|---|---:|")
                fnames = is_reduced ? reduced_feature_names : res.feature_names
                for (fname, beta) in zip(fnames, betas)
                    println(io, "| $(fname) | ", @sprintf("%.8f", beta), " |")
                end
                println(io)
            end
        end
    end
end

function run_gabls3_predictive_baseline(; input_csv::String=DEFAULT_INPUT, md_path::String=DEFAULT_MD, metrics_dir::String=DEFAULT_METRICS_DIR)
    isfile(input_csv) || error("Missing trajectory file: $input_csv")
    mkpath(metrics_dir)

    raw = CSV.read(input_csv, DataFrame)
    df, counts, have_lag = _prepare_dataframe(raw)

    X, feature_names = _build_design_matrix(df; have_lag=have_lag)
    w = Float64.(df.SampleWeight)

    targets = Dict(
        "HeatFluxTruth" => :HeatFluxTruth,
        "MomUFluxTruth" => :MomUFluxTruth,
        "MomVFluxTruth" => :MomVFluxTruth,
    )

    target_results = Dict{String, Any}()

    for (name, col) in targets
        y = Float64.(df[!, col])

        # When lagged features are present, the first row of each continuous run
        # has NaN lags (no prior state). Build a finite-row mask for the global
        # fit so those rows don't poison the Cholesky standardisation.
        if have_lag
            lag_finite = Float64.(df.HasLag)  # 1.0 where all lags are finite
            global_mask = map(i -> isfinite(y[i]) && Bool(df.HasLag[i]), 1:nrow(df))
        else
            global_mask = map(i -> isfinite(y[i]), 1:nrow(df))
        end
        df_g = df[global_mask, :]
        X_g  = X[global_mask, :]
        y_g  = Float64.(df_g[!, col])
        w_g  = Float64.(df_g.SampleWeight)

        # --- Global weighted linear baseline ---
        beta, yhat_g = _fit_weighted_linear(X_g, y_g, w_g)
        # Fill predictions back into the full-length vector (lag-NaN rows stay NaN).
        yhat_global = fill(NaN, nrow(df))
        yhat_global[global_mask] .= yhat_g
        pred_col = Symbol(name * "_Pred")
        df[!, pred_col] = yhat_global

        global_m = _metrics(y_g, yhat_g, w_g)
        # Score only rows that have a valid global prediction.
        by_tag = _metrics_by_tag(df[isfinite.(yhat_global), :], pred_col, col)
        metrics_path = joinpath(metrics_dir, lowercase(name) * "_metrics_by_window.csv")
        CSV.write(metrics_path, by_tag)

        # --- Regime-conditioned sub-model ---
        # Pass the full (df, X) so each regime can apply its own finite-row mask.
        min_samples_full = 15
        yhat_regime, regime_betas, reduced_tags = _fit_regime_conditioned(
            df, X, y; min_samples_full=min_samples_full,
        )
        regime_pred_col = Symbol(name * "_RegimePred")
        df[!, regime_pred_col] = yhat_regime

        # Only score rows where the regime model produced a valid prediction.
        regime_valid = isfinite.(yhat_regime)
        regime_by_tag = _metrics_by_tag(df[regime_valid, :], regime_pred_col, col)
        regime_metrics_path = joinpath(metrics_dir, lowercase(name) * "_regime_metrics_by_window.csv")
        CSV.write(regime_metrics_path, regime_by_tag)

        target_results[name] = (
            beta=beta,
            global_metrics=global_m,
            by_tag=by_tag,
            regime_betas=regime_betas,
            regime_by_tag=regime_by_tag,
            reduced_tags=reduced_tags,
            min_samples_full=min_samples_full,
            feature_names=feature_names,
        )
    end

    preds_path = joinpath(metrics_dir, "predictions.csv")
    CSV.write(preds_path, df)

    _write_markdown_report(md_path, counts, target_results, have_lag)

    println("GABLS3 predictive baseline report written to: " * md_path)
    println("GABLS3 predictive baseline metrics written to: " * metrics_dir)
end

if abspath(PROGRAM_FILE) == @__FILE__
    if isempty(ARGS)
        run_gabls3_predictive_baseline()
    else
        custom_in = ARGS[1]
        base_name = splitext(basename(custom_in))[1]
        custom_md = joinpath("reports", "gabs3", "$(base_name)_flux_skill.md")
        custom_dir = joinpath("reports", "gabs3", "metrics_$(base_name)")
        run_gabls3_predictive_baseline(
            input_csv=custom_in,
            md_path=custom_md,
            metrics_dir=custom_dir,
        )
    end
end
