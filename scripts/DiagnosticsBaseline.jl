# scripts/DiagnosticsBaseline.jl
module DiagnosticsBaseline

using CSV
using DataFrames
using Statistics
using Dates
using Printf

export BASELINE_VERSION, BASELINE_SOURCE, load_master_trajectory, clean_diagnostics_frame,
       d_eff_summary, diagnostics_summary, diagnostics_summary_macros, tex_provenance_comments

const BASELINE_VERSION = "v1.0.0"
const BASELINE_SOURCE = "src/SpectralDiagnostics.jl"

function _find_colkey(df::DataFrame, colname::AbstractString)
    for k in names(df)
        if string(k) == colname
            return k
        end
    end
    return nothing
end

function _to_float(x)
    if ismissing(x)
        return NaN
    elseif x isa Number
        return Float64(x)
    end

    sx = strip(string(x))
    if isempty(sx)
        return NaN
    end

    v = tryparse(Float64, sx)
    return v === nothing ? NaN : v
end

function load_master_trajectory(data_dir::AbstractString="data"; merged_filename::AbstractString="diagnostic_trajectory.csv")
    merged_path = joinpath(data_dir, merged_filename)
    if isfile(merged_path) && filesize(merged_path) > 0
        try
            df = CSV.read(merged_path, DataFrame)
            if nrow(df) > 0
                return df, merged_path
            end
        catch
            # Fall back to shards below.
        end
    end

    shard_files = sort(filter(f -> startswith(f, "trajectory_") && endswith(f, ".csv"), readdir(data_dir)))
    parts = DataFrame[]
    for f in shard_files
        p = joinpath(data_dir, f)
        if !isfile(p) || filesize(p) == 0
            continue
        end
        part = CSV.read(p, DataFrame)
        if nrow(part) > 0
            push!(parts, part)
        end
    end

    if isempty(parts)
        return DataFrame(), merged_path
    end

    merged = reduce((a, b) -> vcat(a, b; cols=:union), parts)
    CSV.write(merged_path, merged)
    return merged, merged_path
end

function clean_diagnostics_frame(df::DataFrame)
    d_key = _find_colkey(df, "D_eff")
    fw_key = _find_colkey(df, "F_W")
    chi_key = _find_colkey(df, "chi_N")
    time_key = _find_colkey(df, "TimeIdx")
    file_key = _find_colkey(df, "FileDate")
    status_key = _find_colkey(df, "RunStatus")
    ri_b_key = _find_colkey(df, "Ri_b")
    ri_g_key = _find_colkey(df, "Ri_g")
    ri_f_key = _find_colkey(df, "Ri_f")
    e_meso_key = _find_colkey(df, "E_meso")
    e_int_key = _find_colkey(df, "E_interaction")

    missing = String[]
    if d_key === nothing; push!(missing, "D_eff"); end
    if fw_key === nothing; push!(missing, "F_W"); end
    if chi_key === nothing; push!(missing, "chi_N"); end

    ri_source = ""
    ri_key = nothing
    if ri_b_key !== nothing
        ri_key = ri_b_key
        ri_source = "Ri_b"
    elseif ri_g_key !== nothing
        ri_key = ri_g_key
        ri_source = "Ri_g"
    elseif ri_f_key !== nothing
        ri_key = ri_f_key
        ri_source = "Ri_f"
    else
        push!(missing, "Ri_b|Ri_g|Ri_f")
    end

    if !isempty(missing)
        return DataFrame(), ri_source, 0, missing
    end

    n = nrow(df)
    d_eff = [_to_float(df[i, d_key]) for i in 1:n]
    f_w = [_to_float(df[i, fw_key]) for i in 1:n]
    chi_n = [_to_float(df[i, chi_key]) for i in 1:n]
    ri = [_to_float(df[i, ri_key]) for i in 1:n]

    valid = isfinite.(d_eff) .& isfinite.(f_w) .& isfinite.(chi_n) .& isfinite.(ri)
    dropped = count(!, valid)

    out = DataFrame(
        D_eff = d_eff[valid],
        F_W = f_w[valid],
        chi_N = chi_n[valid],
        Ri = ri[valid]
    )

    if time_key !== nothing
        out[!, :TimeIdx] = [_to_float(v) for v in df[valid, time_key]]
    else
        out[!, :TimeIdx] = collect(1.0:nrow(out))
    end

    if file_key !== nothing
        out[!, :FileDate] = Int64.(round.([_to_float(v) for v in df[valid, file_key]]))
    end

    if status_key !== nothing
        out[!, :RunStatus] = string.(df[valid, status_key])
    end

    if e_meso_key !== nothing
        e_meso = [_to_float(v) for v in df[valid, e_meso_key]]
        out[!, :E_meso] = e_meso
    end

    if e_int_key !== nothing
        e_int = [_to_float(v) for v in df[valid, e_int_key]]
        out[!, :E_interaction] = e_int
    end

    return out, ri_source, dropped, String[]
end

function d_eff_summary(d_eff_vals::AbstractVector{<:Real}; early_count::Int=3024)
    if isempty(d_eff_vals)
        error("Cannot summarize D_eff: no samples provided.")
    end

    vals = Float64.(d_eff_vals)
    split_idx = min(early_count, length(vals))

    early = mean(@view vals[1:split_idx])
    late = split_idx < length(vals) ? mean(@view vals[(split_idx + 1):end]) : early

    return (
        mean = mean(vals),
        early = early,
        late = late,
        min = minimum(vals),
        max = maximum(vals),
        samples = length(vals)
    )
end

function _metric_summary(vals::AbstractVector{<:Real}; early_count::Int=3024)
    values = Float64.(vals)
    if isempty(values)
        error("Cannot summarize metric: no samples provided.")
    end

    split_idx = min(early_count, length(values))
    early = mean(@view values[1:split_idx])
    late = split_idx < length(values) ? mean(@view values[(split_idx + 1):end]) : early

    return (
        mean = mean(values),
        early = early,
        late = late,
        min = minimum(values),
        max = maximum(values),
        samples = length(values)
    )
end

function diagnostics_summary(clean_df::DataFrame; early_count::Int=3024)
    if nrow(clean_df) == 0
        error("Cannot summarize diagnostics: no valid rows provided.")
    end

    out = (
        d_eff = _metric_summary(clean_df.D_eff; early_count=early_count),
        f_w = _metric_summary(clean_df.F_W; early_count=early_count),
        chi_n = _metric_summary(clean_df.chi_N; early_count=early_count),
        samples = nrow(clean_df)
    )

    return out
end

function diagnostics_summary_macros(clean_df::DataFrame; early_count::Int=3024)
    summary = diagnostics_summary(clean_df; early_count=early_count)
    colkeys = Dict(string(k) => k for k in names(clean_df))

    return Dict{String,String}(
        "DiagnosticsFormulaVersion" => BASELINE_VERSION,
        "DiagnosticsSamples" => string(summary.samples),

        "DefEffMean" => @sprintf("%.2f", summary.d_eff.mean),
        "DefEffEarly" => @sprintf("%.2f", summary.d_eff.early),
        "DefEffLate" => @sprintf("%.2f", summary.d_eff.late),
        "DefEffMin" => @sprintf("%.2f", summary.d_eff.min),
        "DefEffMax" => @sprintf("%.2f", summary.d_eff.max),

        "FwMean" => @sprintf("%.3f", summary.f_w.mean),
        "FwEarly" => @sprintf("%.3f", summary.f_w.early),
        "FwLate" => @sprintf("%.3f", summary.f_w.late),
        "FwMin" => @sprintf("%.3f", summary.f_w.min),
        "FwMax" => @sprintf("%.3f", summary.f_w.max),

        "ChiNMean" => @sprintf("%.3f", summary.chi_n.mean),
        "ChiNEarly" => @sprintf("%.3f", summary.chi_n.early),
        "ChiNLate" => @sprintf("%.3f", summary.chi_n.late),
        "ChiNMin" => @sprintf("%.3f", summary.chi_n.min),
        "ChiNMax" => @sprintf("%.3f", summary.chi_n.max)
    )

    if haskey(colkeys, "E_meso")
        e_meso_col = colkeys["E_meso"]
        e_meso_raw = clean_df[!, e_meso_col]
        e_meso_vals = e_meso_raw[isfinite.(e_meso_raw)]
        if !isempty(e_meso_vals)
            e_meso = _metric_summary(e_meso_vals; early_count=early_count)
            out["EMesoMean"] = @sprintf("%.3f", e_meso.mean)
            out["EMesoEarly"] = @sprintf("%.3f", e_meso.early)
            out["EMesoLate"] = @sprintf("%.3f", e_meso.late)
            out["EMesoMin"] = @sprintf("%.3f", e_meso.min)
            out["EMesoMax"] = @sprintf("%.3f", e_meso.max)
        end
    end

    if haskey(colkeys, "E_interaction")
        e_int_col = colkeys["E_interaction"]
        e_int_raw = clean_df[!, e_int_col]
        e_int_vals = e_int_raw[isfinite.(e_int_raw)]
        if !isempty(e_int_vals)
            e_int = _metric_summary(e_int_vals; early_count=early_count)
            out["EIntMean"] = @sprintf("%.3f", e_int.mean)
            out["EIntEarly"] = @sprintf("%.3f", e_int.early)
            out["EIntLate"] = @sprintf("%.3f", e_int.late)
            out["EIntMin"] = @sprintf("%.3f", e_int.min)
            out["EIntMax"] = @sprintf("%.3f", e_int.max)
        end
    end

    return out
end

function tex_provenance_comments(sample_count::Int, source_path::AbstractString)
    timestamp = Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS")
    return [
        "% Baseline-Version: $(BASELINE_VERSION)",
        "% Baseline-Source: $(BASELINE_SOURCE)",
        "% Source-Data: $(source_path)",
        "% Samples: $(sample_count)",
        "% Generated-UTC: $(timestamp)"
    ]
end

end # module
