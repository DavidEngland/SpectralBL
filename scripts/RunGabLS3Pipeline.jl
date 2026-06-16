# scripts/RunGabLS3Pipeline.jl

const SUMMARY_SRC = joinpath(@__DIR__, "..", "src", "Summary.jl")
const DIAGNOSTICS_SRC = joinpath(@__DIR__, "..", "src", "SpectralDiagnostics.jl")
const CASES_INGEST_SRC = joinpath(@__DIR__, "..", "src", "CasesIngestion.jl")
const GABLS3_INGEST_SRC = joinpath(@__DIR__, "..", "src", "GabLS3Ingestion.jl")

isfile(SUMMARY_SRC) || error("Missing module source: $SUMMARY_SRC")
isfile(DIAGNOSTICS_SRC) || error("Missing module source: $DIAGNOSTICS_SRC")
isfile(CASES_INGEST_SRC) || error("Missing module source: $CASES_INGEST_SRC")
isfile(GABLS3_INGEST_SRC) || error("Missing module source: $GABLS3_INGEST_SRC")

include(SUMMARY_SRC)
include(DIAGNOSTICS_SRC)
include(CASES_INGEST_SRC)
include(GABLS3_INGEST_SRC)

using .AtmosphericDataPipeline
using .SpectralDiagnostics
using .GabLS3Ingestion
using UnifiedManifold
using CSV
using DataFrames
using NCDatasets
using ProgressMeter
using Statistics

const GABLS3_WORKSPACE = UnifiedManifoldWorkspace(32, 0.0, 800.0, 0.01; invert_windows=true)
const DEFAULT_INPUT = joinpath("data", "gabs3", "gabls3_scm_cabauw_obs_v33.nc")
const DEFAULT_OUTPUT = joinpath("data", "trajectory_gabls3.csv")
const DEFAULT_AGG_DIR = joinpath("data", "gabs3", "aggregates")

function _mean_or_nan(v::AbstractVector{<:Real})
    vals = filter(isfinite, Float64.(v))
    return isempty(vals) ? NaN : mean(vals)
end

function _near_surface_or_nan(z::AbstractVector{<:Real}, v::AbstractVector{<:Real})
    if isempty(z) || isempty(v)
        return NaN
    end
    i = argmin(Float64.(z))
    return Float64(v[i])
end

function _mean_abs_shear_or_nan(z::AbstractVector{<:Real}, v::AbstractVector{<:Real})
    if length(z) < 2 || length(v) < 2
        return NaN
    end
    zz = Float64.(z)
    vv = Float64.(v)
    s = Float64[]
    for i in 1:(length(zz) - 1)
        dz = zz[i + 1] - zz[i]
        if abs(dz) > eps(Float64)
            push!(s, abs((vv[i + 1] - vv[i]) / dz))
        end
    end
    return isempty(s) ? NaN : mean(s)
end

function execute_gabls3_sweep(; input_nc::String=DEFAULT_INPUT, output_csv::String=DEFAULT_OUTPUT)
    isfile(input_nc) || error("GABLS3 input file missing: $input_nc")
    mkpath(dirname(output_csv))
    mkpath(DEFAULT_AGG_DIR)

    dataset = NetCDFDataset(input_nc)
    run_pipeline_check(dataset, ["th", "u", "wt", "uw", "vw"])

    rows = DataFrame()

    # Rolling state for lagged inertial-oscillation features.
    # These track the column-mean of u and v from the previous time slice so
    # that Δu/Δt and Δv/Δt can be computed as finite-difference tendency terms.
    # NaN on the first slice is intentional; the predictor filters these rows.
    prev_u_mean   = NaN
    prev_v_mean   = NaN
    prev_t_hours  = NaN

    Dataset(input_nc, "r") do ds
        t_steps = length(ds["time"])

        @showprogress "Processing GABLS3 timeline: " for t_idx in 1:t_steps
            slice = ingest_and_project_gabls3_slice!(input_nc, t_idx, GABLS3_WORKSPACE)
            slice === nothing && continue

            gate = run_validation_gate(
                dataset,
                "th",
                Float64.(slice.z_mean),
                Float64.(slice.theta_profile);
                signal=Float64.(slice.u_profile),
                N=GABLS3_WORKSPACE.N,
                α_stretch=GABLS3_WORKSPACE.alpha_stretch,
                spike_threshold=3.5,
            )

            status = slice.status
            if !gate.physical_gradients_pass
                status = string(status, " | PhysicalGateWarn")
            end
            if !gate.spectral_conditioning_pass
                status = string(status, " | SpectralGateWarn")
            end

            metrics = process_timestamp_metrics(t_idx, slice.c_theta, slice.c_u, GABLS3_WORKSPACE, status)
            f_w_adaptive, peak_m, n_min_eff, in_window, _, compression_factor = calculate_adaptive_wave_fraction(
                GABLS3_WORKSPACE,
                slice.c_u,
                metrics.D_eff;
                alpha_floor=1.5,
            )

            heat_flux_mean = _mean_or_nan(slice.heat_flux)
            mom_u_flux_mean = _mean_or_nan(slice.mom_u_flux)
            mom_v_flux_mean = _mean_or_nan(slice.mom_v_flux)

            heat_flux_surface = _near_surface_or_nan(slice.z_flux, slice.heat_flux)
            mom_u_flux_surface = _near_surface_or_nan(slice.z_flux, slice.mom_u_flux)
            mom_v_flux_surface = _near_surface_or_nan(slice.z_flux, slice.mom_v_flux)

            u_shear_mean = _mean_abs_shear_or_nan(slice.z_mean, slice.u_profile)
            v_shear_mean = _mean_abs_shear_or_nan(slice.z_mean, slice.v_profile)
            shear_mag_mean = isfinite(u_shear_mean) && isfinite(v_shear_mean) ? sqrt(u_shear_mean^2 + v_shear_mean^2) : NaN
            deff_x_shear = isfinite(shear_mag_mean) ? metrics.D_eff * shear_mag_mean : NaN

            # Lagged inertial-oscillation tendency terms (Δu/Δt, Δv/Δt).
            # These provide the phase information that a static snapshot cannot
            # capture for the ageostrophic V-component oscillation.
            cur_u_mean = _mean_or_nan(slice.u_profile)
            cur_v_mean = _mean_or_nan(slice.v_profile)
            dt_hours = isfinite(prev_t_hours) ? (slice.time_hours - prev_t_hours) : NaN
            delta_u_dt = isfinite(prev_u_mean) && isfinite(dt_hours) && abs(dt_hours) > eps(Float64) ?
                         (cur_u_mean - prev_u_mean) / dt_hours : NaN
            delta_v_dt = isfinite(prev_v_mean) && isfinite(dt_hours) && abs(dt_hours) > eps(Float64) ?
                         (cur_v_mean - prev_v_mean) / dt_hours : NaN
            inertial_mag_dt = isfinite(delta_u_dt) && isfinite(delta_v_dt) ?
                              sqrt(delta_u_dt^2 + delta_v_dt^2) : NaN

            # Advance rolling state for next slice.
            prev_u_mean  = cur_u_mean
            prev_v_mean  = cur_v_mean
            prev_t_hours = slice.time_hours

            row = DataFrame(
                SourceCase = ["GABLS3"],
                TimeIdx = [metrics.time_idx],
                SourceTimeHours = [slice.time_hours],
                WindowTag = [slice.window_tag],
                Ri_g = [metrics.Ri_g],
                Ri_b = [metrics.Ri_b],
                R_W = [metrics.R_W],
                F_W = [f_w_adaptive],
                chi_N = [metrics.chi_N],
                D_eff = [metrics.D_eff],
                E_total = [metrics.E_total],
                E_wave = [metrics.E_total * f_w_adaptive],
                E_turb = [metrics.E_total * (1.0 - f_w_adaptive)],
                E_meso = [metrics.E_meso],
                E_interaction = [metrics.E_interaction],
                compression_factor = [compression_factor],
                peak_mode = [peak_m],
                wave_window_min = [n_min_eff],
                wave_window_max = [metrics.wave_window_max],
                peak_in_wave_window = [in_window],
                HeatFluxTruth = [heat_flux_mean],
                MomUFluxTruth = [mom_u_flux_mean],
                MomVFluxTruth = [mom_v_flux_mean],
                HeatFluxSurfaceTruth = [heat_flux_surface],
                MomUFluxSurfaceTruth = [mom_u_flux_surface],
                MomVFluxSurfaceTruth = [mom_v_flux_surface],
                UShearMean = [u_shear_mean],
                VShearMean = [v_shear_mean],
                ShearMagMean = [shear_mag_mean],
                Deff_x_ShearMag = [deff_x_shear],
                DeltaU_Dt = [delta_u_dt],
                DeltaV_Dt = [delta_v_dt],
                InertialMagDt = [inertial_mag_dt],
                MeanLevels = [length(slice.z_mean)],
                FluxLevels = [length(slice.z_flux)],
                RunStatus = [status],
            )
            append!(rows, row)
        end
    end

    CSV.write(output_csv, rows)

    if nrow(rows) > 0
        agg = combine(
            groupby(rows, :WindowTag),
            :D_eff => mean => :D_eff_mean,
            :chi_N => mean => :chi_N_mean,
            :F_W => mean => :F_W_mean,
            :ShearMagMean => mean => :ShearMagMean_mean,
            :InertialMagDt => (x -> mean(filter(isfinite, x))) => :InertialMagDt_mean,
            :HeatFluxTruth => mean => :HeatFluxTruth_mean,
            :MomUFluxTruth => mean => :MomUFluxTruth_mean,
            :MomVFluxTruth => mean => :MomVFluxTruth_mean,
            nrow => :n_samples,
        )
        CSV.write(joinpath(DEFAULT_AGG_DIR, "window_summary.csv"), agg)
    end

    println("GABLS3 trajectory written to: " * output_csv)
end

if abspath(PROGRAM_FILE) == @__FILE__
    input_nc = isempty(ARGS) ? DEFAULT_INPUT : ARGS[1]
    output_csv = length(ARGS) >= 2 ? ARGS[2] : DEFAULT_OUTPUT
    execute_gabls3_sweep(input_nc=input_nc, output_csv=output_csv)
end
