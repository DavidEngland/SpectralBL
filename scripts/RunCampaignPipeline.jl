# scripts/RunCampaignPipeline.jl

# 1. Force-load the module source file explicitly before bringing the namespace into scope
const MODULE_SRC = joinpath(@__DIR__, "..", "src", "Summary.jl")
if isfile(MODULE_SRC)
    include(MODULE_SRC)
else
    error("Pipeline aborted: Core module source file missing at $MODULE_SRC")
end

const DIAGNOSTICS_SRC = joinpath(@__DIR__, "..", "src", "SpectralDiagnostics.jl")
if isfile(DIAGNOSTICS_SRC)
    include(DIAGNOSTICS_SRC)
else
    error("Pipeline aborted: Diagnostics module source file missing at $DIAGNOSTICS_SRC")
end

# 2. Bring your module and the rest of your development environment into scope
using .AtmosphericDataPipeline  # Note the leading dot (.) which indicates a locally included module namespace
using .SpectralDiagnostics
using CasesIngestion
using UnifiedManifold
using ProgressMeter, CSV, DataFrames, NCDatasets, Statistics

# --- FIXED: SBL Compactification Workspace Geometry ---
# 1. Fixed z_max from 50.0 to 55.0 to align exactly with your top 'u_55m' sonic metadata
# 2. Hardcoded your paper's exact tuning choice (alpha_stretch = 0.05)
const CAMPAIGN_WORKSPACE = UnifiedManifoldWorkspace(32, 1.5, 55.0, 0.05)

function execute_campaign_sweep()
    data_dir = "data/ncar_eol_dee0099881"
    output_csv = "data/diagnostic_trajectory.csv"

    # Define the exact physical height array mapping to your NetCDF sonic names
    heights = [1.5, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 55.0]
    target_vars = ["u_1_5m", "u_5m", "u_10m", "u_20m", "u_30m", "u_40m", "u_50m", "u_55m"]
    tc_vars     = ["tc_1_5m", "tc_5m", "tc_10m", "tc_20m", "tc_30m", "tc_40m", "tc_50m", "tc_55m"]

    # If an NC file is passed explicitly, process only that file; else process the full campaign month.
    nc_files = String[]
    if !isempty(ARGS)
        provided = ARGS[1]
        if isfile(provided)
            push!(nc_files, basename(provided))
            data_dir = dirname(provided)
        else
            error("Provided NetCDF path does not exist: $provided")
        end
    else
        nc_files = sort(filter(f -> match(r"^cases\.\d+\.nc$", f) !== nothing, readdir(data_dir)))
    end

    if isempty(nc_files)
        error("No campaign files found in $data_dir. Check paths.")
    end
    
    println("Beginning processing sweep over $(length(nc_files)) daily target logs...")
    master_df = DataFrame()

    valid_sample(x) = !ismissing(x) && x != -1037.0 && !isnan(Float64(x))

    for nc_file in nc_files
        full_nc_path = joinpath(data_dir, nc_file)
        println("\nIngesting Campaign File: $nc_file")
        
        # Instantiate your NetCDF wrapper type for the Multiple Dispatch engine
        pipeline_dataset = NetCDFDataset(full_nc_path)

        # Run a global validation sweep across the target day just to print data quality logs
        run_pipeline_check(pipeline_dataset, target_vars)

        # Open dataset frame to safely extract and slice time-dependent arrays
        Dataset(full_nc_path, "r") do ds
            # Handle both 1D (time) and 2D (sample x time) layouts.
            u55 = ds["u_55m"]
            t_steps = ndims(u55) == 1 ? length(u55) : size(u55, 2)

            slice_mean(var, t_idx) = begin
                if ndims(var) == 1
                    value = var[t_idx]
                    return valid_sample(value) ? Float64(value) : NaN
                end
                samples = filter(valid_sample, var[:, t_idx])
                return isempty(samples) ? NaN : mean(Float64.(samples))
            end

            @showprogress "Slicing Profile Timeline: " for t in 1:t_steps

                # --- FIXED: Dynamically build the vertical profiles for this time slice (t) ---
                # NCAR missing values (-1037.0) are handled safely here via filtering
                theta_profile = Float64[]
                u_profile     = Float64[]

                for i in 1:length(heights)
                    tc_val = slice_mean(ds[tc_vars[i]], t)
                    u_val  = slice_mean(ds[target_vars[i]], t)

                    push!(theta_profile, tc_val)
                    push!(u_profile, u_val)
                end

                # Skip slice entirely if the raw vector profiling contains corrupted NaN fields
                if any(isnan, theta_profile) || any(isnan, u_profile); continue; end

                # --- INTERCEPT: Evaluate via the AtmosphericDataPipeline Quality Gate ---
                # Use the top level stream 'u_55m' slice as our raw sonic spike test input
                u_top_stream = if ndims(ds["u_55m"]) == 1
                    [Float64(coalesce(ds["u_55m"][t], NaN))]
                else
                    Float64.(coalesce.(vec(ds["u_55m"][:, t]), NaN))
                end

                gate_result = run_validation_gate(
                    pipeline_dataset,
                    "u_55m",
                    heights,
                    theta_profile;
                    signal=u_top_stream,
                    N=32,
                    α_stretch=0.05,
                    spike_threshold=3.5
                )

                # Keep spectral conditioning as a hard gate; treat physical-gradient failures as warnings.
                if !gate_result.spectral_conditioning_pass
                    continue
                end

                # --- MANIFOLD PROJECTION ---
                # Invoke your fixed ingestion function passing the valid workspace type
                result = ingest_and_project_slice!(full_nc_path, t, CAMPAIGN_WORKSPACE)
                if result === nothing || result[1] === nothing; continue; end

                c_theta, c_u, run_status = result

                status_with_gate = gate_result.physical_gradients_pass ? run_status : string(run_status, " | PhysicalGateWarn")
                metrics = process_timestamp_metrics(t, c_theta, c_u, CAMPAIGN_WORKSPACE, status_with_gate; theta_ref=293.15)

                file_date_match = match(r"cases\.(\d+)\.nc$", nc_file)
                file_date = file_date_match === nothing ? "unknown" : file_date_match.captures[1]

                # Build individual row entry matching your output schemas
                row_data = DataFrame(
                    FileDate = file_date,
                    TimeIdx = metrics.time_idx,
                    Ri_f = metrics.Ri_f,
                    R_W = metrics.R_W,
                    F_W = metrics.F_W,
                    chi_N = metrics.chi_N,
                    D_eff = metrics.D_eff,
                    E_total = metrics.E_total,
                    E_wave = metrics.E_wave,
                    E_turb = metrics.E_turb,
                    peak_mode = metrics.peak_mode,
                    wave_window_min = metrics.wave_window_min,
                    wave_window_max = metrics.wave_window_max,
                    peak_in_wave_window = metrics.peak_in_wave_window,
                    RunStatus = metrics.Status
                )
                append!(master_df, row_data)
            end
        end
    end
    
    # Save clean, non-stationary continuous profiles to your CSV artifact
    CSV.write(output_csv, master_df)
    println("\n✓ Successfully finalized trajectory tracking. Data compiled in: $output_csv")
end

# Trigger calculation handle when run from command line
if abspath(PROGRAM_FILE) == @__FILE__
    execute_campaign_sweep()
end