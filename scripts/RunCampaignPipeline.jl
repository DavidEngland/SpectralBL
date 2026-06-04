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

            # --- HARDENED SCHEMA EXTRACTION VALVE ---
            # Dynamically look for the highest available sonic variable to determine the time step ceiling
            reference_var = nothing
            for var_name in reverse(target_vars)
                if haskey(ds, var_name)
                    reference_var = ds[var_name]
                    break
                end
            end

            if reference_var === nothing
                println("⚠️ Warning: No valid target velocity variables found in $nc_file. Skipping file.")
                return
            end

            t_steps = ndims(reference_var) == 1 ? length(reference_var) : size(reference_var, 2)

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
                theta_profile = Float64[]
                u_profile     = Float64[]

                for i in 1:length(heights)
                    # If the variable exists, extract its slice mean; otherwise assign NaN safely
                    tc_val = haskey(ds, tc_vars[i])     ? slice_mean(ds[tc_vars[i]], t)     : NaN
                    u_val  = haskey(ds, target_vars[i]) ? slice_mean(ds[target_vars[i]], t) : NaN

                    push!(theta_profile, tc_val)
                    push!(u_profile, u_val)
                end

                # Skip slice entirely if the raw vector profiling contains corrupted NaN fields
                if any(isnan, theta_profile) || any(isnan, u_profile); continue; end

                # --- INTERCEPT: Evaluate via the AtmosphericDataPipeline Quality Gate ---
                # Check for alternative available upper boundary channels if u_55m is entirely missing
                top_stream_var_name = haskey(ds, "u_55m") ? "u_55m" : (haskey(ds, "u_50m") ? "u_50m" : (haskey(ds, "u_30m") ? "u_30m" : nothing))

                if top_stream_var_name === nothing
                    continue # No functional top anemometer available to test signal quality
                end

                top_stream_var = ds[top_stream_var_name]
                u_top_stream = if ndims(top_stream_var) == 1
                    [Float64(coalesce(top_stream_var[t], NaN))]
                else
                    Float64.(coalesce.(vec(top_stream_var[:, t]), NaN))
                end

                gate_result = run_validation_gate(
                    pipeline_dataset,
                    top_stream_var_name,
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
                result = ingest_and_project_slice!(full_nc_path, t, CAMPAIGN_WORKSPACE)
                if result === nothing || result[1] === nothing; continue; end

                c_theta, c_u, run_status = result

                status_with_gate = gate_result.physical_gradients_pass ? run_status : string(run_status, " | PhysicalGateWarn")
                metrics = process_timestamp_metrics(t, c_theta, c_u, CAMPAIGN_WORKSPACE, status_with_gate; theta_ref=293.15)

                # --- ADAPTIVE INTERCEPT HOOK ---
                f_w_adaptive, peak_m, n_min_eff, in_window, run_log = calculate_adaptive_wave_fraction(
                    CAMPAIGN_WORKSPACE,
                    c_u,
                    metrics.D_eff;
                    alpha_floor=1.5
                )

                # Merge the quality gate warning flags back into the updated adaptive status string
                final_status = gate_result.physical_gradients_pass ? run_log : string(run_log, " | PhysicalGateWarn")

                file_date_match = match(r"cases\.(\d+)\.nc$", nc_file)
                file_date = file_date_match === nothing ? "unknown" : file_date_match.captures[1]

                # Build individual row entry matching your output schemas with updated adaptive parameters
                row_data = DataFrame(
                    FileDate = fill(parse(Int, file_date), 1),
                    TimeIdx = fill(metrics.time_idx, 1),
                    Ri_f = fill(metrics.Ri_f, 1),
                    R_W = fill(metrics.R_W, 1),
                    F_W = fill(f_w_adaptive, 1),
                    chi_N = fill(metrics.chi_N, 1),
                    D_eff = fill(metrics.D_eff, 1),
                    E_total = fill(metrics.E_total, 1),
                    E_wave = fill(metrics.E_total * f_w_adaptive, 1),
                    E_turb = fill(metrics.E_total * (1.0 - f_w_adaptive), 1),
                    peak_mode = fill(peak_m, 1),
                    wave_window_min = fill(n_min_eff, 1),
                    wave_window_max = fill(metrics.wave_window_max, 1),
                    peak_in_wave_window = fill(in_window, 1),
                    RunStatus = fill(final_status, 1)
                )
                append!(master_df, row_data)
            end
        end
    end
    
    # Save clean, non-stationary continuous profiles to your CSV artifact
    # Extract the file prefix directly from the loop-level 'nc_file' string
    file_date_match = match(r"cases\.(\d+)\.nc$", basename(nc_files[1]))
    file_date = file_date_match === nothing ? "unknown" : file_date_match.captures[1]

    # Force a unique output CSV name per day so parallel tasks never overlap
    output_csv = joinpath("data", "trajectory_$(file_date).csv")

    # Save clean daily profiles to its own unique shard
    CSV.write(output_csv, master_df)
    println("\n✓ Successfully finalized trajectory tracking. Data compiled in: $output_csv")
end

# Trigger calculation handle when run from command line
if abspath(PROGRAM_FILE) == @__FILE__
    execute_campaign_sweep()
end