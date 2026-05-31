push!(LOAD_PATH, joinpath(@__DIR__, "../src"))

using UnifiedManifold, CasesIngestion, SpectralDiagnostics
using NCDatasets, DataFrames, CSV, Logging

function execute_pipeline(nc_path::String, output_csv::String)
    # Match the domain top precisely to the highest physical instrument (50.0m)
    @info "Initializing Unified Manifold Architecture Layer (N=32)"
    ws = UnifiedManifoldWorkspace(32, 1.5, 50.0, 0.15, K_q=72)

    # Secure total timesteps using a validated, guaranteed NCAR variable from the snapshot
    nt = Dataset(nc_path, "r") do ds
        if haskey(ds, "tc_1_5m")
            return length(ds["tc_1_5m"])
        elseif haskey(ds, "tc_50m")
            return length(ds["tc_50m"])
        else
            # Fallback: find any variable that isn't coordinate metadata
            for key in keys(ds)
                if key != "height" && key != "time" && length(size(ds[key])) == 1
                    return length(ds[key])
                end
            end
            error("Could not dynamically resolve a valid time dimension in the NetCDF schema.")
        end
    end
    
    @info "Resolved dataset dimensions: $nt sequential 5-minute snapshots found."
    
    # Initialize tabular tracking archive
    records = DataFrame(
        TimeIdx = Int[],
        Ri_f = Float64[],
        R_W = Float64[],
        F_W = Float64[],
        chi_N = Float64[],
        D_eff = Float64[],
        RunStatus = String[]
    )

    @info "Looping through campaign timeline steps..."
    for t in 1:nt
        # Call our updated SVD-truncated station unroller
        c_theta, c_u, run_status = ingest_and_project_slice!(nc_path, t, ws)
        
        # Safely step over data dropouts or sensor calibration hours
        if isnothing(c_theta) || isnothing(c_u)
            continue
        end

        # Extract metrics
        metrics = process_timestamp_metrics(t, c_theta, c_u, ws, run_status)

        push!(records, (
            metrics.time_idx,
            metrics.Ri_f,
            metrics.R_W,
            metrics.F_W,
            metrics.chi_N,
            metrics.D_eff,
            metrics.Status
        ))
    end
    
    # Save the output trajectory matrix
    CSV.write(output_csv, records)
    @info "Pipeline Run Terminated. Trajectory written to: $output_csv"
end

if length(ARGS) >= 2
    execute_pipeline(ARGS[1], ARGS[2])
else
    println("Usage: julia RunCampaignPipeline.jl <input_nc> <output_csv>")
end