push!(LOAD_PATH, joinpath(@__DIR__, "../src"))

using UnifiedManifold, CasesIngestion, SpectralDiagnostics
using NCDatasets, DataFrames, CSV, Logging

function execute_pipeline(nc_path::String, output_csv::String)
    @info "Initializing Unified Manifold Architecture Layer (N=32)"
    ws = UnifiedManifoldWorkspace(32, 1.5, 60.0, 0.15, K_q=72)
    
    # Extract structural run time allocations
    nt = Dataset(nc_path, "r") do ds
        size(ds["theta"], size(ds["theta"], 1) == length(ds["height"][:]) ? 2 : 1)
    end
    
    records = DataFrame(TimeIdx=Int[], Ri_f=Float64[], R_W=Float64[], F_W=Float64[], chi_N=Float64[], D_eff=Float64[], RunStatus=String[])
    
    @info "Looping through $nt observational steps..."
    for t in 1:nt
        c_theta, c_u, run_status = ingest_and_project_slice!(nc_path, t, ws)
        if isnothing(c_theta); continue; end
        
        metrics = process_timestamp_metrics(t, c_theta, c_u, ws, run_status)
        push!(records, (metrics.time_idx, metrics.Ri_f, metrics.R_W, metrics.F_W, metrics.chi_N, metrics.D_eff, metrics.Status))
    end
    
    CSV.write(output_csv, records)
    @info "Pipeline Run Terminated. Output written to: $output_csv"
end

if length(ARGS) >= 2
    execute_pipeline(ARGS[1], ARGS[2])
else
    println("Usage: julia RunCampaignPipeline.jl <input_nc> <output_csv>")
end
