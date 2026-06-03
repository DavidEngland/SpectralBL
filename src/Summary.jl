module AtmosphericDataPipeline

using NCDatasets
using Printf
using Statistics
using LinearAlgebra

include("validate_physical_gradients.jl")
include("audit_spectral_conditioning.jl")
include("clean_sonic_spikes.jl")

export AbstractAtmosphericDataset
export NetCDFDataset
export ProfileSummary
export ValidationGateResult
export summarize_conditions
export run_pipeline_check
export validate_physical_gradients
export audit_spectral_conditioning
export clean_sonic_spikes
export run_validation_gate

# ============================================================================
# 1. Type Hierarchy for Extensible Ingestion
# ============================================================================
abstract type AbstractAtmosphericDataset end

struct NetCDFDataset <: AbstractAtmosphericDataset
    filepath::String
end

# Future expansion placeholders:
# struct GRIB2Dataset <: AbstractAtmosphericDataset; filepath::String; end
# struct HDF5Dataset  <: AbstractAtmosphericDataset; filepath::String; end

# Struct to hold unified profile summaries across any format
struct ProfileSummary
    total_records::Int
    missing_count::Int
    missing_percentage::Float64
    mean_value::Float64
    min_value::Float64
    max_value::Float64
end

struct ValidationGateResult
    summary::ProfileSummary
    physical_gradients_pass::Bool
    spectral_conditioning_pass::Bool
    spikes_filtered::Bool
    downstream_allowed::Bool
end

# ============================================================================
# 2. Pipeline Core Functions (Multiple Dispatch Engine)
# ============================================================================

"""
    summarize_conditions(dataset::NetCDFDataset, target_variable::String)

Dispatches NetCDF-specific parsing logic to extract and validate vertical profile metrics.
"""
function summarize_conditions(dataset::NetCDFDataset, target_variable::String)::ProfileSummary
    if !isfile(dataset.filepath)
        throw(SystemError("File not found at target path: $(dataset.filepath)"))
    end

    # Safely scope the NetCDF file handling
    local_summary = Dataset(dataset.filepath, "r") do ds
        if !haskey(ds, target_variable)
            available = join(keys(ds), ", ")
            throw(KeyError("Variable '$target_variable' not found. Available variables: [$available]"))
        end

        # Read into memory efficiently
        raw_data = ds[target_variable][:]

        # SBL Quality Assurance Bounds
        # NCAR ISFF conventions use -1037.0. Standardize with NaN checks.
        ncar_missing_flag = -1037.0

        is_missing(x) = (x == ncar_missing_flag || isnan(x) || ismissing(x))

        total_elements = length(raw_data)
        missing_count = count(is_missing, raw_data)
        missing_pct = (missing_count / total_elements) * 100

        # Extract validated physical vector
        clean_data = filter(!is_missing, raw_data)

        if isempty(clean_data)
            return ProfileSummary(total_elements, missing_count, missing_pct, NaN, NaN, NaN)
        end

        return ProfileSummary(
            total_elements,
            missing_count,
            missing_pct,
            mean(clean_data),
            minimum(clean_data),
            maximum(clean_data)
        )
    end

    return local_summary
end

"""
    run_validation_gate(dataset::AbstractAtmosphericDataset,
                        target_variable::String,
                        heights::Vector{Float64},
                        theta::Vector{Float64};
                        signal::Union{Nothing, Vector{Float64}}=nothing,
                        N::Int=32,
                        α_stretch::Float64=0.3,
                        spike_threshold::Float64=3.5)

Runs the ingestion summary stage plus all validation stages. Downstream execution
is permitted only when both physical and spectral validators pass.
"""
function run_validation_gate(
    dataset::AbstractAtmosphericDataset,
    target_variable::String,
    heights::Vector{Float64},
    theta::Vector{Float64};
    signal::Union{Nothing, Vector{Float64}}=nothing,
    N::Int=32,
    α_stretch::Float64=0.3,
    spike_threshold::Float64=3.5
)::ValidationGateResult
    summary = summarize_conditions(dataset, target_variable)
    physical_ok = validate_physical_gradients(heights, theta)
    spectral_ok = audit_spectral_conditioning(heights, N, α_stretch)

    spikes_filtered = false
    if signal !== nothing
        cleaned = clean_sonic_spikes(signal; threshold=spike_threshold)
        spikes_filtered = any(isnan, cleaned)
    end

    return ValidationGateResult(
        summary,
        physical_ok,
        spectral_ok,
        spikes_filtered,
        physical_ok && spectral_ok
    )
end

# Future Dispatch Implementations would look like this:
# function summarize_conditions(dataset::GRIB2Dataset, target_variable::String)::ProfileSummary
#     # GRIB2 parsing logic using GRIB.jl
# end


# ============================================================================
# 3. Execution & Validation Entrypoint
# ============================================================================

"""
    run_pipeline_check(dataset::AbstractAtmosphericDataset, variables::Vector{String})

Top-level orchestration tool for processing automated batch sweeps across your paper's variables.
"""
function run_pipeline_check(dataset::AbstractAtmosphericDataset, variables::Vector{String})
    println("="^70)
    println("     METRIC-CONSISTENT SBL PIPELINE INTEGRITY SWEEP")
    println("="^70)
    println("Target Source File : $(dataset.filepath)")
    println("Data Context Type  : $(typeof(dataset))")
    println("-"^70)

    for var in variables
        try
            stats = summarize_conditions(dataset, var)

            @printf("Variable: %-12s | Quality: %6.2f%% Valid | Range: [%6.2f, %6.2f] | μ: %6.2f\n",
                var,
                100.0 - stats.missing_percentage,
                stats.min_value,
                stats.max_value,
                stats.mean_value
            )

            # Critical Alert: Low-Level Jet Shear Warning
            if var == "u_55m" && stats.max_value > 15.0
                println("  [ALERT] High velocity jet structure detected at upper boundary (max: $(round(stats.max_value, digits=2)) m/s). Check CFL stability criteria.")
            end

        catch e
            println("  [ERROR processing $var]: $(e)")
        end
    end
    println("="^70)
end

end # module