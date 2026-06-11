using Test
using NCDatasets
using DataFrames

include("../src/Summary.jl")
include("../scripts/DiagnosticsBaseline.jl")
using .AtmosphericDataPipeline
using .DiagnosticsBaseline

function create_dummy_netcdf(values::Vector{Float64}, varname::String)
    path = tempname() * ".nc"

    NCDataset(path, "c") do ds
        defDim(ds, "time", length(values))
        v = defVar(ds, varname, Float64, ("time",))
        v[:] = values
    end

    return path
end

@testset "clean_sonic_spikes" begin
    signal = [1.0, 1.2, 0.8, 60.0, 1.1, 0.9, 1.0]
    cleaned = clean_sonic_spikes(signal; threshold=3.5)
    @test isnan(cleaned[4])
    @test cleaned[1] == signal[1]
    @test cleaned[2] == signal[2]
    @test length(cleaned) == length(signal)

    flat_signal = fill(5.0, 8)
    @test clean_sonic_spikes(flat_signal; threshold=3.5) == flat_signal

    clean_signal = [1.0, 1.1, 0.9, 1.0, 1.05]
    result = clean_sonic_spikes(clean_signal; threshold=3.5)
    @test result == clean_signal
    @test !any(isnan, result)
end

@testset "validate_physical_gradients" begin
    heights = [0.0, 5.0, 10.0, 20.0]
    theta_stable = [290.0, 291.0, 292.2, 294.0]
    @test validate_physical_gradients(heights, theta_stable)

    theta_unstable = [290.0, 289.0, 288.0, 287.0]
    @test !validate_physical_gradients(heights, theta_unstable)

    @test_throws AssertionError validate_physical_gradients([0.0, 5.0], [290.0])
end

@testset "audit_spectral_conditioning" begin
    heights_good = collect(0.0:5.0:60.0)
    @test audit_spectral_conditioning(heights_good, 6, 0.3)

    heights_degenerate = [10.0, 10.0, 10.0]
    @test !audit_spectral_conditioning(heights_degenerate, 6, 0.3)
end

@testset "synthetic ingestion integration" begin
    values = [280.0, 281.0, 282.0, -1037.0, 283.0, 284.0]
    nc_path = create_dummy_netcdf(values, "theta")
    dataset = NetCDFDataset(nc_path)

    heights = [0.0, 5.0, 10.0, 20.0, 40.0]
    theta_ok = [290.0, 291.0, 292.0, 293.5, 295.0]
    theta_bad = [290.0, 289.0, 288.0, 287.0, 286.0]
    sonic_signal = [0.2, 0.21, 0.22, 6.5, 0.19, 0.2]

    pass_result = run_validation_gate(
        dataset,
        "theta",
        heights,
        theta_ok;
        signal=sonic_signal,
        N=6,
        α_stretch=0.3,
        spike_threshold=3.5
    )

    @test pass_result.physical_gradients_pass
    @test pass_result.spectral_conditioning_pass
    @test pass_result.spikes_filtered
    @test pass_result.downstream_allowed
    @test pass_result.summary.total_records == length(values)
    @test pass_result.summary.missing_count == 1

    fail_result = run_validation_gate(
        dataset,
        "theta",
        heights,
        theta_bad;
        signal=sonic_signal,
        N=6,
        α_stretch=0.3,
        spike_threshold=3.5
    )

    @test !fail_result.physical_gradients_pass
    @test fail_result.spectral_conditioning_pass
    @test !fail_result.downstream_allowed

    rm(nc_path; force=true)
end

@testset "diagnostics baseline cleaning and summaries" begin
    raw = DataFrame(
        D_eff = ["1.0", "2.0", "bad", 4.0],
        F_W = [0.1, 0.2, 0.3, "0.4"],
        chi_N = [0.01, 0.02, 0.03, 0.04],
        Ri_g = [0.5, "0.6", 0.7, NaN],
        TimeIdx = [1, 2, 3, 4]
    )

    clean, ri_name, dropped, missing_cols = clean_diagnostics_frame(raw)
    @test isempty(missing_cols)
    @test ri_name == "Ri_g"
    @test dropped == 2
    @test nrow(clean) == 2
    @test clean.D_eff == [1.0, 2.0]

    dsum = d_eff_summary(clean.D_eff; early_count=1)
    @test isapprox(dsum.mean, 1.5; atol=1e-12)
    @test isapprox(dsum.early, 1.0; atol=1e-12)
    @test isapprox(dsum.late, 2.0; atol=1e-12)

    macros = diagnostics_summary_macros(clean; early_count=1)
    @test macros["DefEffMean"] == "1.50"
    @test macros["DiagnosticsFormulaVersion"] == BASELINE_VERSION
end
