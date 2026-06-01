using LinearAlgebra
using Statistics
using DataFrames
using CSV
using Plots

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using UnifiedManifold

# Reuse validated scenario utilities without re-running the script entrypoint.
include(joinpath(@__DIR__, "test_wave_reflection.jl"))

const CAMPAIGN_GEOMETRY = Dict(
    "CASES_99" => (1.5, 50.0, 0.05),
    "FLOSS_II" => (2.0, 30.0, 0.08),
    "SHEBA" => (0.5, 20.0, 0.02),
    "BLLAST" => (1.0, 60.0, 0.06),
)

function campaign_output_root(campaign_name::String)
    return joinpath("data", "universal_sponge", lowercase(campaign_name))
end

function run_campaign_sponge_validation(campaign_name::String, z_min::Float64, z_max::Float64, alpha::Float64;
                                        N::Int = 32, K_q::Int = 72,
                                        c_wave::Float64 = 2.0,
                                        dt_nominal::Float64 = 0.02,
                                        total_time::Float64 = 30.0,
                                        damping_rate::Float64 = 0.05,
                                        sample_windows::Int = 100)
    println("Initializing Verification Suite for Campaign: [", campaign_name, "]")
    println("Geometry: z_min=", z_min, " m, z_max=", z_max, " m, alpha=", alpha)

    ws = UnifiedManifoldWorkspace(N, z_min, z_max, alpha, K_q = K_q)
    B = chebyshev_basis_matrix(Vector{Float64}(ws.xi_target), ws.N)

    z_center = z_min + 0.73 * (z_max - z_min)
    sigma_z = 0.04 * (z_max - z_min)
    wavelength = 0.16 * (z_max - z_min)

    u_initial = [exp(-0.5 * ((z - z_center) / sigma_z)^2) * sin(2.0 * pi * (z - z_center) / wavelength) for z in ws.z_atm]
    u_initial = Vector{Float64}(u_initial)
    v_initial = -c_wave .* (Matrix{Float64}(ws.Dz_atm) * u_initial)

    D2_raw = Matrix{Float64}(ws.Dz_atm) * Matrix{Float64}(ws.Dz_atm)
    D2_sym = 0.5 .* (D2_raw .+ transpose(D2_raw))
    eig = eigen(Symmetric(D2_sym))
    λ_stable = min.(eig.values, 0.0)
    D2 = Matrix(eig.vectors * Diagonal(λ_stable) * transpose(eig.vectors))

    dt = stable_timestep(D2, c_wave, dt_nominal)
    n_steps = max(sample_windows, Int(round(total_time / dt)))
    sample_stride = max(1, n_steps ÷ sample_windows)

    println("Using stable dt = ", round(dt, sigdigits = 5), " with ", n_steps, " steps")

    energy_with = run_scenario!(u_initial, v_initial, D2, B, ws, c_wave, dt, n_steps;
                                apply_sponge = true, damping_rate = damping_rate, sample_stride = sample_stride)
    energy_without = run_scenario!(u_initial, v_initial, D2, B, ws, c_wave, dt, n_steps;
                                   apply_sponge = false, damping_rate = 0.0, sample_stride = sample_stride)

    n_tail = min(20, length(energy_with), length(energy_without))
    E_masked = mean(energy_with[end - n_tail + 1:end])
    E_unmasked = mean(energy_without[end - n_tail + 1:end])
    suppression_ratio = E_unmasked / (E_masked + 1e-12)

    out_dir = campaign_output_root(campaign_name)
    mkpath(out_dir)

    n_samples = min(length(energy_with), length(energy_without))
    time_axis = collect(1:n_samples)

    diag = DataFrame(
        sample = time_axis,
        energy_top_with_sponge = energy_with[1:n_samples],
        energy_top_without_sponge = energy_without[1:n_samples],
    )
    csv_path = joinpath(out_dir, "wave_reflection_metrics.csv")
    CSV.write(csv_path, diag)

    p = plot(time_axis, energy_with[1:n_samples], label = "WITH Sponge (Modal Masking)", linewidth = 2,
             xlabel = "Sample Window (Time)", ylabel = "Energy at Upper Boundary",
             title = "$(campaign_name): Sponge Layer Validation", legend = :topright,
             left_margin = 14Plots.mm, bottom_margin = 8Plots.mm, size = (1100, 550))
    plot!(p, time_axis, energy_without[1:n_samples], label = "WITHOUT Sponge (Free Reflection)", linewidth = 2)
    png_path = joinpath(out_dir, "wave_reflection_test.png")
    pdf_path = joinpath(out_dir, "wave_reflection_test.pdf")
    savefig(p, png_path)
    savefig(p, pdf_path)

    summary_path = joinpath(out_dir, "wave_reflection_summary.md")
    open(summary_path, "w") do io
        write(io, "# Universal Sponge Validation: $(campaign_name)\n\n")
        write(io, "## Geometry\n")
        write(io, "- z_min: $(z_min) m\n")
        write(io, "- z_max: $(z_max) m\n")
        write(io, "- alpha: $(alpha)\n")
        write(io, "- modes N: $(N)\n\n")
        write(io, "## Metrics\n")
        write(io, "- Mean boundary energy without sponge (tail-20): $(round(E_unmasked, digits = 6))\n")
        write(io, "- Mean boundary energy with sponge (tail-20): $(round(E_masked, digits = 6))\n")
        write(io, "- Reflection suppression ratio: $(round(suppression_ratio, digits = 2))x\n\n")
        write(io, "## Artifacts\n")
        write(io, "- wave_reflection_metrics.csv\n")
        write(io, "- wave_reflection_test.png\n")
    end

    println("=== UNIVERSAL SPONGE VALIDATION RESULTS ===")
    println("Campaign: ", campaign_name)
    println("Mean Boundary Energy (WITHOUT Sponge): ", round(E_unmasked, digits = 6))
    println("Mean Boundary Energy (WITH Sponge):    ", round(E_masked, digits = 6))
    println("Reflection Suppression Performance:    ", round(suppression_ratio, digits = 2), "x Ratio")
    println("✓ Metrics saved to: ", csv_path)
    println("✓ Plot saved to: ", png_path)
    println("✓ Plot saved to: ", pdf_path)
    println("✓ Summary saved to: ", summary_path)

    return (; campaign_name, E_unmasked, E_masked, suppression_ratio, csv_path, png_path, summary_path)
end

function run_from_args(args::Vector{String})
    if isempty(args)
        z_min, z_max, alpha = CAMPAIGN_GEOMETRY["CASES_99"]
        run_campaign_sponge_validation("CASES_99", z_min, z_max, alpha)
        return
    end

    if length(args) == 1 && args[1] == "--all"
        for campaign_name in sort(collect(keys(CAMPAIGN_GEOMETRY)))
            z_min, z_max, alpha = CAMPAIGN_GEOMETRY[campaign_name]
            run_campaign_sponge_validation(campaign_name, z_min, z_max, alpha)
        end
        return
    end

    if length(args) == 1
        campaign_name = uppercase(args[1])
        if !haskey(CAMPAIGN_GEOMETRY, campaign_name)
            error("Unknown campaign: $(args[1]). Use one of $(join(sort(collect(keys(CAMPAIGN_GEOMETRY))), ", ")) or pass custom geometry args.")
        end
        z_min, z_max, alpha = CAMPAIGN_GEOMETRY[campaign_name]
        run_campaign_sponge_validation(campaign_name, z_min, z_max, alpha)
        return
    end

    if length(args) == 4
        campaign_name = args[1]
        z_min = parse(Float64, args[2])
        z_max = parse(Float64, args[3])
        alpha = parse(Float64, args[4])
        run_campaign_sponge_validation(campaign_name, z_min, z_max, alpha)
        return
    end

    error("Usage: julia --project=\".\" scripts/run_universal_sponge_test.jl [CAMPAIGN_NAME|--all|campaign z_min z_max alpha]")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_from_args(ARGS)
end
