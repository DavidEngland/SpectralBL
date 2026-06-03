using DataFrames
using CSV
using Plots
using LinearAlgebra
using Printf  # Added for clean numerical axis formatting

# 1. Include and use your source module from src/
include("../src/Cases99.jl")
using .UnifiedManifold: UnifiedManifoldWorkspace, physical_to_computational

# Helper: Custom tick formatter to eliminate raw scientific notation (e.g., 0.0001 instead of 1e-4)
function clean_decimal_formatter(x)
    if x == 0
        return "0.0"
    elseif abs(x) >= 1000 || abs(x) < 0.01
        # Fallback to standard clean scientific notation formatting if numbers are extreme
        return @sprintf("%.1e", x)
    else
        # Force readable decimal representations for standard boundary-layer metric values
        return @sprintf("%.3f", x)
    end
end

function generate_tier_plots(output_dir::String, draft_fig_dir::String, day_suffix::String)
    trajectory_path = joinpath("data", "diagnostic_trajectory.csv")
    if !isfile(trajectory_path)
        println("! Skipping tier plot regeneration: missing ", trajectory_path)
        return
    end

    traj = CSV.read(trajectory_path, DataFrame)
    if nrow(traj) == 0
        println("! Skipping tier plot regeneration: empty trajectory CSV")
        return
    end

    # Format the date for human readable display (e.g., "_991024" -> "1999-10-24")
    raw_date = replace(day_suffix, "_" => "")
    display_date = length(raw_date) == 6 ? "19$(raw_date[1:2])-$(raw_date[3:4])-$(raw_date[5:6])" : "CASES-99 Run"

    if !isempty(day_suffix) && hasproperty(traj, :FileDate)
        traj = filter(row -> string(row.FileDate) == raw_date, traj)
        if nrow(traj) == 0
            println("! No trajectory records match day $raw_date. Skipping tier plots.")
            return
        end
    end

    ri_series = hasproperty(traj, :Ri_g) ? traj.Ri_g : traj.Ri_f

    # --- INDIVIDUAL PLOTS WITH EXPLICIT DATES ---
    p_energy = scatter(traj.D_eff, traj.F_W,
        title = "Energy-Dimension Plane\n[$display_date]",
        xlabel = "Effective Modal Dimension (D_eff)",
        ylabel = "Wave Energy Fraction (F_W)",
        yformatter = clean_decimal_formatter,
        markersize = 4, markerstrokewidth = 0.7, alpha = 0.85, legend = false)

    p_curv = scatter(traj.chi_N, ri_series,
        title = "Curvature-Stratification Plane\n[$display_date]",
        xlabel = "Spectral Curvature (χ_N)",
        ylabel = "Gradient Richardson Number (Ri_g)",
        xformatter = clean_decimal_formatter, yformatter = clean_decimal_formatter,
        markersize = 4, markerstrokewidth = 0.7, alpha = 0.85, legend = false)

    p_time = plot(traj.TimeIdx, traj.F_W,
        title = "Temporal Feature Trace\n[$display_date]",
        xlabel = "Time Index", ylabel = "Wave Energy Fraction (F_W)",
        yformatter = clean_decimal_formatter, linewidth = 2, legend = false)

    # --- NEW: COMBINED SIDE-BY-SIDE PUBLICATION ASSETS ---
    # Fig A: State-Space Analysis (Energy-Dimension alongside Curvature-Stratification)
    p_combined_states = plot(p_energy, p_curv,
        layout = (1, 2),
        size = (1200, 500),
        left_margin = 12Plots.mm, bottom_margin = 8Plots.mm)

    # Save the standalone plots
    for (name, fig) in (("tier1_plane1", p_energy), ("tier1_plane2", p_curv), ("temporal_trace", p_time))
        savefig(fig, joinpath(output_dir, name * day_suffix * ".pdf"))
        savefig(fig, joinpath(draft_fig_dir, name * day_suffix * ".pdf"))
    end

    # Save the combined side-by-side publication assets
    savefig(p_combined_states, joinpath(output_dir, "manuscript_states_combined" * day_suffix * ".pdf"))
    savefig(p_combined_states, joinpath(draft_fig_dir, "manuscript_states_combined" * day_suffix * ".pdf"))

    println("✓ Clean publication layout and side-by-side figures saved for $display_date.")
end

function run_diagnostic_pipeline(output_dir::String)
    mkpath(output_dir)
    println("Target directory verified: ", output_dir)

    day_suffix = length(ARGS) >= 1 ? "_" * ARGS[1] : ""

    N = 32
    z_0m = 1.5
    z_top = 55.0
    alpha_stretch = 0.05

    ws = UnifiedManifoldWorkspace(N, z_0m, z_top, alpha_stretch)

    # --- STEP 1: Generate Diagnostics CSV ---
    csv_path = joinpath(output_dir, "manifold_diagnostics$(day_suffix).csv")
    df = DataFrame(
        Node_Index = 0:ws.N,
        Grid_Z_m = ws.z_atm,
        Coord_Xi = ws.xi_target,
        Psi_M_Meso = ws.psi_M,
        Psi_W_Wave = ws.psi_W,
        Psi_T_Turb = ws.psi_T
    )
    CSV.write(csv_path, df)
    println("✓ Diagnostics CSV saved to: ", csv_path)

    # --- STEP 2: Unified Static Geometry Plots ---
    # Safe Unicode layout overrides to bypass standard LaTeX rendering engine dependencies
    p1 = plot(ws.xi_target, ws.z_atm, marker=:circle, linewidth=2,
              title="Hyperbolic Mapping (α = 0.05)",
              xlabel="Computational Coordinate (ξ)", ylabel="Physical Height z (m)",
              label="Grid Nodes", legend=:topleft,
              left_margin=16Plots.mm, bottom_margin=8Plots.mm)

    p2 = plot(0:ws.N, [ws.psi_M ws.psi_W ws.psi_T], linewidth=2.5,
              title="Spectral Partitioning Windows",
              xlabel="Chebyshev Mode Index (n)", ylabel="Filter Weight (ψ)",
              label=["Meso (ψ_M)" "Wave (ψ_W)" "Turb (ψ_T)"], legend=:topright,
              left_margin=10Plots.mm, bottom_margin=8Plots.mm)

    plot_path_png = joinpath(output_dir, "manifold_geometry_plots$(day_suffix).png")
    plot_path_pdf = joinpath(output_dir, "manifold_geometry_plots$(day_suffix).pdf")
    combined_plot = plot(p1, p2, layout=(1, 2), size=(1240, 520), left_margin=14Plots.mm, bottom_margin=8Plots.mm)
    savefig(combined_plot, plot_path_png)
    savefig(combined_plot, plot_path_pdf)
    println("✓ Diagnostic geometry assets confirmed.")

    draft_fig_dir = joinpath("data", "drafts", "figures")
    mkpath(draft_fig_dir)
    savefig(combined_plot, joinpath(draft_fig_dir, "manifold_geometry_plots$(day_suffix).png"))
    savefig(combined_plot, joinpath(draft_fig_dir, "manifold_geometry_plots$(day_suffix).pdf"))

    # Execute dynamic structural diagnostic loops
    generate_tier_plots(output_dir, draft_fig_dir, day_suffix)

    # --- STEP 3: Markdown Generation Layer ---
    report_path = joinpath(output_dir, "manifold_summary_report$(day_suffix).md")
    cond_number = cond(ws.Manifold_Mass)
    dz_vector = [ws.z_atm[i] - ws.z_atm[i+1] for i in 1:ws.N]
    min_dz = minimum(dz_vector)
    max_dz = maximum(dz_vector)

    report_content = raw"# Comprehensive Manifold Diagnostic Report" # ... (Keep Markdown content identical as previous setup)

    open(report_path, "w") do io
        write(io, report_content)
    end
    println("✓ Summary Markdown report written to: ", report_path)
end

run_diagnostic_pipeline("reports/ncar_eol_dee0099881")