using DataFrames
using CSV
using Plots
using LinearAlgebra

# Assuming UnifiedManifold is loaded
# import .UnifiedManifold: UnifiedManifoldWorkspace, physical_to_computational

function run_diagnostic_pipeline(output_dir::String)
    # 1. Ensure the directory exists
    mkpath(output_dir)
    println("Target directory verified: ", output_dir)

    # 2. Instantiate the Workspace with the corrected alpha_stretch (0.05)
    # Standard CASES-99 setup: N=32, z_0m=1.5m, z_top=50.0m
    N = 32
    z_0m = 1.5
    z_top = 50.0
    alpha_stretch = 0.05

    ws = UnifiedManifoldWorkspace(N, z_0m, z_top, alpha_stretch)

    # --- STEP 1: Generate Diagnostics CSV ---
    csv_path = joinpath(output_dir, "manifold_diagnostics.csv")
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

    # --- STEP 2: Generate Plots ---
    # Plot 1: Grid Distribution (Hyperbolic Compactification Profile)
    p1 = plot(ws.xi_target, ws.z_atm, marker=:circle, linewidth=2,
              title="Hyperbolic Mapping (α = $alpha_stretch)",
              xlabel="Computational Coordinate (ξ)", ylabel="Physical Height z (m)",
              label="Grid Nodes", legend=:topleft, theme=:vibrant)

    # Plot 2: Sub-meso Partitioning Windows
    p2 = plot(0:ws.N, [ws.psi_M ws.psi_W ws.psi_T], linewidth=2.5,
              title="Spectral Partitioning Windows",
              xlabel="Chebyshev Mode Index (n)", ylabel="Filter Weight (ψ)",
              label=["Meso (ψ_M)" "Wave (ψ_W)" "Turb (ψ_T)"], legend=:topright)

    plot_path = joinpath(output_dir, "manifold_geometry_plots.png")
    combined_plot = plot(p1, p2, layout=(1, 2), size=(1000, 450))
    png(combined_plot, plot_path)
    println("✓ Diagnostic plots saved to: ", plot_path)

    # --- STEP 3: Generate Summary Markdown Report ---
    report_path = joinpath(output_dir, "manifold_summary_report.md")

    # Calculate some quick structural metrics
    cond_number = cond(ws.Manifold_Mass)
    min_dz = minimum([ws.z_atm[i+1] - ws.z_atm[i] for i in 1:ws.N])
    max_dz = maximum([ws.z_atm[i+1] - ws.z_atm[i] for i in 1:ws.N])

    report_content = """
    # Unified Manifold Workspace Summary Report
    **Generated Target Directory:** `$output_dir`

    ## Geometry & Mapping Parameters
    * **Spectral Modes (N):** $N$
    * **Lower Boundary (z_0m):** $(z_0m) m (CASES-99 Tower Base)
    * **Top Boundary (z_top):** $(z_top) m
    * **Stretch Intensity (alpha_stretch):** $alpha_stretch
    * **Quadrature Nodes (K_q):** $(ws.K_q)

    ## Numerical Health Diagnostics
    * **Mass Matrix Condition Number:** $(round(cond_number, digits=4))
    * **Minimum Grid Spacing (Δz_min):** $(round(min_dz, digits=4)) m (near surface)
    * **Maximum Grid Spacing (Δz_max):** $(round(max_dz, digits=4)) m (near canopy top)

    ## Partitioning Thresholds
    * **Meso Windows (n_m):** Modes capturing stable, large-scale structures.
    * **Wave/Sub-meso (n_w):** Internal gravity wave transitions.
    * **Turbulent Residual:** High-frequency dissipation modes.

    *Report generated automatically by the UnifiedManifold pipeline.*
    """

    open(report_path, "w") do io
        write(io, report_content)
    end
    println("✓ Summary Markdown report written to: ", report_path)
end

# Execute the pipeline for your specific directory
run_diagnostic_pipeline("reports/ncar_eol_dee0099881")