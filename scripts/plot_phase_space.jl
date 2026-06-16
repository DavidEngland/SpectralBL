# scripts/plot_phase_space.jl
using CSV, DataFrames, Plots

df = CSV.read("data/diagnostic_trajectory.csv", DataFrame)
# Filter out non-finite artifacts
clean_df = filter(row -> !isnan(row.D_eff) && !isnan(row.F_W) && row.Regime in [1, 2, 3], df)

# Allocate categorical colors based on your newly calibrated thresholds
c_map = [row.Regime == 1 ? :red : (row.Regime == 2 ? :blue : :gold) for row in eachrow(clean_df)]
labels = [row.Regime == 1 ? "Turbulent" : (row.Regime == 2 ? "Wave-Dominated" : "Intermittent") for row in eachrow(clean_df)]

scatter(clean_df.F_W, clean_df.D_eff, group=labels,
        color=c_map, alpha=0.4, markersize=3, markerstrokewidth=0,
        xlabel="Adaptive Wave Fraction (F_W)",
        ylabel="Effective Modal Dimension (D_eff)",
        title="CASES-99 Atmospheric Regime Phase-Space Attractors",
        legend=:topright, size=(800, 600))

savefig("reports/ncar_eol_dee0099881/regime_attractor_scatter.png")