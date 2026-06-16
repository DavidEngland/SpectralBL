# scripts/plot_iop_case_study.jl
using CSV, DataFrames, Plots

df = CSV.read("data/diagnostic_trajectory.csv", DataFrame)
# Isolate a pristine IOP day
day_df = filter(row -> row.FileDate == 991026, df)

p1 = plot(day_df.TimeIdx, day_df.D_eff, color=:blue, linewidth=2, ylabel="D_eff", label="Modal Dimension")
hline!([8.0], color=:black, linestyle=:dash, label="Wave Collapse Threshold")

p2 = plot(day_df.TimeIdx, day_df.Ri_f, color=:green, linewidth=2, ylabel="Ri_f", label="Flux Richardson")
hline!([0.25], color=:red, linestyle=:dot, label="Critical Ri")

plot(p1, p2, layout=(2,1), xlabel="Time Step (Diurnal Sequence)",
     title="Structural Manifold Response During IOP Inversion (Oct 26)", size=(900, 500))

savefig("reports/ncar_eol_dee0099881/case_study_oct26.pdf")