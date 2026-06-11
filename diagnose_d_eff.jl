using LinearAlgebra

include("src/UnifiedManifold.jl")
include("src/SpectralDiagnostics.jl")

using .UnifiedManifold
using .SpectralDiagnostics

# Setup workspace
ws = UnifiedManifoldWorkspace(32, 0.0, 60.0, 0.15)

# Create a test profile (wave-dominated) in full (N+1) modal vectors
c_theta_test = zeros(33)
c_theta_test[1:4] .= [5.0, 2.0, 0.8, 0.3]
c_u_test = zeros(33)
c_u_test[1:3] .= [1.0, 0.5, 0.2]

# Compute using NEW code path
record = process_timestamp_metrics(1, c_theta_test, c_u_test, ws, "TEST")

println("=== NEW D_eff from process_timestamp_metrics ===")
println("D_eff = ", record.D_eff)
println("F_W = ", record.F_W)
println("chi_N = ", record.chi_N)
println("Ri_g = ", record.Ri_g)
println("")

# Reconstruct OLD D_eff using full modal coefficients
energies_old = c_theta_test .^ 2
sum_e_old = sum(energies_old)
p_old = energies_old ./ sum_e_old
entropy_old = sum(-p_old[p_old .> 0.0] .* log.(p_old[p_old .> 0.0]))
D_eff_old = exp(entropy_old)

println("=== OLD D_eff (reconstructed) ===")
println("D_eff = ", D_eff_old)
println("")
println("Difference: ", abs(record.D_eff - D_eff_old))
