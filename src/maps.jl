# src/maps.jl

# 1. Explicitly include the local source file path
include("transforms.jl")

# 2. Use the compiled module
using .Transforms
import Plots
using Plots.Measures  # <--- CRITICAL: Enables unit measurements like mm, cm, inch

# Setup grid matching your spectral execution
N = 64
ξ = cos.(range(0, stop=pi, length=N)) # Chebyshev Gauss-Lobatto nodes

# Instantiate the structures
m_tanh = TanhMap(0.0, 1000.0, 2.5)
m_hyp  = HyperbolicMap(0.0, 1000.0, 3.0)

# Compute inverted physical node heights
z_tanh = inverse.(Ref(m_tanh), ξ)
z_hyp  = inverse.(Ref(m_hyp), ξ)

# Compute the Metric Jacobians
dzdξ_tanh = dzdξ.(Ref(m_tanh), ξ)
dzdξ_hyp  = dzdξ.(Ref(m_hyp), ξ)

# --- Plot 1: Resolution Mapping Distribution ---
# Added left_margin and bottom_margin padding to prevent label clipping
p1 = Plots.plot(ξ, z_tanh, label="TanhMap", color=:blue, lw=2,
            left_margin=10mm, bottom_margin=8mm)
Plots.plot!(p1, ξ, z_hyp, label="HyperbolicMap", color=:red, lw=2,
            xlabel="Computational Coordinate (ξ)",
            ylabel="Physical Height z (m)",
            title="Geometry Manifold Profiles",
            legend=:topleft)

# --- Plot 2: Metric Jacobian Weighting ---
# Added left_margin, bottom_margin, and right_margin for balanced side-by-side padding
p2 = Plots.plot(z_tanh, dzdξ_tanh, label="TanhMap", color=:blue, lw=2,
            left_margin=10mm, bottom_margin=8mm, right_margin=5mm)
Plots.plot!(p2, z_hyp, dzdξ_hyp, label="HyperbolicMap", color=:red, lw=2,
            xlabel="Physical Height z (m)",
            ylabel="Jacobian Metric (dz/dξ)",
            title="Resolution Density Distribution",
            legend=:topright)

# Combine them side-by-side with an all-around safety margin buffer
manifold_figure = Plots.plot(p1, p2, layout=(1, 2), size=(1000, 450), margin=5mm)

# Save directly to your workspace
Plots.savefig(manifold_figure, "manifold_comparison.pdf")
println("Successfully generated un-chopped manifold_comparison.pdf")