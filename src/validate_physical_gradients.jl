"""
    validate_physical_gradients(heights::Vector{Float64}, theta::Vector{Float64})

Checks SBL physical profiles for strong vertical anomalies. Flags extreme super-adiabatic
layers and warns if an intense surface inversion layer (> 5 K/m) could challenge the grid.
"""
function validate_physical_gradients(heights::Vector{Float64}, theta::Vector{Float64})
    @assert length(heights) == length(theta) "Height and Theta vectors must match in size."

    passed = true
    for i in 1:(length(heights)-1)
        dz = heights[i+1] - heights[i]
        dtheta = theta[i+1] - theta[i]
        gradient = dtheta / dz

        # 1. Super-adiabatic threshold violation check for nocturnal SBL
        if gradient < -0.05
            println("  ❌ [PHYSICAL VIOLATION] Unphysical super-adiabatic drop at z = $(heights[i])m to $(heights[i+1])m (G: $(round(gradient, digits=3)) K/m).")
            passed = false
        end

        # 2. Extreme inversion warning (Useful for checking alpha_stretch performance)
        if gradient > 5.0
            println("  💡 [STRUCTURE WARNING] Extreme thermal gradient detected at z = $(heights[i])m ($(round(gradient, digits=2)) K/m). Ensure compactification is active.")
        end
    end
    return passed
end