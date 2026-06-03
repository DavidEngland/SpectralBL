using LinearAlgebra

"""
    audit_spectral_conditioning(heights::Vector{Float64}, N::Int, α_stretch::Float64)

Constructs a test observation matrix H and audits its condition number.
Returns a boolean indicating if the grid layout is numerically stable.
"""
function audit_spectral_conditioning(heights::Vector{Float64}, N::Int, α_stretch::Float64)
    if length(heights) < 2
        println("--- [NUMERICAL AUDIT] ---")
        println("  [CRITICAL] Need at least two vertical levels for conditioning audit.")
        return false
    end

    # Replicate paper's coordinate compactification mapping physical z to computational ξ
    z_min, z_max = minimum(heights), maximum(heights)
    if z_max == z_min
        println("--- [NUMERICAL AUDIT] ---")
        println("  [CRITICAL] Degenerate grid detected: all heights are identical.")
        return false
    end

    σ = (z_max - z_min) * α_stretch / 2.0

    # Map heights to ξ ∈ [-1, 1]
    ξ = map(z -> 1.0 - (2.0 * σ) / (z - z_min + σ / α_stretch), heights)
    if any(x -> !isfinite(x), ξ)
        println("--- [NUMERICAL AUDIT] ---")
        println("  [CRITICAL] Non-finite compactified coordinates detected.")
        return false
    end

    # Construct evaluation matrix H (M observations x N+1 modes)
    M = length(heights)
    H = zeros(Float64, M, N + 1)

    for i in 1:M
        for j in 1:(N + 1)
            # Standard Chebyshev polynomial evaluation via recurrence relation
            n = j - 1
            if n == 0
                H[i, j] = 1.0
            elseif n == 1
                H[i, j] = ξ[i]
            else
                # T_n(x) = 2x * T_{n-1}(x) - T_{n-2}(x)
                t_m2 = 1.0
                t_m1 = ξ[i]
                val = 0.0
                for k in 2:n
                    val = 2.0 * ξ[i] * t_m1 - t_m2
                    t_m2 = t_m1
                    t_m1 = val
                end
                H[i, j] = val
            end
        end
    end

    # Compute condition number using Singular Values
    s = svdvals(H)
    cond_H = s[1] / s[end]

    println("--- [NUMERICAL AUDIT] ---")
    println("  Matrix H Grid Dimension : $(M)x$(N+1)")
    println("  Calculated Condition No : $(cond_H)")

    # Arbitrary regularization ceiling (Adjust based on double precision limits)
    if !isfinite(cond_H) || s[end] == 0.0 || cond_H > 1e12
        println("  ❌ [CRITICAL] Grid layout is severely ill-conditioned. Spectral tracking will fail.")
        return false
    else
        println("  ✅ [PASS] Matrix condition number is safe for SVD rank-truncation mapping.")
        return true
    end
end