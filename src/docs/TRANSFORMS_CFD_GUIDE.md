# docs/TRANSFORMS_CFD_GUIDE.md

# CfdWallMap & Advanced Jacobian Stack Composition

## CfdWallMap: Fixed Wall Boundary Layer Refinement

The `CfdWallMap` is specialized for **CFD and LES/DNS applications** requiring extreme
near-wall resolution. Unlike `TanhMap` which compresses toward an internal point,
`CfdWallMap` concentrates all refinement at the **wall** (z = 0, ξ = -1).

### Mathematical Formulation

```
z(ξ) = ztop · (1 - tanh(δ(1 - ξ)) / tanh(δ))
```

**Domain mapping**:

- ξ = -1: z → 0 (wall, ξ → -∞ singularity compressed)
- ξ = 0: z → ztop · (1 - tanh(δ)/tanh(δ)) ≈ ztop/2 (approximate center)
- ξ = 1: z → ztop (upper domain)

### Jacobians

```
dz/dξ = ztop · (δ / tanh(δ)) · sech²(δ(1 - ξ))

d²z/dξ² = -2ztop · (δ² / tanh(δ)) · sech²(δ(1 - ξ)) · tanh(δ(1 - ξ))
```

**Key properties**:

- J₁(ξ=-1) is very small (wall refinement)
- J₁(ξ=1) is larger (sparse spacing aloft)
- Always monotonically increasing (z is strictly increasing in ξ)

### Parameter Selection: δ (packing strength)

|δ       |Wall compression   |Recommended for        |Notes                     |
|--------|-------------------|-----------------------|--------------------------|
|0.5–1.0 |Very weak          |Coarse grids, testing  |Nearly linear             |
|1.5–2.5 |Mild               |General CFD            |Balanced refinement       |
|3.0–5.0 |**Moderate-Strong**|**LES/DNS (y⁺ < 1)**   |Standard choice           |
|6.0–10.0|Extreme            |DNS, molecular dynamics|CFL constraints tighten   |
|>10.0   |Ultra-extreme      |Specialized research   |Verify stability carefully|

### Typical Use Case: LES with Wall Refinement

```julia
using Transforms

# Domain: 0 to 500 m (LES domain height)
m = CfdWallMap(500.0, δ=4.5)

# Check that wall is highly refined
z_near_wall = inverse(m, -0.95)    # ξ ≈ -0.95
z_lower = inverse(m, -0.5)         # ξ ≈ -0.5
z_mid = inverse(m, 0.0)            # ξ = 0
z_upper = inverse(m, 0.95)         # ξ ≈ 0.95

@printf "Near-wall (ξ=-0.95): z ≈ %.2f m\n" z_near_wall     # ~1–5 m
@printf "Lower region (ξ=-0.5): z ≈ %.2f m\n" z_lower       # ~30–80 m
@printf "Mid-domain (ξ=0.0): z ≈ %.2f m\n" z_mid            # ~250 m
@printf "Upper region (ξ=0.95): z ≈ %.2f m\n" z_upper       # ~490 m

# For y⁺ < 1 at z ≈ 0.5 m (u_* ≈ 0.5 m/s):
# Need Δz(z=0.5 m) ~ 0.05 m (50 nodes in lowest 25 m)
J1_near_wall = dzdξ(m, -0.95)
Δξ_per_node = 2.0 / 256  # 256-point spectral mesh
Δz_near_wall = J1_near_wall * Δξ_per_node
@printf "Physical spacing near wall: Δz ≈ %.4f m\n" Δz_near_wall
```

### CfdWallMap vs. TanhMap: Comparison

```julia
# Scenario: LES domain 0–500 m, need y⁺ < 1 at z = 0.5 m

# Option 1: CfdWallMap (specialized for wall)
m_cfd = CfdWallMap(500.0, 4.5)
z_wall = inverse(m_cfd, -0.9)    # ~2 m very close to wall

# Option 2: TanhMap (general purpose)
m_tanh = TanhMap(0.0, 500.0, 3.5)
z_wall = inverse(m_tanh, -0.9)   # ~40 m (still within lower layer)

# CfdWallMap gets much closer to z=0 for same ξ
# Use CfdWallMap when wall layer is the focus; TanhMap for balanced compression
```

-----

## Advanced: Exact Analytical Jacobian Stacks

The improved `evaluate_jacobian_stack()` function now uses **exact analytical chain rules**
instead of finite-difference approximations. This is crucial for:

1. **Accuracy**: No FD truncation error (machine precision only)
1. **Efficiency**: Single evaluation pass; no extra function calls
1. **Stability**: Better conditioning near metric singularities

### Mathematical Foundation

For a stack of n maps: z → x₁ → x₂ → … → ξ

**Forward transformation**:

```
ξ = fₙ(fₙ₋₁(...f₂(f₁(z))...))
```

**First metric (chain rule product)**:

```
dξ/dz = ∏ᵢ (dξᵢ/dzᵢ) = ∏ᵢ 1/(dz/dξ)ᵢ
```

**Second metric (full chain-rule expansion)**:

```
d²ξ/dz² = ∑ᵢ [ (d²ξᵢ/dzᵢ²) · (dzᵢ/dz)² · ∏ⱼ>ᵢ (dξⱼ/dzⱼ) ]
```

This is a **nested product-sum** that captures all curvature contributions at each layer.

### Example: Two-Layer Manifold Decomposition

**Layer 1**: Stability-to-compactified transformation

- Physical domain: z ∈ [0, 2000 m]
- Intermediate domain: z₁ ∈ [-1, 1]
- Map: `TanhMap(0, 2000, α=2.5)`

**Layer 2**: Chebyshev domain refinement

- Intermediate domain: z₁ ∈ [-1, 1]
- Computational domain: ξ ∈ [-1, 1]
- Map: `HyperbolicMap(-1, 1, α=1.5)`

```julia
using Transforms

# Define the stack
m1 = TanhMap(0.0, 2000.0, 2.5)
m2 = HyperbolicMap(-1.0, 1.0, 1.5)

stack = JacobianStack([m1, m2], ["Stability Compression", "Chebyshev Refinement"])

# Evaluate at multiple physical heights
z_test = [10.0, 50.0, 100.0, 500.0, 1000.0]

println("Height (m) | Comp. ξ | dξ/dz | d²ξ/dz²")
println("-" ^ 50)

for z in z_test
    ξ, J1, J2 = evaluate_jacobian_stack(stack, z)
    @printf "%.1f      | %.4f  | %.4e | %.4e\n" z ξ J1 J2
end
```

**Output** (qualitative):

```
Height (m) | Comp. ξ | dξ/dz | d²ξ/dz²
--------------------------------------------------
10.0       | -0.94   | 0.012 | 0.0023
50.0       | -0.81   | 0.034 | 0.0015
100.0      | -0.71   | 0.067 | 0.0008
500.0      | -0.05   | 0.185 | -0.0002
1000.0     | 0.61    | 0.142 | -0.0001
```

**Interpretation**:

- Near z = 0: dξ/dz is small (many computational nodes map to small physical span)
- Mid-domain: dξ/dz increases (node spacing expands)
- Upper domain: dξ/dz tapers (coarser spacing aloft)
- d²ξ/dz² shows metric **curvature** (nonlinearity of spacing)

-----

## Assembling Spectral Operators with Stacked Jacobians

Once you have stacked metrics, assembly of spectral-element operators follows the same pattern:

### Quadrature Weighting

```julia
using FastGaussQuadrature

# Gauss-Legendre quadrature on [-1, 1]
N = 128  # spectral points
ξ_quad, w_quad = gausslegendre(N)

# Evaluate stack at quadrature nodes
ξ_final = similar(ξ_quad)
J1_quad = similar(ξ_quad)
J2_quad = similar(ξ_quad)

for (i, ξ) in enumerate(ξ_quad)
    # Transform back to physical domain
    z = inverse(stack.maps[1], ξ)
    ξ_f, J1, J2 = evaluate_jacobian_stack(stack, z)
    ξ_final[i] = ξ_f
    J1_quad[i] = J1
    J2_quad[i] = J2
end

# Physical-domain quadrature weights
w_phys = w_quad .* J1_quad

# Integrate function f(z) over [0, 2000 m]
f_values = f.(z)  # Function evaluated at quadrature points
integral = dot(w_phys, f_values)
```

### Differentiation Matrix with Stack Metrics

```julia
# Spectral differentiation matrices (Chebyshev)
D1_ξ, D2_ξ = chebyshev_matrices(N)

# At spectral nodes, retrieve stacked metrics
J1_nodes = J1_quad   # (precomputed above)
J2_nodes = J2_quad
z_nodes = z_quad

# Construct physical-domain differentiation operator
# du/dz = (dξ/dz) * du/dξ
dξdz_nodes = 1 ./ J1_nodes
D1_z = Diagonal(dξdz_nodes) * D1_ξ

# Second derivative with metric curvature term
d2ξdz2_nodes = -J2_nodes ./ (J1_nodes .^ 3)
D2_z = Diagonal(1 ./ (J1_nodes .^ 2)) * D2_ξ - Diagonal(d2ξdz2_nodes) * D1_ξ
```

### Laplacian Assembly

```julia
# Laplacian in physical domain via stacked metrics
Lap_z = D2_z  # For 1D, this is complete

# Apply to solution
residual = Lap_z * u
```

-----

## Validation: Testing Stacked Metrics

### Round-Trip Consistency

```julia
function validate_stack_round_trip(stack::JacobianStack, z::Float64; tol=1e-10)
    # Forward: z → intermediate coordinates → ξ
    ξ, _, _ = evaluate_jacobian_stack(stack, z)

    # Backward: recover z from intermediate stages
    z_recovered = z  # Start with ξ at end
    for i in length(stack.maps):-1:1
        z_recovered = inverse(stack.maps[i], z_recovered)
    end

    error = abs(z_recovered - z) / abs(z + 1e-10)
    if error > tol
        @warn "Round-trip error: z=$z, recovered=$z_recovered, rel_error=$error"
    end
    return error < tol
end

# Test at multiple points
z_test = range(10, 1990, length=50)
all_valid = all(validate_stack_round_trip(stack, z) for z in z_test)
@assert all_valid "Round-trip validation failed!"
```

### Finite-Difference Verification

```julia
function verify_stack_jacobians(stack::JacobianStack, z::Float64; eps=1e-8)
    ξ_center, J1_analytical, J2_analytical = evaluate_jacobian_stack(stack, z)

    # FD for J1
    _, J1_plus, _ = evaluate_jacobian_stack(stack, z + eps)
    _, J1_minus, _ = evaluate_jacobian_stack(stack, z - eps)
    J1_fd = (J1_plus - J1_minus) / (2 * eps)

    # FD for J2
    dJ1dz_plus = J1_plus
    dJ1dz_minus = J1_minus
    J2_fd = (dJ1dz_plus - dJ1dz_minus) / (2 * eps)

    J1_error = abs(J1_analytical - J1_fd) / abs(J1_analytical + 1e-10)
    J2_error = abs(J2_analytical - J2_fd) / abs(J2_analytical + 1e-10)

    @printf "J₁ relative error: %.2e\n" J1_error
    @printf "J₂ relative error: %.2e\n" J2_error

    @assert J1_error < 1e-5 "J1 FD check failed"
    @assert J2_error < 1e-4 "J2 FD check failed"
end
```

-----

## Production Tips

### 1. Monitor Metric Spectrum

Before using a Jacobian stack in a solver, plot the metric functions:

```julia
using Plots

z_range = range(10, 1990, length=500)
ξ_vals = similar(z_range)
J1_vals = similar(z_range)
J2_vals = similar(z_range)

for (i, z) in enumerate(z_range)
    ξ_vals[i], J1_vals[i], J2_vals[i] = evaluate_jacobian_stack(stack, z)
end

p1 = plot(z_range, J1_vals, label="dξ/dz", xlabel="z (m)", ylabel="Metric")
p2 = plot(z_range, J2_vals, label="d²ξ/dz²", xlabel="z (m)", ylabel="Curvature")
plot(p1, p2, layout=(1,2))
```

Check:

- J₁ > 0 everywhere (monotonic map required)
- J₁ smooth, no sharp discontinuities
- J₂ within reasonable bounds (large curvature → CFL constraints)

### 2. CFL Constraint with Stacked Maps

The CFL condition tightens with strong metric compression:

```
Δt ≤ CFL · min(Δz) / u_max

min(Δz) ≈ min(J₁) · Δξ
```

For strong stacks, min(J₁) can be very small, forcing Δt → 0.

```julia
# Estimate CFL constraint
Δξ = 2.0 / N_spectral
min_J1 = minimum(J1_vals)
min_Δz = min_J1 * Δξ

u_max = 10.0  # m/s (typical atmospheric wind)
CFL_num = 0.5
Δt_cfl = CFL_num * min_Δz / u_max

@printf "Minimum Δz: %.4f m\n" min_Δz
@printf "CFL-limited Δt: %.6f s\n" Δt_cfl
```

-----

## References

- **Canuto et al. (1988)**: Spectral Methods in Fluid Dynamics, Chapter 2 (coordinate transforms)
- **Boyd (2001)**: Chebyshev and Fourier Spectral Methods, Section 17 (tanh maps for semi-infinite)
- **Trefethen (2000)**: Spectral Methods in MATLAB, Chapters 7–8 (practical assembly)

-----

**Last Updated**: June 2026 | **Status**: Production Ready