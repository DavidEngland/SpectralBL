# docs/TRANSFORMS_GUIDE.md

# Transforms Module: Usage Guide & Examples

## Overview

The `Transforms` module implements five coordinate map strategies for transforming physical
height profiles (z ∈ [0, z_max]) into a bounded computational domain (ξ ∈ [-1, 1]). This is
essential for:

1. **Spectral methods**: Chebyshev polynomials are orthogonal on [-1, 1]; physical domains must be mapped.
1. **Metric consistency**: Exact Jacobians (dz/dξ, d²z/dξ²) enable correct differentiation & integration.
1. **Adaptive refinement**: Different maps compress degrees-of-freedom toward regions of interest.

## Quick Start

### Example 1: Simple Linear Transformation

```julia
using Transforms

# Define domain: surface (z=0) to 2 km altitude
m = LinearMap(0.0, 2000.0)

# Transform a height
z = 500.0
ξ = forward(m, z)        # ξ ≈ -0.5 (computational domain)
z_back = inverse(m, ξ)   # z ≈ 500.0 (round-trip)

# Get metric Jacobians
J1 = dzdξ(m, ξ)          # ≈ 1000.0 (constant for linear map)
J2 = d2zdξ2(m, ξ)        # = 0.0 (no curvature)
```

### Example 2: Production TanhMap with Profile

```julia
# Moderate compression toward surface (boundary layer focus)
m = TanhMap(0.0, 2000.0, 2.5)

# Observational heights (e.g., CASES-99 tower)
z_obs = [5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0, 1000.0, 2000.0]

# Transform to computational domain
ξ_obs = profile_transform(m, z_obs)

# Inspect metric behavior
for (z, ξ) in zip(z_obs, ξ_obs)
    J1 = dzdξ(m, ξ)
    J2 = d2zdξ2(m, ξ)
    println("z=$z m → ξ=$ξ, J₁=$J1, J₂=$J2")
end
```

**Output** (qualitative):

```
z=5.0 m → ξ≈-0.90, J₁≈80, J₂≈0.1      [compressed region]
z=500.0 m → ξ≈0.0, J₁≈480, J₂≈-0.05   [center, reduced refinement]
z=2000.0 m → ξ≈1.0, J₁≈320, J₂≈0.02   [upper boundary]
```

**Interpretation**: Near z=0, the metric J₁ is small (many computational nodes compressed into
small physical space). Higher up, J₁ expands, spacing nodes further apart. The second derivative
J₂ indicates metric curvature (nonlinearity of stretching).

-----

## Map Comparison

|Map               |Purpose                   |Strength                 |Caution                         |
|------------------|--------------------------|-------------------------|--------------------------------|
|**LinearMap**     |Testing, verification     |Simple, predictable      |No refinement                   |
|**HyperbolicMap** |Extreme boundary focus    |O(1/n²) convergence      |Can cause ill-conditioning      |
|**LogarithmicMap**|Exponential stratification|Natural for stable layers|Requires z > 0; unbounded growth|
|**TanhMap**       |General production use    |Smooth, stable, flexible |Choose α carefully              |
|**CustomMap**     |Specialized domains       |Maximum control          |User validates metric           |

-----

## Domain Constraints & Parameter Selection

### LinearMap & HyperbolicMap & TanhMap

**Domain**: zmin ≥ 0, zmin < zmax

**When to use**:

- LinearMap: Unit tests, method verification, debugging
- HyperbolicMap/TanhMap: Physical simulations where lower atmosphere requires high resolution

**Example with realistic SBL**:

```julia
# Nocturnal stable boundary layer (SBL): 0–500 m
m = TanhMap(0.0, 500.0, α=3.5)

# This allocates ~60% of spectral points to z < 100 m (inversion layer)
# while covering the entire domain [0, 500 m]
```

### LogarithmicMap

**Domain**: zmin > 0, zmin < zmax (strictly positive heights)

**When to use**:

- Exponentially decaying variance profiles (e.g., SMEAR-II forest canopy)
- Thick stratification with weak surface gradients
- Tall domains (e.g., 10 m to 2 km) where log-scaling is physically natural

**Parameter guidance**:

```julia
# Small domain (10 m to 100 m) — use TanhMap
m1 = TanhMap(10.0, 100.0, 2.5)

# Medium domain (10 m to 1000 m) — LogarithmicMap competitive
m2 = LogarithmicMap(10.0, 1000.0)

# Tall domain (10 m to 2000 m) — TanhMap preferred; LogarithmicMap strained
m3 = TanhMap(10.0, 2000.0, 3.0)
```

### TanhMap: Choosing α

The parameter α controls compression strength. Guidance:

|α      |Compression|Recommendation                        |
|-------|-----------|--------------------------------------|
|0.5–1.0|Very weak  |Exploratory; near-linear behavior     |
|1.5–2.5|Moderate   |**Default choice**; stable, robust    |
|3.0–4.0|Strong     |High-resolution lower layer; CFL tight|
|5.0+   |Extreme    |Specialty use; verify stability       |

**Automatic selection heuristic**:

```julia
# If lower 10% of domain contains most energy variance, use stronger α
z_critical = 0.1 * (zmax - zmin)  # e.g., 20 m in 200 m domain

# Estimate variance concentration
if energy_fraction_below(z_critical) > 0.6
    α = 3.5  # Strong compression
else
    α = 2.0  # Moderate (default)
end

m = TanhMap(zmin, zmax, α)
```

-----

## Validation & Robustness

### Checking Map Validity

All maps perform domain validation:

```julia
using Transforms

try
    m = LogarithmicMap(0.0, 1000.0)  # ERROR: zmin = 0 (log singularity)
catch e
    println("Caught error: $e")
end

try
    m_ok = LogarithmicMap(10.0, 1000.0)  # OK: zmin > 0
    is_valid(m_ok)  # Returns true
catch e
    println("Unexpected error: $e")
end
```

### Round-Trip Consistency Test

Always verify forward–inverse composition:

```julia
m = TanhMap(0.0, 1000.0, 2.5)

z_test = 250.0
ξ = forward(m, z_test)
z_recovered = inverse(m, ξ)

error = abs(z_recovered - z_test)
println("Round-trip error: $error (should be < 1e-14)")

# For profiles
z_profile = range(0.0, 1000.0, length=50)
ξ_profile = profile_transform(m, z_profile)
z_profile_recovered = inverse.(Ref(m), ξ_profile)
max_error = maximum(abs.(z_profile_recovered .- z_profile))
println("Max profile error: $max_error")
```

### Metric Jacobian Verification

Check dz/dξ via finite differences (expected error ~ 1e-6 for ε = 1e-8):

```julia
function verify_dzdxi(m::CoordinateMap, ξ_test; eps=1e-8)
    # Analytical
    dz_analytical = dzdξ(m, ξ_test)

    # Finite difference
    z_plus = inverse(m, ξ_test + eps)
    z_minus = inverse(m, ξ_test - eps)
    dz_fd = (z_plus - z_minus) / (2 * eps)

    error = abs(dz_analytical - dz_fd) / abs(dz_analytical)
    println("J₁ relative error: $error (analytical vs. FD)")
    return error < 1e-5
end

m = TanhMap(0.0, 1000.0, 2.5)
verify_dzdxi(m, 0.0)   # Check at center
verify_dzdxi(m, 0.5)   # Check at boundary region
```

-----

## Integration with Spectral Methods

### Quadrature Weighting

When integrating a function f(z) over physical domain [zmin, zmax]:

```
∫_zmin^zmax f(z) dz = ∫_{-1}^{1} f(z(ξ)) * (dz/dξ) dξ
```

The metric Jacobian **dz/dξ** is the integration weight. Implementation:

```julia
using FastGaussQuadrature  # or similar quadrature library

m = TanhMap(0.0, 1000.0, 2.5)
N = 50  # Number of spectral points

# Get Gauss-Legendre quadrature nodes & weights on [-1, 1]
ξ_quad, w_quad = gausslegendre(N)

# Map to physical domain
z_quad = inverse.(Ref(m), ξ_quad)

# Apply metric weighting
dz_quad = dzdξ.(Ref(m), ξ_quad)
w_phys = w_quad .* dz_quad  # Physical-domain weights

# Integrate f(z) over [zmin, zmax]
f_values = f.(z_quad)
integral = dot(w_phys, f_values)
```

### Differentiation in Mapped Coordinates

To compute df/dz using spectral derivatives:

```
d/dz = (dξ/dz) * d/dξ = (1 / (dz/dξ)) * d/dξ
```

Implementation with Chebyshev differentiation matrix D_ξ:

```julia
# Compute du/dξ (spectral derivative matrix applied to coefficients)
du_dxi = D_ξ * u_coeff

# Scale to physical domain
dξdz_vals = dξdz.(Ref(m), ξ_nodes)
du_dz = dξdz_vals .* du_dxi
```

### Laplacian in Mapped Coordinates

The Laplacian transforms as:

```
d²u/dz² = (1/(dz/dξ)²) * d²u/dξ² - (d²z/dξ²)/(dz/dξ)³ * du/dξ
```

**Operator assembly**:

```julia
# First and second derivatives (spectral matrices)
D1_ξ, D2_ξ = chebyshev_matrices(N)  # [user-supplied]

# Metric arrays at nodes
J1 = dzdξ.(Ref(m), ξ_nodes)
J2 = d2zdξ2.(Ref(m), ξ_nodes)

# Laplacian operator in physical domain
# L_z = diag(1/J1²) * D2_ξ - diag(J2/J1³) * D1_ξ
L_z = Diagonal(1 ./ (J1 .^ 2)) * D2_ξ - Diagonal(J2 ./ (J1 .^ 3)) * D1_ξ

# Apply to solution vector
residual = L_z * u
```

-----

## Advanced: Jacobian Stack (Multi-Layer Mapping)

For the full manifold decomposition (Tier 2 theory), a Jacobian stack composes
multiple transformations:

```
z_physical → z_intermediate1 → z_intermediate2 → ξ_computational
```

**Use case**: Stability-to-compactified + Chebyshev domain mapping.

```julia
# Layer 1: Compress stable stratification
m1 = TanhMap(0.0, 2000.0, 2.5)

# Layer 2: Further refinement in Chebyshev domain (maps [-1,1] → [-1,1])
# (Example: emphasize wavy regions)
m2 = HyperbolicMap(-1.0, 1.0, 1.5)

# Compose stack
stack = JacobianStack([m1, m2], ["Stability", "Chebyshev"])

# Evaluate at physical height
z = 150.0
ξ, J1, J2 = evaluate_jacobian_stack(stack, z)

println("z=$z m → ξ=$ξ (computational)")
println("Composite metric: J₁=$J1, J₂=$J2")
```

**Caution**: Jacobian stacks require careful validation. The chain-rule composition
can amplify metric curvature (J₂) dramatically. Check spectral conditioning before deployment.

-----

## Testing & Troubleshooting

### Common Errors

**Error 1: DomainError for LogarithmicMap**

```julia
m = LogarithmicMap(0.0, 1000.0)  # ERROR: zmin ≤ 0
```

**Fix**: Use zmin > 0, e.g., `LogarithmicMap(1.0, 1000.0)`.

**Error 2: atanh domain overflow in TanhMap.forward**

```julia
m = TanhMap(0.0, 1000.0, 10.0)  # Very strong α
z = 10.0  # Near boundary
ξ = forward(m, z)  # Possible numerical issue
```

**Fix**: The forward transform includes clamping to prevent overflow. If you see
warnings, reduce α or increase domain size.

**Error 3: Singular Jacobian (dz/dξ ≈ 0)**

```julia
J1 = dzdξ(m, ξ)  # Returns very small number
J_inv = dξdz(m, z)  # Computes 1/J1; numerical division error
```

**Fix**: Check that the map is valid and monotonic. Plot dz/dξ(ξ) to visualize problem regions.

### Diagnostics

**Plot metric Jacobians**:

```julia
using Plots

m = TanhMap(0.0, 1000.0, 2.5)
ξ = range(-1, 1, length=200)
z = inverse.(Ref(m), ξ)
J1_vals = dzdξ.(Ref(m), ξ)
J2_vals = d2zdξ2.(Ref(m), ξ)

p1 = plot(z, J1_vals, xlabel="z (m)", ylabel="dz/dξ", title="First Metric")
p2 = plot(z, J2_vals, xlabel="z (m)", ylabel="d²z/dξ²", title="Second Metric Curvature")
plot(p1, p2)
```

**Check metric monotonicity**:

```julia
m = TanhMap(0.0, 1000.0, 2.5)
ξ = range(-1, 1, length=100)
J1 = dzdξ.(Ref(m), ξ)

# J1 must be strictly positive for monotonic z(ξ)
if all(J1 .> 0)
    println("✓ Monotonic map (good)")
else
    println("✗ Non-monotonic region detected!")
    idx = findfirst(J1 .<= 0)
    println("  Problem at ξ = $(ξ[idx])")
end
```

-----

## References & Further Reading

1. **Canuto et al. (1988)**: *Spectral Methods in Fluid Dynamics*. Chapter 2 covers coordinate transforms.
1. **Boyd (2001)**: *Chebyshev and Fourier Spectral Methods*. Excellent treatment of tanh maps (Section 17.4).
1. **Stull (1988)**: *An Introduction to Boundary Layer Meteorology*. Physical motivation for MOST + log profiles.
1. **Trefethen (2000)**: *Spectral Methods in MATLAB*. Practical implementation guidance.

-----

## Appendix: Test Suite

```julia
# test_transforms.jl
using Transforms
using Test

@testset "LinearMap" begin
    m = LinearMap(0.0, 1000.0)

    @test forward(m, 0.0) ≈ -1.0
    @test forward(m, 1000.0) ≈ 1.0
    @test forward(m, 500.0) ≈ 0.0

    @test inverse(m, -1.0) ≈ 0.0
    @test inverse(m, 1.0) ≈ 1000.0
    @test inverse(m, 0.0) ≈ 500.0

    @test dzdξ(m, 0.0) ≈ 500.0
    @test d2zdξ2(m, 0.0) ≈ 0.0
end

@testset "TanhMap" begin
    m = TanhMap(0.0, 2000.0, 2.5)

    # Boundary conditions
    @test forward(m, 0.0) ≈ -1.0 atol=1e-10
    @test forward(m, 2000.0) ≈ 1.0 atol=1e-10

    # Round-trip
    z_test = 500.0
    ξ_test = forward(m, z_test)
    z_recovered = inverse(m, ξ_test)
    @test z_recovered ≈ z_test atol=1e-10

    # Metric positivity
    ξ_range = range(-0.9999, 0.9999, length=50)
    J1_vals = dzdξ.(Ref(m), ξ_range)
    @test all(J1_vals .> 0)  # Monotonic
end

@testset "LogarithmicMap" begin
    # Valid: zmin > 0
    m = LogarithmicMap(1.0, 1000.0)

    # Boundary conditions
    @test forward(m, 1.0) ≈ -1.0 atol=1e-10
    @test forward(m, 1000.0) ≈ 1.0 atol=1e-10

    # Invalid: zmin ≤ 0
    @test_throws DomainError LogarithmicMap(0.0, 1000.0)
    @test_throws DomainError LogarithmicMap(-10.0, 1000.0)
end

@testset "profile_transform" begin
    m = TanhMap(0.0, 1000.0, 2.0)
    z_profile = [10.0, 100.0, 500.0, 1000.0]
    ξ_profile = profile_transform(m, z_profile)

    @test ξ_profile[1] < ξ_profile[2] < ξ_profile[3] < ξ_profile[4]  # Monotonic
    @test ξ_profile[end] ≈ 1.0 atol=1e-10  # Upper boundary
end

@testset "JacobianStack" begin
    m1 = TanhMap(0.0, 1000.0, 2.0)
    m2 = HyperbolicMap(-1.0, 1.0, 1.0)
    stack = JacobianStack([m1, m2], ["TanhMap", "HyperbolicMap"])

    z = 500.0
    ξ, J1, J2 = evaluate_jacobian_stack(stack, z)

    @test -1.0 < ξ < 1.0  # Within computational domain
    @test J1 > 0  # Monotonic
end
```

-----

**End of guide.** For module API, see the docstrings in `src/transforms.jl`.
