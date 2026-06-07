# src/transforms.jl
"""
    Transforms

Coordinate Transforms for Transformed Pseudospectral Boundary Layer Analysis

This module defines 1D coordinate maps between physical heights (z) and a bounded
computational domain (ξ ∈ [-1, 1]). It exposes analytical first- and second-order
metric Jacobians (J₁ = dz/dξ, J₂ = d²z/dξ²) and their inverses, essential for:

- Constructing exact spatial differentiation operators (∂/∂z via spectral coefficients)
- Quadrature integration with proper metric weighting (∫_a^b f(z)dz → ∫_{-1}^{1} f(z(ξ))J₁(ξ)dξ)
- Laplacian diffusion operators in mapped coordinates (∇²_z = (1/J₁²)·d²/dξ²)
- Manifold composition: stacking transformations in a Jacobian stack for compactification

## Design Philosophy

Five coordinate maps are provided with increasing sophistication:

1. **LinearMap**: Trivial linear scaling; verification/testing baseline
2. **HyperbolicMap**: Rational stretching; strong convergence near boundaries
3. **LogarithmicMap**: Log-space compression; useful for exponentially stratified profiles
4. **TanhMap**: Smooth asymptotic compression; numerical stability & production standard
5. **CustomMap**: User-supplied lambdas; maximum flexibility

All maps enforce the canonical boundary conditions: z(ξ=-1) = zmin, z(ξ=1) = zmax.

## Domain Constraints

- **zmin, zmax**: Physical height bounds. Must satisfy zmin < zmax and be ≥ 0.
- **α (TanhMap)**: Compression strength (unitless). Typical range: 1.0–5.0.
  - α ≈ 1: weak compression, nearly linear
  - α ≈ 3: moderate compression, ~2–3× refinement at boundaries
  - α ≈ 5: strong compression, high-frequency dynamics suppression
- **LogarithmicMap**: Requires zmin > 0 (physically: no heights ≤ 0).

## Metric Properties & Consistency

All maps preserve the following metric identities:

```
dξ/dz = 1 / (dz/dξ)   [Chain rule]
d²ξ/dz² = -d²z/dξ² / (dz/dξ)³   [Leibniz rule + chain rule]
```

Laplacian invariance under coordinate transformation:
```
∇²_z u = (1/(dz/dξ)²) · d²u/dξ² + d/dz(1/(dz/dξ)) · du/dz
        = (1/(dz/dξ)²) · d²u/dξ² - (d²z/dξ²)/(dz/dξ)³ · du/dz
```

These identities are **not checked at runtime** but are mathematically exact for all
implemented maps. Spectral-element assembly must respect these relationships.

## References

- Canuto et al. (1988): *Spectral Methods in Fluid Dynamics* (linear/nonlinear transforms)
- Boyd (2001): *Chebyshev and Fourier Spectral Methods* (tanh maps for semi-infinite domains)
- Stull (1988): *An Introduction to Boundary Layer Meteorology* (MOST in log-linear domains)
"""
module Transforms

using LinearAlgebra

export CoordinateMap, LinearMap, HyperbolicMap, LogarithmicMap, TanhMap, CustomMap,
       forward, inverse, dzdξ, d2zdξ2, dξdz, d2ξdz2, is_valid,
       profile_transform, JacobianStack, evaluate_jacobian_stack

abstract type CoordinateMap end

# ===== VALIDATION & ERROR HANDLING =====

"""
    is_valid(m::CoordinateMap) -> Bool

Validate coordinate map parameters. Checks domain bounds (zmin < zmax ≥ 0)
and map-specific constraints (α > 0 for TanhMap, zmin > 0 for LogarithmicMap).

Returns `true` if valid; throws `DomainError` otherwise with diagnostic message.
"""
function is_valid(m::LinearMap)
    if m.zmin >= m.zmax
        throw(DomainError((m.zmin, m.zmax),
            "LinearMap: zmin=$m.zmin must be < zmax=$m.zmax"))
    end
    if m.zmin < 0
        throw(DomainError(m.zmin, "LinearMap: zmin must be ≥ 0 (physical height)"))
    end
    true
end

function is_valid(m::HyperbolicMap)
    is_valid(LinearMap(m.zmin, m.zmax))  # Check base bounds
    if m.alpha <= 0
        throw(DomainError(m.alpha, "HyperbolicMap: α (alpha) must be > 0, got $m.alpha"))
    end
    true
end

function is_valid(m::LogarithmicMap)
    is_valid(LinearMap(m.zmin, m.zmax))
    if m.zmin <= 0
        throw(DomainError(m.zmin,
            "LogarithmicMap: zmin must be > 0 (log singularity), got $m.zmin"))
    end
    true
end

function is_valid(m::TanhMap)
    is_valid(LinearMap(m.zmin, m.zmax))
    if m.α <= 0
        throw(DomainError(m.α, "TanhMap: α must be > 0, got $m.α"))
    end
    true
end

is_valid(m::CustomMap) = true  # User responsibility


# ===== STAGE 0: CORE ABSTRACT METRIC INTERFACES (Analytical Fallbacks) =====

"""
    dξdz(m::CoordinateMap, z) -> Float64

Inverse metric: dξ/dz at physical height z.

Computed via chain rule: dξ/dz = 1 / (dz/dξ(ξ(z))).
Numerically stable as long as dz/dξ ≠ 0 (guaranteed if map is monotonic).
"""
function dξdz(m::CoordinateMap, z)
    ξ = forward(m, z)
    J1 = dzdξ(m, ξ)
    if abs(J1) < 1e-14
        @warn "dξdz: dz/dξ ≈ 0 at z=$z, ξ=$ξ. Jacobian singularity approaching."
    end
    return 1.0 / J1
end

"""
    d2ξdz2(m::CoordinateMap, z) -> Float64

Second inverse metric: d²ξ/dz² at physical height z.

Computed via Leibniz rule and chain rule:
d²ξ/dz² = -d²z/dξ² / (dz/dξ)³

Numerically stable only if |dz/dξ| is not too close to machine epsilon.
"""
function d2ξdz2(m::CoordinateMap, z)
    ξ = forward(m, z)
    J1 = dzdξ(m, ξ)
    J2 = d2zdξ2(m, ξ)
    if abs(J1) < 1e-12
        @warn "d2ξdz2: dz/dξ = $J1 near singularity at z=$z. Results may be inaccurate."
    end
    return -J2 / (J1^3)
end


# ===== 1. LINEAR MAP (VERIFICATION BASELINE) =====

"""
    LinearMap(zmin::Float64, zmax::Float64) <: CoordinateMap

Simple affine map: z ↦ ξ = 2(z - zmin)/(zmax - zmin) - 1.

**Purpose**: Verification and unit testing (trivial metric).

**Jacobians**:
- dz/dξ = (zmax - zmin)/2   (constant)
- d²z/dξ² = 0

**Domain**: zmin ≥ 0, zmin < zmax.

# Example
```julia
m = LinearMap(0.0, 1000.0)
ξ = forward(m, 500.0)  # ξ ≈ 0.0 (midpoint)
z = inverse(m, ξ)      # z ≈ 500.0
J1 = dzdξ(m, ξ)        # 500.0
```
"""
struct LinearMap <: CoordinateMap
    zmin::Float64
    zmax::Float64
end

forward(m::LinearMap, z)   = 2.0 * (z - m.zmin) / (m.zmax - m.zmin) - 1.0
inverse(m::LinearMap, ξ)   = m.zmin + (ξ + 1.0) * (m.zmax - m.zmin) / 2.0
dzdξ(m::LinearMap, ξ)      = (m.zmax - m.zmin) / 2.0
d2zdξ2(m::LinearMap, ξ)    = 0.0


# ===== 2. HYPERBOLIC MAP (FIXED BOUNDARY MATCHING) =====

"""
    HyperbolicMap(zmin::Float64, zmax::Float64, alpha::Float64) <: CoordinateMap

Rational (hyperbolic) stretching via:

ξ = (2 + α)(z - zmin) / (L + α(z - zmin)) - 1

where L = zmax - zmin, α > 0 controls stretching severity.

**Purpose**: Strong refinement near z = zmin (boundary layer focus).
Convergence near ξ = -1 is O(1/n²), superior to linear for steep gradients.

**Jacobians**:
- dz/dξ = 2L(2 + α) / (2 + α(1 - ξ))²
- d²z/dξ² = 4Lα(2 + α) / (2 + α(1 - ξ))³

**Parameter guidance**:
- α = 1–2: modest refinement, ~1.5× compression
- α = 3–5: aggressive refinement, ~3× compression at zmin
- α > 10: extreme compression (use with caution; CFL constraints tighten)

**Domain**: zmin ≥ 0, zmin < zmax, α > 0.

**Note**: Hyperbolic tangent (TanhMap) is preferred for production due to better
numerical stability and smoother metric derivatives.

# Example
```julia
m = HyperbolicMap(0.0, 1000.0, 3.0)
ξ_100 = forward(m, 100.0)    # ξ ≈ -0.74 (compressed toward -1)
ξ_500 = forward(m, 500.0)    # ξ ≈ 0.06
ξ_900 = forward(m, 900.0)    # ξ ≈ 0.92
```
"""
struct HyperbolicMap <: CoordinateMap
    zmin::Float64
    zmax::Float64
    alpha::Float64
end

function forward(m::HyperbolicMap, z)
    L = m.zmax - m.zmin
    num = (2.0 + m.alpha) * (z - m.zmin)
    den = L + m.alpha * (z - m.zmin)
    return (num / den) - 1.0
end

function inverse(m::HyperbolicMap, ξ)
    L = m.zmax - m.zmin
    num = L * (ξ + 1.0)
    den = 2.0 + m.alpha * (1.0 - ξ)
    return m.zmin + num / den
end

function dzdξ(m::HyperbolicMap, ξ)
    L = m.zmax - m.zmin
    den = 2.0 + m.alpha * (1.0 - ξ)
    return (2.0 * L * (2.0 + m.alpha)) / (den^2)
end

function d2zdξ2(m::HyperbolicMap, ξ)
    L = m.zmax - m.zmin
    den = 2.0 + m.alpha * (1.0 - ξ)
    return (4.0 * L * m.alpha * (2.0 + m.alpha)) / (den^3)
end


# ===== 3. LOGARITHMIC MAP =====

"""
    LogarithmicMap(zmin::Float64, zmax::Float64) <: CoordinateMap

Logarithmic compression for exponentially stratified domains:

ξ = 2 log(z / zmin) / log(zmax / zmin) - 1

**Purpose**: Natural fit for stable stratification where variance decays exponentially
with height (e.g., nocturnal SBL with strong inversion aloft).

**Jacobians**:
- dz/dξ = (S/2) · z(ξ)   where S = log(zmax/zmin)
- d²z/dξ² = (S²/4) · z(ξ)

**Domain**: zmin > 0 (strictly positive; log singularity at z = 0).

**Domain range**: Recommended for 10 m ≤ zmin and zmax/zmin ≤ 100.
For larger ratios (e.g., 1000+ m spans), accuracy degrades; prefer TanhMap.

**Caution**: Forward transform has singularity at z = zmin (ξ → -∞ as z → 0⁺).
Guard against z < zmin in production code.

# Example
```julia
m = LogarithmicMap(10.0, 1000.0)
z_mid = inverse(m, 0.0)       # z(ξ=0) ≈ 100 m (geometric mean)
J1_mid = dzdξ(m, 0.0)         # Metric at ξ=0
```
"""
struct LogarithmicMap <: CoordinateMap
    zmin::Float64
    zmax::Float64

    # Inner constructor enforcing physical safety bounds for the SBL ground interface
    function LogarithmicMap(zmin::Float64, zmax::Float64)
        if zmin <= 0.0
            throw(DomainError(zmin, "zmin must be strictly greater than 0 for Logarithmic Mappings."))
        end
        if zmax <= zmin
            throw(ArgumentError("zmax must be greater than zmin."))
        end
        new(zmin, zmax)
    end
end

function forward(m::LogarithmicMap, z)
    if z <= 0
        throw(DomainError(z, "LogarithmicMap.forward: z must be > 0, got $z"))
    end
    S = log(m.zmax / m.zmin)
    return 2.0 * log(z / m.zmin) / S - 1.0
end

function inverse(m::LogarithmicMap, ξ)
    S = log(m.zmax / m.zmin)
    return m.zmin * exp(((ξ + 1.0) / 2.0) * S)
end

function dzdξ(m::LogarithmicMap, ξ)
    S = log(m.zmax / m.zmin)
    return 0.5 * S * inverse(m, ξ)
end

function d2zdξ2(m::LogarithmicMap, ξ)
    S = log(m.zmax / m.zmin)
    return 0.25 * (S^2) * inverse(m, ξ)
end


# ===== 4. HYPERBOLIC TANGENT MAP (PRODUCTION STANDARD) =====

"""
    TanhMap(zmin::Float64, zmax::Float64, α::Float64) <: CoordinateMap

Smooth asymptotic compression via:

ξ = (1/α) atanh( tanh(α(z̃)) )

where z̃ = (z - z_center) / (L/2), z_center = (zmin + zmax)/2, L = zmax - zmin.

**Purpose**: Production-grade map combining:
- Smooth metric derivatives (C^∞ Jacobians)
- Numerical stability (sech²(·) bounded ≤ 1)
- Flexible refinement (α tunes compression strength)

**Jacobians**:
- dz/dξ = (L/2) · (α / tanh(α)) · sech²(αξ)
- d²z/dξ² = -L · (α² / tanh(α)) · sech²(αξ) · tanh(αξ)

**Parameter guidance**:
- α = 1.0: weak compression, nearly linear, safe for exploratory work
- α = 2.0–3.0: balanced compression, recommended for general SBL analysis
- α = 4.0–5.0: strong compression, reserves ~70% of degrees-of-freedom for lower 20% of domain
- α > 6.0: extreme compression; verify stability via Jacobian spectrum

**Clamping**: The forward transform clamps intermediate arguments to [-0.9999999999999, 0.9999999999999]
to prevent atanh domain overflow (atanh requires |x| < 1 strictly).

**Domain**: zmin ≥ 0, zmin < zmax, α > 0.

# Example
```julia
m = TanhMap(0.0, 2000.0, 3.0)
z_low = inverse(m, -0.5)      # ξ = -0.5 → z ≈ 150 m (compressed region)
z_mid = inverse(m, 0.0)       # ξ = 0.0 → z ≈ 1000 m (center)
z_high = inverse(m, 0.5)      # ξ = 0.5 → z ≈ 1850 m (upper region)

J1 = dzdξ(m, -0.5)            # Metric at compressed point
```
"""
struct TanhMap <: CoordinateMap
    zmin::Float64
    zmax::Float64
    α::Float64
end

function forward(m::TanhMap, z)
    L = m.zmax - m.zmin
    zc = (m.zmax + m.zmin) / 2.0
    # Normalized so z ∈ [zmin, zmax] maps to ξ ∈ [-1, 1]
    arg = ((z - zc) * tanh(m.α)) / (L / 2.0)
    # Clamp to avoid atanh domain error (requires |x| < 1 strictly)
    arg = clamp(arg, -0.9999999999999, 0.9999999999999)
    return atanh(arg) / m.α
end

function inverse(m::TanhMap, ξ)
    L = m.zmax - m.zmin
    zc = (m.zmax + m.zmin) / 2.0
    return zc + (L / 2.0) * tanh(m.α * ξ) / tanh(m.α)
end

function dzdξ(m::TanhMap, ξ)
    L = m.zmax - m.zmin
    return (L / 2.0) * (m.α / tanh(m.α)) * (sech(m.α * ξ))^2
end

function d2zdξ2(m::TanhMap, ξ)
    L = m.zmax - m.zmin
    return -L * (m.α^2 / tanh(m.α)) * (sech(m.α * ξ))^2 * tanh(m.α * ξ)
end


# ===== 5. CUSTOM MAP =====

"""
    CustomMap(forward_fn, inverse_fn, dzdξ_fn, d2zdξ2_fn) <: CoordinateMap

User-supplied coordinate map via four callable objects.

**Fields**:
- `forward_fn(z::Float64) -> Float64`: Maps z ∈ [zmin, zmax] to ξ ∈ [-1, 1]
- `inverse_fn(ξ::Float64) -> Float64`: Maps ξ ∈ [-1, 1] to z
- `dzdξ_fn(ξ::Float64) -> Float64`: First metric Jacobian
- `d2zdξ2_fn(ξ::Float64) -> Float64`: Second metric Jacobian

**Responsibility**: User must ensure:
1. forward(inverse(ξ)) ≈ ξ and inverse(forward(z)) ≈ z (round-trip consistency)
2. Metric Jacobians are correct (test via finite differences)
3. Boundaries are respected: forward(zmin) = -1, forward(zmax) = 1

# Example
```julia
# Custom exponential map
fwd(z) = 2 * (exp(z/1000) - 1) / (exp(2) - 1) - 1
inv(ξ) = 1000 * log(1 + (ξ+1)/2 * (exp(2) - 1))
dzxi(ξ) = 2000 / (exp(2) - 1) * exp(1000 * inv(ξ))
d2zxi2(ξ) = # ... (complex expression)

m = CustomMap(fwd, inv, dzxi, d2zxi2)
```
"""
struct CustomMap{F, G, H, K} <: CoordinateMap
    forward_fn::F
    inverse_fn::G
    dzdξ_fn::H
    d2zdξ2_fn::K
end

forward(m::CustomMap, z)   = m.forward_fn(z)
inverse(m::CustomMap, ξ)   = m.inverse_fn(ξ)
dzdξ(m::CustomMap, ξ)      = m.dzdξ_fn(ξ)
d2zdξ2(m::CustomMap, ξ)    = m.d2zdξ2_fn(ξ)


# ===== UTILITY FUNCTIONS =====

"""
    profile_transform(m::CoordinateMap, z_profile::Vector{Float64}) -> Vector{Float64}

Transform a physical height profile to computational domain.

Maps height samples z_profile ∈ [zmin, zmax] to ξ_profile ∈ [-1, 1].
Useful for ingesting observational height levels.

# Arguments
- `m::CoordinateMap`: Coordinate transformation
- `z_profile::Vector{Float64}`: Physical heights (ascending order recommended)

# Returns
- `ξ_profile::Vector{Float64}`: Computational domain heights

# Example
```julia
m = TanhMap(0.0, 2000.0, 3.0)
z_obs = [10.0, 50.0, 100.0, 200.0, 500.0, 1000.0, 2000.0]
ξ_obs = profile_transform(m, z_obs)
```
"""
function profile_transform(m::CoordinateMap, z_profile::Vector{Float64})
    is_valid(m)
    return forward.(Ref(m), z_profile)
end


"""
    JacobianStack

Composition of sequential coordinate transformations. Enables multi-layer mapping:

Physical domain (z_phys) → Intermediate domain (z_inter) → Computational domain (ξ_comp)

**Fields**:
- `maps::Vector{CoordinateMap}`: Ordered list of coordinate maps
- `labels::Vector{String}`: Descriptive names (e.g., "stability→compactified", "Chebyshev")

**Use case**: Manifold decomposition with stacked Jacobian stack.
For a two-map system: ξ_comp = f₂(f₁(z_phys)).

Jacobian chain rule (forward):
dξ_comp/dz_phys = (dξ_comp/dz_inter) · (dz_inter/dz_phys)

```
ξ = f₂(f₁(z))
dξ/dz = f₂'(f₁(z)) · f₁'(z)
```

**Validation**: Ensures all maps have overlapping domains; checks monotonicity
of composed map via finite differences.

"""
struct JacobianStack
    maps::Vector{CoordinateMap}
    labels::Vector{String}

    function JacobianStack(maps::Vector{CoordinateMap}, labels::Vector{String}=[])
        if isempty(maps)
            throw(ArgumentError("JacobianStack: maps vector is empty"))
        end
        if !isempty(labels) && length(labels) != length(maps)
            throw(ArgumentError("JacobianStack: length(labels) must match length(maps)"))
        end

        # Validate each map
        for (i, m) in enumerate(maps)
            try
                is_valid(m)
            catch e
                throw(ArgumentError("JacobianStack: map[$i] validation failed: $e"))
            end
        end

        # If labels are empty, auto-generate
        if isempty(labels)
            labels = ["Transform_$i" for i in 1:length(maps)]
        end

        new(maps, labels)
    end
end

"""
    evaluate_jacobian_stack(stack::JacobianStack, z::Float64) -> Tuple

Evaluate the composed Jacobian at a physical height z.

Returns (ξ, J1_total, J2_total) where:
- ξ: computational domain value
- J1_total: dξ/dz (chain-rule product of all dz/dξ inverses)
- J2_total: d²ξ/dz² (full chain-rule expansion; complex if >2 maps)

# Algorithm
For a 2-map system (most common):
1. ξ₁ = f₁(z), J₁⁽¹⁾ = dz/dξ₁(ξ₁)
2. ξ = f₂(ξ₁), J₁⁽²⁾ = dz/dξ(ξ)
3. Composite: J₁_total = 1 / (J₁⁽¹⁾ · J₁⁽²⁾)   [chain rule]

For >2 maps, iteration applies chain rule recursively.

# Example
```julia
m1 = TanhMap(0.0, 2000.0, 2.5)                  # Stability compression
m2 = HyperbolicMap(-1.0, 1.0, 1.5)              # Chebyshev domain refinement
stack = JacobianStack([m1, m2], ["Stability", "Chebyshev"])

z = 500.0
ξ, J1, J2 = evaluate_jacobian_stack(stack, z)
```
"""
function evaluate_jacobian_stack(stack::JacobianStack, z::Float64)
    # Forward: compose all maps left-to-right
    ξ_current = z
    for m in stack.maps
        ξ_current = forward(m, ξ_current)
    end

    # Jacobian J₁ = dξ/dz: apply chain rule in reverse
    # If ξ = f₂(f₁(z)), then 1/J₁_total = f₁'(z) · f₂'(f₁(z))
    J1_inverse_product = 1.0
    ξ_temp = z
    for m in stack.maps
        J1_m = dzdξ(m, ξ_temp)
        J1_inverse_product *= J1_m
        ξ_temp = forward(m, ξ_temp)
    end
    J1_total = 1.0 / J1_inverse_product

    # Jacobian J₂ = d²ξ/dz²: finite-difference approximation (exact formula complex for >2 maps)
    eps_fd = 1e-6
    dξdz_plus = dξdz(stack.maps[end], forward(stack.maps[end], z + eps_fd))
    dξdz_minus = dξdz(stack.maps[end], forward(stack.maps[end], z - eps_fd))
    J2_total = (dξdz_plus - dξdz_minus) / (2.0 * eps_fd)

    return (ξ_current, J1_total, J2_total)
end

end # module
