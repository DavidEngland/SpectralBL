# docs/TRANSFORMS_QUICK_REFERENCE.md

# Transforms Module: Quick Reference Card

## 1. Import & Validation

```julia
using Transforms

# Every map should be validated before use
m = TanhMap(0.0, 2000.0, 2.5)
is_valid(m)  # Throws DomainError if invalid
```

## 2. Map Selection Flowchart

```
Do you need refinement near z=0?
    ├─ YES, extreme → HyperbolicMap(zmin, zmax, α)
    ├─ YES, moderate → TanhMap(zmin, zmax, α)  ← **RECOMMENDED**
    └─ NO
        └─ Is z_min > 0 and log-spacing natural?
            ├─ YES → LogarithmicMap(zmin, zmax)
            └─ NO → LinearMap(zmin, zmax)  [testing only]
```

## 3. One-Liners

|Task                    |Code                                 |
|------------------------|-------------------------------------|
|Create map              |`m = TanhMap(0.0, 1000.0, 2.5)`      |
|Transform height        |`ξ = forward(m, z)`                  |
|Recover height          |`z = inverse(m, ξ)`                  |
|Get metric J₁ = dz/dξ   |`J1 = dzdξ(m, ξ)`                    |
|Get metric J₂ = d²z/dξ² |`J2 = d2zdξ2(m, ξ)`                  |
|Inverse metric          |`dξdz(m, z)`                         |
|Transform profile vector|`ξ_vec = profile_transform(m, z_vec)`|

## 4. Parameter Defaults by Use Case

### Nocturnal Stable Boundary Layer (0–500 m)

```julia
m = TanhMap(0.0, 500.0, 3.0)  # Strong compression for inversion
```

### Convective Daytime (0–2000 m)

```julia
m = TanhMap(0.0, 2000.0, 2.0)  # Moderate compression
```

### Forest Canopy with Log Decay (10–100 m)

```julia
m = LogarithmicMap(10.0, 100.0)  # Natural log scaling
```

### Research/Testing (Any domain)

```julia
m = LinearMap(zmin, zmax)  # Trivial, no surprises
```

## 5. Integration Weights for Quadrature

```julia
# Given Gauss-Legendre nodes ξ_quad and weights w_quad on [-1, 1]:
J1 = dzdξ.(Ref(m), ξ_quad)
w_physical = w_quad .* J1

# Then ∫_zmin^zmax f(z) dz ≈ sum(w_physical .* f(z_quad))
```

## 6. Spectral Differentiation

```julia
# Physical derivative from spectral derivative
du_dξ = D * u  # D = Chebyshev differentiation matrix
du_dz = (1 ./ dzdξ.(Ref(m), ξ)) .* du_dξ
```

## 7. Laplacian Assembly

```julia
J1 = dzdξ.(Ref(m), ξ_nodes)
J2 = d2zdξ2.(Ref(m), ξ_nodes)

Lap_z = Diagonal(1 ./ (J1 .^ 2)) * D2 - Diagonal(J2 ./ (J1 .^ 3)) * D1
```

## 8. Finite-Difference Verification of Metrics

```julia
# Check J₁ against finite difference
ε = 1e-8
z_plus = inverse(m, ξ + ε)
z_minus = inverse(m, ξ - ε)
dz_fd = (z_plus - z_minus) / (2ε)
dz_analytical = dzdξ(m, ξ)
error = abs(dz_analytical - dz_fd) / abs(dz_analytical)
@assert error < 1e-6  # Should be <1e-5
```

## 9. TanhMap α Selection

|Domain  |Lower BL structure      |Recommended α|
|--------|------------------------|-------------|
|0–200 m |Strong inversion        |3.5–4.5      |
|0–500 m |Moderate gradient       |2.5–3.5      |
|0–2000 m|Weak gradient           |1.5–2.5      |
|>5000 m |Tall/weak stratification|1.0–2.0      |

**Heuristic**: If >50% of energy/variance in lowest 10% of domain, increase α by 1–2.

## 10. Common Errors & Fixes

|Error                           |Cause                   |Fix                          |
|--------------------------------|------------------------|-----------------------------|
|`DomainError: zmin ≤ 0`         |LogarithmicMap          |Use `zmin > 0`               |
|`DomainError: zmin ≥ zmax`      |Any map                 |Ensure `zmin < zmax`         |
|`DomainError: α ≤ 0`            |TanhMap/HyperbolicMap   |Use `α > 0`                  |
|`atanh domain error`            |TanhMap.forward overflow|Reduce α or increase domain  |
|Small J₁ causing division errors|Ill-conditioned map     |Reduce α; reconsider map type|

## 11. Debugging: Metric Plots

```julia
using Plots
m = TanhMap(0.0, 1000.0, 2.5)

ξ = range(-0.99, 0.99, length=200)
z = inverse.(Ref(m), ξ)
J1 = dzdξ.(Ref(m), ξ)
J2 = d2zdξ2.(Ref(m), ξ)

plot(layout=(1,2), size=(1000, 400))
plot!(z, J1; xlabel="z (m)", ylabel="dz/dξ", subplot=1)
plot!(z, J2; xlabel="z (m)", ylabel="d²z/dξ²", subplot=2)
```

## 12. Jacobian Stack (Advanced)

```julia
# Multi-layer mapping: z_phys → ξ_comp
m1 = TanhMap(0.0, 2000.0, 2.5)           # Stability compress
m2 = HyperbolicMap(-1.0, 1.0, 1.5)       # Chebyshev refine
stack = JacobianStack([m1, m2])

z = 150.0
ξ, J1, J2 = evaluate_jacobian_stack(stack, z)
```

-----

## Function Signatures

```julia
# Core operations
forward(m::CoordinateMap, z::Float64) -> Float64
inverse(m::CoordinateMap, ξ::Float64) -> Float64
dzdξ(m::CoordinateMap, ξ::Float64) -> Float64
d2zdξ2(m::CoordinateMap, ξ::Float64) -> Float64

# Inverse metrics (chain rule)
dξdz(m::CoordinateMap, z::Float64) -> Float64
d2ξdz2(m::CoordinateMap, z::Float64) -> Float64

# Utilities
is_valid(m::CoordinateMap) -> Bool
profile_transform(m::CoordinateMap, z_vec::Vector) -> Vector
evaluate_jacobian_stack(stack::JacobianStack, z::Float64) -> Tuple(ξ, J1, J2)
```

-----

**Last Updated**: June 2026 | **Version**: 2.0 (Fully Validated)
