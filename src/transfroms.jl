# src/transforms.jl
"""
    Transforms

Coordinate Transforms for Transformed Pseudospectral Boundary Layer Analysis

This module defines 1D coordinate maps between physical heights (z) and a bounded
computational domain (ξ ∈ [-1, 1]). It exposes analytical first- and second-order
metric Jacobians (J₁ = dz/dξ, J₂ = d²z/dξ²) and their inverses.
"""
module Transforms

using LinearAlgebra

export CoordinateMap, LinearMap, HyperbolicMap, LogarithmicMap, TanhMap, CfdWallMap, CustomMap,
       forward, inverse, dzdξ, d2zdξ2, dξdz, d2ξdz2, is_valid,
       profile_transform, JacobianStack, evaluate_jacobian_stack

abstract type CoordinateMap end

# ===== VALIDATION & ERROR HANDLING =====

function is_valid(m::LinearMap)
    if m.zmin >= m.zmax
        throw(DomainError((m.zmin, m.zmax), "LinearMap: zmin=$m.zmin must be < zmax=$m.zmax"))
    end
    if m.zmin < 0
        throw(DomainError(m.zmin, "LinearMap: zmin must be ≥ 0 (physical height)"))
    end
    true
end

function is_valid(m::HyperbolicMap)
    is_valid(LinearMap(m.zmin, m.zmax))
    if m.alpha <= 0
        throw(DomainError(m.alpha, "HyperbolicMap: α (alpha) must be > 0, got $m.alpha"))
    end
    true
end

function is_valid(m::LogarithmicMap)
    is_valid(LinearMap(m.zmin, m.zmax))
    if m.zmin <= 0
        throw(DomainError(m.zmin, "LogarithmicMap: zmin must be > 0 (log singularity), got $m.zmin"))
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

function is_valid(m::CfdWallMap)
    if m.ztop <= 0
        throw(DomainError(m.ztop, "CfdWallMap: ztop must be > 0"))
    end
    if m.δ <= 0
        throw(DomainError(m.δ, "CfdWallMap: δ packing parameter must be > 0"))
    end
    true
end

is_valid(m::CustomMap) = true


# ===== STAGE 0: CORE ABSTRACT METRIC INTERFACES =====

function dξdz(m::CoordinateMap, z)
    ξ = forward(m, z)
    return 1.0 / dzdξ(m, ξ)
end

function d2ξdz2(m::CoordinateMap, z)
    ξ = forward(m, z)
    J1 = dzdξ(m, ξ)
    J2 = d2zdξ2(m, ξ)
    return -J2 / (J1^3)
end


# ===== 1. LINEAR MAP =====

struct LinearMap <: CoordinateMap
    zmin::Float64
    zmax::Float64
end

forward(m::LinearMap, z)   = 2.0 * (z - m.zmin) / (m.zmax - m.zmin) - 1.0
inverse(m::LinearMap, ξ)   = m.zmin + (ξ + 1.0) * (m.zmax - m.zmin) / 2.0
dzdξ(m::LinearMap, ξ)      = (m.zmax - m.zmin) / 2.0
d2zdξ2(m::LinearMap, ξ)    = 0.0


# ===== 2. HYPERBOLIC MAP =====

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

struct LogarithmicMap <: CoordinateMap
    zmin::Float64
    zmax::Float64

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


# ===== 4. HYPERBOLIC TANGENT MAP =====

struct TanhMap <: CoordinateMap
    zmin::Float64
    zmax::Float64
    α::Float64
end

function forward(m::TanhMap, z)
    L = m.zmax - m.zmin
    zc = (m.zmax + m.zmin) / 2.0
    arg = ((z - zc) * tanh(m.α)) / (L / 2.0)
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


# ===== 5. CFD WALL BOUNDARY LAYER MAP (FIXED) =====

struct CfdWallMap <: CoordinateMap
    ztop::Float64
    δ::Float64
end

forward(m::CfdWallMap, z)  = 1.0 - (1.0 / m.δ) * atanh(clamp((1.0 - z / m.ztop) * tanh(m.δ), -0.9999999999999, 0.9999999999999))
inverse(m::CfdWallMap, ξ)  = m.ztop * (1.0 - tanh(m.δ * (1.0 - ξ)) / tanh(m.δ))
dzdξ(m::CfdWallMap, ξ)     = m.ztop * (m.δ / tanh(m.δ)) * sech(m.δ * (1.0 - ξ))^2
d2zdξ2(m::CfdWallMap, ξ)    = -2.0 * m.ztop * (m.δ^2 / tanh(m.δ)) * sech(m.δ * (1.0 - ξ))^2 * tanh(m.δ * (1.0 - ξ))


# ===== 6. CUSTOM MAP =====

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


# ===== UTILITY FUNCTIONS & STACK COMPOSITION =====

function profile_transform(m::CoordinateMap, z_profile::Vector{Float64})
    is_valid(m)
    return forward.(Ref(m), z_profile)
end

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
        for (i, m) in enumerate(maps)
            is_valid(m)
        end
        if isempty(labels)
            labels = ["Transform_$i" for i in 1:length(maps)]
        end
        new(maps, labels)
    end
end

"""
    evaluate_jacobian_stack(stack::JacobianStack, z::Float64) -> Tuple{Float64, Float64, Float64}

Evaluates a multi-layer stacked coordinate manifold at physical location z.
Returns (ξ_final, dξ/dz, d²ξ/dz²) calculated entirely via exact analytical chain rules,
preserving scale invariance without finite-difference precision loss.
"""
function evaluate_jacobian_stack(stack::JacobianStack, z::Float64)
    # Track the intermediate coordinates x_i and computational positions ξ_i
    coords = Vector{Float64}(undef, length(stack.maps) + 1)
    coords[1] = z
    for i in 1:length(stack.maps)
        coords[i+1] = forward(stack.maps[i], coords[i])
    end
    ξ_final = coords[end]

    # Compute dξ/dz via sequential chain rule multiplication
    # dξ/dz = ∏ (dξ_i / dx_i) where dξ_i/dx_i = 1 / dzdξ_i(ξ_i)
    dξdx = Vector{Float64}(undef, length(stack.maps))
    for i in 1:length(stack.maps)
        dξdx[i] = 1.0 / dzdξ(stack.maps[i], coords[i+1])
    end
    dξdz_total = prod(dξdx)

    # Compute d²ξ/dz² using the exact analytical derivative chain rule expansion:
    # d²ξ/dz² = ∑_i [ (d²ξ_i / dx_i²) * (dx_i / dz)² * ∏_{j>i} (dξ_j / dx_j) ]
    d2ξdz2_total = 0.0
    dxdz_current = 1.0 # dx_1/dz is trivially 1.0

    for i in 1:length(stack.maps)
        m = stack.maps[i]
        ξ_i = coords[i+1]

        # Pull core component metrics
        J1 = dzdξ(m, ξ_i)
        J2 = d2zdξ2(m, ξ_i)

        # Calculate localized spatial mapping curvature
        d2ξdx2 = -J2 / (J1^3)

        # Product of all downstream forward derivatives: ∏_{j>i} (dξ_j / dx_j)
        downstream_prod = 1.0
        for j in (i+1):length(stack.maps)
            downstream_prod *= dξdx[j]
        end

        d2ξdz2_total += d2ξdx2 * (dxdz_current^2) * downstream_prod

        # Update current step derivative mapping for next nested level
        dxdz_current *= dξdx[i]
    end

    return (ξ_final, dξdz_total, d2ξdz2_total)
end

end # module