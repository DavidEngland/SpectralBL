# src/transforms.jl
"""
Coordinate Transforms for Transformed Pseudospectral Boundary Layer Analysis

This module defines 1D coordinate maps between physical heights (z) and the
bounded computational domain (ξ ∈ [-1, 1]). It exposes analytical first- and
second-order metric Jacobians necessary for constructing exact spatial
differentiation, quadrature integration, and Laplacian diffusion operators.
"""
module Transforms

using LinearAlgebra

export CoordinateMap, LinearMap, HyperbolicMap, LogarithmicMap, TanhMap, CustomMap,
       forward, inverse, dzdξ, d2zdξ2, dξdz, d2ξdz2

abstract type CoordinateMap end

# --- STAGE 0: CORE ABSTRACT METRIC INTERFACES (Analytical Fallbacks) ---
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

# --- 1. LINEAR MAP (Verification Baseline) ---
struct LinearMap <: CoordinateMap
    zmin::Float64
    zmax::Float64
end

forward(m::LinearMap, z)   = 2.0 * (z - m.zmin) / (m.zmax - m.zmin) - 1.0
inverse(m::LinearMap, ξ)   = m.zmin + (ξ + 1.0) * (m.zmax - m.zmin) / 2.0
dzdξ(m::LinearMap, ξ)      = (m.zmax - m.zmin) / 2.0
d2zdξ2(m::LinearMap, ξ)    = 0.0

# --- 2. HYPERBOLIC MAP (FIXED BOUNDARY MATCHING) ---
struct HyperbolicMap <: CoordinateMap
    zmin::Float64
    zmax::Float64
    alpha::Float64 # Stretching severity control
end

function forward(m::HyperbolicMap, z)
    L = m.zmax - m.zmin
    # Explicitly normalized so that z = zmax maps strictly to ξ = 1
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

# --- 3. LOGARITHMIC MAP ---
struct LogarithmicMap <: CoordinateMap
    zmin::Float64
    zmax::Float64
end

forward(m::LogarithmicMap, z) = 2.0 * log(z / m.zmin) / log(m.zmax / m.zmin) - 1.0
inverse(m::LogarithmicMap, ξ) = m.zmin * exp(((ξ + 1.0) / 2.0) * log(m.zmax / m.zmin))

function dzdξ(m::LogarithmicMap, ξ)
    S = log(m.zmax / m.zmin)
    return 0.5 * S * inverse(m, ξ)
end

function d2zdξ2(m::LogarithmicMap, ξ)
    S = log(m.zmax / m.zmin)
    return 0.25 * (S^2) * inverse(m, ξ)
end

# --- 4. HYPERBOLIC TANGENT MAP (Production Standard) ---
struct TanhMap <: CoordinateMap
    zmin::Float64
    zmax::Float64
    α::Float64 # Resolution compression parameters
end

function forward(m::TanhMap, z)
    L = m.zmax - m.zmin
    zc = (m.zmax + m.zmin) / 2.0
    # Guard against float drift overflow before entering inverse hyperbolic domain
    arg = clamp((z - zc) * tanh(m.α) / (L / 2.0), -0.9999999999999, 0.9999999999999)
    return atanh(arg) / m.α
end

function inverse(m::TanhMap, ξ)
    L = m.zmax - m.zmin
    zc = (m.zmax + m.zmin) / 2.0
    return zc + (L / 2.0) * tanh(m.α * ξ) / tanh(m.α)
end

function dzdξ(m::TanhMap, ξ)
    L = m.zmax - m.zmin
    return (L / 2.0) * (m.α / tanh(m.α)) * sech(m.α * ξ)^2
end

function d2zdξ2(m::TanhMap, ξ)
    L = m.zmax - m.zmin
    return -L * (m.α^2 / tanh(m.α)) * sech(m.α * ξ)^2 * tanh(m.α * ξ)
end

# --- 5. CUSTOM MAP ---
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

end # module
