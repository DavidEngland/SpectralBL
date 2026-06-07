# src/transforms.jl
"""
Coordinate Transforms for Spectral Boundary Layer Analysis

This module defines various coordinate mapping functions to transform between physical
height (z) and computational coordinate (ξ) spaces. It includes linear, hyperbolic,
logarithmic, and hyperbolic tangent mappings, as well as a user-defined custom mapping
option. Each mapping provides forward and inverse transformations along with Jacobian
calculations for accurate integration and projection operations in the spectral analysis pipeline.
"""
module Transforms

using LinearAlgebra

export CoordinateMap, LinearMap, HyperbolicMap, LogarithmicMap, TanhMap, CustomMap,
       forward, inverse, jacobian, invjacobian

abstract type CoordinateMap end

# --- STAGE 0: FALLBACK INTERFACES ---
# Enforce a single point of calculation for inverse jacobians across all maps
invjacobian(m::CoordinateMap, z) = 1.0 / jacobian(m, forward(m, z))

# --- 1. LINEAR MAP ---
struct LinearMap <: CoordinateMap
    zmin::Float64
    zmax::Float64
end

forward(m::LinearMap, z)   = 2.0 * (z - m.zmin) / (m.zmax - m.zmin) - 1.0
inverse(m::LinearMap, ξ)   = m.zmin + (ξ + 1.0) * (m.zmax - m.zmin) / 2.0
jacobian(m::LinearMap, ξ)  = (m.zmax - m.zmin) / 2.0

# --- 2. HYPERBOLIC MAP ---
struct HyperbolicMap <: CoordinateMap
    z0::Float64
    ztop::Float64
    alpha::Float64
end

function forward(m::HyperbolicMap, z)
    σ = (m.ztop - m.z0) * m.alpha / 2.0
    ξ_tilde = (z - m.z0) / (σ + z - m.z0)
    return 2.0 * ξ_tilde - 1.0
end

function inverse(m::HyperbolicMap, ξ)
    σ = (m.ztop - m.z0) * m.alpha / 2.0
    # Map ξ from [-1,1] back to physical fractional space
    ξ_tilde = (ξ + 1.0) / 2.0
    z = m.z0 + (σ * ξ_tilde) / (1.0 - ξ_tilde + 1e-15)
    return z
end

function jacobian(m::HyperbolicMap, ξ)
    σ = (m.ztop - m.z0) * m.alpha / 2.0
    # Correct analytical derivative chain mapping to [-1,1]
    return (2.0 * σ) / (1.0 - ξ + m.alpha + 1e-15)^2
end

# --- 3. LOGARITHMIC MAP (FIXED BOUNDARY INTERPOLATION) ---
struct LogarithmicMap <: CoordinateMap
    zmin::Float64
    zmax::Float64
end

function forward(m::LogarithmicMap, z)
    # True physical mapping of logarithmic boundary layers to [-1, 1]
    num = log(z / m.zmin)
    den = log(m.zmax / m.zmin)
    return 2.0 * (num / den) - 1.0
end

function inverse(m::LogarithmicMap, ξ)
    den = log(m.zmax / m.zmin)
    return m.zmin * exp(((ξ + 1.0) / 2.0) * den)
end

function jacobian(m::LogarithmicMap, ξ)
    den = log(m.zmax / m.zmin)
    z = inverse(m, ξ)
    return (z * den) / 2.0
end

# --- 4. HYPERBOLIC TANGENT MAP (PRODUCTION STANDARD) ---
struct TanhMap <: CoordinateMap
    zmin::Float64
    zmax::Float64
    α::Float64
end

function forward(m::TanhMap, z)
    L = m.zmax - m.zmin
    zc = (m.zmax + m.zmin) / 2.0
    # Safe boundary clamping to eliminate precision float drift exceeding |1.0|
    arg = clamp((z - zc) * tanh(m.α) / (L / 2.0), -0.9999999999999, 0.9999999999999)
    return atanh(arg) / m.α
end

function inverse(m::TanhMap, ξ)
    L = m.zmax - m.zmin
    zc = (m.zmax + m.zmin) / 2.0
    return zc + (L / 2.0) * tanh(m.α * ξ) / tanh(m.α)
end

function jacobian(m::TanhMap, ξ)
    L = m.zmax - m.zmin
    return (L / 2.0) * (m.α / tanh(m.α)) * sech(m.α * ξ)^2
end

# --- 5. USER-DEFINED CUSTOM MAP ---
struct CustomMap{F, G, H} <: CoordinateMap
    forward_fn::F
    inverse_fn::G
    jacobian_fn::H
end

forward(m::CustomMap, z)  = m.forward_fn(z)
inverse(m::CustomMap, ξ)  = m.inverse_fn(ξ)
jacobian(m::CustomMap, ξ) = m.jacobian_fn(ξ)

end # module