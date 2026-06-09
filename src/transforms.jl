# src/transforms.jl
module Transforms

using LinearAlgebra

export CoordinateMap, LinearMap, HyperbolicMap, LogarithmicMap, TanhMap, CfdWallMap, CustomMap,
       forward, inverse, dzdξ, d2zdξ2, dξdz, d2ξdz2, is_valid,
    profile_transform, JacobianStack, evaluate_jacobian_stack,
    map_name, map_parameter_pairs, describe_map

abstract type CoordinateMap end

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

    function LinearMap(zmin::Float64, zmax::Float64)
        zmin >= zmax && throw(DomainError((zmin, zmax), "LinearMap: zmin must be < zmax"))
        zmin < 0.0 && throw(DomainError(zmin, "LinearMap: zmin must be ≥ 0"))
        new(zmin, zmax)
    end
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

    function HyperbolicMap(zmin::Float64, zmax::Float64, alpha::Float64)
        zmin >= zmax && throw(DomainError((zmin, zmax), "HyperbolicMap: zmin must be < zmax"))
        zmin < 0.0 && throw(DomainError(zmin, "HyperbolicMap: zmin must be ≥ 0"))
        alpha <= 0.0 && throw(DomainError(alpha, "HyperbolicMap: α must be > 0"))
        new(zmin, zmax, alpha)
    end
end

function forward(m::HyperbolicMap, z)
    L = m.zmax - m.zmin
    num = 2.0 * (2.0 + m.alpha) * (z - m.zmin)
    den = 2.0 * L + m.alpha * (z - m.zmin)
    return (num / den) - 1.0
end

function inverse(m::HyperbolicMap, ξ)
    L = m.zmax - m.zmin
    num = 2.0 * L * (ξ + 1.0)
    den = 2.0 * (2.0 + m.alpha) - m.alpha * (ξ + 1.0)
    return m.zmin + num / den
end

function dzdξ(m::HyperbolicMap, ξ)
    L = m.zmax - m.zmin
    term = 2.0 * (2.0 + m.alpha) - m.alpha * (ξ + 1.0)
    return (4.0 * L * (2.0 + m.alpha)) / (term^2)
end

function d2zdξ2(m::HyperbolicMap, ξ)
    L = m.zmax - m.zmin
    term = 2.0 * (2.0 + m.alpha) - m.alpha * (ξ + 1.0)
    return (8.0 * L * m.alpha * (2.0 + m.alpha)) / (term^3)
end

# ===== 3. LOGARITHMIC MAP =====
struct LogarithmicMap <: CoordinateMap
    zmin::Float64
    zmax::Float64
    function LogarithmicMap(zmin::Float64, zmax::Float64)
        zmin <= 0.0 && throw(DomainError(zmin, "LogarithmicMap: zmin must be > 0"))
        zmax <= zmin && throw(ArgumentError("LogarithmicMap: zmax must be > zmin"))
        new(zmin, zmax)
    end
end
function forward(m::LogarithmicMap, z)
    z <= 0 && throw(DomainError(z, "LogarithmicMap.forward: z must be > 0"))
    return 2.0 * log(z / m.zmin) / log(m.zmax / m.zmin) - 1.0
end
function inverse(m::LogarithmicMap, ξ)
    return m.zmin * exp(((ξ + 1.0) / 2.0) * log(m.zmax / m.zmin))
end
function dzdξ(m::LogarithmicMap, ξ)
    return 0.5 * log(m.zmax / m.zmin) * inverse(m, ξ)
end
function d2zdξ2(m::LogarithmicMap, ξ)
    return 0.25 * (log(m.zmax / m.zmin)^2) * inverse(m, ξ)
end

# ===== 4. HYPERBOLIC TANGENT MAP =====
struct TanhMap <: CoordinateMap
    zmin::Float64
    zmax::Float64
    α::Float64

    function TanhMap(zmin::Float64, zmax::Float64, α::Float64)
        zmin >= zmax && throw(DomainError((zmin, zmax), "TanhMap: zmin must be < zmax"))
        zmin < 0.0 && throw(DomainError(zmin, "TanhMap: zmin must be ≥ 0"))
        α <= 0.0 && throw(DomainError(α, "TanhMap: α must be > 0"))
        new(zmin, zmax, α)
    end
end
function forward(m::TanhMap, z)
    L = m.zmax - m.zmin
    zc = (m.zmax + m.zmin) / 2.0
    arg = clamp(((z - zc) * tanh(m.α)) / (L / 2.0), -0.9999999999999, 0.9999999999999)
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

# ===== 5. CFD WALL BOUNDARY LAYER MAP =====
struct CfdWallMap <: CoordinateMap
    ztop::Float64
    δ::Float64

    # Direct validation upon object instantiation
    function CfdWallMap(ztop::Float64, δ::Float64)
        if ztop <= 0.0
            throw(DomainError(ztop, "CfdWallMap: ztop must be > 0"))
        end
        if δ <= 0.0
            throw(DomainError(δ, "CfdWallMap: δ must be > 0"))
        end
        new(ztop, δ)
    end
end

function forward(m::CfdWallMap, z)
    arg = clamp((1.0 - z / m.ztop) * tanh(m.δ), -0.9999999999999, 0.9999999999999)
    return 1.0 - (2.0 / m.δ) * atanh(arg)
end

function inverse(m::CfdWallMap, ξ)
    return m.ztop * (1.0 - tanh(0.5 * m.δ * (1.0 - ξ)) / tanh(m.δ))
end

function dzdξ(m::CfdWallMap, ξ)
    return m.ztop * (0.5 * m.δ / tanh(m.δ)) * (sech(0.5 * m.δ * (1.0 - ξ)))^2
end

function d2zdξ2(m::CfdWallMap, ξ)
    return -m.ztop * (0.5 * m.δ^2 / tanh(m.δ)) * (sech(0.5 * m.δ * (1.0 - ξ)))^2 * tanh(0.5 * m.δ * (1.0 - ξ))
end

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

# ===== VALIDATION & ERROR HANDLING =====
# Fallback checks hook cleanly into type constructors where possible
is_valid(m::CoordinateMap) = true

function map_name(m::CoordinateMap)
    return string(nameof(typeof(m)))
end

function map_parameter_pairs(m::CoordinateMap)
    pairs = Vector{Pair{String,String}}()
    for field in fieldnames(typeof(m))
        value = getfield(m, field)
        value_str = value isa AbstractFloat ? string(round(value, sigdigits=8)) : string(value)
        push!(pairs, string(field) => value_str)
    end
    return pairs
end

function describe_map(m::CoordinateMap)
    return (
        map_type = map_name(m),
        parameters = map_parameter_pairs(m),
    )
end

# ===== UTILITY FUNCTIONS & STACK COMPOSITION =====
function profile_transform(m::CoordinateMap, z_profile::Vector{Float64})
    is_valid(m)
    return forward.(Ref(m), z_profile)
end

struct JacobianStack
    maps::Vector{CoordinateMap}
    labels::Vector{String}
    function JacobianStack(maps::Vector{CoordinateMap}, labels::Vector{String}=[])
        isempty(maps) && throw(ArgumentError("JacobianStack: maps vector is empty"))
        !isempty(labels) && length(labels) != length(maps) && throw(ArgumentError("JacobianStack: length mismatch"))
        for m in maps; is_valid(m); end
        isempty(labels) && (labels = ["Transform_$i" for i in 1:length(maps)])
        new(maps, labels)
    end
end

function evaluate_jacobian_stack(stack::JacobianStack; z::Float64)
    coords = Vector{Float64}(undef, length(stack.maps) + 1)
    coords[1] = z
    for i in 1:length(stack.maps)
        coords[i+1] = forward(stack.maps[i], coords[i])
    end
    ξ_final = coords[end]

    dξdz_chain = Vector{Float64}(undef, length(stack.maps))
    for i in 1:length(stack.maps)
        dξdz_chain[i] = 1.0 / dzdξ(stack.maps[i], coords[i+1])
    end
    dξdz_total = prod(dξdz_chain)

    d2ξdz2_total = 0.0
    dzsdz_current = 1.0
    for i in 1:length(stack.maps)
        m = stack.maps[i]
        J1 = dzdξ(m, coords[i+1])
        J2 = d2zdξ2(m, coords[i+1])
        d2ξdz2_local = -J2 / (J1^3)

        downstream_prod = 1.0
        for j in (i+1):length(stack.maps)
            downstream_prod *= dξdz_chain[j]
        end
        d2ξdz2_total += d2ξdz2_local * (dzsdz_current^2) * downstream_prod
        dzsdz_current *= dξdz_chain[i]
    end
    return (ξ_final, dξdz_total, d2ξdz2_total)
end

end # module