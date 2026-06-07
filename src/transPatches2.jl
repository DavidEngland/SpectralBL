# src/transforms.jl
module Transforms

using LinearAlgebra

export CoordinateMap, LinearMap, HyperbolicMap, LogarithmicMap, TanhMap, CfdWallMap, CustomMap,
       forward, inverse, dzdξ, d2zdξ2, dξdz, d2ξdz2, is_valid,
       profile_transform, JacobianStack, evaluate_jacobian_stack

abstract type CoordinateMap end

# ===== STAGE 0: CORE ABSTRACT METRIC INTERFACES =====
function dξdz(m::CoordinateMap, z) ... end
function d2ξdz2(m::CoordinateMap, z) ... end

# ===== 1. LINEAR MAP =====
struct LinearMap <: CoordinateMap ... end
# functions...

# ===== 2. HYPERBOLIC MAP =====
struct HyperbolicMap <: CoordinateMap ... end
# functions...

# ===== 3. LOGARITHMIC MAP =====
struct LogarithmicMap <: CoordinateMap ... end
# functions...

# ===== 4. HYPERBOLIC TANGENT MAP =====
struct TanhMap <: CoordinateMap ... end
# functions...

# ===== 5. CFD WALL BOUNDARY LAYER MAP =====
struct CfdWallMap <: CoordinateMap ... end
# functions...

# ===== 6. CUSTOM MAP =====
struct CustomMap ... end


# ===== ===== VALIDATION & ERROR HANDLING (MOVED HERE) ===== =====
# Moving this below the structs allows the compiler to resolve type definitions!

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
    is_valid(LinearMap(m.zmin, m.zmax)) # Compiler now perfectly knows what LinearMap is!
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


# ===== UTILITY FUNCTIONS & STACK COMPOSITION =====
# profile_transform, JacobianStack, evaluate_jacobian_stack ...

end # module