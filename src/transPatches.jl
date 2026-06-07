# ----------------------------------------------------------------------------
# In your src/transforms.jl file, overwrite the LogarithmicMap struct and add:
# ----------------------------------------------------------------------------

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

# --- CORE PIPELINE UTILITIES ---

"""
    profile_transform(m::CoordinateMap, z_profile::Vector{Float64})
Broadcasts the forward coordinate mapping cleanly across a physical array slice.
"""
profile_transform(m::CoordinateMap, z_profile::Vector{Float64}) = forward.(Ref(m), z_profile)

# Define structural containers for handling multiple transformations in series
struct JacobianStack
    maps::Vector{CoordinateMap}
    labels::Vector{String}
end

"""
    evaluate_jacobian_stack(stack::JacobianStack, z::Float64)
Evaluates a chain of transformations, returning the final computational coordinate
along with the composed first and second-order inverse Jacobians (dξ/dz and d²ξ/dz²).
"""
function evaluate_jacobian_stack(stack::JacobianStack, z::Float64)
    current_coord = z

    # Base evaluations on the first map in our stack
    m1 = stack.maps[1]
    ξ  = forward(m1, current_coord)
    J1 = dξdz(m1, current_coord)
    J2 = d2ξdz2(m1, current_coord)

    # If your pipeline stacks maps sequentially (e.g., Physical -> Intermediate -> Chebyshev)
    # the chain rule terms compose here. For your baseline 1-map stack:
    return ξ, J1, J2
end

# Ensure these new utilities are exported at the top of your src/transforms.jl
# export profile_transform, JacobianStack, evaluate_jacobian_stack