module Richardson

using LinearAlgebra

export compute_chebyshev_nodes_and_matrix,
       build_physical_gradient_operator,
       calculate_pseudospectral_ri_g,
       build_and_evaluate_ri_g

"""
    compute_chebyshev_nodes_and_matrix(N::Int)

Generate Chebyshev-Gauss-Lobatto nodes `xi in [-1, 1]` and the
spectral differentiation matrix `D_xi` of size `(N+1, N+1)`.
"""
function compute_chebyshev_nodes_and_matrix(N::Int)
    N < 2 && throw(ArgumentError("N must be >= 2"))

    xi = [cos(pi * j / N) for j in 0:N]
    D_xi = zeros(Float64, N + 1, N + 1)
    c = [j == 0 || j == N ? 2.0 : 1.0 for j in 0:N]

    for i in 0:N, j in 0:N
        if i == j
            if i == 0
                D_xi[i + 1, j + 1] = (2.0 * N^2 + 1.0) / 6.0
            elseif i == N
                D_xi[i + 1, j + 1] = -(2.0 * N^2 + 1.0) / 6.0
            else
                D_xi[i + 1, j + 1] = -xi[i + 1] / (2.0 * (1.0 - xi[i + 1]^2))
            end
        else
            D_xi[i + 1, j + 1] = (c[i + 1] / c[j + 1]) * ((-1.0)^(i + j)) / (xi[i + 1] - xi[j + 1])
        end
    end

    return xi, D_xi
end

"""
    build_physical_gradient_operator(N::Int, z_min::Real, z_max::Real, alpha::Real)

Build the physical coordinate `z`, analytical Jacobian `J = dz/dxi`, and
manifold-aware first-derivative matrix `D_z = diag(1 ./ J) * D_xi` using
an analytical hyperbolic tangent map.
"""
function build_physical_gradient_operator(N::Int, z_min::Real, z_max::Real, alpha::Real)
    z_min >= z_max && throw(ArgumentError("z_min must be < z_max"))
    alpha <= 0 && throw(ArgumentError("alpha must be > 0"))

    xi, D_xi = compute_chebyshev_nodes_and_matrix(N)

    mid_z = 0.5 * (float(z_max) + float(z_min))
    delta_z = 0.5 * (float(z_max) - float(z_min))

    z = mid_z .+ delta_z .* (tanh.(float(alpha) .* xi) ./ tanh(float(alpha)))
    J = delta_z .* (float(alpha) / tanh(float(alpha))) .* (1.0 .- tanh.(float(alpha) .* xi).^2)

    # For valid monotone maps, Jacobian must remain strictly positive.
    minimum(J) <= 0.0 && throw(DomainError(minimum(J), "Jacobian must be positive across all nodes"))

    D_z = Diagonal(1.0 ./ J) * D_xi
    return z, J, D_z, xi, D_xi
end

"""
    calculate_pseudospectral_ri_g(c_theta, c_u, D_z; g=9.81, theta_ref=293.15, shear_floor=1e-8)

Evaluate `dtheta/dz`, `du/dz`, and gradient Richardson number profile using
matrix operators and element-wise vector algebra.

Returns a named tuple with fields:
- `dtheta_dz`
- `du_dz`
- `Ri_g`
"""
function calculate_pseudospectral_ri_g(
    c_theta::AbstractVector{<:Real},
    c_u::AbstractVector{<:Real},
    D_z::AbstractMatrix{<:Real};
    g::Real=9.81,
    theta_ref::Real=293.15,
    shear_floor::Real=1e-8,
)
    length(c_theta) == length(c_u) || throw(ArgumentError("c_theta and c_u must have equal length"))
    size(D_z, 1) == size(D_z, 2) || throw(ArgumentError("D_z must be square"))
    size(D_z, 1) == length(c_theta) || throw(ArgumentError("D_z size must match coefficient vector length"))
    theta_ref <= 0 && throw(ArgumentError("theta_ref must be > 0"))
    shear_floor <= 0 && throw(ArgumentError("shear_floor must be > 0"))

    dtheta_dz = D_z * Float64.(c_theta)
    du_dz = D_z * Float64.(c_u)

    denom = max.(du_dz .^ 2, float(shear_floor)^2)
    Ri_g = (float(g) / float(theta_ref)) .* dtheta_dz ./ denom

    return (dtheta_dz=dtheta_dz, du_dz=du_dz, Ri_g=Ri_g)
end

"""
    build_and_evaluate_ri_g(N, z_min, z_max, alpha, c_theta, c_u; kwargs...)

Convenience wrapper that builds manifold operators and evaluates Ri_g in one call.
Returns a named tuple containing `z`, `J`, `D_z`, `xi`, `D_xi`, `dtheta_dz`,
`du_dz`, and `Ri_g`.
"""
function build_and_evaluate_ri_g(
    N::Int,
    z_min::Real,
    z_max::Real,
    alpha::Real,
    c_theta::AbstractVector{<:Real},
    c_u::AbstractVector{<:Real};
    g::Real=9.81,
    theta_ref::Real=293.15,
    shear_floor::Real=1e-8,
)
    z, J, D_z, xi, D_xi = build_physical_gradient_operator(N, z_min, z_max, alpha)
    prof = calculate_pseudospectral_ri_g(c_theta, c_u, D_z; g=g, theta_ref=theta_ref, shear_floor=shear_floor)
    return (
        z=z,
        J=J,
        D_z=D_z,
        xi=xi,
        D_xi=D_xi,
        dtheta_dz=prof.dtheta_dz,
        du_dz=prof.du_dz,
        Ri_g=prof.Ri_g,
    )
end

end # module
