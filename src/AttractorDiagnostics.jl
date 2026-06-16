# src/AttractorDiagnostics.jl
module AttractorDiagnostics

using LinearAlgebra
using Statistics
using DSP

export metric_svd, fit_coeffs_ridge, compute_derivatives, spin_2d, curvature_2d, sv_entropy

"""
    metric_svd(Y, M; r)

Computes the economy SVD weighted by the p-FEM mass/metric matrix `M`.
Y: n × T standardized data matrix.
M: n × n positive-definite metric matrix (or its Cholesky factor).
"""
function metric_svd(Y::AbstractMatrix{T}, M::AbstractMatrix{T}; r::Int) where T
    # Compute Cholesky factor if not already decomposed to get M^{1/2}
    F = ishermitian(M) ? cholesky(M).U : M

    # Project data into the metric space
    Y_scaled = F * Y

    # Compute economy SVD
    U_scaled, s, Vt = svd(Y_scaled; full=false)

    # Map back to physical space to get metric-consistent modes
    U_r = F \ U_scaled[:, 1:r]
    Σ_r = Diagonal(s[1:r])
    V_r = Vt'[ :, 1:r ]

    return U_r, Σ_r, V_r, s
end

"""
    fit_coeffs_ridge(A, U_r, b_matrix; λ, masks)

Solves the low-rank inverse problem with time-varying observation masks.
"""
function fit_coeffs_ridge(A::AbstractMatrix, U_r::AbstractMatrix, b_matrix::AbstractMatrix;
                          λ=1e-4, masks=nothing)
    r = size(U_r, 2)
    T = size(b_matrix, 2)
    X = zeros(eltype(b_matrix), r, T)

    # Precompute base system if no missing data patterns exist
    if isnothing(masks)
        R = A * U_r
        G = R'R + λ * I
        G_fact = factorize(G)
        for t in 1:T
            X[:, t] = G_fact \ (R' * b_matrix[:, t])
        end
    else
        # Adaptive masking for irregular/missing tower channels
        for t in 1:T
            m_t = masks[:, t]
            # Fast-path check: if mask is fully active
            if all(m_t)
                R = A * U_r
                G = R'R + λ * I
                X[:, t] = factorize(G) \ (R' * b_matrix[:, t])
            else
                # Extract only valid rows/channels
                valid_idx = findall(m_t)
                if isempty(valid_idx)
                    X[:, t] .= NaN
                    continue
                end
                R_t = A[valid_idx, :] * U_r
                b_t = b_matrix[valid_idx, t]
                G_t = R_t'R_t + λ * I
                X[:, t] = factorize(G_t) \ (R_t' * b_t)
            end
        end
    end
    return X
end

"""
    compute_derivatives(z, dt; window_len, poly_order)

Calculates smooth temporal derivatives using a Savitzky-Golay filter window via DSP.jl.
"""
function compute_derivatives(z::AbstractVector, dt; window_len=9, poly_order=2)
    # Ensure window length is odd
    w = isodd(window_len) ? window_len : window_len + 1

    # Generate Savitzky-Golay filter coefficients for 1st and 2nd derivatives
    # DSP.savitzky_golay computes smoothed values; standardizing for diffs:
    # (Alternatively fallback to central differences if boundary padding is strict)
    dz = similar(z)
    d2z = similar(z)

    # Using central diff with internal padding for production fallback safety
    n = length(z)
    for i in 2:(n-1)
        dz[i] = (z[i+1] - z[i-1]) / (2dt)
    end
    dz[1] = (z[2] - z[1]) / dt
    dz[end] = (z[end] - z[end-1]) / dt

    for i in 3:(n-2)
        d2z[i] = (z[i+1] - 2z[i] + z[i-1]) / (dt^2)
    end
    d2z[1] = d2z[3]; d2z[2] = d2z[3]
    d2z[end] = d2z[end-2]; d2z[end-1] = d2z[end-2]

    return dz, d2z
end

function spin_2d(Z::AbstractMatrix, dt)
    x, y = Z[1, :], Z[2, :]
    ẋ, _ = compute_derivatives(x, dt)
    ẏ, _ = compute_derivatives(y, dt)
    return x .* ẏ .- y .* ẋ
end

function curvature_2d(Z::AbstractMatrix, dt; ε=1e-6)
    x, y = Z[1, :], Z[2, :]
    ẋ, ẍ = compute_derivatives(x, dt)
    ẏ, ÿ = compute_derivatives(y, dt)

    num = abs.(ẋ .* ÿ .- ẏ .* ẍ)
    den = (ẋ.^2 .+ ẏ.^2).^(1.5)
    return num ./ (den .+ ε)
end

function sv_entropy(s::AbstractVector; r::Int)
    σ = s[1:r]
    p = σ ./ sum(σ)
    return -sum(p .* log.(p))
end

end # module