# src/MeshProjection.jl
module MeshProjection

using NCDatasets
using LinearAlgebra

export compute_mesh_mass_matrix, compute_metric_weights, project_snapshots_to_metric_space

const DEFAULT_FILL_VALUE = Float32(-9999.0)

"""
    compute_mesh_mass_matrix(heights)::Matrix{Float32}

Compute the analytical 1D hat-function mass matrix on a nonuniform vertical mesh.
The input heights must be strictly increasing.
"""
function compute_mesh_mass_matrix(heights::AbstractVector{<:Real})::Matrix{Float32}
    z = Float32.(collect(heights))
    m = length(z)

    @assert m >= 2 "Tower profile must contain at least 2 active sensor heights."
    @assert issorted(z) "Heights must be strictly sorted in ascending order."

    h = diff(z)
    @assert all(>(0.0f0), h) "Heights must be strictly increasing with no duplicates."

    mass_matrix = zeros(Float32, m, m)

    for i in 1:m
        if i == 1
            mass_matrix[1, 1] = h[1] / 3.0f0
            mass_matrix[1, 2] = h[1] / 6.0f0
        elseif i == m
            mass_matrix[m, m] = h[end] / 3.0f0
            mass_matrix[m, m - 1] = h[end] / 6.0f0
        else
            mass_matrix[i, i] = (h[i - 1] + h[i]) / 3.0f0
            mass_matrix[i, i - 1] = h[i - 1] / 6.0f0
            mass_matrix[i, i + 1] = h[i] / 6.0f0
        end
    end

    return mass_matrix
end

"""
    compute_metric_weights(mass_matrix)::Matrix{Float32}

Compute the symmetric square root M^(1/2) from a symmetric positive semidefinite
mass matrix using eigendecomposition.
"""
function compute_metric_weights(mass_matrix::AbstractMatrix{<:Real})::Matrix{Float32}
    m = Float64.(mass_matrix)
    @assert size(m, 1) == size(m, 2) "Mass matrix must be square."

    symmetry_error = norm(m - m', Inf)
    @assert symmetry_error <= 1e-8 "Mass matrix must be symmetric."

    eig = eigen(Hermitian(m))
    min_eig = minimum(eig.values)
    @assert min_eig >= -1e-8 "Mass matrix must be positive semidefinite within tolerance."

    clamped_values = max.(eig.values, 0.0)
    m_half = eig.vectors * Diagonal(sqrt.(clamped_values)) * eig.vectors'

    return Float32.(m_half)
end

"""
    project_snapshots_to_metric_space(nc_filename, variable_name)

Load tower snapshots from NetCDF, drop invalid frames containing fill values or NaNs,
and return raw plus metric-weighted snapshots.
"""
function project_snapshots_to_metric_space(nc_filename::AbstractString, variable_name::AbstractString)
    NCDataset(nc_filename, "r") do ds
        haskey(ds, "sensor_height") || error("Missing required variable: sensor_height")
        haskey(ds, variable_name) || error("Missing required variable: $(variable_name)")

        raw_heights = Float32.(coalesce.(ds["sensor_height"][:], NaN32))
        raw_data = Float32.(coalesce.(ds[variable_name][:, :], NaN32))

        fill_val_attr = get(ds[variable_name].attrib, "_FillValue", DEFAULT_FILL_VALUE)
        fill_val = Float32(fill_val_attr)

        valid_indices = Int[]
        num_times = size(raw_data, 2)

        for t in 1:num_times
            column = raw_data[:, t]
            if all(x -> isfinite(x) && x != fill_val, column)
                push!(valid_indices, t)
            end
        end

        y_raw = raw_data[:, valid_indices]

        mass_matrix = compute_mesh_mass_matrix(raw_heights)
        m_half = compute_metric_weights(mass_matrix)
        y_tilde = m_half * y_raw

        return y_raw, y_tilde, m_half
    end
end

end # module
