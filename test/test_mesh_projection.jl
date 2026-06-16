using Test
using NCDatasets

include("../src/MeshProjection.jl")
using .MeshProjection

function create_dummy_mesh_projection_netcdf()
    path = tempname() * ".nc"
    heights = Float32[16.8, 33.6, 67.2, 125.0]
    fill_val = Float32(-9999.0)

    NCDataset(path, "c") do ds
        defDim(ds, "height", length(heights))
        defDim(ds, "time", 3)

        z = defVar(ds, "sensor_height", Float32, ("height",))
        z[:] = heights

        w = defVar(ds, "wind_speed", Float32, ("height", "time"); fillvalue=fill_val)

        data = Float32[
            3.2  fill_val  4.0;
            4.1  fill_val  4.8;
            5.0  fill_val  5.5;
            5.6  fill_val  6.2
        ]
        w[:, :] = data
        w.attrib["units"] = "m s-1"
        w.attrib["standard_name"] = "wind_speed"
    end

    return path
end

@testset "Mesh Projection & Metric Integrity Tests" begin
    heights = Float32[16.8, 33.6, 67.2, 125.0]
    total_span = heights[end] - heights[1]

    m = compute_mesh_mass_matrix(heights)
    @test size(m) == (4, 4)
    @test issymmetric(m)

    u_uniform = ones(Float32, 4)
    continuous_integral = (u_uniform' * m * u_uniform)[1]
    @test continuous_integral ≈ total_span rtol=1e-5

    m_half = compute_metric_weights(m)
    @test m_half * m_half ≈ m rtol=1e-4
end

@testset "Metric Projection Snapshot Filtering" begin
    nc_path = create_dummy_mesh_projection_netcdf()

    y_raw, y_tilde, m_half = project_snapshots_to_metric_space(nc_path, "wind_speed")

    @test size(y_raw) == (4, 2)
    @test size(y_tilde) == (4, 2)
    @test size(m_half) == (4, 4)
    @test all(isfinite, y_raw)
    @test all(isfinite, y_tilde)

    rm(nc_path; force=true)
end
