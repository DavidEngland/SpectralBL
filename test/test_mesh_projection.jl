using Test
using NCDatasets
using DataFrames
using LinearAlgebra
using Statistics

include("../src/MeshProjection.jl")
using .MeshProjection

@testset "Mesh Projection & Metric Integrity Verification" begin
    smear_heights = Float32[16.8, 33.6, 67.2, 125.0]
    total_span = smear_heights[end] - smear_heights[1]

    @testset "Analytical Mass Matrix Structure" begin
        m = compute_mesh_mass_matrix(smear_heights)
        @test size(m) == (4, 4)
        @test issymmetric(m)

        # For a 4x4 tridiagonal mass matrix, these extreme off-diagonals must be zero.
        @test m[1, 3] == 0.0f0
        @test m[1, 4] == 0.0f0
        @test m[2, 4] == 0.0f0

        @test_throws AssertionError compute_mesh_mass_matrix(Float32[10.0])
        @test_throws AssertionError compute_mesh_mass_matrix(Float32[20.0, 10.0])
        @test_throws AssertionError compute_mesh_mass_matrix(Float32[10.0, 10.0, 30.0])
    end

    @testset "Physical Volume Preservation" begin
        m = compute_mesh_mass_matrix(smear_heights)
        u_uniform = ones(Float32, 4)
        continuous_integral = (u_uniform' * m * u_uniform)[1]
        @test continuous_integral ≈ total_span rtol=1e-5
    end

    @testset "Metric Operator Spectral Decomposition" begin
        m = compute_mesh_mass_matrix(smear_heights)
        m_half = compute_metric_weights(m)

        @test size(m_half) == (4, 4)
        @test issymmetric(m_half)
        @test m_half * m_half ≈ m rtol=1e-5
    end

    @testset "NetCDF Ingestion & Column Filtering" begin
        mktempdir() do tmpdir
            test_nc = joinpath(tmpdir, "test_synthetic_snapshots.nc")

            num_heights = length(smear_heights)
            num_times = 3

            synthetic_data = Float32[
                3.0   4.5   -9999.0;
                4.0   NaN32  6.0;
                5.0   5.2    7.0;
                6.0   6.1    8.0
            ]

            NCDataset(test_nc, "c") do ds
                defDim(ds, "height", num_heights)
                defDim(ds, "time", num_times)

                sh = defVar(ds, "sensor_height", Float32, ("height",))
                sh[:] = smear_heights

                v = defVar(ds, "wind_speed", Float32, ("height", "time"); fillvalue=Float32(-9999.0))
                v[:] = synthetic_data
            end

            y_raw, y_tilde, m_half = project_snapshots_to_metric_space(test_nc, "wind_speed")

            @test size(y_raw, 2) == 1
            @test size(y_tilde, 2) == 1
            @test y_raw[:, 1] == Float32[3.0, 4.0, 5.0, 6.0]
            @test y_tilde ≈ m_half * y_raw[:, 1:1]
        end
    end
end
