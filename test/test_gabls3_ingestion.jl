#test/test_gabls3_ingestion.jl
using Test
using NCDatasets

include("../src/CasesIngestion.jl")
include("../src/GabLS3Ingestion.jl")

using UnifiedManifold
using .GabLS3Ingestion

function create_dummy_gabls3_netcdf()
    path = tempname() * ".nc"

    NCDataset(path, "c") do ds
        defDim(ds, "time", 3)
        defDim(ds, "zf", 6)
        defDim(ds, "zh", 4)

        t = defVar(ds, "time", Float64, ("time",))
        zf = defVar(ds, "zf", Float64, ("zf", "time"))
        zh = defVar(ds, "zh", Float64, ("zh", "time"))

        u = defVar(ds, "u", Float64, ("zf", "time"))
        v = defVar(ds, "v", Float64, ("zf", "time"))
        th = defVar(ds, "th", Float64, ("zf", "time"))
        wt = defVar(ds, "wt", Float64, ("zh", "time"))
        uw = defVar(ds, "uw", Float64, ("zh", "time"))
        vw = defVar(ds, "vw", Float64, ("zh", "time"))

        # Epoch hours matching sim-relative hours 7.25, 8.25, 10.0
        # (simulation starts at epoch hour 12.0).
        t[:] = [19.25, 20.25, 22.0]

        zf_base = [2.0, 10.0, 20.0, 40.0, 80.0, 120.0]
        zh_base = [5.0, 20.0, 60.0, 100.0]

        for ti in 1:3
            zf[:, ti] = zf_base
            zh[:, ti] = zh_base

            for zi in 1:length(zf_base)
                z = zf_base[zi]
                th[zi, ti] = 289.0 + 0.03 * z + 0.2 * ti
                u[zi, ti] = 2.0 + 0.02 * z + 0.1 * ti
                v[zi, ti] = -1.0 + 0.01 * z
            end

            for zi in 1:length(zh_base)
                z = zh_base[zi]
                wt[zi, ti] = 0.02 + 0.0002 * z
                uw[zi, ti] = -0.10 + 0.0005 * z
                vw[zi, ti] = 0.01 + 0.0003 * z
            end
        end
    end

    return path
end

@testset "GABLS3 window tagging" begin
    # Epoch hours 19-20 and 20-21 correspond to sim hours 7-8 and 8-9.
    @test infer_hour_window(19.1) == "HOUR8"
    @test infer_hour_window(20.5) == "HOUR9"
    @test infer_hour_window(12.0) == "OTHER"
    @test infer_hour_window(35.0) == "OTHER"
end

@testset "GABLS3 ingestion and projection" begin
    nc_path = create_dummy_gabls3_netcdf()
    ws = UnifiedManifoldWorkspace(16, 0.0, 150.0, 0.05; invert_windows=true)

    out = ingest_and_project_gabls3_slice!(nc_path, 1, ws)

    @test out !== nothing
    @test length(out.c_theta) == ws.N + 1
    @test length(out.c_u) == ws.N + 1
    @test length(out.z_mean) == 6
    @test length(out.z_flux) == 4
    @test out.window_tag == "HOUR8"   # t[1]=19.25 → sim_hour=7.25 → HOUR8
    @test all(isfinite, out.c_theta)
    @test all(isfinite, out.c_u)

    rm(nc_path; force=true)
end
