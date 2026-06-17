# test/test_smear_fetch.jl
using Test
using DataFrames
using NCDatasets

include("../src/smear_fetch.jl")

function create_dummy_smear_dataframe()
    samptime = [
        "2026-01-01T00:00:00.000Z",
        "2026-01-01T00:00:00.000Z",
        "2026-01-01T00:00:00.000Z",
        "2026-01-01T00:00:00.000Z",
        "2026-01-01T00:00:00.000Z",
        "2026-01-01T00:00:00.000Z",
        "2026-01-01T00:00:00.000Z",
        "2026-01-01T00:00:00.000Z",
        "2026-01-01T00:00:00.000Z",
        "2026-01-01T00:10:00.000Z",
        "2026-01-01T00:10:00.000Z",
        "2026-01-01T00:10:00.000Z",
        "2026-01-01T00:10:00.000Z",
        "2026-01-01T00:10:00.000Z",
        "2026-01-01T00:10:00.000Z",
        "2026-01-01T00:10:00.000Z",
        "2026-01-01T00:10:00.000Z",
    ]

    tablevariable = [
        "HYY_META.WS168",
        "HYY_META.WS336",
        "HYY_META.WS672",
        "HYY_META.WS125",
        "HYY_META.T168",
        "HYY_META.T336",
        "HYY_META.T672",
        "HYY_META.T125",
        "HYY_META.P_R",
        "HYY_META.WS168",
        "HYY_META.WS336",
        "HYY_META.WS125",
        "HYY_META.T168",
        "HYY_META.T336",
        "HYY_META.T672",
        "HYY_META.T125",
        "HYY_META.P_R",
    ]

    value = [
        3.2,
        4.1,
        5.0,
        5.6,
        12.5,
        12.0,
        11.5,
        11.0,
        1008.0,
        3.4,
        4.0,
        5.7,
        12.7,
        12.2,
        11.8,
        11.1,
        1007.5,
    ]

    return DataFrame(samptime=samptime, tablevariable=tablevariable, value=value)
end

@testset "smear ingestion pivot and export" begin
    raw_df = create_dummy_smear_dataframe()
    mktempdir() do dir
        cd(dir) do
            process_and_save_profiles(raw_df)

            out_path = joinpath(dir, "smear_ii_processed_profiles.nc")
            @test isfile(out_path)

            NCDataset(out_path, "r") do ds
                @test haskey(ds, "wind_speed")
                @test haskey(ds, "potential_temperature")
                @test haskey(ds, "air_pressure")
                @test haskey(ds, "sensor_height")

                @test size(ds["wind_speed"]) == (4, 2)
                @test size(ds["potential_temperature"]) == (4, 2)
                @test size(ds["air_pressure"]) == (4, 2)
                @test ds["sensor_height"][:] == SMARTSMEAR_HEIGHTS
                @test ds["wind_speed"].attrib["_FillValue"] == SMARTSMEAR_FILL_VALUE
                @test ds.attrib["observer_type"] == "tower_profile"
                @test ds.attrib["coordinate_system"] == "physical_height"

                @test size(ds["wind_speed"]) == (4, 2)
                @test size(ds["potential_temperature"]) == (4, 2)
                @test size(ds["air_pressure"]) == (4, 2)
            end
        end
    end
end

@testset "smear partial temperature dropout handling" begin
    raw_df = create_dummy_smear_dataframe()
    allowmissing!(raw_df, :value)

    # Simulate a single-level temperature dropout at the first timestamp.
    idx = findfirst(i -> raw_df.samptime[i] == "2026-01-01T00:00:00.000Z" && raw_df.tablevariable[i] == "HYY_META.T125", eachindex(raw_df.samptime))
    @test idx !== nothing
    raw_df.value[idx] = missing

    mktempdir() do dir
        cd(dir) do
            process_and_save_profiles(raw_df)

            out_path = joinpath(dir, "smear_ii_processed_profiles.nc")
            @test isfile(out_path)

            NCDataset(out_path, "r") do ds
                pot_raw = ds["potential_temperature"].var[:, :]
                wind_raw = ds["wind_speed"].var[:, :]

                # Find a timestamp where all wind levels are valid (no fill sentinels).
                full_wind_col = findfirst(t -> count(==(SMARTSMEAR_FILL_VALUE), wind_raw[:, t]) == 0, 1:size(wind_raw, 2))
                @test full_wind_col !== nothing

                # One missing temperature level should not force all levels to fill value.
                @test count(==(SMARTSMEAR_FILL_VALUE), pot_raw[:, full_wind_col]) == 1
                @test count(==(SMARTSMEAR_FILL_VALUE), pot_raw[:, full_wind_col]) < size(pot_raw, 1)
            end
        end
    end
end
