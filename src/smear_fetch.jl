# src/smear_fetch.jl
using HTTP
using JSON3
using DataFrames
using NCDatasets
using Dates
using LinearAlgebra
using Statistics

const SMARTSMEAR_HEIGHTS = Float32[16.8, 33.6, 67.2, 125.0]
const SMARTSMEAR_FILL_VALUE = Float32(-9999.0)
const GRAVITY = 9.80665f0
const R_D = 287.05f0

function _finite_float32(value)
    if ismissing(value)
        return NaN32
    end
    numeric = Float32(value)
    return isfinite(numeric) ? numeric : NaN32
end

function _group_value(gdf::AbstractDataFrame, variable_name::AbstractString)
    if !(:tablevariable in names(gdf)) || !(:value in names(gdf))
        return NaN32
    end

    idx = findfirst(==(variable_name), gdf.tablevariable)
    idx === nothing && return NaN32
    return _finite_float32(gdf.value[idx])
end

function _hydrostatic_pressure_profile(p_station::Float32, heights::AbstractVector{<:Real}, temperature_k::AbstractVector{<:Real})
    finite_t = filter(isfinite, Float32.(temperature_k))
    temp_ref = isempty(finite_t) ? 288.15f0 : Float32(mean(finite_t))
    return Float32[p_station * exp(-(GRAVITY * Float32(z)) / (R_D * temp_ref)) for z in heights]
end

function _cf_fill_values(array::AbstractArray{Float32})
    return replace(array, NaN32 => SMARTSMEAR_FILL_VALUE)
end

"""
    fetch_smear_data(start_date::DateTime, end_date::DateTime)

Fetches raw boundary layer tower profiles from the SmartSMEAR API with network error handling.
"""
function fetch_smear_data(start_date::DateTime, end_date::DateTime)
    endpoint = "https://smear-backend.rahtiapp.fi/search/timeseries"

    # Format datetimes to match API expectations
    from_str = Dates.format(start_date, "yyyy-mm-ddTHH:MM:SS.000")
    to_str   = Dates.format(end_date, "yyyy-mm-ddTHH:MM:SS.000")

    payload = Dict(
        "tablevariable" => [
            "HYY_META.WS125", "HYY_META.WS672", "HYY_META.WS336", "HYY_META.WS168",
            "HYY_META.T125",  "HYY_META.T672",  "HYY_META.T336",  "HYY_META.T168",
            "HYY_META.P_R" # Adding reference pressure for true potential temperature scaling
        ],
        "from" => from_str,
        "to" => to_str,
        "quality" => "ANY",
        "aggregation" => "NONE"
    )

    headers = ["Content-Type" => "application/json"]

    try
        println("Requesting SMEAR II profile data from $from_str to $to_str...")
        response = HTTP.post(endpoint, headers, JSON3.write(payload), retry=true, redirect=true)

        if response.status == 200
            parsed_json = JSON3.read(response.body)
            return DataFrame(parsed_json[:data])
        else
            error("SMEAR API returned non-200 status code: $(response.status)")
        end
    catch e
        println("Critical Network or Parse Failure in SMEAR Ingestion Loop.")
        rethrow(e)
    end
end

"""
    process_and_save_profiles(raw_df::DataFrame)

Pivots long-format records, handles NaN missing states, corrects to potential temperature,
and writes a CF-compliant NetCDF array.
"""
function process_and_save_profiles(raw_df::DataFrame)
    if isempty(raw_df)
        error("Input DataFrame is empty. Ingestion aborted.")
    end

    @assert issorted(SMARTSMEAR_HEIGHTS)

    # 1. Clean timestamps and normalize string keys
    raw_df.timestamp = DateTime.(replace.(raw_df.samptime, "Z" => ""), dateformat"yyyy-mm-ddTHH:MM:SS.s")

    # 2. Split-Apply-Combine to pivot from API long format to analytical wide format
    println("Executing split-apply-combine to resolve long-format rows...")
    df_wide = combine(groupby(raw_df, :timestamp)) do gdf
        (
            WS168 = _group_value(gdf, "HYY_META.WS168"),
            WS336 = _group_value(gdf, "HYY_META.WS336"),
            WS672 = _group_value(gdf, "HYY_META.WS672"),
            WS125 = _group_value(gdf, "HYY_META.WS125"),
            T168  = _group_value(gdf, "HYY_META.T168"),
            T336  = _group_value(gdf, "HYY_META.T336"),
            T672  = _group_value(gdf, "HYY_META.T672"),
            T125  = _group_value(gdf, "HYY_META.T125"),
            P_surf = _group_value(gdf, "HYY_META.P_R")
        )
    end

    sort!(df_wide, :timestamp)

    # 3. Establish Dimensions
    heights = copy(SMARTSMEAR_HEIGHTS)
    num_heights = length(heights)
    num_times = nrow(df_wide)
    time_coords = [datetime2unix(t) for t in df_wide.timestamp]

    # 4. Allocate Tensors (Height × Time for column-major efficiency)
    wind_speed = fill(NaN32, num_heights, num_times)
    pot_temp   = fill(NaN32, num_heights, num_times)
    pressure_profile = fill(NaN32, num_heights, num_times)

    p_0 = 1000.0f0  # Standard reference pressure (hPa)
    R_Cp = 0.286f0 # Poisson constant for dry air

    println("Computing potential temperature corrections and mapping NaNs...")
    for t in 1:num_times
        row = df_wide[t, :]

        # Assign wind speed profile directly (preserving NaNs)
        wind_speed[1, t] = row.WS168
        wind_speed[2, t] = row.WS336
        wind_speed[3, t] = row.WS672
        wind_speed[4, t] = row.WS125

        # Extract temperatures in Kelvin
        tk1 = row.T168 + 273.15f0
        tk2 = row.T336 + 273.15f0
        tk3 = row.T672 + 273.15f0
        tk4 = row.T125 + 273.15f0

        if any(isnan, (tk1, tk2, tk3, tk4))
            continue
        end

        # Use a height-aware hydrostatic pressure profile so θ remains thermodynamically consistent.
        p_station = isnan(row.P_surf) ? 1013.25f0 : row.P_surf
        temp_profile = Float32[tk1, tk2, tk3, tk4]
        p_profile = _hydrostatic_pressure_profile(p_station, heights, temp_profile)
        pressure_profile[:, t] .= p_profile

        for i in 1:num_heights
            pot_temp[i, t] = temp_profile[i] * (p_0 / p_profile[i])^R_Cp
        end
    end

    # 5. Build NetCDF Storage
    filename = "smear_ii_processed_profiles.nc"
    if isfile(filename) rm(filename) end

    println("Writing validated NetCDF data array...")

    NCDataset(filename, "c") do ds
        defDim(ds, "height", num_heights)
        defDim(ds, "time", num_times)

        sensor_height = defVar(ds, "sensor_height", Float32, ("height",))
        sensor_height[:] = heights
        sensor_height.attrib["units"] = "m"
        sensor_height.attrib["standard_name"] = "height"

        time_var = defVar(ds, "time", Float64, ("time",))
        time_var[:] = time_coords
        time_var.attrib["units"] = "seconds since 1970-01-01 00:00:00"

        wind_var = defVar(ds, "wind_speed", Float32, ("height", "time"); fillvalue=SMARTSMEAR_FILL_VALUE)
        wind_var[:] = _cf_fill_values(wind_speed)
        wind_var.attrib["units"] = "m s-1"
        wind_var.attrib["standard_name"] = "wind_speed"

        theta_var = defVar(ds, "potential_temperature", Float32, ("height", "time"); fillvalue=SMARTSMEAR_FILL_VALUE)
        theta_var[:] = _cf_fill_values(pot_temp)
        theta_var.attrib["units"] = "K"
        theta_var.attrib["standard_name"] = "air_potential_temperature"

        pressure_var = defVar(ds, "air_pressure", Float32, ("height", "time"); fillvalue=SMARTSMEAR_FILL_VALUE)
        pressure_var[:] = _cf_fill_values(pressure_profile)
        pressure_var.attrib["units"] = "hPa"
        pressure_var.attrib["standard_name"] = "air_pressure"

        ds.attrib["title"] = "SMEAR II Processed Boundary Layer Profile Space"
        ds.attrib["institution"] = "University of Helsinki / SpectralBL-Analytics"
        ds.attrib["comment"] = "Pivoted long-format array with height-aware potential temperature and CF fill values."
        ds.attrib["observer_type"] = "tower_profile"
        ds.attrib["coordinate_system"] = "physical_height"
        ds.attrib["recommended_projection"] = "mass_weighted_pFEM"
    end

    println("File compiled: ", filename)
end

function main()
    # 6. Operational Execution (Using validated past windows relative to June 2026)
    start_dt = DateTime(2026, 01, 01, 0, 0, 0)
    end_dt   = DateTime(2026, 01, 07, 23, 59, 59)

    raw_data = fetch_smear_data(start_dt, end_dt)
    process_and_save_profiles(raw_data)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end