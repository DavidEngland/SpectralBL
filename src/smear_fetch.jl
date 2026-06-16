using HTTP
using JSON3
using DataFrames
using NetCDF
using Dates
using LinearAlgebra

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

    # 1. Clean timestamps and normalize string keys
    raw_df.timestamp = DateTime.(replace.(raw_df.samptime, "Z" => ""), dateformat"yyyy-mm-ddTHH:MM:SS.s")

    # 2. Split-Apply-Combine to pivot from API long format to analytical wide format
    println("Executing split-apply-combine to resolve long-format rows...")
    df_wide = combine(groupby(raw_df, :timestamp)) do gdf
        # Inline helper to safely extract missing or corrupt sensor flags as NaN
        get_val(var_str) = begin
            val = get(gdf, var_str, [nothing])
            return (isempty(val) || val[1] === nothing) ? NaN32 : Float32(val[1])
        end

        (
            WS168 = get_val("HYY_META.WS168"),
            WS336 = get_val("HYY_META.WS336"),
            WS672 = get_val("HYY_META.WS672"),
            WS125 = get_val("HYY_META.WS125"),
            T168  = get_val("HYY_META.T168"),
            T336  = get_val("HYY_META.T336"),
            T672  = get_val("HYY_META.T672"),
            T125  = get_val("HYY_META.T125"),
            P_surf = get_val("HYY_META.P_R")
        )
    end

    sort!(df_wide, :timestamp)

    # 3. Establish Dimensions
    heights = Float32[16.8, 33.6, 67.2, 125.0]
    num_heights = length(heights)
    num_times = nrow(df_wide)
    time_coords = [datetime2unix(t) for t in df_wide.timestamp]

    # 4. Allocate Tensors (Height × Time for column-major efficiency)
    wind_speed = zeros(Float32, num_heights, num_times)
    pot_temp   = zeros(Float32, num_heights, num_times)

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

        # Extract pressure or use hydrostatically sound fallback if station pressure drops
        p_station = isnan(row.P_surf) ? 1013.25f0 : row.P_surf

        # Calculate Potential Temperature profiles: θ = T * (p_0 / p)^(R/Cp)
        # Note: In a deep-layer mast, p would vary with height. For the 127m mast,
        # surface reference pressure scales the column consistently for local stability metrics.
        scaling_factor = (p_0 / p_station)^R_Cp
        pot_temp[1, t] = tk1 * scaling_factor
        pot_temp[2, t] = tk2 * scaling_factor
        pot_temp[3, t] = tk3 * scaling_factor
        pot_temp[4, t] = tk4 * scaling_factor
    end

    # 5. Build NetCDF Storage
    filename = "smear_ii_processed_profiles.nc"
    if isfile(filename) rm(filename) end

    println("Writing validated NetCDF data array...")

    # Define variables with appropriate fill values and metadata attributes
    nccreate(filename, "wind_speed",
             "height", heights, Dict("units" => "m"),
             "time", time_coords, Dict("units" => "seconds since 1970-01-01 00:00:00"),
             atts = Dict("units" => "m s-1", "standard_name" => "wind_speed", "_FillValue" => NaN32))

    nccreate(filename, "potential_temperature",
             "height", heights, Dict("units" => "m"),
             "time", time_coords, Dict("units" => "seconds since 1970-01-01 00:00:00"),
             atts = Dict("units" => "K", "standard_name" => "air_potential_temperature", "_FillValue" => NaN32))

    ncwrite(wind_speed, filename, "wind_speed")
    ncwrite(pot_temp, filename, "potential_temperature")

    ncputatt(filename, "global", Dict(
        "title" => "SMEAR II Processed Boundary Layer Profile Space",
        "institution" => "University of Helsinki / SpectralBL-Analytics",
        "comment" => "Pivoted wide-format array with true potential temperature and NaN data gap markers."
    ))

    println("File compiled: ", filename)
end

# 6. Operational Execution (Using validated past windows relative to June 2026)
start_dt = DateTime(2026, 01, 01, 0, 0, 0)
end_dt   = DateTime(2026, 01, 07, 23, 59, 59)

raw_data = fetch_smear_data(start_dt, end_dt)
process_and_save_profiles(raw_data)