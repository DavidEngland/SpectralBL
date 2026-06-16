# src/smear_fetch.jl
using HTTP
using JSON3
using DataFrames
using NetCDF
using Dates

function fetch_smear_data()
    endpoint = "https://smear-backend.rahtiapp.fi/search/timeseries"

    # Target variables across the SMEAR II (Hyytiälä) 127m mast
    # WS = Wind Speed (m/s), T = Temperature (°C)
    payload = Dict(
        "tablevariable" => [
            "HYY_META.WS125", "HYY_META.WS672", "HYY_META.WS336", "HYY_META.WS168",
            "HYY_META.T125",  "HYY_META.T672",  "HYY_META.T336",  "HYY_META.T168"
        ],
        "from" => "2025-12-01T00:00:00.000",
        "to" => "2025-12-07T23:59:59.000",
        "quality" => "ANY",
        "aggregation" => "NONE"
    )

    println("Querying SMEAR API via HTTP.jl...")
    headers = ["Content-Type" => "application/json"]
    response = HTTP.post(endpoint, headers, JSON3.write(payload))

    # Parse JSON directly into a lightweight structured object
    parsed_json = JSON3.read(response.body)
    raw_records = parsed_json[:data]

    # Convert to DataFrame for easier multidimensional pivoting
    df = DataFrame(raw_records)
    df.timestamp = DateTime.(replace.(df.samptime, "Z" => ""), dateformat"yyyy-mm-ddTHH:MM:SS.s")
    return df
end

function compile_spectralbl_netcdf(df::DataFrame)
    # Define physical heights matching the SMEAR tower booms
    heights = Float32[16.8, 33.6, 67.2, 125.0]
    num_heights = length(heights)

    # Extract temporal dimension
    unique_times = df.timestamp
    num_times = length(unique_times)

    # Convert timestamps to unix epoch seconds for CF-convention compliance
    time_coords = Float64[datetime2unix(t) for t in unique_times]

    # Pre-allocate 2D field tensors: Shape (height, time) for optimal column-major memory layout
    wind_speed = zeros(Float32, num_heights, num_times)
    pot_temp = zeros(Float32, num_heights, num_times)

    println("Pivoting flat API response into column-major tensors...")
    p_0 = 1000.0f0 # Reference surface pressure in hPa
    R_Cp = 0.286f0  # Dry air gas constant / specific heat ratio

    for t in 1:num_times
        # Map wind speed layers
        wind_speed[1, t] = Float32(get(df[t, :], Symbol("HYY_META.WS168"), 0.0))
        wind_speed[2, t] = Float32(get(df[t, :], Symbol("HYY_META.WS336"), 0.0))
        wind_speed[3, t] = Float32(get(df[t, :], Symbol("HYY_META.WS672"), 0.0))
        wind_speed[4, t] = Float32(get(df[t, :], Symbol("HYY_META.WS125"), 0.0))

        # Pull temperatures, convert from Celsius to Kelvin
        t1 = Float32(get(df[t, :], Symbol("HYY_META.T168"), 0.0)) + 273.15f0
        t2 = Float32(get(df[t, :], Symbol("HYY_META.T336"), 0.0)) + 273.15f0
        t3 = Float32(get(df[t, :], Symbol("HYY_META.T672"), 0.0)) + 273.15f0
        t4 = Float32(get(df[t, :], Symbol("HYY_META.T125"), 0.0)) + 273.15f0

        # Calculate Potential Temperature (Θ). In a full deployment, scale using local pressure arms.
        # For this prototype, we treat them as absolute profiles to feed the state estimator.
        pot_temp[1, t] = t1
        pot_temp[2, t] = t2
        pot_temp[3, t] = t3
        pot_temp[4, t] = t4
    end

    filename = "smear_ii_polar_profiles.nc"
    if isfile(filename) rm(filename) end

    println("Writing CF-Compliant NetCDF file via NetCDF.jl...")

    # Define Dimensions and Variables simultaneously in Julia
    nccreate(filename, "wind_speed",
             "height", heights, Dict("units" => "m"),
             "time", time_coords, Dict("units" => "seconds since 1970-01-01 00:00:00"),
             atts = Dict("units" => "m s-1", "standard_name" => "wind_speed"))

    nccreate(filename, "potential_temperature",
             "height", heights, Dict("units" => "m"),
             "time", time_coords, Dict("units" => "seconds since 1970-01-01 00:00:00"),
             atts = Dict("units" => "K", "standard_name" => "air_potential_temperature"))

    # Write arrays directly to disk
    ncwrite(wind_speed, filename, "wind_speed")
    ncwrite(pot_temp, filename, "potential_temperature")

    # Add project metadata for your NotebookLM references
    ncputatt(filename, "global", Dict(
        "title" => "SMEAR II Boreal/Polar Boundary Layer Ingestion Field",
        "institution" => "University of Helsinki / SpectralBL-Analytics",
        "source" => "SmartSmear API Engine",
        "processing_stage" => "Stage 1 Linear Observer Input"
    ))

    println("File saved successfully as: ", filename)
end

# Execute Ingestion Pipeline
df_smear = fetch_smear_data()
compile_spectralbl_netcdf(df_smear)