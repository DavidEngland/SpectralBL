# scripts/GenerateCampaignReport.jl
using NCDatasets, LinearAlgebra, Plots, Dates

# Ensure output directories exist
mkpath("reports")
mkpath("reports/plots")

function generate_report()
    data_dir = "data/ncar_eol_dee0099881"
    # Filter and sort daily campaign files
    nc_files = sort(filter(f -> match(r"^cases\.\d+\.nc$", f) !== nothing, readdir(data_dir)))

    println("Found $(length(nc_files)) campaign files to process.")

    # Initialize the master Markdown report file
    open("reports/CASES99_Campaign_Report.md", "w") do report
        write(report, "# CASES-99 Field Campaign Diagnostic Report\n\n")
        write(report, "Generated on: $(Dates.today())  \n")
        write(report, "Data Source: NCAR EOL Asset `ncar_eol_dee0099881`\n\n")
        write(report, "---\n\n")

        write(report, "## Executive Summary\n")
        write(report, "Automated processing sweep of $(length(nc_files)) sequential daily boundary-layer profiles from the 1999 intensive operational period.\n\n")

        # Main campaign iteration sweep
        for (idx, nc_file) in enumerate(nc_files)
            full_path = joinpath(data_dir, nc_file)
            # Pull date stamp from filename (e.g., cases.990919.nc -> Sep 19, 1999)
            date_str = "19" * nc_file[7:12]
            parsed_date = Dates.DateTime(date_str, "yyyymmdd")
            formatted_date = Dates.format(parsed_date, "u d, yyyy")

            println("Processing: $nc_file ($formatted_date)...")

            # Extract profile variance statistics for the day
            mean_wind, mean_temp, times = process_daily_profiles(full_path)

            if mean_wind === nothing
                continue
            end

            # Plot 1: Daily Boundary Layer Evolution Contour
            plot_name = "daily_evolution_day_$(idx).png"
            plot_path = "reports/plots/\$plot_name"

            # Simple visualization using Plots.jl
            p = heatmap(times, [1.5, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0], mean_wind,
                        title="Wind Profile Evolution - \$formatted_date",
                        xlabel="Time Index / Dimension", ylabel="Height (z) [m]",
                        color=:viridis, clabel="Velocity [m/s]")
            savefig(p, plot_path)

            # Append findings to the Markdown report
            write(report, "### Day $idx: $formatted_date (`$nc_file`)\n\n")
            write(report, "| Metric | Value |\n")
            write(report, "| :--- | :--- |\n")
            write(report, "| Profile Snapshots | \$(size(mean_wind, 2)) |\n")
            write(report, "| Max Recorded Wind | \$(round(maximum(mean_wind), sigdigits=4)) m/s |\n")
            write(report, "| Core Temp Range | \$(round(minimum(mean_temp), sigdigits=4)) to \$(round(maximum(mean_temp), sigdigits=4)) °C/K |\n\n")

            write(report, "#### Boundary Layer Structure\n")
            write(report, "![Wind Profile Evolution for \$formatted_date](./plots/\$plot_name)\n\n")
            write(report, "---\n\n")
        end
    end
    println("Campaign report generated successfully at: `reports/CASES99_Campaign_Report.md`")
end

function process_daily_profiles(nc_path)
    Dataset(nc_path, "r") do ds
        z_obs = [1.5, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0]
        M_obs = length(z_obs)

        # Safely extract temporal dimensions
        t_len = haskey(ds, "time") ? length(ds["time"]) : 100 # Fallback sizing

        w_matrix = zeros(Float64, M_obs, t_len)
        t_matrix = zeros(Float64, M_obs, t_len)

        try
            for i in 1:M_obs
                h = z_obs[i]
                h_str = h == 1.5 ? "1_5m" : string(Int(h)) * "m"

                # Check for schema variant fallbacks
                u_key = haskey(ds, "u_\$h_str") ? "u_\$h_str" : (haskey(ds, "U_\$h_str") ? "U_\$h_str" : return nothing, nothing, nothing)
                t_key = haskey(ds, "tc_\$h_str") ? "tc_\$h_str" : (haskey(ds, "T_\$h_str") ? "T_\$h_str" : return nothing, nothing, nothing)

                # Fill profiles across the timeline
                w_matrix[i, :] .= [ismissing(x) ? 0.0 : Float64(x) for x in ds[u_key][:]]
                t_matrix[i, :] .= [ismissing(x) ? 273.15 : Float64(x) for x in ds[t_key][:]]
            end
            return w_matrix, t_matrix, 1:t_len
        catch e
            return nothing, nothing, nothing
        end
    end
end

generate_report()