using NCDatasets, Logging

function validate_and_dump_schema(nc_path::String, output_txt_path::String)
    @info "Opening NetCDF for validation check: $nc_path"
    if !isfile(nc_path)
        error("Target NetCDF path invalid or missing.")
    end
    open(output_txt_path, "w") do io
        Dataset(nc_path, "r") do ds
            println(io, "==================================================================")
            println(io, "CASES-99 NETCDF FILE METADATA SNAPSHOT FOR PIPELINE CONFIGURATION")
            println(io, "==================================================================")
            println(io, "Variables present:\n  ", join(keys(ds), ", "), "\n")
            for var_name in ["height", "theta", "u", "v"]
                if haskey(ds, var_name)
                    println(io, "Variable '$var_name' Info -> Size: ", size(ds[var_name]))
                end
            end
        end
    end
    @info "Schema configuration written to: $output_txt_path"
end

if length(ARGS) >= 2
    validate_and_dump_schema(ARGS[1], ARGS[2])
else
    println("Usage: julia validate_netcdf_schema.jl <input_nc> <output_txt>")
end
