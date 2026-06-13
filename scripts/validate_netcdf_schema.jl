using NCDatasets, Logging

function validate_and_dump_schema(nc_path::String, output_txt_path::String)
    @info "Opening NetCDF for validation check: $nc_path"
    if !isfile(nc_path)
        error("Target NetCDF path invalid or missing.")
    end
    open(output_txt_path, "w") do io
        Dataset(nc_path, "r") do ds
            println(io, "==================================================================")
            println(io, "NETCDF FILE METADATA SNAPSHOT FOR PIPELINE CONFIGURATION")
            println(io, "SOURCE FILE: ", nc_path)
            println(io, "==================================================================")
            println(io, "\nDimensions:")
            for (dim_name, dim_size) in ds.dim
                println(io, "  ", dim_name, " = ", dim_size)
            end

            println(io, "\nVariables present:\n  ", join(keys(ds), ", "), "\n")

            println(io, "Variable metadata:")
            for var_name in sort!(String.(collect(keys(ds))))
                var = ds[var_name]
                dims = join(dimnames(var), ",")
                sz = join(size(var), "x")
                units = haskey(var.attrib, "units") ? string(var.attrib["units"]) : ""
                println(io, "  ", var_name, " | dims=", dims, " | size=", sz, " | units=", units)
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
