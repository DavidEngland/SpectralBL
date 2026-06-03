using NCDatasets
using Statistics

function verify_cases99_netcdf(filepath::String)
    println("="^60)
    println("CASES-99 NETCDF MANIFEST & SANITY CHECK")
    println("="^60)

    # Open dataset in read-only mode
    Dataset(filepath, "r") do ds
        # 1. Dimensions Check
        println("\n--- [1] DIMENSIONS ---")
        for (dimname, dimlen) in ds.dim
            println("  Dimension: $dimname => Size: $dimlen")
        end

        # 2. Variable & Height Level Verification
        println("\n--- [2] HEIGHT LEVEL CROSS-CHECK ---")
        # Extract sonic variable names to isolate your 8 targeted levels
        all_vars = keys(ds)
        u_vars = filter(v -> occursin(r"^u_\d+m(\d+cm)?", v), all_vars)

        println("Found $(length(u_vars)) horizontal wind 'u' variable levels:")
        for uv in sort(u_vars)
            # Fetch height attribute if available, or parse from variable name
            if haskey(ds[uv].attrib, "height")
                println("  Variable: $uv | Height Attrib: $(ds[uv].attrib["height"])")
            else
                println("  Variable: $uv (Verify height manually from token)")
            end
        end

        # 3. Missing Value & Out-of-Bounds Summary (October 22-31 IOP)
        println("\n--- [3] DATA INTEGRITY & MISSING VALUE STATISTICS ---")
        # Let's inspect a typical target level, e.g., the contested top sonic level
        # (Change to match your variable schema, e.g., "u_55m" or "u_60m")
        target_var = haskey(ds, "u_55m") ? "u_55m" : (haskey(ds, "u_60m") ? "u_60m" : nothing)

        if target_var !== nothing
            u_data = ds[target_var][:]

            # NCAR missing value conventions can use -1037.0
            ncar_missing_val = -1037.0
            missing_count = count(x -> x == ncar_missing_val || isnan(x), u_data)
            total_elements = length(u_data)
            missing_pct = (missing_count / total_elements) * 100

            # Filter clean data to find physical min/max overrides
            clean_data = filter(x -> x != ncar_missing_val && !isnan(x), u_data)

            println("Variable Selected for Review: $target_var")
            println("  Total Datapoints: $total_elements")
            println("  Missing/Flagged Datapoints: $missing_count ($(round(missing_pct, digits=2))%)")
            if !isEmpty(clean_data)
                println("  Physical Min Range: $(minimum(clean_data)) m/s")
                println("  Physical Max Range: $(maximum(clean_data)) m/s")
                println("  Mean Parameter Flow: $(mean(clean_data)) m/s")
            end
        else
            println("  ⚠️ Target top level variable (u_55m / u_60m) not detected directly. Check schema naming.")
        end
    end
    println("="^60)
end

# To execute, replace with your local file path:
# verify_cases99_netcdf("data/cases99_tower_hr.nc")