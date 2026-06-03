As currently written, **no, it is not locked exclusively to CASES-99.** The pipeline architecture is intentionally engineered to be **campaign-agnostic**, meaning you can reuse this exact framework for other field datasets (such as MATERHORN, SGP, or Arctic tower networks) with minimal configuration changes.

However, because you are pulling data from a different project, there are three explicit places where campaign-specific variables are defined. To pivot this code to a new dataset, you only need to modify these three configuration blocks.

---

### The 3 Points of Adaptation for Other Datasets

#### 1. The File Naming Regex (`scripts/RunCampaignPipeline.jl`)

Currently, the pipeline looks for files matching the CASES-99 naming convention:

```julia
# Looks for "cases.991024.nc"
nc_files = sort(filter(f -> match(r"^cases\.\d+\.nc$", f) !== nothing, readdir(data_dir)))

```

If you switch to a dataset where the files are named `materhorn_2012_v1.nc`, you simply update this regular expression string to match your new file names.

#### 2. The Instrument Metadata Arrays (`scripts/RunCampaignPipeline.jl`)

This is the most critical change. Different campaigns mount sensors at different heights and name their NetCDF variables according to their own conventions. You just update these three array maps at the top of your execution function:

```julia
# Example: Shifting from CASES-99 to a hypothetical 4-level Arctic Tower
heights     = [2.0, 10.0, 20.0, 40.0]
target_vars = ["wnd_spd_2m", "wnd_spd_10m", "wnd_spd_20m", "wnd_spd_40m"]
tc_vars     = ["temp_2m", "temp_10m", "temp_20m", "temp_40m"]

```

#### 3. Hardware-Specific Quality Flags (`src/Summary.jl`)

In the `clean_by_ncar_quality_flags` function inside your module, the data-cleaning bitmasks are tied to the specific hardware used in 1999 (Campbell CSAT3 diagnostics and ATI sampling limits):

```julia
if height_string in ["1_5m", "5m", "30m", "50m"]
    diag_var = "diag_" * height_string
    # ... bitwise cleaning ...

```

If your new dataset does not use these exact variables for quality assurance (or uses standard native NetCDF attributes like `_FillValue` and `valid_range`), this function will safely fall through to a pass-through without crashing. If the new campaign uses a unique quality flag variable (e.g., a simple QC flag matrix where `0 = good`, `1 = bad`), you can use Julia's multiple dispatch or add a new conditional branch to parse those specific flags.

---

### Summary of What Makes It Modular

Because your core mathematical tools inside `src/Summary.jl` (the Gradient Test, the Matrix Conditioner, and the Spike Filter) ingest raw, standard `Vector{Float64}` arrays, **they do not care what campaign the data came from.** As long as your execution script opens the target file format, extracts the variables into a 1D vector, and passes them along, your pseudospectral pipeline will analyze the underlying fluid mechanics flawlessly, regardless of the campaign.