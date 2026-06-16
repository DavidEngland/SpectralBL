#  <#Title#>

Here is an updated, highly structured version of your `README.md` for the **SpectralBL** repository.

This update reflects your new architectural direction: it explicitly documents the status of **GABLS3**, clarifies the workflows for both campaign types, organizes the generated artifacts, and integrates the reporting/manuscript separation you discussed.

---

# SpectralBL (Formerly CASES99-SpectralBL)

A high-order Spectral Finite Element ($p$-FEM) toolkit written in Julia, built specifically to ingest sparse, irregular tower profiles and map them onto a metric-consistent geometry for stable atmospheric boundary layer simulations.

---

## 📋 Campaign Execution Status

This framework supports multiple observational campaigns. The table below outlines the current implementation and validation status for each framework target:

| Campaign ID | Focus / Target | Status | Validation Scripts |
| --- | --- | --- | --- |
| **CASES_99** | Stable Boundary Layer (SBL) | **Production-Ready** | `RunCampaignPipeline.jl`, `run_universal_sponge_test.jl` |
| **GABLS3** | Diurnal Cycle & Low-Level Jet | ⚙️ **Active Integration** | Use `GABLS3` target flag with the universal framework |
| **FLOSS_II** | Shallow Cold Pools | **Supported** | Universal framework verification |
| **SHEBA** | Arctic Ice Boundary Layer | **Supported** | Universal framework verification |
| **BLLAST** | Boundary Layer Late Afternoon Transition | **Supported** | Universal framework verification |

---

## 🚀 Quick Start Configuration

### 1. Environment Setup

Initialize dependencies inside Julia's active runtime package space:

```bash
julia --project="." -e 'using Pkg; Pkg.instantiate()'

```

### 2. Pre-Flight Header Validation

Run validation on your custom NetCDF header tracking scheme:

```bash
julia --project="." scripts/validate_netcdf_schema.jl data/sample.nc data/cases99_netcdf_schema.txt

```

### 3. Production Pipeline Execution

Execute the full feature-extraction production matrix loop across all data steps:

```bash
julia --project="." scripts/RunCampaignPipeline.jl data/sample.nc data/diagnostic_trajectory.csv

```

> ⚠️ **Note:** After running `make clean`, remember to regenerate your schema and trajectory artifacts before initiating reporting or manuscript builds.

---

## 🛠️ Unified Campaign Workflow

Use this sequence when rebuilding a manuscript-ready draft from a clean workspace. The pipeline isolates data processing from presentation layouts.

```bash
make setup
make validate
make run
julia --project="." scripts/Report.jl
julia --project="." scripts/run_synoptic_analysis.jl
julia --project="." scripts/run_universal_sponge_test.jl CASES_99
make ms

```

### Generated Assets File Tree

Running the production workflow produces the figure and table assets consumed by your manuscript draft (`drafts/main.tex`).

```text
├── data/
│   ├── drafts/figures/             # Generated pipeline figure assets
│   └── universal_sponge/
│       ├── cases_99/               # CASES-99 specific outputs
│       └── gabls3/                 # GABLS3 specific outputs
└── drafts/
    └── figures/                    # Hard links/copies decoupled for manuscript compilation

```

---

## 🌊 Synthetic Wave Sponge Validation

Run standalone synthetic gravity-wave reflection experiments using native scripts or shorthand `make` targets.

### Direct Execution

```bash
julia --project="." scripts/test_wave_reflection.jl

```

### Automation via Make

```bash
make wave_test             # Synthetic wave-only validation
make full                  # Full Pipeline: setup + validate + run + report + wave_test
make universal_wave_test   # Campaign-geometry replication run (Defaults to CASES_99)

```

### Synthetic Artifact Pipeline Outputs

* `data/wave_reflection_test.png`
* `data/wave_reflection_test.pdf`
* `data/wave_reflection_metrics.csv`

---

## 🌍 Universal Campaign Replication Framework

Run the universal boundary layer tracking pipeline using predefined or custom experimental campaign geometries.

### Predefined Campaigns

Execute the pipeline against standard experimental geometries (e.g., `CASES_99`, `GABLS3`, `SHEBA`):

```bash
julia --project="." scripts/run_universal_sponge_test.jl CASES_99
julia --project="." scripts/run_universal_sponge_test.jl GABLS3

```

### Batch Execution

Run all predefined campaigns sequentially through the pipeline matrix:

```bash
julia --project="." scripts/run_universal_sponge_test.jl --all

```

### Custom Campaign Geometry

Define your own domain properties on the fly by passing explicit configurations:

```bash
julia --project="." scripts/run_universal_sponge_test.jl MY_CAMPAIGN [Domain_Scale] [Sponge_Depth] [Damping_Coeff]
# Example:
julia --project="." scripts/run_universal_sponge_test.jl MY_CAMPAIGN 1.5 50.0 0.05

```

### Campaign Output Architecture

Per-campaign artifacts are written under `data/universal_sponge/<campaign_name>/` in strictly **lowercase** directories (e.g., `GABLS3` outputs map directly to `data/universal_sponge/gabls3/`):

```text
data/universal_sponge/<campaign_name>/
├── wave_reflection_metrics.csv     # Raw tabular error metrics
├── wave_reflection_test.png        # Diagnostic quick-look visualization
├── wave_reflection_test.pdf        # Publication-ready vector plot
└── wave_reflection_summary.md      # Auto-generated markdown metric card

```
