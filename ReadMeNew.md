# CASES99-SpectralBL

A high-order Spectral Finite Element ($p$-FEM) toolkit written in Julia, built specifically to ingest sparse, irregular tower profiles and map them onto a metric-consistent geometry.

## Recent Updates & Changelog

### [v1.2.0] - Stable Boundary Layer (SBL) Geometry Expansion
- **NetCDF Structural Data Curation:** Enhanced `scripts/validate_netcdf_schema.jl` to natively parse, validate, and standardize binary self-describing NetCDF files from sparse tower streams before manifold ingestion.
- **Tri-Modal GMM Segmentation Matrix:** Deployed a GMM segmentation schema over the multi-regime information-theoretic array ($D_{\mathrm{eff}}$, $F_W$, $\chi_N$).
- **Silhouette Validation Optimization:** Integrated rigorous silhouette width evaluations confirming a definitive global peak at $K=3$ ($\bar{S} = 0.582$). This formally establishes the *Intermittent Shear Bursting State* as a statistically independent, long-lived physical regime rather than a transient tracking artifact.
- **Virtual Tower Interfacing Operator ($\mathcal{P}_{\text{VT}}$):** Formulated a generalized, regularized least-squares mapping scheme to bridge the grid-topology mismatch by projecting coarse numerical weather prediction (NWP) levels ($\eta$-coordinates) directly onto the continuous continuous computational manifold.
- **Spectral Anti-Aliasing Array:** Implemented a metric-consistent $2/3$-rule exponential spectral filter ($n_{\mathrm{crit}} = 21$ for $N=32$) to control non-linear energy accumulation in forward-predictive implementations.

---

## Quick Start Configuration

1. Initialize dependencies inside Julia's active runtime package space:

```bash
julia --project="." -e 'using Pkg; Pkg.instantiate()'

```

2. Run pre-flight validation on your custom NetCDF header tracking scheme:

```bash
julia --project="." scripts/validate_netcdf_schema.jl data/sample.nc data/cases99_netcdf_schema.txt

```

3. Execute the full feature-extraction production matrix loop across all data steps:

```bash
julia --project="." scripts/RunCampaignPipeline.jl data/sample.nc data/diagnostic_trajectory.csv

```

---

## Synthetic Wave Sponge Validation

Run the standalone synthetic gravity-wave reflection experiment:

```bash
julia --project="." scripts/test_wave_reflection.jl

```

Or use Make targets:

```bash
make wave_test   # synthetic wave-only validation
make full        # setup + validate + run + report + wave_test
make universal_wave_test   # campaign-geometry replication run (default CASES_99)

```

Artifacts are written to:

* `data/wave_reflection_test.png`
* `data/wave_reflection_metrics.csv`

---

## Universal Campaign Replication Framework

Run with predefined campaign geometry:

```bash
julia --project="." scripts/run_universal_sponge_test.jl CASES_99
julia --project="." scripts/run_universal_sponge_test.jl FLOSS_II
julia --project="." scripts/run_universal_sponge_test.jl SHEBA
julia --project="." scripts/run_universal_sponge_test.jl BLLAST

```

Run all predefined campaigns:

```bash
julia --project="." scripts/run_universal_sponge_test.jl --all

```

Run custom geometry (Arguments: `NAME z_min z_max alpha_stretch`):

```bash
julia --project="." scripts/run_universal_sponge_test.jl MY_CAMPAIGN 1.5 50.0 0.05

```

Per-campaign artifacts are written under `data/universal_sponge/<campaign_name>/`:

* `wave_reflection_metrics.csv`
* `wave_reflection_test.png`
* `wave_reflection_summary.md`
