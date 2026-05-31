# CASES99-SpectralBL

A high-order Spectral Finite Element ($p$-FEM) toolkit written in Julia, built specifically to ingest sparse, irregular tower profiles and map them onto a metric-consistent geometry.

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
