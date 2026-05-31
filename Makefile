# Makefile for CASES99-SpectralBL repo automation

.PHONY: env validate run quicktest clean

env:
	julia --project="." -e 'using Pkg; Pkg.instantiate(); Pkg.status()'

validate:
	julia --project="." scripts/validate_netcdf_schema.jl data/sample.nc data/cases99_netcdf_schema.txt

run:
	julia --project="." scripts/RunCampaignPipeline.jl data/sample.nc data/diagnostic_trajectory.csv

quicktest:
	julia --project="." -e 'include("scripts/RunCampaignPipeline.jl"); execute_pipeline("data/sample.nc", "data/test_trajectory.csv", nt_limit=100)'

clean:
	rm -f data/diagnostic_trajectory.csv data/test_trajectory.csv data/cases99_netcdf_schema.txt