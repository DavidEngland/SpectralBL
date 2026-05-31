.PHONY: validate run quicktest clean

validate:
	julia --project="." scripts/validate_netcdf_schema.jl data/sample.nc data/cases99_netcdf_schema.txt

run:
	julia --project="." scripts/RunCampaignPipeline.jl data/sample.nc data/diagnostic_trajectory.csv

quicktest:
	julia --project="." -e 'include("scripts/RunCampaignPipeline.jl"); execute_pipeline("data/sample.nc", "data/test_trajectory.csv")'

clean:
	rm -f data/diagnostic_trajectory.csv data/test_trajectory.csv data/cases99_netcdf_schema.txt