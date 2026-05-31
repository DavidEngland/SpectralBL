# ==============================================================================
# CASES-99 Pipeline Orchestration Engine
# ==============================================================================

.PHONY: all validate run report quicktest clean

# Default target runs validation, the core pipeline, and generates the report
all: validate run report

# 1. Validation Layer: Verifies NCAR ISFS file structures against target definitions
validate:
	julia --project="." scripts/validate_netcdf_schema.jl data/sample.nc data/cases99_netcdf_schema.txt

# 2. Production Execution Layer: Sweeps and processes the primary campaign targets
run:
	julia --project="." scripts/RunCampaignPipeline.jl

# 3. Diagnostics & Business Intelligence: Generates analytical summaries and plots
report:
	julia --project="." -e 'using Plots; include("scripts/GenerateCampaignReport.jl")'

# 4. CI/CD Sandbox: Rapid regression test against a local sample payload
quicktest:
	julia --project="." -e 'include("scripts/RunCampaignPipeline.jl"); execute_pipeline("data/sample.nc", "data/test_trajectory.csv")'

# 5. Idempotency Maintenance: Wipes transient artifacts and resets workspace paths
clean:
	rm -rf reports/
	rm -f data/diagnostic_trajectory.csv data/test_trajectory.csv data/cases99_netcdf_schema.txt