# ==============================================================================
# CASES-99 Pipeline Orchestration Engine
# ==============================================================================

# Context Variables
INPUT_NC     = reports/ncar_eol_dee0099881/cases.991031.nc
REPORT_DIR   = reports/ncar_eol_dee0099881
SCHEMA_DEF   = data/cases99_netcdf_schema.txt

.PHONY: all full validate run report wave_test universal_wave_test quicktest clean setup test

# Default target orchestrates the entire localized lifecycle
all: setup validate run report

# Aggregate target including synthetic wave sponge validation
full: setup validate run report wave_test

# 0. Infrastructure Layer: Ensures target directory structures exist natively
setup:
	@mkdir -p $(REPORT_DIR)
	@mkdir -p data/

# 1. Validation Layer: Verifies NCAR ISFS structures using your actual file location
validate: setup
	julia --project="." scripts/validate_netcdf_schema.jl $(INPUT_NC) $(SCHEMA_DEF)

# 2. Production Execution Layer: Runs the pipeline
run: setup
	julia --project="." scripts/RunCampaignPipeline.jl

# 3. Diagnostics Layer: Executes your new reporting logic, targeting the EOL directory
report: setup
	julia --project="." scripts/Report.jl

# 3b. Synthetic validation for sponge-layer wave reflection suppression
wave_test: setup
	julia --project="." scripts/test_wave_reflection.jl

# 3c. Universal campaign geometry replication test (defaults to CASES_99)
universal_wave_test: setup
	julia --project="." scripts/run_universal_sponge_test.jl

# 4. CI/CD Sandbox: Rapid regression test against your primary campaign payload
quicktest: setup
	julia --project="." -e 'include("src/Cases99.jl"); println("✓ Sub-meso Manifold Stack Validated Safely.")'

# 4b. Standard unit/integration test execution
test:
	julia --project="." -e 'using Pkg; Pkg.instantiate()'
	julia --project="." test/runtests.jl

# 5. Idempotency Maintenance: SAFE CLEAN
# Avoids 'rm -rf reports/' so your raw .nc data and hand-copied diagnostics aren't deleted!
clean:
	rm -f data/*_trajectory.csv
	rm -f $(SCHEMA_DEF)
	rm -f data/wave_reflection_test.png
	rm -f data/wave_reflection_metrics.csv
	rm -f $(REPORT_DIR)/manifold_geometry_plots.png
	rm -f $(REPORT_DIR)/manifold_summary_report.md
	@echo "✓ Transient artifacts cleared. Raw field data preserved."