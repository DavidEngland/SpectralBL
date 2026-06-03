# ==============================================================================
# CASES-99 Pipeline Orchestration Engine (Dynamic Data Ingestion)
# ==============================================================================

# Core Paths - FIXED: Point to your actual data directory
REPORT_DIR   = data/ncar_eol_dee0099881
SCHEMA_DEF   = data/cases99_netcdf_schema.txt

# --- DYNAMIC TARGET LOGIC ---
ifdef DAY
    INPUT_NC := $(REPORT_DIR)/cases.$(DAY).nc
else
    INPUT_NC := $(REPORT_DIR)/cases.991031.nc
endif

# Find ALL netcdf day files available in your folder for bulk runs
ALL_NC_FILES := $(wildcard $(REPORT_DIR)/cases.9910*.nc)


# 2. Dynamic metaprogramming rule mapping each day to an isolated runtime pipeline
$(CAMPAIGN_DAYS): setup
	@echo "--- [BATCH RUN] Launching Pipeline for Campaign Day: 19$$@ ---"
	JULIA_LOAD_PATH="src:@:@v#.#" julia --project="." scripts/RunCampaignPipeline.jl $(REPORT_DIR)/cases.$@.nc
	@echo "--- [BATCH REPORT] Generating Diagnostics for Campaign Day: 19$$@ ---"
	julia --project="." scripts/Report.jl $@

.PHONY: all full validate run report wave_test universal_wave_test quicktest clean setup test run-all-parallel $(CAMPAIGN_DAYS)

# 1. Main orchestration hook for multi-core scaling
# Run with 'make run-all-parallel -j4' to use 4 CPU cores simultaneously!
run-all-parallel: $(CAMPAIGN_DAYS)
	@echo "✓ All available campaign days processed in parallel clusters."

# 2. Dynamic metaprogramming rule mapping each day to an isolated runtime pipeline
$(CAMPAIGN_DAYS): setup
	@echo "--- [BATCH RUN] Launching Pipeline for Campaign Day: 19$$@ ---"
	julia --project="." scripts/RunCampaignPipeline.jl $(REPORT_DIR)/cases.$@.nc

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

# 2. Production Execution Layer: Runs the pipeline (FIXED: Passes INPUT_NC dynamically)
run: setup
	julia --project="." scripts/RunCampaignPipeline.jl $(INPUT_NC)

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
clean:
	rm -f data/*_trajectory.csv
	rm -f $(SCHEMA_DEF)
	rm -f data/wave_reflection_test.png
	rm -f data/wave_reflection_metrics.csv
	rm -f $(REPORT_DIR)/manifold_geometry_plots*.png
	rm -f $(REPORT_DIR)/manifold_summary_report*.md
	@echo "✓ All transient artifacts cleared. Raw field data preserved."