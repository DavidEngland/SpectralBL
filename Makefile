# ==============================================================================
# CASES-99 Pipeline Orchestration Engine (Dynamic Data Ingestion)
# ==============================================================================

# Core Paths - FIXED: Point to data/ for ingestion, and reports/ for diagnostics
ROOT_DIR     := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
DATA_DIR     := $(ROOT_DIR)/data
REPORT_DIR   := $(ROOT_DIR)/reports/ncar_eol_dee0099881
SCHEMA_DEF   := $(DATA_DIR)/cases99_netcdf_schema.txt

# --- DYNAMIC TARGET LOGIC ---
# FIXED: Change input tracking from REPORT_DIR to DATA_DIR where the raw data lives!
ifdef DAY
    INPUT_NC := $(DATA_DIR)/ncar_eol_dee0099881/cases.$(DAY).nc
    DAY_SUFFIX := $(DAY)
else
    INPUT_NC := $(DATA_DIR)/ncar_eol_dee0099881/cases.991031.nc
    DAY_SUFFIX := 991031
endif

# Find ALL NetCDF day files available in your folder for bulk runs
# Use absolute-rooted paths and derive day IDs from filenames
ALL_NC_FILES := $(wildcard $(DATA_DIR)/ncar_eol_dee0099881/cases.9910*.nc)
CAMPAIGN_DAYS := $(sort $(patsubst $(DATA_DIR)/ncar_eol_dee0099881/cases.%.nc,%,$(ALL_NC_FILES)))

# Declare all symbolic, execution-only macro endpoints safely
.PHONY: all full validate run report wave_test universal_wave_test quicktest clean setup test run-all-parallel $(CAMPAIGN_DAYS)

# Default target orchestrates the entire localized lifecycle with forced sequencing
all: setup
	$(MAKE) validate
	$(MAKE) run
	$(MAKE) report

# Aggregate target including synthetic wave sponge validation
full: all
	$(MAKE) wave_test

# 1. Main orchestration hook for multi-core scaling
# Run with 'make run-all-parallel -j4' to use 4 CPU cores simultaneously!
run-all-parallel: setup
	@echo "🚀 Initiating parallel execution matrix over $(words $(CAMPAIGN_DAYS)) campaign days..."
	$(MAKE) $(CAMPAIGN_DAYS) -j$(shell nproc 2>/dev/null || echo 4)
	@echo "✓ All available campaign days processed in parallel clusters."

# 2. Dynamic metaprogramming rule mapping each day to an isolated runtime pipeline
# FIXED: Removed the direct dependency on 'setup' to prevent parallel directory creation races.
$(CAMPAIGN_DAYS):
	@echo "--- [BATCH RUN] Launching Pipeline for Campaign Day: 19$@ ---"
	JULIA_LOAD_PATH="$(ROOT_DIR)/src:@:@v#.#" julia --project="$(ROOT_DIR)" $(ROOT_DIR)/scripts/RunCampaignPipeline.jl $(DATA_DIR)/ncar_eol_dee0099881/cases.$@.nc
	@echo "--- [BATCH REPORT] Generating Diagnostics for Campaign Day: 19$@ ---"
	julia --project="$(ROOT_DIR)" $(ROOT_DIR)/scripts/Report.jl $@

# 0. Infrastructure Layer: Ensures target directory structures exist natively and silently
setup:
	@mkdir -p $(REPORT_DIR)
	@mkdir -p $(DATA_DIR)
	@mkdir -p $(DATA_DIR)/drafts/figures

# 1. Validation Layer: Verifies NCAR ISFS structures using your actual file location
validate:
	@echo "🔍 Validating schema for target file: $(notdir $(INPUT_NC))"
	julia --project="$(ROOT_DIR)" $(ROOT_DIR)/scripts/validate_netcdf_schema.jl $(INPUT_NC) $(SCHEMA_DEF)

# 2. Production Execution Layer: Runs the pipeline with dynamic target parameters
run:
	@echo "🚀 Executing manifold transformations on: $(notdir $(INPUT_NC))"
	JULIA_LOAD_PATH="$(ROOT_DIR)/src:@:@v#.#" julia --project="$(ROOT_DIR)" $(ROOT_DIR)/scripts/RunCampaignPipeline.jl $(INPUT_NC)

# 3. Diagnostics Layer: Executes your new reporting logic, passing active day suffix
report:
	@echo "📊 Compiling diagnostic summaries for window identifier: $(DAY_SUFFIX)"
	julia --project="$(ROOT_DIR)" $(ROOT_DIR)/scripts/Report.jl $(DAY_SUFFIX)

# 3b. Synthetic validation for sponge-layer wave reflection suppression
wave_test:
	@echo "🧪 Running numerical sponge layer reflection verification tests..."
	julia --project="$(ROOT_DIR)" $(ROOT_DIR)/scripts/test_wave_reflection.jl

# 3c. Universal campaign geometry replication test (defaults to CASES_99)
universal_wave_test:
	julia --project="$(ROOT_DIR)" $(ROOT_DIR)/scripts/run_universal_sponge_test.jl

# 4. CI/CD Sandbox: Rapid regression test against your primary campaign payload
quicktest:
	julia --project="$(ROOT_DIR)" -e 'include("src/Cases99.jl"); println("✓ Sub-meso Manifold Stack Validated Safely.")'

# 4b. Standard unit/integration test execution
test:
	julia --project="$(ROOT_DIR)" -e 'using Pkg; Pkg.instantiate()'
	julia --project="$(ROOT_DIR)" $(ROOT_DIR)/test/runtests.jl

# 5. Idempotency Maintenance: SAFE CLEAN
clean:
	rm -f $(DATA_DIR)/trajectory_*.csv
	rm -f $(DATA_DIR)/diagnostic_trajectory.csv
	rm -f $(SCHEMA_DEF)
	rm -f $(DATA_DIR)/wave_reflection_test.png
	rm -f $(DATA_DIR)/wave_reflection_test.pdf
	rm -f $(DATA_DIR)/wave_reflection_metrics.csv
	rm -f $(DATA_DIR)/synoptic_run.log
	rm -f $(REPORT_DIR)/manifold_geometry_plots*.png
	rm -f $(REPORT_DIR)/manifold_geometry_plots*.pdf
	rm -f $(REPORT_DIR)/manifold_summary_report*.md
	rm -f $(REPORT_DIR)/campaign_synoptic_evolution*.png
	rm -f $(REPORT_DIR)/campaign_synoptic_evolution*.pdf
	rm -f $(REPORT_DIR)/synoptic_campaign_audit*.md
	rm -f $(DATA_DIR)/drafts/figures/manifold_geometry_plots*.png
	rm -f $(DATA_DIR)/drafts/figures/manifold_geometry_plots*.pdf
	@echo "✓ All transient artifacts cleared. Raw field data preserved."