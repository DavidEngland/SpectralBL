# ==============================================================================
# CASES-99 Unified Pipeline & Manuscript Orchestration Engine
# ==============================================================================

# Core Paths - FIXED: Point to data/ for ingestion, and reports/ for diagnostics
ROOT_DIR     := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
DATA_DIR     := $(ROOT_DIR)/data
REPORT_DIR   := $(ROOT_DIR)/reports/ncar_eol_dee0099881
SCHEMA_DEF   := $(DATA_DIR)/cases99_netcdf_schema.txt

# Manuscript Paths
DRAFT_DIR    := $(ROOT_DIR)/drafts
FIG_DIR      := $(DRAFT_DIR)/figures
GENERATED_DIR := $(DRAFT_DIR)/sections/generated
DOC          := main

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
ALL_NC_FILES := $(wildcard $(DATA_DIR)/ncar_eol_dee0099881/cases.9910*.nc)
CAMPAIGN_DAYS := $(sort $(patsubst $(DATA_DIR)/ncar_eol_dee0099881/cases.%.nc,%,$(ALL_NC_FILES)))

# Declare all symbolic, execution-only macro endpoints safely
.PHONY: all full validate run report wave_test universal_wave_test quicktest clean setup test run-all-parallel ms clear_ms_artifacts $(CAMPAIGN_DAYS)

# Default target orchestrates the entire localized lifecycle with forced sequencing
all: setup
	$(MAKE) validate
	$(MAKE) run
	$(MAKE) report

# Aggregate target including synthetic wave sponge validation
full: all
	$(MAKE) wave_test

# ==============================================================================
# 1. Parallel Multi-Core Scaling Data Pipeline
# ==============================================================================
run-all-parallel: setup
	@echo "🚀 Initiating parallel execution matrix over $(words $(CAMPAIGN_DAYS)) campaign days..."
	$(MAKE) $(CAMPAIGN_DAYS) -j$(shell nproc 2>/dev/null || echo 4)
	@echo "✓ All available campaign days processed in parallel clusters."

$(CAMPAIGN_DAYS):
	@echo "--- [BATCH RUN] Launching Pipeline for Campaign Day: 19$@ ---"
	JULIA_LOAD_PATH="$(ROOT_DIR)/src:@:@v#.#" julia --project="$(ROOT_DIR)" $(ROOT_DIR)/scripts/RunCampaignPipeline.jl $(DATA_DIR)/ncar_eol_dee0099881/cases.$@.nc
	@echo "--- [BATCH REPORT] Generating Diagnostics for Campaign Day: 19$@ ---"
	julia --project="$(ROOT_DIR)" $(ROOT_DIR)/scripts/Report.jl $@

# ==============================================================================
# 2. Production Execution & Validation Infrastructure
# ==============================================================================
setup:
	@mkdir -p $(REPORT_DIR)
	@mkdir -p $(DATA_DIR)
	@mkdir -p $(FIG_DIR)
	@mkdir -p $(GENERATED_DIR)

validate:
	@echo "🔍 Validating schema for target file: $(notdir $(INPUT_NC))"
	julia --project="$(ROOT_DIR)" $(ROOT_DIR)/scripts/validate_netcdf_schema.jl $(INPUT_NC) $(SCHEMA_DEF)

run:
	@echo "🚀 Executing manifold transformations on: $(notdir $(INPUT_NC))"
	JULIA_LOAD_PATH="$(ROOT_DIR)/src:@:@v#.#" julia --project="$(ROOT_DIR)" $(ROOT_DIR)/scripts/RunCampaignPipeline.jl $(INPUT_NC)

report:
	@echo "📊 Compiling diagnostic summaries for window identifier: $(DAY_SUFFIX)"
	julia --project="$(ROOT_DIR)" $(ROOT_DIR)/scripts/Report.jl $(DAY_SUFFIX)
	@echo "🧾 Exporting generated TeX snippet overlays from synoptic analysis..."
	julia --project="$(ROOT_DIR)" $(ROOT_DIR)/scripts/run_synoptic_analysis.jl
	@echo "✅ Generated TeX snippets refreshed under $(GENERATED_DIR)"

wave_test:
	@echo "🧪 Running numerical sponge layer reflection verification tests..."
	julia --project="$(ROOT_DIR)" $(ROOT_DIR)/scripts/test_wave_reflection.jl

universal_wave_test:
	julia --project="$(ROOT_DIR)" $(ROOT_DIR)/scripts/run_universal_sponge_test.jl

quicktest:
	julia --project="$(ROOT_DIR)" -e 'include("src/Cases99.jl"); println("✓ Sub-meso Manifold Stack Validated Safely.")'

test:
	julia --project="$(ROOT_DIR)" -e 'using Pkg; Pkg.instantiate()'
	julia --project="$(ROOT_DIR)" $(ROOT_DIR)/test/runtests.jl

# ==============================================================================
# 3. Manuscript Compilation & Assembly Layer (AMS ametsoc2015 Sequence)
# ==============================================================================
ms: setup
	@echo "📝 Compiling text manuscript: $(DRAFT_DIR)/$(DOC).tex"
	cd $(DRAFT_DIR) && pdflatex -synctex=1 -interaction=nonstopmode -file-line-error $(DOC).tex
	@echo "📚 Parsing ametsoc2014.bst references database..."
	cd $(DRAFT_DIR) && bibtex $(DOC)
	@echo "🔄 Updating layout cross-references (Pass 2)..."
	cd $(DRAFT_DIR) && pdflatex -interaction=nonstopmode $(DOC).tex
	@echo "🔏 Finalizing document state tracking (Pass 3)..."
	cd $(DRAFT_DIR) && pdflatex -interaction=nonstopmode $(DOC).tex
	@echo "💾 Mirroring compiled layout output to target asset folder..."
	cp $(DRAFT_DIR)/$(DOC).pdf $(FIG_DIR)/
	@echo "✓ Manuscript successfully compiled and preserved in $(FIG_DIR)/$(DOC).pdf"

clear_ms_artifacts:
	@echo "🧹 Purging transient auxiliary files inside $(DRAFT_DIR)/..."
	cd $(DRAFT_DIR) && rm -f *.aux *.log *.bbl *.blg *.out *.toc *.synctex.gz *.run.xml *.bcf

# ==============================================================================
# 4. Idempotency Maintenance
# ==============================================================================
clean: clear_ms_artifacts
	@echo "🧹 Sweeping transient diagnostic logs, metrics, and report figures..."
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
	rm -f $(FIG_DIR)/manifold_geometry_plots*.png
	rm -f $(FIG_DIR)/manifold_geometry_plots*.pdf
	rm -f $(DRAFT_DIR)/$(DOC).pdf
	rm -f $(FIG_DIR)/$(DOC).pdf
	rm -f $(GENERATED_DIR)/*.tex
	@echo "✓ All transient data structures and manuscript outputs cleared. Raw field payloads preserved."