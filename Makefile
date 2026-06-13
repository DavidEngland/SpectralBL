# ==============================================================================
# CASES-99 Unified Pipeline & Manuscript Orchestration Engine
# ==============================================================================

# Core Paths - FIXED: Point to data/ for ingestion, and reports/ for diagnostics
ROOT_DIR     := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
DATA_DIR     := $(ROOT_DIR)/data
REPORT_DIR   := $(ROOT_DIR)/reports/ncar_eol_dee0099881
SCHEMA_DEF   := $(DATA_DIR)/cases99_netcdf_schema.txt
GABLS3_NC    := $(DATA_DIR)/gabs3/gabls3_scm_cabauw_obs_v33.nc
GABLS3_SCHEMA_DEF := $(DATA_DIR)/gabs3/gabls3_netcdf_schema.txt

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

# Versioned manuscript export filename, e.g. SpectralBL_manuscript_CASES99_v991031
EXPORT_NAME := SpectralBL_manuscript_CASES99_v$(DAY_SUFFIX)

# Find ALL NetCDF day files available in your folder for bulk runs
ALL_NC_FILES := $(wildcard $(DATA_DIR)/ncar_eol_dee0099881/cases.9910*.nc)
CAMPAIGN_DAYS := $(sort $(patsubst $(DATA_DIR)/ncar_eol_dee0099881/cases.%.nc,%,$(ALL_NC_FILES)))

# Declare all symbolic, execution-only macro endpoints safely
.PHONY: all full validate run report transform-report attractor-figure wave_test universal_wave_test quicktest clean setup test run-all-parallel ms clear_ms_artifacts purge-generated verify-manuscript gabls3-validate gabls3-run gabls3-predict gabls3-plots gabls3 $(CAMPAIGN_DAYS)

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

gabls3-validate:
	@echo "🔍 Validating schema for GABLS3 target file: $(notdir $(GABLS3_NC))"
	julia --project="$(ROOT_DIR)" $(ROOT_DIR)/scripts/validate_netcdf_schema.jl $(GABLS3_NC) $(GABLS3_SCHEMA_DEF)

gabls3-run:
	@echo "🚀 Executing GABLS3 LES manifold pipeline on: $(notdir $(GABLS3_NC))"
	JULIA_LOAD_PATH="$(ROOT_DIR)/src:@:@v#.#" julia --project="$(ROOT_DIR)" $(ROOT_DIR)/scripts/RunGabLS3Pipeline.jl $(GABLS3_NC) $(DATA_DIR)/trajectory_gabls3.csv

gabls3-predict:
	@echo "🧠 Training weighted GABLS3 predictive baseline from trajectory features..."
	JULIA_LOAD_PATH="$(ROOT_DIR)/src:@:@v#.#" julia --project="$(ROOT_DIR)" $(ROOT_DIR)/scripts/GabLS3PredictiveBaseline.jl

gabls3-plots:
	@echo "📈 Generating GABLS3 diagnostic plots via CairoMakie..."
	JULIA_LOAD_PATH="$(ROOT_DIR)/src:@:@v#.#" julia --project="$(ROOT_DIR)" $(ROOT_DIR)/scripts/PlotGabLS3Diagnostics.jl

gabls3: setup gabls3-validate gabls3-run gabls3-predict gabls3-plots
	@echo "✓ GABLS3 ingestion, diagnostics, and predictive baseline run complete."

report:
	@echo "📊 Compiling diagnostic summaries for window identifier: $(DAY_SUFFIX)"
	julia --project="$(ROOT_DIR)" $(ROOT_DIR)/scripts/Report.jl $(DAY_SUFFIX)
	@echo "🧾 Exporting generated TeX snippet overlays from synoptic analysis..."
	julia --project="$(ROOT_DIR)" $(ROOT_DIR)/scripts/run_synoptic_analysis.jl
	@echo "🌀 Rendering physical-to-spectral attractor collapse figure..."
	julia --project="$(ROOT_DIR)" $(ROOT_DIR)/scripts/generate_coefficient_trajectory_figure.jl
	@echo "✅ Generated TeX snippets refreshed under $(GENERATED_DIR)"

transform-report:
	@echo "🧭 Exporting transform and diagnostic snippets for window identifier: $(DAY_SUFFIX)"
	julia --project="$(ROOT_DIR)" $(ROOT_DIR)/scripts/Report.jl $(DAY_SUFFIX)
	@echo "🧾 Refreshing synoptic diagnostics macros when trajectory data is available..."
	julia --project="$(ROOT_DIR)" $(ROOT_DIR)/scripts/run_synoptic_analysis.jl
	@echo "🌀 Rendering physical-to-spectral attractor collapse figure..."
	julia --project="$(ROOT_DIR)" $(ROOT_DIR)/scripts/generate_coefficient_trajectory_figure.jl

attractor-figure:
	@echo "🌀 Rendering physical-to-spectral attractor collapse figure..."
	julia --project="$(ROOT_DIR)" $(ROOT_DIR)/scripts/generate_coefficient_trajectory_figure.jl

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
	@echo "♻️ Regenerating synoptic diagnostics tables/macros for manuscript consistency..."
	julia --project="$(ROOT_DIR)" $(ROOT_DIR)/scripts/run_synoptic_analysis.jl
	@echo "♻️ Regenerating D_eff diagnostics macros for manuscript consistency..."
	julia --project="$(ROOT_DIR)" $(ROOT_DIR)/regenerate_tex_exports.jl
	@echo "📝 Compiling text manuscript: $(DRAFT_DIR)/$(DOC).tex"
	cd $(DRAFT_DIR) && pdflatex -synctex=1 -interaction=nonstopmode -file-line-error $(DOC).tex
	@echo "📚 Parsing ametsoc2014.bst references database..."
	cd $(DRAFT_DIR) && bibtex $(DOC)
	@echo "🔄 Updating layout cross-references (Pass 2)..."
	cd $(DRAFT_DIR) && pdflatex -interaction=nonstopmode $(DOC).tex
	@echo "🔏 Finalizing document state tracking (Pass 3)..."
	cd $(DRAFT_DIR) && pdflatex -interaction=nonstopmode $(DOC).tex
	@echo "💾 Exporting compiled manuscript with versioned tag..."
	cp $(DRAFT_DIR)/$(DOC).pdf $(DRAFT_DIR)/$(EXPORT_NAME).pdf
	cp $(DRAFT_DIR)/$(DOC).pdf $(FIG_DIR)/
	@echo "✓ Manuscript successfully compiled and preserved as $(DRAFT_DIR)/$(EXPORT_NAME).pdf"

purge-generated: setup
	@echo "🧹 Purging generated manuscript fragments in $(GENERATED_DIR)..."
	rm -f $(GENERATED_DIR)/*.tex
	@echo "✓ Generated fragment directory reset."

verify-manuscript: purge-generated report ms
	@echo "🔍 Verification complete: generated fragments refreshed and manuscript compiled."

clear_ms_artifacts:
	@echo "🧹 Purging transient auxiliary files inside $(DRAFT_DIR)/..."
	cd $(DRAFT_DIR) && rm -f *.aux *.log *.bbl *.blg *.out *.toc *.synctex.gz *.run.xml *.bcf

# ==============================================================================
# 4. Idempotency Maintenance
# ==============================================================================
clean: clear_ms_artifacts
	@echo "🧹 Sweeping transient diagnostic logs, metrics, and report figures..."
	@echo "ℹ Preserving trajectory shards and merged trajectory cache to keep manuscript diagnostics reproducible."
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
	rm -f $(DRAFT_DIR)/SpectralBL_manuscript_*.pdf
	rm -f $(GENERATED_DIR)/*.tex
	@echo "✓ All transient data structures and manuscript outputs cleared. Raw field payloads preserved."