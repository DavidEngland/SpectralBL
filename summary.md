# Project Summary: CASES-99 Boundary Layer Manifold Projection Pipeline

## 1. Project Objective & Scope

The goal of this project is to process, clean, and project micro-meteorological tower data from the **CASES-99 (Cooperative Atmosphere-Surface Exchange Study 1999)** intensive field campaign into a stable mathematical manifold.

Using high-resolution instrumentation snapshots, the pipeline strips out sensor noise and drops, reconstructing clean vertical profiles for wind velocity ($u$) and potential temperature ($\theta$). These continuous profiles allow the system to extract precise secondary fluid-dynamics metrics—such as the Flux Richardson Number ($Ri_f$), spectral entropy diversity ($D_{\text{eff}}$), and scalar destruction rates ($\chi_N$)—to analyze atmospheric boundary layer stability and turbulence transitions.

---

## 2. System Architecture & Component Mapping

The pipeline is organized into a modular framework in Julia, divided into ingestion, deconstruction, and reporting layers:

```
                      [ Raw CASES-99 NetCDF Stream ]
                                    │
                                    ▼
       ┌─────────────────────────────────────────────────────────┐
       │                  module CasesIngestion                  │
       ├─────────────────────────────────────────────────────────┤
       │  • 7-Point Vertical Tower Extraction (1.5m to 50m)      │
       │  • Dynamic Schema Fallbacks ("tc_zm" vs "T_zm")         │
       │  • Type-Safety Fence & Out-of-Bounds Gating             │
       └────────────────────────────┬────────────────────────────┘
                                    │ (Design Matrix H & Vectors)
                                    ▼
       ┌─────────────────────────────────────────────────────────┐
       │              project_with_svd_truncation                │
       ├─────────────────────────────────────────────────────────┤
       │  • High-Order Chebyshev Column Normalization            │
       │  • Economy Singular Value Decomposition (SVD)           │
       │  • Fixed Machine Precision Rank Preservation Floor      │
       └────────────────────────────┬────────────────────────────┘
                                    │ (Spectral Coefficients c)
                                    ▼
       ┌─────────────────────────────────────────────────────────┐
       │             Downstream Processing & Reporting           │
       ├─────────────────────────────────────────────────────────┤
       │  • Secondary Physical Diagnostic Calculations           │
       │  • Campaign-Oriented Daily Processing and Aggregation   │
       │  • Automated Markdown and manuscript figure rendering   │
       └─────────────────────────────────────────────────────────┘

```

### Core Architecture Components:

* **Data Ingestion Layer (`CasesIngestion.jl`):** Dynamically tracks and maps EOL instrument heights ($z = [1.5, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0]\text{ m}$) across varying NetCDF variable conventions. Intercepts missing payload values and filters out non-physical sensor dropouts.
* **Spectral Projection Layer (`project_with_svd_truncation`):** Maps discrete physical elevations into a continuous spectral space using high-degree Chebyshev polynomials ($N = 32$). It scales columns to prevent Vandermonde-style conditioning degradation, applies an economy Singular Value Decomposition, and scales coefficients back into matching physical units.
* **Reporting Layer (`scripts/Report.jl`):** Generates diagnostic geometry assets, state-space planes, temporal traces, and manuscript-facing figures from the trajectory outputs.
* **Synoptic Layer (`scripts/run_synoptic_analysis.jl`):** Aggregates campaign windows and writes timeline diagnostics including `campaign_synoptic_evolution.pdf` and synoptic audit outputs.

---

## 3. Mathematical Resolution & Boundary Conditions

A primary mathematical focus of the project is resolving the **Rank Collapse** caused by over-parameterization. Because the system projects a continuous space using a 33-column polynomial basis matrix ($H$) sampled at only 7 physical tower locations, columns beyond the first few are highly collinear.

To prevent the SVD solver from aggressively truncating these modes down to an uninformative `Rank=1` state, the mathematical constraints were hardened:

1. **Column Regularization:** Matrix columns are scaled by their $L_2$ norm ($\|H_{:, j}\|_2$) before factorization to balance out the extreme numerical ranges of high-order polynomials.
2. **Precision Floor Override:** The relative singular value truncation threshold is hardcoded to a precise floor ($1e-12 \times \sigma_1$), overriding loose upstream parameter fallbacks. This forces the solver to preserve all 7 available physical degrees of freedom.

---

## 4. Current Operational Status & Deliverables

* **Ingestion Pipeline:** Operational and production-ready. Successfully isolates valid data from NaN drops across the complete 42-day campaign sub-directory.
* **Data Artifacts:** Generates `data/diagnostic_trajectory.csv` tracking fluid states, advection tendencies, and turbulence signatures over the timeline.
* **Reporting Layer:** Integrated into the project `Makefile`. Executing `make report` writes diagnostics under `reports/ncar_eol_dee0099881/`, and manuscript-ready figure copies under `data/drafts/figures/`.
* **Synoptic Diagnostics:** `scripts/run_synoptic_analysis.jl` writes `campaign_synoptic_evolution.pdf` to both report and manuscript-figure directories.
* **Synthetic Sponge Validation:** `scripts/test_wave_reflection.jl` writes `data/wave_reflection_test.png`, `data/wave_reflection_test.pdf`, and `data/wave_reflection_metrics.csv`.
* **Universal Sponge Replication:** `scripts/run_universal_sponge_test.jl` writes per-campaign artifacts under `data/universal_sponge/<campaign_lowercase>/`.