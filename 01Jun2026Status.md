# 01 Jun 2026 Status: Synthetic Wave Sponge Validation

## Objective
Validate the spectral sponge-layer hypothesis in a controlled synthetic setting before campaign-level coupling.

Hypothesis: high-frequency modal damping weighted by `psi_T` reduces upper-boundary wave reflection energy relative to an unmasked run.

## Implementation Summary

### Added Script
- `scripts/test_wave_reflection.jl`

### Experiment Design
- Constructs `UnifiedManifoldWorkspace(32, 1.5, 50.0, 0.05, K_q=72)`.
- Initializes a Gaussian-modulated sinusoidal wave packet centered near 35 m.
- Evolves the linear wave equation via RK4 using `Dz_atm * Dz_atm` as the second-derivative operator.
- Runs two scenarios:
  1. With sponge: modal damping applied each step using `exp(-damping_rate * psi_T[n])`.
  2. Without sponge: no spectral damping.
- Tracks top-boundary energy time series and computes a late-window suppression ratio.

### Diagnostics and Artifacts
- `data/wave_reflection_test.png`
- `data/wave_reflection_metrics.csv`

### Build Integration
- Added `make wave_test` for standalone execution.
- Added `make universal_wave_test` for campaign-geometry replication execution.
- Added `make full` for aggregate execution (`setup + validate + run + report + wave_test`).
- Kept `make all` unchanged to avoid default pipeline regression.

## Replication Framework
- Added `scripts/run_universal_sponge_test.jl` to mirror the same sponge verification strategy across campaign geometries.
- Predefined campaign profiles include CASES_99, FLOSS_II, SHEBA, and BLLAST, with optional custom geometry arguments.
- Per-campaign outputs are written to `data/universal_sponge/<campaign>/`.

## Acceptance Criteria
- Script executes under project environment.
- Wave diagnostics artifacts are generated and non-empty.
- Suppression ratio > 1 indicates reflected-energy reduction with sponge masking.

## Notes
This pass is intentionally synthetic-only and does not modify campaign ingestion logic (`CasesIngestion.jl`) or NetCDF schema flow.
