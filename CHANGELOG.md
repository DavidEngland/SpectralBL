# Changelog

All notable changes to this repository are documented in this file.

## 2026-06-01

### Added
- Added synthetic sponge-layer validation script at `scripts/test_wave_reflection.jl`.
- Added universal replication script at `scripts/run_universal_sponge_test.jl` for CASES-99, FLOSS-II, SHEBA, and BLLAST geometry testing.
- Added `make wave_test` target for opt-in wave reflection suppression testing.
- Added `make universal_wave_test` target for campaign-geometry replication testing.
- Added `make full` aggregate target to include the wave validation run.
- Added output artifacts for wave experiment diagnostics:
  - `data/wave_reflection_test.png`
  - `data/wave_reflection_metrics.csv`

### Changed
- Updated `README.md` with instructions for synthetic wave validation and new Make targets.
- Updated `Makefile` clean target to remove generated wave-test artifacts.
