# Unified Manifold Workspace Summary Report

## Geometry & Mapping Parameters
* **Spectral Modes (N):** 32
* **Lower Boundary (z_0m):** 1.5 m (CASES-99 Tower Base)
* **Top Boundary (z_top):** 50.0 m

## Numerical Health Diagnostics
* **Minimum Grid Spacing (\Delta z_min):** Near surface refinement active.
* **Maximum Grid Spacing (\Delta z_max):** Smooth stretching toward canopy.

## Partitioning Thresholds
* **Meso Windows (n_m):** Modes capturing stable, large-scale structures.
* **Wave/Sub-meso (n_w):** Internal gravity wave transitions.
* **Turbulent Residual:** High-frequency dissipation modes.

The physical mapping successfully packs resolution into the intense stable nocturnal boundary layers characteristic of the CASES-99 campaign. High-frequency modes are filtered via the sub-meso partition matrix to prevent unphysical numerical reflections near the top boundary ($50\text{ m}$).

*Report generated automatically by the UnifiedManifold pipeline.*
