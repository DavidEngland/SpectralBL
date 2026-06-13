# GABLS3 NetCDF Schema Map (Set A / Set C)

Source file: `data/gabs3/gabls3_scm_cabauw_obs_v33.nc`

## Dimensions
- `time = 144`
- `zf = 38` (mean-profile levels)
- `zh = 5` (flux-profile levels)
- `zt = 38`
- `zs = 8`

## Coordinates
- `time(time)` units: `hours since 2006-07-01 00:00:00 0:00`
- `zf(zf,time)` units: `m` (mean-profile heights)
- `zh(zh,time)` units: `m` (flux-profile heights)

## Set A (Predictive Inputs X)
- `u(zf,time)` units: `m/s`
- `v(zf,time)` units: `m/s`
- `th(zf,time)` units: `K`

Current ingestion uses `u` and `th` directly for projection/diagnostics.
`v` is retained for future multicomponent feature engineering.

## Set C (Targets Y)
- `wt(zh,time)` units: `K m/s` (resolved heat flux)
- `uw(zh,time)` units: `m2/s2` (resolved momentum flux u-component)
- `vw(zh,time)` units: `m2/s2` (resolved momentum flux v-component)

## Window Tags
- `HOUR8`: `7.0 <= time < 8.0`
- `HOUR9`: `8.0 <= time < 9.0`
- `OTHER`: all remaining times

## Notes
- Height coordinates are time-indexed in this file (`zf(:,t)`, `zh(:,t)`), so ingestion must read geometry per time slice.
- Projection path is geometry-aware and uses per-slice physical heights before Chebyshev mapping.
