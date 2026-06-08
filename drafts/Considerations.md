## 1. Diagnostics to Replace MOST in the HSNBL

Monin-Obukhov Similarity Theory (MOST) breaks down in the Highly Stable Nocturnal Boundary Layer (HSNBL) primarily because the continuous, body-force driven turbulence ceases, giving way to global intermittency, gravity waves, and sub-mesoscale drainage flows (Steeneveld, 2014). To track and parameterize these conditions without running full, high-overhead multi-order closure schemes (like Mellor-Yamada-Nakanishi-Niino), several low-overhead shortcuts and diagnostic framework replacements exist.

### A. Local Scaling & The Flux-Richardson Framework (Nieuwstadt Scaling)

Instead of surface-based parameters ($u_*$, $\theta_*$, $L$), the easiest diagnostic pivot is to drop MOST entirely and switch to **Local Scaling** (Nieuwstadt scaling) (Steeneveld, 2014). This replaces the Monin-Obukhov length $L$ with a local length scale $\Lambda$:

$$\Lambda = -\frac{\tau^{3/2}}{\kappa \frac{g}{\theta_0} \overline{w'\theta'}}$$

where $\tau$ and $\overline{w'\theta'}$ are the local momentum and heat fluxes extracted directly at height $z$ from your numerical grid, rather than assuming constant flux layers.

* **The Shortcut:** Pair local scaling with the **Beljaars-Holtslag formulation**, which enforces a non-critical Flux Richardson number ($R_f \rightarrow 1$) (Casasanta et al., 2021). This prevents the numerical model from completely artificial "shut-off" behavior where vertical mixing drops to zero during extreme stability.

### B. Machine Learning (ML) Shortcuts: Physics-Informed ML Surrogates

Rather than trying to capture complex non-local wave-turbulence interactions analytically, modern workflows deploy small neural networks or random forests as structural flux/profile surrogates.

* **The Input State:** Feed the model local grid diagnostics: the non-local bulk Richardson number $R_{b}$, the vertical grid resolution $\Delta z$, the local buoyancy frequency $N = \sqrt{\frac{g}{\theta}\frac{\partial\theta}{\partial z}}$, and the geostrophic wind shear.
* **The Output:** Instead of solving transcendental MOST equations iteratively, the ML surrogate outputs the stability functions ($\phi_m, \phi_h$) or directly yields the effective eddy diffusivities ($K_m, K_h$). This reduces computational overhead to a simple vector-matrix multiplication per grid column.

### C. The Multi-Limit Height Diagnostic

If your ultimate goal is diagnosing the actual depth ($h$) of the HSNBL without tracking complex budget equations, a highly effective dimensional shortcut ignores MOST entirely and leverages the **Steeneveld-van de Wiel-Holtslag** diagnostic regime (Steeneveld et al., 2007). For moderate to highly stable regimes, the boundary layer depth scales entirely based on the surface buoyancy flux ($B_s$) and the free-flow stability ($N$):

$$h \sim \left(\frac{|B_s|}{N^3}\right)^{1/2}$$

This offers a massive structural shortcut because it bypasses calculating friction velocity or dealing with iterative $z/L$ convergence loops entirely (Steeneveld et al., 2007).

---

## 2. Complementary NetCDF Datasets (Order of Priority)

When branching out from CASES-99 to evaluate or train your structural models, datasets should be prioritized based on how well they resolve the vertical structure of the stable boundary layer, their data clean-room quality, and how straightforward they are to ingest into NetCDF processing pipelines (`xarray`, `netCDF4`).

The recommended sequence of data ingestion is structured below:

| Order | Dataset / Campaign | Primary SBL Regime Covered | Key Advantage for Numerical Work |
| --- | --- | --- | --- |
| **1** | **GABLS Profiles** (GABLS1 / GABLS2 / GABLS3) | IDEALized SBL to Moderately Stable | The global standard for single-column model benchmarking. Clean, pre-formatted baseline NetCDF datasets designed specifically to test where MOST fails (JAX-SCM v1.0, 2026). |
| **2** | **SHEBA** (Surface Heat Budget of the Arctic Ocean) | Highly Stable, Perennial Cryosphere | Critical for testing HSNBL without diurnal transitions. Features massive data logs under extreme stability conditions ($\zeta > 1$) over Arctic ice (Casasanta et al., 2021). |
| **3** | **Cabauw Observatory Data** (CESAR - Netherlands) | Weakly Stable to Moderately Intermittent | Decades of continuous tall-tower (200m) profiling. Perfect for validating bulk/gradient curvature diagnostics because of the long vertical array of sensor baselines (Casasanta et al., 2021). |
| **4** | **ISOBAR** (Innovative Strategies for Observations in the Arctic Atmospheric Boundary Layer) | Extreme Stable, Coastal Sea-Ice | Modern campaign utilizing Unmanned Aerial Vehicles (UAVs) alongside surface masts (Wenta et al., 2020). Provides high vertical resolution profile sheets over sea ice in NetCDF format via repositories like PANGAEA. |
| **5** | **MARCUS** (Measurements of Aerosols, Radiation, and Clouds over the Southern Ocean) | Marine / High-Latitude Stable Boundary Layer | Shipborne profiling over varying sea-ice concentrations (Knight et al., 2024). Excellent if your curvature and Richardson corrections need to account for high-latitude marine/ice boundary transitions. |

### Rationale for this Ingestion Order:

1. Start with **GABLS** because the files are already perfectly formatted for numerical intercomparisons and synthetic column setups (JAX-SCM v1.0, 2026).
2. Move to **SHEBA** and **Cabauw** to get real-world, high-stability data over flat terrain to stress-test your gradient curvature corrections (Casasanta et al., 2021).
3. Conclude with **ISOBAR** and **MARCUS** to introduce complex transitions (such as coastal boundaries and sea-ice interactions) where traditional local similarity diagnostics completely fall apart (Knight et al., 2024; Wenta et al., 2020).

## References

* Casasanta, G., Sozzi, R., Petenko, I., & Argentini, S. (2021). Flux–Profile Relationships in the Stable Boundary Layer—A Critical Discussion. *Atmosphere*, *12*(9), 1197. [https://doi.org/10.3390/atmos12091197](https://www.google.com/search?q=https://doi.org/10.3390/atmos12091197)
Cited by: 6
* JAX-SCM v1.0: a modern atmospheric single-column model for boundary layer research. (2026). *arXiv preprint arXiv:2605.24544*. [https://arxiv.org/abs/2605.24544](https://www.google.com/search?q=https://arxiv.org/abs/2605.24544)
* Knight, C. L., Mallet, M. D., Alexander, S. P., Fraser, A. D., Protat, A., & McFarquhar, G. M. (2024). Cloud Properties and Boundary Layer Stability Above Southern Ocean Sea Ice and Coastal Antarctica. *Journal of Geophysical Research: Atmospheres*, *129*(10). [https://doi.org/10.1029/2022jd038280](https://www.google.com/search?q=https://doi.org/10.1029/2022jd038280)
Cited by: 4
* Steeneveld, G. J. (2014). Current challenges in understanding and forecasting stable boundary layers over land and ice. *Frontiers in Environmental Science*, *2*, 41. [https://doi.org/10.3389/fenvs.2014.00041](https://www.google.com/search?q=https://doi.org/10.3389/fenvs.2014.00041)
Cited by: 69
* Steeneveld, G. J., van de Wiel, B. J. H., & Holtslag, A. A. M. (2007). Diagnostic Equations for the Stable Boundary Layer Height: Evaluation and Dimensional Analysis. *Journal of Applied Meteorology and Climatology*, *46*(2), 212–225. [https://doi.org/10.1175/jam2454.1](https://www.google.com/search?q=https://doi.org/10.1175/jam2454.1)
Cited by: 85
* Wenta, M., Brus, D., Doulgeris, K., Vakkari, V., & Herman, A. (2020). Winter atmospheric boundary layer observations over sea ice in the coastal zone of the Bothnian Bay (Baltic Sea). *Earth System Science Data*, *13*(1), 33–49. [https://doi.org/10.5194/essd-2020-153](https://www.google.com/search?q=https://doi.org/10.5194/essd-2020-153)
Cited by: 6