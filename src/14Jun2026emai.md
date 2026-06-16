Dear Dick and Arastoo,

I wanted to provide an update on the status of the `SpectralBL` toolkit as my contract transition approaches, particularly regarding several developments that emerged while running the analysis pipeline against the CASES-99 and GABLS3 (I found a NetCDF for GABLS3) datasets.

One note before diving into the technical details: the public repository currently lags substantially behind the active development branch. I have intentionally delayed merging several of these changes because the CASES-99 manuscript analysis is now producing stable and reproducible results, and I wanted to avoid introducing GABLS3-related modifications that might inadvertently affect the existing workflow. The core framework is considerably more mature than what is currently reflected in the repository.

A second practical consideration is computational cost. Throughout this work I have been operating at approximately five times my allocated AI-compute budget in order to accelerate testing, diagnostics, manuscript preparation, and code development. Consequently, a significant emphasis has been placed on developing methods that remain computationally lean while still capturing the underlying nonlinear dynamics.

### 1. The Multi-Regime Nature of MOST

A recurring issue in Stable Boundary Layer (SBL) analysis is that many diagnostic approaches implicitly rely on near-neutral linearizations or piecewise-linear similarity assumptions. These approximations often perform reasonably well within a single regime but become problematic during transitions between stable and convective states.

The CASES-99 profiles reveal that the boundary layer does not evolve along a static equilibrium path. Instead, it repeatedly expands and contracts through the diurnal cycle—effectively “breathing”—and produces a coherent trajectory in state space. Traditional profile fitting tends to average over this motion, obscuring the underlying dynamical structure.

The boundary layer behaves less like a collection of independent profiles and more like a low-dimensional attractor that is continuously traversed as radiative cooling, surface heating, shear production, and turbulent mixing compete throughout the day.

### 2. Metric-Consistent Low-Rank Optimization

To isolate this structure without introducing significant computational overhead, I implemented a Thin Singular Value Decomposition directly within the $p$-FEM mass metric using $M^{1/2}Y$. This formulation preserves the physical inner product induced by the finite-element mass matrix, ensuring that modal energy rankings remain physically meaningful and are not distorted by grid spacing or numerical weighting choices.

Sparse tower observations are incorporated through an adaptive observation operator $A$, and the reduced-state coordinates $\eta$ are obtained by solving:

$$\min_{\eta} \|A U_r \eta - b\|_2^2 + \lambda \|\eta\|_2^2$$

For low truncation ranks ($r \le 3$), the optimization executes in microseconds per time step while naturally accommodating missing observations and irregular sampling without gap-filling.

### 3. Geometric Diagnostics of Regime Transitions

Once projected into a reduced coordinate space $(\eta_1, \eta_2, \eta_3)$, three diagnostic quantities emerge naturally:

* **Singular Value Entropy ($H$):** The Shannon entropy computed from the normalized singular value spectrum, measuring the effective dimensionality and structural complexity of the evolving boundary layer.
* **Phase Curvature ($\kappa$):** Quantifies abrupt directional changes in attractor motion and highlights rapid stability transitions, frontal passages, and low-level jet decay.
* **Attractor Spin ($\Omega$):**
$$\Omega = \eta_1 \dot{\eta}_2 - \eta_2 \dot{\eta}_1$$

Using robust central differencing, the sign of $\Omega$ acts as a signed rotational phase metric, identifying the direction of attractor rotation within the reduced phase space. Preliminary results suggest that the sign of this spin separates two distinct transition pathways:

1. Shear-driven breakout events associated with turbulence regeneration.
2. Radiatively driven stabilization events associated with collapse toward strongly stable conditions.

If this interpretation continues to hold across both CASES-99 and GABLS3, $\Omega$ may provide a compact scalar diagnostic of transition directionality, potentially identifying when the boundary layer is moving toward or away from regimes where classical MOST assumptions become unreliable.

### Next Steps and Repository Structure

The reporting layer has now been fully separated from the simulation and transformation framework. Diagnostic extraction, figure generation, and manuscript assets are being migrated into an independent reporting repository, allowing campaign-specific analyses to be generated through simple command-line workflows such as:

```bash
julia extract_attractor_diagnostics.jl GABLS3

```

Generated figures are written directly into the active manuscript directories, insulating the scientific narrative from changes occurring within the computational framework.

I would be interested in your thoughts regarding the phase-space interpretation of these SBL transitions, particularly whether the attractor spin diagnostic ($\Omega$) might provide a useful pathway for constraining sub-grid parameterizations in regimes where classical MOST assumptions begin to break down.

Cheers,

Dave

SpectralBL Development Team

University of Alabama in Huntsville (UAH)

**Repository:** [https://github.com/DavidEngland/SpectralBL](https://github.com/DavidEngland/SpectralBL)