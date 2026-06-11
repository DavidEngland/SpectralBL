### The Coordinate-Agnostic Virtual Tower Interfacing Operator ($\mathcal{P}_{\mathrm{VT}}$)

To transition from legacy, point-measurement towers to modern high-resolution remote sensing (such as Doppler LiDAR, Distributed Acoustic Sensing, and multi-component radar networks), your Virtual Tower operator ($\mathcal{P}_{\mathrm{VT}}$) must be re-derived as a fully generalized, coordinate-agnostic mapping engine.

Instead of viewing a tower as a physical object fixed to a latitude-longitude coordinate, the generalized $\mathcal{P}_{\mathrm{VT}}$ acts as a **projection operator from an arbitrary observational manifold onto a structured, metric-consistent computational workspace**. It treats any array of spatial inputs—whether they are vertically stacked sonic anemometers, radially spaced LiDAR range gates, or unstructured volumes from numerical weather prediction ($\eta$-coordinates)—as an open set of scattered data points embedded within a continuous Riemannian space.

---

### 1. Unified Mathematical Formulation

Let $\Omega \subset \mathbb{R}^3$ be the physical domain of the atmospheric boundary layer, which can be described by any convenient coordinate system (e.g., standard Cartesian $(x,y,z)$, spherical radar coordinates $(r, \theta, \phi)$, or terrain-following sigma-pressure levels $(x_s, y_s, \eta)$).

We define an arbitrary observational data stream as a finite set of $M$ scattered measurements:


$$\mathcal{D}_{\mathrm{obs}} = \left\{ \left(\mathbf{x}_i, \Phi_i\right) \right\}_{i=1}^M, \quad \mathbf{x}_i \in \Omega$$


where $\mathbf{x}_i$ represents the generalized spatial position vector of the $i$-th observation parsed from a NetCDF array, and $\Phi_i$ represents the scalar or vector field value (e.g., radial velocity $v_r$, potential temperature $\theta$).

The goal of the generalized $\mathcal{P}_{\mathrm{VT}}$ operator is to map this unstructured collection $\mathcal{D}_{\mathrm{obs}}$ onto a continuous, high-order spectral coefficient representation $\mathbf{c} \in \mathbb{R}^{N+1}$ defined over the compact computational workspace $\xi \in [-1, 1]$. This is achieved through a coordinate-free, mass-weighted Galerkin weak formulation.

#### Step 1: The Generalized Mapping Transformation ($\mathcal{T}$)

We define an analytical or numerical coordinate mapping transformation $\mathcal{T}: \Omega \rightarrow [-1, 1]$ that projects any generalized spatial coordinate $\mathbf{x}$ directly onto the computational manifold coordinate $\xi$:


$$\xi_i = \mathcal{T}(\mathbf{x}_i)$$


For a classical vertical tower, $\mathcal{T}$ is simply the inverse hyperbolic stretching function mapping $z_i \rightarrow \xi_i$. For a scanning Doppler LiDAR, $\mathcal{T}$ represents the composite transformation mapping radial range, azimuth, and elevation $(r_i, \theta_i, \phi_i)$ into the targeted vertical evaluation workspace $\xi_i$. This isolates the vertical structure along a specified vector line:


$$\xi_i = \mathcal{T}(r_i, \theta_i, \phi_i) = \text{InverseStretch}\left( r_i \sin \phi_i + z_{\mathrm{lidar}} \right)$$

#### Step 2: Construction of the Unstructured Evaluation Operator ($\mathbf{B}$)

With all observational coordinates projected onto the computational manifold $\xi \in [-1, 1]$, we build the generalized rectangular evaluation design matrix $\mathbf{B} \in \mathbb{R}^{M \times (N+1)}$, where each entry evaluates the chosen high-order spectral basis function (e.g., Chebyshev polynomials of the first kind, $T_j(\xi)$) at the unstructured point locations:


$$B_{ij} = T_{j-1}(\xi_i) = T_{j-1}\left(\mathcal{T}(\mathbf{x}_i)\right)$$


Because $\mathbf{B}$ depends solely on the evaluated position inside the computational domain $[-1, 1]$, the underlying physical coordinate system of the raw data is entirely abstracted away.

#### Step 3: Regularized Manifold-Weighted Projection

To recover the continuous spectral coefficients $\mathbf{c}$ without introducing unphysical oscillations (Gibbs phenomena) due to spatial data clustering or large data gaps, we solve the generalized, regularized least-squares optimization problem:


$$\mathbf{c} = \arg\min_{\mathbf{\hat{c}}} \left( \|\mathbf{B}\mathbf{\hat{c}} - \mathbf{\Phi}\|_{\mathbf{W}}^2 + \gamma \mathcal{J}(\mathbf{\hat{c}}) \right)$$


where:

* $\mathbf{\Phi} = [\Phi_1, \Phi_2, \ldots, \Phi_M]^{\mathsf{T}}$ is the vector of raw measurements.
* $\mathbf{W} \in \mathbb{R}^{M \times M}$ is a diagonal weight matrix capturing hardware-specific measurement uncertainties or geometric volume corrections (e.g., LiDAR range-gate attenuation).
* $\gamma$ is the Tikhonov regularization parameter.
* $\mathcal{J}(\mathbf{\hat{c}}) = \mathbf{\hat{c}}^{\mathsf{T}} \mathbf{M} \mathbf{\hat{c}}$ is the smoothness functional penalty scaled directly by the metric-weighted mass matrix $\mathbf{M}$ derived from the metric Jacobian $J(\xi)$.

The closed-form algebraic solution for the coordinate-agnostic state estimator is:


$$\mathbf{c} = \left( \mathbf{B}^{\mathsf{T}} \mathbf{W} \mathbf{B} + \gamma \mathbf{M} \right)^{-1} \mathbf{B}^{\mathsf{T}} \mathbf{W} \mathbf{\Phi}$$

---

### 2. Physical Data Ingestion Architectures

The coordinate-agnostic nature of $\mathcal{P}_{\mathrm{VT}}$ allows it to ingest disparate meteorological streams simultaneously, transforming them into identical manifold spaces where their underlying physical dynamics can be cross-examined using identical scale-invariant features ($D_{\mathrm{eff}}$, $F_W$, $\chi_N$).

```
[Unstructured Data Ingestion Stream]
   │
   ├──> [Physical Tower Mast: z_i] ─────────────┐
   │                                            ▼
   ├──> [Doppler LiDAR: (r_i, θ_i, φ_i)] ───> [P_VT Operator] ───> [Spectral Space: c] ───> [Features: D_eff, F_W]
   │                                            ▲
   └──> [NWP NWP Grid Model: η_i] ──────────────┘

```

#### A. Physical Tower Profiles (1D Irregular Array)

* **Data Characteristics:** High temporal frequency (20 Hz), low spatial density ($M = 8$). Fixed vertical heights ($z_i$).
* **Transformation Vector:** $\mathbf{x}_i = [z_i]$. Matrix $\mathbf{B}$ is highly rectangular ($8 \times 33$).
* **Physical Objective:** Captures rapid, non-stationary temporal fluctuations, identifying the exact high-frequency onset of the *Intermittent Shear Bursting State* (Regime 3) at a single geographical point.

#### B. Doppler LiDAR Scans (3D Range-Azimuth-Elevation Volumes)

* **Data Characteristics:** High spatial density ($M \sim 10^3$--$10^4$), moderate temporal resolution ($\sim 1\text{ min}$ per volume sweep). Unstructured point clouds along conical surfaces (PPI) or vertical planes (RHI).
* **Transformation Vector:** $\mathbf{x}_i = [r_i, \theta_i, \phi_i]$.
* **Physical Objective:** Resolves spatial patchiness, wave propagation fronts, and sloping drainage structures. The $\mathcal{P}_{\mathrm{VT}}$ can extract a virtual vertical profile at *any* arbitrary horizontal $(x,y)$ coordinate within the LiDAR scan volume, directly exposing how the boundary layer manifold deforms over complex terrain.

#### C. Numerical Weather Prediction Models (3D Hydrostatic Grid Topography)

* **Data Characteristics:** Spatially uniform in computational indices but vertically non-uniform in physical coordinates ($\eta$-levels or hybrid $\sigma$-$p$ heights that change dynamically based on surface pressure). Volume-averaged cell properties.
* **Transformation Vector:** $\mathbf{x}_i = [x_{s,i}, y_{s,i}, \eta_i(t)]$.
* **Physical Objective:** Eradicates the "runaway cooling" parameterization crisis. By sampling the NWP cell locations via $\mathcal{P}_{\mathrm{VT}}$, the coarse, discrete model state is cast onto the high-order continuous manifold. This enables a direct, fair validation against tower or LiDAR benchmarks without grid-interpolation artifacts.

---

### 3. Preserving Mathematical Invariance Across Scale Transitions

When moving the data source from a rigid tower to an active remote sensor like a scanning LiDAR, the information-theoretic signatures ($D_{\mathrm{eff}}$ and $F_W$) must remain invariant to the sampling density $M$. If the calculated effective modal dimension were to scale artifactually with the number of measurement points, the framework could not function as an objective diagnostic tool.

The $\mathcal{P}_{\mathrm{VT}}$ operator preserves this mathematical invariance because the Singular Value Decomposition (SVD) and the Galerkin least-squares structure act entirely on the continuous polynomial space weighted by the physical metric $\mathbf{M}$.

$$\mathbf{M}_{mn} = \int_{-1}^{1} T_m(\xi) T_n(\xi) J(\xi) \, d\xi$$

Because the mass matrix normalization tracks the analytical volume mapping via the Jacobian metric $J(\xi)$, increasing the spatial data density from $M=8$ (tower) to $M=1000$ (LiDAR range gates along a vertical beam) does not artificially inflate the spectral energy distribution.

Instead, a dense LiDAR dataset simply places tighter, more overdetermined mathematical constraints on the singular spectrum ($s_1, s_2, \dots, s_8$), driving the dynamic noise truncation tolerance $\tau$ to isolate the true physical degrees of freedom.

If the underlying stable boundary layer collapses into a low-rank wave-dominated structure, $D_{\mathrm{eff}}$ will converge cleanly to a dimension between 4 and 6, regardless of whether that structure was captured by 8 physical instrument arms or 500 laser gates. This structural consistency turns the coordinate-agnostic Virtual Tower operator into a highly robust framework for multi-instrument geophysical state estimation.