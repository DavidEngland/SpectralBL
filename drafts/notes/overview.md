# The `SpectralBL-Analytics` Core Architecture

The entire pipeline can be formalized as an information-preserving reduction mapping sequence:

$$\mathbf{b}(t) \xrightarrow{\text{Section 1}} \mathbf{u}(z,t) \xrightarrow{\text{Section 2}} \mathbf{\tilde{Y}} \xrightarrow{\text{Section 3}} \mathbf{U}_r \xrightarrow{\text{Section 4}} \boldsymbol{\eta}(t) \xrightarrow{\text{Section 5}} (R, \Omega) \xrightarrow{\text{Section 6}} \dot{\boldsymbol{\eta}} = \mathbf{F}(\boldsymbol{\eta})$$

```
[Raw Sensors b(t)] ──► [p-FEM Expansion u(z,t)] ──► [Mass-Weighted tSVD Ỹ] ──► [Attractor Projection η(t)] ──► [Polar Metrics (R, Ω)] ──► [SINDy Forecasting]

```

---

## Section 1: Front-End Data Ingestion and Continuous Field Expansion

### 1.1 The Domain Problem

Raw meteorological tower measurements arrive as a discrete, sparse, and often irregularly spaced vector of observations $\mathbf{b}(t) = [U(z_m, t), \Theta(z_m, t)]^T$ sampled at $M$ discrete sensor heights $z_m$.

### 1.2 Mathematical Formulation

To eliminate dependence on the arbitrary physical placement of tower booms, the front-end project maps these discrete measurements into a continuous vector-valued Hilbert space $\mathcal{H}(\Omega_z)$ spanning the vertical domain $\Omega_z = [z_0, z_{\text{top}}]$. This is accomplished via an analytic partition-of-unity piece-wise finite-element ($p$-FEM) expansion:

$$\mathbf{u}(z, t) = \begin{bmatrix} U(z, t) \\ \Theta(z, t) \end{bmatrix} = \sum_{i=1}^{N} \mathbf{c}_i(t) \phi_i(z)$$

where $\phi_i(z)$ are deterministic, localized shape functions with overlapping support satisfying the partition-of-unity constraint $\sum_{i} \phi_i(z) = 1$. This interpolation step shields downstream calculations from data gaps and localized grid irregularities (*Li et al., 2023; Hui et al., 2024*).

---

## Section 2: Metric-Consistent Manifold Discretization

### 2.1 The Domain Problem

If a standard algebraic Principal Component Analysis (PCA) or Euclidean Singular Value Decomposition (SVD) is performed on a field snapshot matrix $\mathbf{Y}$, the resulting spatial modes become direct functions of the mesh resolution and local node clustering. If sensors are closely grouped near the canopy, the Euclidean norm artificially over-weights that region.

### 2.2 Mathematical Formulation

To ensure that the learned coordinates reflect continuous physical energy invariants rather than discretization geometry, the pipeline constructs a symmetric, positive-definite $p$-FEM mass matrix $\mathbf{M}$ (*Yadalam & Feeny, 2011; Davis & Constantinescu, 2023*):

$$\mathbf{M}_{ij} = \int_{z_0}^{z_{\text{top}}} \phi_i(z) \phi_j(z) \, dz$$

The raw snapshot matrix $\mathbf{Y}$ is transformed into a metric-consistent snapshot array $\mathbf{\tilde{Y}}$ by applying the symmetric square root of this mass matrix (*Van Schie et al., 2025*):

$$\mathbf{\tilde{Y}} = \mathbf{M}^{1/2} \mathbf{Y}$$

This mass-matrix weighting guarantees that modal energy is measured with respect to the continuous physical inner product $\langle u,v \rangle_{\mathbf{M}} = u^T \mathbf{M} v$ rather than the Euclidean norm of the observation vector, eliminating dependence on sensor spacing and mesh resolution.

---

## Section 3: Empirical Attractor Filtering ($\boldsymbol{\psi}$ Masks)

### 3.1 The Domain Problem

Isolating the underlying geometric surface—the **Nocturnal Manifold** $\mathcal{M}$—from complex, non-equilibrium atmospheric flows requires finding an optimal, low-rank basis that remains valid across highly stratified and intermittent regimes (*Poveda-Jaramillo & Puente, 1993; Howell & Mahrt, 1997*).

### 3.2 Mathematical Formulation

A Thin Singular Value Decomposition (tSVD) is executed on the metric-consistent snapshot array:

$$\mathbf{\tilde{Y}} = \mathbf{\tilde{U}} \mathbf{\Sigma} \mathbf{V}^T$$

Because $\mathbf{\tilde{U}}$ lives within the artificially weighted metric space, the true physical, unweighted empirical orthogonal mode matrix $\mathbf{U}_r$ (truncated at rank $r=3$) is recovered by inverting the mass tensor:

$$\mathbf{U}_r = \mathbf{M}^{-1/2} \mathbf{\tilde{U}}_r = \begin{bmatrix} \vert & \vert & \vert \\ \boldsymbol{\psi}_1(z) & \boldsymbol{\psi}_2(z) & \boldsymbol{\psi}_3(z) \\ \vert & \vert & \vert \end{bmatrix}$$

These **$\boldsymbol{\psi}$ masks** are globally orthogonal with respect to the continuous column mass: $\boldsymbol{\psi}_i^T \mathbf{M} \boldsymbol{\psi}_j = \delta_{ij}$. Across CASES-99 and GABLS3 benchmark runs, these profiles appear as emergent physical structures matching the continuous governing fluid constraints:

* **$\boldsymbol{\psi}_1(z)$ (Emergent Bulk Inversion Mask):** Monotonic structure ($N_{\text{zeros}} = 0$) tracking macro-scale radiative cooling and boundary layer inversion depth.
* **$\boldsymbol{\psi}_2(z)$ (Emergent Shear/LLJ Mask):** Biphasic structure ($N_{\text{zeros}} = 1$) isolating the mechanical wind shear production and Low-Level Jet (LLJ) altitude migrations aloft.
* **$\boldsymbol{\psi}_3(z)$ (Emergent Structural Curvature Mask):** Triphasic structure ($N_{\text{zeros}} = 2$) capturing submesoscale perturbations, internal gravity waves, and localized pre-concursive gradient distortions prior to turbulence collapse (*Basu, 2006; Shimizu & Kawahara, 2017*).

---

## Section 4: Reduced Latent Space Estimation (The Observation Operator)

### 4.1 The Domain Problem

In operational settings, field sensors may experience dropouts, vary in altitude, or provide highly sparse profiles. The framework must estimate the global latent manifold coordinates safely from imperfect, sparse data vectors $\mathbf{b}(t)$ (*Sondak & Smith, 2024; Cyr et al., 2024*).

### 4.2 Mathematical Formulation

Let $\mathbf{A}$ be a sparse linear **observation operator** that samples continuous profiles specifically at active sensor heights. The time-dependent coordinate vector $\boldsymbol{\eta}(t) = [\eta_1(t), \eta_2(t), \eta_3(t)]^T$ is isolated by solving a regularized adaptive optimization problem:

$$\boldsymbol{\eta}(t) = \arg\min_{\boldsymbol{\eta}} \|\mathbf{A} \mathbf{U}_r \boldsymbol{\eta}(t) - \mathbf{b}(t)\|_2^2 + \lambda \|\boldsymbol{\eta}(t)\|_2^2$$

This formulation uses a Tikhonov regularization term ($\lambda \|\boldsymbol{\eta}\|_2^2$) to stabilize the projection against high-frequency instrumental noise. The full continuous profiles across all heights can then be analytically recovered at any instant via $\mathbf{u}_{\text{recon}}(z,t) = \sum_{k=1}^3 \eta_k(t)\boldsymbol{\psi}_k(z)$.

---

## Section 5: Topological Dynamics and Non-Equilibrium Polar Metrics

### 5.1 The Domain Problem

Traditional Monin-Obukhov similarity parameters are fundamentally static and local, assuming a structural equilibrium that completely collapses during intermittent nocturnal shear-burst events. The pipeline requires diagnostics that track the non-equilibrium velocity of the system as it moves between regimes.

### 5.2 Mathematical Formulation

Because the bulk inversion coordinate ($\eta_1$) and mechanical shear coordinate ($\eta_2$) map the primary energetic plane of the system, we define a complete polar description of the manifold path via the radius $R(t)$ and the attractor spin $\Omega(t)$:

$$R(t) = \sqrt{\eta_1^2(t) + \eta_2^2(t)}$$

$$\Omega(t) = \boldsymbol{\eta} \wedge \dot{\boldsymbol{\eta}} = \det \begin{pmatrix} \eta_1(t) & \eta_2(t) \\ \dot{\eta}_1(t) & \dot{\eta}_2(t) \end{pmatrix} = \eta_1(t)\dot{\eta}_2(t) - \eta_2(t)\dot{\eta}_1(t) = \frac{dA}{dt}$$

* **The Manifold Radius $R(t)$** measures the total structural departure from the unstratified, near-neutral equilibrium core.
* **The Attractor Spin $\Omega(t)$** represents a differential two-form mapping the **signed areal sweep rate** $\frac{dA}{dt}$ in the reduced phase plane.

Instead of taking a static snapshot, $\Omega(t)$ tracks phase acceleration, giving a precise geometric meaning to the diurnal loop:

$$\Omega(t) = \frac{dA}{dt} \quad \Longrightarrow \quad \begin{cases}

> 0 & \text{\textbf{Clockwise Loop}: Radiative cooling driving the system toward stratified configurations.} \
> < 0 & \text{\textbf{Counter-Clockwise Loop}: Mechanical breakout and downward mixing during shear-burst events.} \
> \approx 0 & \text{\textbf{Quasi-Equilibrium Orbit}: Adherence to classical local MOST parameters.}
> \end{cases}$$

---

## Section 6: Nonintrusive System Identification and Operator Inference

### 6.1 The Domain Problem

To transform `SpectralBL` from a purely descriptive diagnostic framework into a predictive engine, the underlying vector field governing the trajectory must be learned directly from data.

### 6.2 Mathematical Formulation

Because the system is constrained to a low-dimensional manifold, its chronological evolution can be modeled as an autonomous dynamical system driven by a continuous vector field $\mathbf{F}(\boldsymbol{\eta})$ (*Peherstorfer & Willcox, 2016; Viola et al., 2025*):

$$\dot{\boldsymbol{\eta}}(t) = \mathbf{F}(\boldsymbol{\eta}(t)) \quad \Longrightarrow \quad \begin{bmatrix} \dot{\eta}_1 \\ \dot{\eta}_2 \\ \dot{\eta}_3 \end{bmatrix} = \begin{bmatrix} f_1(\eta_1, \eta_2, \eta_3) \\ f_2(\eta_1, \eta_2, \eta_3) \\ f_3(\eta_1, \eta_2, \eta_3) \end{bmatrix}$$

Using the coordinates extracted by your observation operator, a symbolic library $\boldsymbol{\Theta}(\boldsymbol{\eta})$ of nonlinear polynomial and trigonometric candidate functions is constructed (*Brunton et al., 2016; Rudy et al., 2017*):

$$\boldsymbol{\Theta}(\boldsymbol{\eta}) = \begin{bmatrix} 1 & \eta_1 & \eta_2 & \eta_1^2 & \eta_1\eta_2 & \eta_2^3 & \dots \end{bmatrix}$$

Applying Sequential Thresholded Ridge Regression solves for a sparse, parsimonious coefficient matrix $\boldsymbol{\Xi}$ (*McQuarrie et al., 2021; Aretz & Willcox, 2025*):

$$\min_{\boldsymbol{\Xi}} \|\dot{\boldsymbol{\eta}} - \boldsymbol{\Theta}(\boldsymbol{\eta})\boldsymbol{\Xi}\|_2^2 + \alpha \|\boldsymbol{\Xi}\|_2^2 \quad \text{subject to} \quad \|\boldsymbol{\Xi}_i\|_\infty > \gamma$$

To bypass numerical differentiation noise amplification from noisy field sensors, the pipeline leverages implicit parallel or decoder-enhanced formulations (*Champion et al., 2020; Meng et al., 2025*), ensuring the discovered governing equations preserve physical invariants. Once $\boldsymbol{\Xi}$ is identified, the non-linear limit cycles and sudden regime transitions of the stable boundary layer become completely deterministic and forecastable.