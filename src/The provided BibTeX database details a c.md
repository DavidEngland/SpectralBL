The provided BibTeX database details a comprehensive research ecosystem centered on data-driven modeling of complex, non-linear dynamical systems through the lens of **Koopman Operator Theory** and **Dynamic Mode Decomposition (DMD)**.

To help frame a literature review or track academic lineage, the bibliography can be systematically categorized into foundational theory, data-driven approximations, deep/online adaptations, and practical validation across scientific fields.

---

## 1. Foundational Theory & Spectral Analysis

This group traces the roots of using infinite-dimensional linear operators to represent finite-dimensional nonlinear systems.

* **Classical Frameworks:** The initial mathematical formulation by Koopman (1931) and its connection to ergodic theory by von Neumann (1932) established that nonlinear Hamiltonian systems could be analyzed via linear transformations in a Hilbert space.
* **Modern Revival:** Mezić (2005) bridged these historical foundations into modern engineering, demonstrating how spectral decompositions relate to model reduction. Lan and Mezić (2013) extended this by matching the Koopman spectrum to the global linearization of attraction basins.

## 2. Dynamic Mode Decomposition (DMD) & Linear Approximations

These works focus on the actual numerical algorithms used to extract linear structures from empirical spatio-temporal data matrices.

* **Core Algorithms:** Rowley et al. (2009) and Schmid (2010) formalized DMD for fluid flows, turning Koopman concepts into a practical data-driven computational tool. Kutz et al. (2016) provided the definitive textbook guide on standard DMD practices.
* **Extensions for Non-linear Observables:** Williams et al. (2015) introduced Extended DMD (EDMD) to lift raw coordinates into user-defined non-linear feature dictionaries. To circumvent hand-crafting these libraries, Li et al. (2017) integrated adaptive dictionary learning directly into the optimization pipeline.

## 3. Deep Koopman Embeddings & Neural Networks

When system trajectories become exceptionally complex, researchers rely on deep neural architectures to learn coordinate transformations (the "lifting function") automatically.

* **Time-Lagged Autoencoders:** Early paradigms mapping deep architectures to slow kinetics were defined by Wehmeyer and Noé (2018).
* **Universal Encoders:** Lusch et al. (2018) popularized using deep autoencoders to force continuous non-linear trajectories into a latent space where they are globally governed by a linear time-invariant system.
* **State-of-the-Art Adaptations:** Recent refinements like *ResKoopNet* (Xu et al., 2025) handle spectral residuals, while *Koopa* (Liu et al., 2023) focuses on learning highly non-stationary sequences with localized predictors.

## 4. Online, Streaming, & Time-Varying Systems

Real-world deployments face streaming data, drifting distributions, and shifting operational regimes, necessitating real-time model updates.

* **Recursive & Spare Approximations:** Sinha et al. (2019, 2023) pioneered operator-theoretic approaches for streaming, noisy environments. Hou et al. (2024) integrated sparse online learning in Reproducing Kernel Hilbert Spaces (RKHS) to optimize nonparametric data limits.
* **Regime Switching:** Loya and Tallapragada (2024) developed frameworks to update tracking matrices as systems cycle through contrasting dynamical environments using Grassmannian distances.

## 5. Conformal Prediction & Uncertainty Quantification

A recent and powerful trend is matching online data-driven predictors with distribution-free uncertainty bounds.

* **Foundations:** Based on algorithmic randomness principles laid down by Vovk et al. (1999) and Vovk et al. (2005), Angelopoulos and Bates (2021) formalized a general paradigm for distribution-free uncertainty bounds.
* **Conformal Online Learning:** Moving beyond static boundaries, Gibbs and Candès (2021) enabled adaptive tracking under distribution shifts. This directly catalyzed architectures like **COLoKe** (Conformal Online Learning of Koopman Embeddings) by Gao et al. (2025), which leverages conformal-style mechanisms to adaptively update model states only when prediction tolerances are statistically violated (Gao, 2025).

## 6. Real-World Applications & Physical Benchmarks

The efficacy of these methodologies is evaluated against real physical systems and noisy tracking problems:

* **Robotics & Fluid Dynamics:** Tracking soft robotic controls (Bruder et al., 2021) and defining magnetized plasma states (Kaptanoglu et al., 2020).
* **Macroenvironmental Forecasters:** Evaluated on non-stationary, multi-scale industrial sets like the *ETDataset* (Zhou, 2021) for power grid transformer temperatures and Beijing air quality matrices (Zhang et al., 2017).

---

### Structural Summary of Interconnectivity

```
[Historical Operator Foundations]
    │  (Koopman 1931; von Neumann 1932; Mezić 2005)
    ▼
[Numerical Implementations (DMD / EDMD)]
    │  (Rowley 2009; Schmid 2010; Williams 2015)
    ▼
[Deep Latent Embeddings (Autoencoders)] ──► [Conformal / Online Tracking]
    │  (Lusch 2018; Otto 2019)                   (Sinha 2023; Gao 2025)
    ▼
[Real-World Deployment & Closed-Loop Control]
       (Bruder 2021; Korda & Mezić 2018)

```

*(Note: The lone entry `@article{greenwade93}` regarding the Comprehensive TeX Archive Network stands completely isolated from this dense control theory and dynamical systems graph, serving exclusively as a utility citation for LaTeX software ecosystem infrastructure.)*

---

## References

* Angelopoulos, A. N., & Bates, S. (2021). A gentle introduction to conformal prediction and distribution-free uncertainty quantification. *arXiv preprint arXiv:2107.07511*.
* Gao, B. (2025). Conformal Online Learning of Deep Koopman Linear Embeddings. *Advances in Neural Information Processing Systems*. [URL](https://arxiv.org/abs/2511.12760).
Cited by: 1
* Gao, B., Patracone, J., Alata, O., & Chrétien, S. (2025). Apprentissage incrémental de l'opérateur de Koopman pour systèmes non-autonomes via prédiction conforme. *GRETSI'25 -- XXXe Colloque Francophone de Traitement du Signal et des Images*.