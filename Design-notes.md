## Architecture Design Note: Model Grid Extensibility & Roadmap

The current stable release of this framework optimizes the `UnifiedManifoldWorkspace` for sparse, non-uniform vertical physical observations ($M < N$) using an economy-size thin Singular Value Decomposition (SVD).

Expanding the framework to ingest dense, fully resolved numerical model outputs (e.g., NetCDF arrays from LES or WRF) requires decoupling the grid configuration from the core spectral solvers. Users planning to adapt this repository for model validation should consult the architecture roadmap sketched in `/docs/model_extension_roadmap.md`.

### Core Development Requirements for Model Ingestion:
1. **Grid Transformation Interface:** Direct projection via matrix multiplication is only valid if the model vertical coordinates are mapped via an analytical Jacobian to a standard domain $\xi \in [-1, 1]$. For arbitrary model vertical spacing, a preprocessing interpolation layer must be implemented.
2. **Dynamic Mass Matrix Allocation:** When evaluating full-domain models, the metric-weighted mass matrix $\mathbf{M}$ must be computed using quadrature weights that precisely match the model's native vertical grid layers to prevent spectral energy estimation errors.
3. **Horizontal Parallelization Vectorization:** Processing a 3D model grid requires a 2D spatial loop over all horizontal columns $(x, y)$ at every simulated timestep. The `svd_project` function should be parallelized using multi-threading (`@threads`) or GPU kernels to efficiently manage the throughput of high-resolution datasets.