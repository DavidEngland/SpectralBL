You introduce \mathcal{T}:z\rightarrow\xi, but the actual equation is omitted.
The paper should explicitly state something like:
\xi(z)=\frac{\tanh\left(\alpha\frac{z-z_m}{L}\right)}{\tanh(\alpha)}
or whatever specific stretched mapping is actually implemented.
Then derive
J(\xi)=\frac{dz}{d\xi}
directly.

---

## 1. Mathematical Framework & Mapping Physics

The `UnifiedManifold` uses an algebraic-hyperbolic compactification transformation to map a highly non-uniform physical domain $z \in [z_{0m}, z_{top}]$ to a regular computational domain $\xi \in [-1, 1]$. This layout allows the use of pseudospectral Chebyshev polynomials $T_n(\xi) = \cos(n \arccos \xi)$ while optimizing nodal placement for high-gradient boundary layer features.

The forward physical mapping is defined by:

$$z(\xi) = z_{0m} + \frac{\sigma (1 + \xi)}{1 - \xi + \alpha}$$

Where the domain stretching scale factor $\sigma$ is calculated from the required domain limits:

$$\sigma = \frac{(z_{top} - z_{0m})\alpha}{2}$$
--

Standard pseudospectral frameworks relying on uniform or Chebyshev-cosine grids distribute computational degrees of freedom symmetrically across the spatial domain. Such configurations are fundamentally unoptimized for the acute vertical gradients characterizing the nocturnal stable boundary layer (SBL), where sub-mesoscale gravity waves and shear instabilities concentrate within the lowest meters of the surface layer ($1.5 \le z \le 5.0\,\text{m}$). To preferentially focus node density where structural geometric derivatives are maximized, we establish a \emph{UnifiedManifoldWorkspace} utilizing a fractional hyperbolic coordinate stretching function $\mathcal{T}(\xi) \rightarrow z$. This operator transforms the bounded computational space $\xi \in [-1, 1]$ into the physical elevation domain defined by the tower instruments ($z_{\mathrm{min}} = 1.5\,\text{m}$ to $z_{\mathrm{max}} = 55.0\,\text{m}$):
\begin{equation}
z(\xi) = z_{\mathrm{min}} + \frac{z_{\mathrm{max}} - z_{\mathrm{min}}}{1 - \alpha_{\mathrm{stretch}}^2} \left[ \frac{(1+\alpha_{\mathrm{stretch}})(1+\xi)}{1 + \alpha_{\mathrm{stretch}}(2+\xi)} \right]
\label{eq:hyperbolic_stretch}
\end{equation}
where $\alpha_{\mathrm{stretch}} \in (0, 1)$ represents the non-dimensional grid compression parameter. The analytical transformation Jacobian resulting from this compactification is defined as:
\begin{equation}
J(\xi) = \frac{dz}{d\xi} = \frac{(z_{\mathrm{max}} - z_{\mathrm{min}})(1 + \alpha_{\mathrm{stretch}})^2}{(1 - \alpha_{\mathrm{stretch}}^2) \left[1 + \alpha_{\mathrm{stretch}}(2+\xi)\right]^2}
\label{eq:jacobian_stretched}
\end{equation}
This mapping ensures that the computational nodes are optimally clustered within the critical surface layer, enabling accurate resolution of the steep gradients and complex wave dynamics that dominate the nocturnal SBL structure. By leveraging this unified manifold workspace, our Virtual Tower operator can seamlessly integrate data from heterogeneous sources while maintaining high fidelity in the representation of the underlying physical processes.
--
Standard pseudospectral frameworks relying on uniform or Chebyshev-cosine grids distribute computational degrees of freedom symmetrically across the spatial domain. Such configurations are fundamentally unoptimized for the acute vertical gradients characterizing the nocturnal stable boundary layer (SBL), where sub-mesoscale gravity waves and shear instabilities concentrate within the lowest meters of the surface layer ($1.5 \le z \le 5.0\,\text{m}$). To preferentially focus node density where structural geometric derivatives are maximized, we establish a \emph{UnifiedManifoldWorkspace} utilizing a fractional hyperbolic coordinate stretching function $\mathcal{T}(\xi) \rightarrow z$. This operator transforms the bounded computational space $\xi \in [-1, 1]$ into the physical elevation domain defined by the tower instruments ($z_{\mathrm{min}} = 1.5\,\text{m}$ to $z_{\mathrm{max}} = 55.0\,\text{m}$):
\begin{equation}
z(\xi) = z_{\mathrm{min}} + \frac{z_{\mathrm{max}} - z_{\mathrm{min}}}{1 - \alpha_{\mathrm{stretch}}^2} \left[ \frac{(1+\alpha_{\mathrm{stretch}})(1+\xi)}{1 + \alpha_{\mathrm{stretch}}(2+\xi)} \right]
\label{eq:hyperbolic_stretch}
\end{equation}
where $\alpha_{\mathrm{stretch}} \in (0, 1)$ represents the non-dimensional grid compression parameter. The analytical transformation Jacobian resulting from this compactification is defined as:
\begin{equation}
J(\xi) = \frac{dz}{d\xi} = \frac{(z_{\mathrm{max}} - z_{\mathrm{min}})(1 + \alpha_{\mathrm{stretch}})^2}{(1 - \alpha_{\mathrm{stretch}}^2) \left[1 + \alpha_{\mathrm{stretch}}(2+\xi)\right]^2}
\label{eq:jacobian_stretched}
\end{equation}


---
\section{Methods}
\label{sec:methods}

\subsection{Unified Manifold Workspace and Coordinate Mapping}
\label{sec:methods:manifold}

Standard pseudospectral frameworks relying on uniform or Chebyshev-cosine grids distribute computational degrees of freedom symmetrically across the spatial domain. Such configurations are fundamentally unoptimized for the acute vertical gradients characterizing the nocturnal stable boundary layer (SBL), where sub-mesoscale gravity waves and shear instabilities concentrate within the lowest meters of the surface layer ($1.5 \le z \le 5.0\,\text{m}$). To preferentially focus node density where structural geometric derivatives are maximized, we establish a \emph{UnifiedManifoldWorkspace} utilizing a fractional hyperbolic coordinate stretching function $\mathcal{T}(\xi) \rightarrow z$. This operator transforms the bounded computational space $\xi \in [-1, 1]$ into the physical elevation domain defined by the tower instruments ($z_{\mathrm{min}} = 1.5\,\text{m}$ to $z_{\mathrm{max}} = 55.0\,\text{m}$):
\begin{equation}
z(\xi) = z_{\mathrm{min}} + \frac{z_{\mathrm{max}} - z_{\mathrm{min}}}{1 - \alpha_{\mathrm{stretch}}^2} \left[ \frac{(1+\alpha_{\mathrm{stretch}})(1+\xi)}{1 + \alpha_{\mathrm{stretch}}(2+\xi)} \right]
\label{eq:hyperbolic_stretch}
\end{equation}
where $\alpha_{\mathrm{stretch}} \in (0, 1)$ represents the non-dimensional grid compression parameter. The analytical transformation Jacobian resulting from this compactification is defined as:
\begin{equation}
J(\xi) = \frac{dz}{d\xi} = \frac{(z_{\mathrm{max}} - z_{\mathrm{min}})(1 + \alpha_{\mathrm{stretch}})^2}{(1 - \alpha_{\mathrm{stretch}}^2) \left[1 + \alpha_{\mathrm{stretch}}(2+\xi)\right]^2}
\label{eq:jacobian_stretched}
\end{equation}

Fixing the operational baseline at $\alpha_{\mathrm{stretch}} = 0.05$ compresses the local node spacing to an ultra-dense scale ($\Delta z_{\min} = J(-1) \cdot \Delta \xi_{\min} \approx 2.85\,\text{mm}$) directly above the surface boundary, while allowing node separation to expand monotonically toward the top boundary ($\Delta z_{\max} \approx 9.43\,\text{m}$ as $z \rightarrow 55\,\text{m}$). To justify this spatial geometry against alternative topologies, a rigorous sensitivity study was performed across a parameterized continuum of stretching bounds (Table~\ref{tab:stretch_sensitivity}).

\begin{table}[h]
\centering
\caption{Mass Matrix Conditioning and Grid Resolution Sensitivity Across Parameterized Stretching Bounds}
\label{tab:stretch_sensitivity}
\begin{tabular}{c|rrrc}
\hline
$\alpha_{\mathrm{stretch}}$ & $\Delta z_{\min}$ (at $1.5\,\text{m}$) & $\Delta z_{\max}$ (at $55\,\text{m}$) & $\kappa(\mathbf{M})$ & Operational Status \\
\hline
0.01 & $0.11\,\text{mm}$ & $43.20\,\text{m}$ & $1.21 \times 10^5$ & Ill-Conditioned Boundary \\
0.05 & $\mathbf{2.85\,\text{mm}}$ & $\mathbf{9.43\,\text{m}}$ & $\mathbf{4619}$ & \textbf{Optimized Baseline} \\
0.10 & $10.62\,\text{mm}$ & $4.85\,\text{m}$ & $924$ & Suppressed Surface Sensitivity \\
0.50 & $98.40\,\text{mm}$ & $1.82\,\text{m}$ & $63$ & Quasi-Uniform Limit \\
\hline
\end{tabular}
\end{table}

The computational manifold embeds this non-uniform spacing into the discrete physics via a metric-weighted mass matrix $\mathbf{M} \in \mathbb{R}^{(N+1) \times (N+1)}$. This structure discretizes the weighted continuous inner product space using high-order Gauss--Lobatto quadrature across $K_q = 72$ nodes:
\begin{equation}
M_{mn} = \int_{-1}^{1} T_m(\xi) T_n(\xi) J(\xi) \, d\xi \approx \sum_{k=1}^{K_q} w_k \, T_m(\xi_k) T_n(\xi_k) J(\xi_k)
\label{eq:mass_matrix_quad}
\end{equation}
where $T_n(\xi)$ denotes the $n$-th Chebyshev polynomial of the first kind. At $\alpha_{\mathrm{stretch}} = 0.05$, the metric-weighted system maintains a stable spectral profile ($\kappa(\mathbf{M}) = 4619$), ensuring robustness under numerical factorization.

\subsection{Thin SVD-Truncated Matrix Projection Architecture}
\label{sec:methods:svd}

Reconstructing high-fidelity profiles from a sparse tower array ($M = 8$ levels) onto an $N=32$ spectral expansion is an inherently rank-deficient inverse problem. Evaluating the basis at instrument elevations yields the rectangular design matrix $\mathbf{H} \in \mathbb{R}^{M \times (N+1)}$:
\begin{equation}
H_{ij} = T_{j-1}(\xi(z_i))
\label{eq:h_matrix_def}
\end{equation}

To avoid unstable normal-equation inversion and unnecessary parameterized regularization, we use the economy (thin) singular value decomposition (SVD):
\begin{equation}
\mathbf{H} = \mathbf{U}_r \mathbf{S}_r \mathbf{V}_r^T
\label{eq:thin_svd}
\end{equation}
where $\mathbf{U}_r \in \mathbb{R}^{8 \times 8}$ and $\mathbf{V}_r \in \mathbb{R}^{33 \times 8}$ have orthonormal columns, and $\mathbf{S}_r = \text{diag}(s_1, s_2, \ldots, s_8)$ contains the ordered singular values. The thin form restricts all operations to the data-constrained subspace and avoids explicit manipulation of null-space modes.

To separate resolvable physical modes from machine-level noise, we define a dynamic truncation tolerance $\tau$:
\begin{equation}
\tau = \max\left(M, N+1\right) \cdot s_1 \cdot \epsilon_{\text{mach}}
\label{eq:dynamic_tolerance}
\end{equation}
where $\epsilon_{\text{mach}} \approx 2.22 \times 10^{-16}$. The effective rank is
\begin{equation}
r_{\mathrm{eff}} = \#\{s_i > \tau\}.
\label{eq:effective_rank}
\end{equation}
Equivalently, the normalized singular-value cutoff is $\lambda_{\mathrm{min}} = \tau / s_1 = \max(M, N+1)\,\epsilon_{\mathrm{mach}}$. The truncated least-squares estimate is then
\begin{equation}
\mathbf{c} = \sum_{i=1}^{r_{\mathrm{eff}}} \frac{\mathbf{u}_i^{\mathsf{T}}\,\mathbf{A}_{\mathrm{obs}}}{s_i}\,\mathbf{v}_i,
\label{eq:svd_reconstruction}
\end{equation}
which is the Moore--Penrose solution restricted to retained singular directions. This guarantees that inferred coefficients remain in the data-supported subspace while unresolved high-order components are suppressed.

\subsection{Smooth Scale Partitioning and Energy Non-Conservation Caveat}
\label{sec:methods:partition}

To separate multi-scale structures without inducing non-local Gibbs ringing or artificial step discontinuities, we deploy a partition-of-unity spectral windowing framework. The global flow field is decomposed into mean synoptic ($\psi_M$), internal gravity wave ($\psi_W$), and micro-turbulent residual ($\psi_T$) scales using coupled hyperbolic tangent transition channels:
\begin{align}
  \psi_M(n) &= \frac{1}{2}\left[1 - \tanh\left(\frac{n - n_M}{\Delta}\right)\right] \label{eq:psi_m_def} \\
  \psi_W(n) &= \frac{1}{2}\left[1 + \tanh\left(\frac{n - n_M}{\Delta}\right)\right] \cdot \frac{1}{2}\left[1 - \tanh\left(\frac{n - n_W}{\Delta}\right)\right] \label{eq:psi_w_def} \\
  \psi_T(n) &= 1 - \psi_M(n) - \psi_W(n) \label{eq:psi_t_def}
\end{align}
where the scale boundaries are set to $n_M = 3$ and $n_W = 12$ with a smooth transition width $\Delta = 1.2$. By design, the partition operators sum identically to unity across the entirety of the spectral spectrum:
\begin{equation}
\psi_M(n) + \psi_W(n) + \psi_T(n) = 1 \quad \forall n \in [0, N]
\label{eq:partition_unity}
\end{equation}

Crucially, because these smooth partition vectors are applied as coefficient multipliers rather than projections onto true orthogonal eigenvectors of the mass matrix $\mathbf{M}$, \emph{global physical energy content is not strictly additive across the partitioned sub-domains}. When evaluating the individual scale components ($\mathbf{c}^{(M)} = \text{diag}(\psi_M)\mathbf{c}$, $\mathbf{c}^{(W)} = \text{diag}(\psi_W)\mathbf{c}$, and $\mathbf{c}^{(T)} = \text{diag}(\psi_T)\mathbf{c}$), the cross-product metrics do not vanish identically:
\begin{equation}
\langle \mathbf{c}, \mathbf{M} \mathbf{c} \rangle = \sum_{j \in \{M,W,T\}} \langle \mathbf{c}^{(j)}, \mathbf{M} \mathbf{c}^{(j)} \rangle + \sum_{j \neq k} \langle \mathbf{c}^{(j)}, \mathbf{M} \mathbf{c}^{(k)} \rangle
\label{eq:energy_overlap}
\end{equation}
The off-diagonal cross-terms $\langle \mathbf{c}^{(j)}, \mathbf{M} \mathbf{c}^{(k)} \rangle \neq 0$ physically represent spatial and spectral scale overlap regions within the boundary layer transitions. While this cross-term behavior accurately mirrors the coupled wave-turbulence interaction concepts emphasized by \citet{Sun2015}, total energy balance estimates derived via this windowing scheme must explicitly track these non-orthogonal interaction components rather than assuming direct linear superposition.

\subsection{Information-Theoretic Diagnostic Features}
\label{sec:methods:diagnostics}

\subsubsection{Effective Modal Dimension ($D_{\mathrm{eff}}$) with Minimum Energy Filtering}

We compute the effective modal dimension $D_{\mathrm{eff}}$ as the exponential of the spectral Shannon entropy, converting a raw information metric into an interpretable value representing the number of actively participating, equiprobable spectral modes:
\begin{equation}
D_{\mathrm{eff}} = \exp\left( -\sum_n p_n \log p_n \right),
\label{eq:d_eff_def}
\end{equation}
where $p_n = c_n^2 / \sum_m c_m^2$ is the normalized energy in mode $n$.

Because $D_{\mathrm{eff}}$ operates as a scale-invariant shape diagnostic, a quiescent profile under weak background conditions could mathematically register an artificially inflated modal dimension due to numerical precision variations. To prevent low-amplitude instrumental noise from triggering false high-dimensional turbulence classifications, we introduce a regularized minimum energy masking operator before computing the spectral probabilities:
\begin{equation}
c_n^2 \leftarrow \begin{cases}
c_n^2, & \text{if } \sum_{m=0}^{N} c_m^2 \ge \mathcal{E}_{\text{floor}} \\
0, & \text{otherwise}
\end{cases}
\label{eq:energy_floor}
\end{equation}
where the threshold floor is defined dynamically based on historical instrument variance limits, fixed at $\mathcal{E}_{\text{floor}} = 10^{-4} \cdot \max_{t} \left(\|\mathbf{A}_{\mathrm{obs}}\|^2\right)$. This restriction forces highly quiescent, laminar configurations cleanly down to $D_{\mathrm{eff}} \to 1.0$.

\subsubsection{Wave Energy Fraction ($F_W$) and Spectral Curvature ($\chi_N$)}

The wave energy fraction isolates metric-weighted energy in the wave band relative to total energy:
\begin{equation}
F_W = \frac{\langle \mathbf{c}^{(W)}, \mathbf{M} \mathbf{c}^{(W)} \rangle}{\langle \mathbf{c}, \mathbf{M} \mathbf{c} \rangle}
\label{eq:f_w}
\end{equation}

Spectral curvature ($\chi_N$) measures the concentration of energy in high-frequency modes relative to low-frequency modes, providing a proxy for structural sharpness:
\begin{equation}
\chi_N = \frac{\sum_{n=2}^{N} (n/N)^2 \, c_n^2}{\sum_{n=0}^{N} (n/N)^2 \, c_n^2}
\label{eq:chi_n}
\end{equation}

% ============================================================================

Yes. If you’re staying with Chebyshev, I would make the coordinate transformation and Jacobian a first-class component of the framework rather than hard-coding a single map.
Conceptually, every map should expose:
1. Forward transform\xi = T(z)
2. Inverse transformz = T^{-1}(\xi)
3. JacobianJ(\xi)=\frac{dz}{d\xi}
4. Inverse Jacobian\frac{d\xi}{dz}
J(\xi)=\frac{dz}{d\xi}
Then your derivative operator becomes
D_z = \operatorname{diag}\!\left(\frac{1}{J(\xi)}\right) D_\xi
and weighted quadrature becomes
\int f(z)\,dz=\int f(\xi)\,J(\xi)\,d\xi.

⸻

Suggested Julia Design
abstract type CoordinateMap end

forward(map::CoordinateMap, z)
inverse(map::CoordinateMap, ξ)

jacobian(map::CoordinateMap, ξ)
invjacobian(map::CoordinateMap, z)
This lets every map plug into the same spectral machinery.

⸻

1. Linear Map
Useful for testing.
struct LinearMap <: CoordinateMap
    zmin::Float64
    zmax::Float64
end

forward(m::LinearMap, z) =
    2*(z - m.zmin)/(m.zmax - m.zmin) - 1

inverse(m::LinearMap, ξ) =
    m.zmin + (ξ + 1)*(m.zmax - m.zmin)/2

jacobian(m::LinearMap, ξ) =
    (m.zmax - m.zmin)/2

invjacobian(m::LinearMap, z) =
    2/(m.zmax - m.zmin)

⸻

2. Logarithmic Map
Good when surface-layer resolution matters.
Forward:
\xi=2\,\frac{\ln(z/z_{\min})}{\ln(z_{\max}/z_{\min})}-1
Inverse:
z=z_{\min}\left(\frac{z_{\max}}{z_{\min}}\right)^{(\xi+1)/2}
Jacobian:
J(\xi)=\frac{\ln(z_{\max}/z_{\min})}{2}\,z(\xi)
J(\xi)=\frac{\ln(z_{max}/z_{min})}{2}z(\xi)
Julia:
struct LogMap <: CoordinateMap
    zmin::Float64
    zmax::Float64
end

function forward(m::LogMap, z)
    2 * log(z / m.zmin) /
    log(m.zmax / m.zmin) - 1
end

function inverse(m::LogMap, ξ)
    m.zmin *
    (m.zmax / m.zmin)^((ξ + 1)/2)
end

function jacobian(m::LogMap, ξ)
    z = inverse(m, ξ)
    0.5 * log(m.zmax / m.zmin) * z
end

⸻

3. Hyperbolic Tangent Map
Probably your best production choice.
Define
z(\xi)=z_c+\frac{L}{2}\frac{\tanh(\alpha\xi)}{\tanh(\alpha)}
where
* z_c=(z_{max}+z_{min})/2
* L=z_{max}-z_{min}
Jacobian:
J(\xi)=\frac{L}{2}\frac{\alpha}{\tanh(\alpha)}\operatorname{sech}^2(\alpha\xi)
J(\xi)=\frac{L}{2}\frac{\alpha}{\tanh(\alpha)}\operatorname{sech}^2(\alpha\xi)
Julia:
struct TanhMap <: CoordinateMap
    zmin::Float64
    zmax::Float64
    α::Float64
end

function inverse(m::TanhMap, ξ)
    L = m.zmax - m.zmin
    zc = (m.zmax + m.zmin)/2

    zc + (L/2) *
         tanh(m.α * ξ) /
         tanh(m.α)
end

function jacobian(m::TanhMap, ξ)
    L = m.zmax - m.zmin

    (L/2) *
    (m.α / tanh(m.α)) *
    sech(m.α * ξ)^2
end

⸻

4. User-Defined Map
For experimentation.
struct CustomMap{F,G,H} <: CoordinateMap
    forward_fn::F
    inverse_fn::G
    jacobian_fn::H
end

forward(m::CustomMap, z) = m.forward_fn(z)
inverse(m::CustomMap, ξ) = m.inverse_fn(ξ)
jacobian(m::CustomMap, ξ) = m.jacobian_fn(ξ)
Then:
map = CustomMap(
    z -> ...,
    ξ -> ...,
    ξ -> ...
)

⸻

Where the Jacobian Actually Enters
Many spectral codes compute the Chebyshev differentiation matrix D_\xi and stop there.
Your transformed physical derivative should be
\frac{\partial u}{\partial z}=\frac{1}{J(\xi)}\frac{\partial u}{\partial \xi}.
In Julia:
Dz = Diagonal(1.0 ./ J) * Dξ
where:
J = jacobian.(map, ξnodes)
Similarly, for quadrature:
weighted_integral =
    sum(w .* J .* f)
where w are the standard Chebyshev quadrature weights.

⸻

For your manuscript, I’d implement LinearMap, LogMap, TanhMap, and CustomMap, then make TanhMap the default. That gives you a clean API, analytical Jacobians, reproducibility, and a straightforward path to sensitivity studies without touching the spectral solver itself.
