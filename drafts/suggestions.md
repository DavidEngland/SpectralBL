### Section 5: The Mathematics of Coordinate-Agnostic Manifold Interfacing

Building upon the metric-consistent Riemannian foundation established in Section~\ref{sec:methods}, we now formalize the complete algebraic and algorithmic architecture of the generalized \emph{Virtual Tower Interfacing Operator}, $\mathcal{P}_{\mathrm{VT}}$. By shifting the computational paradigm away from localized interpolation schemas, $\mathcal{P}_{\mathrm{VT}}$ serves as a coordinate-free projection mapping from an arbitrary, unstructured observational manifold directly onto the continuous polynomial workspace $\xi \in [-1, 1]$.

This mathematical construction enables identical diagnostic features ($D_{\mathrm{eff}}$, $F_W$, $\chi_N$) to be extracted across disparate datasets—ranging from sparse, discrete tower masts to dense remote sensing arrays (Doppler LiDAR) and staggered numerical meshes (DNS/CFD)—without introducing grid-topology artifacts or scale-dependent biases.

\subsection{Generalized Formulation and the Abstract Design Matrix}
\label{sec:extensions:p_vt_math}

Let $\Omega \subset \mathbb{R}^d$ represent a physical atmospheric or laboratory domain described by an arbitrary, user-defined coordinate system (e.g., Cartesian $(x,y,z)$, spherical radar coordinates $(r, \theta, \phi)$, or terrain-following hydrostatic pressure coordinates $\eta$). We define an arbitrary observational ingestion stream as an unstructured set of $M$ scattered measurements:
\begin{equation}
\mathcal{D}_{\mathrm{obs}} = \left{ \left(\mathbf{x}*i, \Phi_i\right) \right}*{i=1}^M, \quad \mathbf{x}_i \in \Omega
\label{eq:scattered_data_set}
\end{equation}
where $\mathbf{x}_i$ is the physical space position vector of the $i$-th observation, and $\Phi_i$ is the target field scalar (e.g., potential temperature $\theta$, local velocity component $u$, or radial Doppler velocity $v_r$).

The coordinate-agnostic projection pipeline proceeds through a three-stage mathematical transformation:

\paragraph{Stage 1: The Coordinate Normalization Vector Map ($\mathcal{T}$)}
We construct a specialized, geometric transformation operator $\mathcal{T}: \Omega \rightarrow [-1, 1]$ that projects any arbitrary spatial position vector $\mathbf{x}_i$ directly onto the bounded computational coordinate $\xi_i$. The mapping function is explicitly configured to isolate the target vertical structure relative to the instrument's sensing orientation:
\begin{equation}
\xi_i = \mathcal{T}(\mathbf{x}_i)
\label{eq:generalized_tau_map}
\end{equation}

Depending on the underlying physical ingestion target, $\mathcal{T}$ alters its internal coordinate parsing vector:
\begin{itemize}
\item \textbf{Physical Meteorological Towers (1D Vertical Mast):} The position vector reduces to the scalar elevation ($\mathbf{x}_i = [z_i]$), and the mapping applies the analytical inverse of the hyperbolic stretching function defined in Eq.~\ref{eq:inverse_stretch}.
\item \textbf{Scanning Doppler LiDAR (3D Range-Azimuth-Elevation Sweeps):} The position vector is tracked in local spherical coordinates ($\mathbf{x}_i = [r_i, \theta_i, \phi_i]$), where $r_i$ is the range gate, $\theta_i$ is the azimuth angle, and $\phi_i$ is the elevation angle. The transformation computes the true physical height relative to the lens center ($z_i = r_i \sin \phi_i + z_{\mathrm{lidar}}$) and maps $z_i \rightarrow \xi_i$ via Eq.~\ref{eq:inverse_stretch}.
\item \textbf{Numerical Weather Prediction Models (Discrete Hydrostatic Cells):} The position vector tracks the time-varying physical altitude of the model's terrain-following coordinates ($\mathbf{x}_i = [x_{s,i}, y_{s,i}, \eta_i(t)]$), mapping the instantaneous layer elevation $\mathbf{z}(\eta_i) \rightarrow \xi_i$.
\end{itemize}

\paragraph{Stage 2: Construction of the Unstructured Evaluation Operator ($\mathbf{B}$)}
Once the scattered physical coordinates are projected onto the uniform computational domain, we construct the rectangular evaluation design matrix $\mathbf{B} \in \mathbb{R}^{M \times (N+1)}$. Each row of the operator evaluates the chosen high-order Chebyshev basis functions at the mapped coordinates:
\begin{equation}
B_{ij} = T_{j-1}(\xi_i) = T_{j-1}\left(\mathcal{T}(\mathbf{x}_i)\right)
\label{eq:b_matrix_unstructured}
\end{equation}
Crucially, because the rows of $\mathbf{B}$ depend exclusively on the normalized manifold coordinate $\xi_i \in [-1, 1]$, the underlying physical coordinate system of the raw instrument or mesh is entirely abstracted away. The design matrix acts as a universal mathematical interface, rendering all downstream matrix arithmetic completely independent of the original data layout.

\paragraph{Stage 3: Regularized Manifold-Weighted Galerkin Projection}
Reconstructing continuous profiles from scattered data streams requires protecting the spectral solution against numerical overfitting and localized data clustering, which would otherwise trigger unphysical Gibbs oscillations. We define an optimization problem that minimizes the mass-weighted residual while penalizing high-frequency spectral curvature:
\begin{equation}
\mathbf{c} = \arg\min_{\mathbf{\hat{c}}} \left( \left| \mathbf{B}\mathbf{\hat{c}} - \mathbf{\Phi} \right|_{\mathbf{W}}^2 + \gamma , \mathbf{\hat{c}}^{\mathsf{T}} \mathbf{M} \mathbf{\hat{c}} \right)
\label{eq:regularized_least_squares}
\end{equation}
where $\mathbf{\Phi} = [\Phi_1, \Phi_2, \dots, \Phi_M]^{\mathsf{T}}$ is the observation vector, $\mathbf{W} \in \mathbb{R}^{M \times M}$ is a diagonal weighting matrix mapping hardware-specific measurement variances or range-gate attenuation profiles, $\gamma$ is the regularization scalar, and $\mathbf{M}$ is the metric-weighted mass matrix derived via Gauss--Lobatto quadrature (Eq.~\ref{eq:mass_matrix_quad}).

Taking the gradient with respect to $\mathbf{\hat{c}}$ and equating to zero yields the closed-form, coordinate-agnostic algebraic state estimator:
\begin{equation}
\mathbf{c} = \left( \mathbf{B}^{\mathsf{T}} \mathbf{W} \mathbf{B} + \gamma \mathbf{M} \right)^{-1} \mathbf{B}^{\mathsf{T}} \mathbf{W} \mathbf{\Phi}
\label{eq:p_vt_closed_form}
\end{equation}

\subsection{Universal Ingestion Adapter: DNS, CFD, and PIV Laboratory Scaling}
\label{sec:extensions:cfd_ingest}

From a computational fluid dynamics (CFD) and laboratory validation perspective, the coordinate-free abstraction of Eq.~\ref{eq:p_vt_closed_form} allows high-fidelity numerical or experimental duct datasets to be piped straight into the boundary layer diagnostic pipeline. In stably stratified laboratory configurations—such as the Stratified Inclined Duct (SID) experiments reviewed by \citet{Lefauve2025}—flow properties are mapped using dense Particle Image Velocimetry (PIV) pixel matrices or resolved over high-resolution Direct Numerical Simulation (DNS) grids.

Standard comparison workflows require interpolating these precise datasets onto an arbitrary regular grid, introducing numerical diffusion. The $\mathcal{P}_{\mathrm{VT}}$ operator bypasses this step entirely by treating the unstructured DNS mesh vertices or experimental PIV laser interrogation points as direct inputs to the design matrix $\mathbf{B}$.

To accommodate the enclosed symmetric geometry of laboratory channels or pipe flows, the asymmetrical atmospheric stretching function of Eq.~\ref{eq:hyperbolic_stretch} is replaced by a symmetric tangent-hyperbolic ($\tanh$) coordinate mapping function:
\begin{equation}
z(\xi) = \frac{H}{2} \left[ 1 + \frac{\tanh(\delta \xi)}{\tanh(\delta)} \right], \quad \xi \in [-1, 1]
\label{eq:tanh_stretch}
\end{equation}
where $H$ is the total channel height, and $\delta$ is the user-defined wall-clustering parameter. The corresponding analytical Jacobian automatically updates the metric-weighted mass matrix $\mathbf{M}$:
\begin{equation}
J(\xi) = \frac{dz}{d\xi} = \frac{H \delta}{2 \tanh(\delta)} \operatorname{sech}^2(\delta \xi)
\label{eq:jacobian_tanh}
\end{equation}

Because the downstream GMM clustering algorithms and scale-partitioning masks (Eq.~\ref{eq:psi_m_def}--\eq{eq:psi_t_def}) interface exclusively with the abstracted spectral coefficients $\mathbf{c}$, this single modification to the mapping Jacobian allows the software toolkit to transition effortlessly between real-world boundary layer meteorology and controlled laboratory fluid mechanics.

\subsection{Preservation of Mathematical Invariance Across Data Densities}
\label{sec:extensions:invariance}

A critical requirement of geophysical state estimation is that the extracted information-theoretic signatures must remain invariant to the sampling density $M$. If the calculated effective modal dimension ($D_{\mathrm{eff}}$) scaled artifactually with the number of operational sensors, it would fail as an objective physical diagnostic tool.

The $\mathcal{P}_{\mathrm{VT}}$ framework guarantees this mathematical invariance across scale transitions because the Singular Value Decomposition (SVD) projection step isolates the dominant singular spectrum ($s_1, s_2, \dots, s_8$) relative to a dynamic noise floor $\tau$ (Eq.~\ref{eq:dynamic_tolerance}). Increasing the spatial sampling density from $M = 8$ (a standard atmospheric tower mast) to $M = 1000$ (a vertically pointing Doppler LiDAR beam or a dense DNS profile line) places tighter, highly overdetermined mathematical constraints on the underlying Chebyshev coefficients without changing the physical volume integration defined by the mass matrix $\mathbf{M}$.

Consequently, if the underlying physical boundary layer collapses into a strongly stratified, gravity-wave dominated state, the energy distribution within the spectral coefficient vector $\mathbf{c}$ remains tightly constrained to the lowest-order modes. The Shannon entropy calculation (Eq.~\ref{eq:d_eff_def}) will converge cleanly to an identical, low-rank value ($D_{\mathrm{eff}} \sim 4$--$6$), completely independent of whether the field was captured by a handful of physical instruments or an intensive laser array. This topological invariance proves that the extracted features isolate true physical properties of the stratified flow, establishing a robust framework for multi-instrument geophysical validation.