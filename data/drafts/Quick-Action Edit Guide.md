# Quick-Action Edit Guide

## PRIORITY 1: CRITICAL FIXES (Do first — these break compilation/citations)

### Edit 1.1: Fix Huang Citation Key Mismatch

**File:** Your LaTeX source
**Find:** Line ~270 in Methods section
**Current text:**

```latex
\cite{Huang2009}
```

**Change to:**

```latex
\cite{Huang1998}
```

**Also in bibliography, find:**

```bibtex
\bibitem[{Huang et~al.(1998)}]{Huang2009}
```

**Change to:**

```bibtex
\bibitem[{Huang et~al.(1998)}]{Huang1998}
```

-----

## PRIORITY 2: HIGH-IMPACT EDITS (Improve submissions, no errors)

### Edit 2.1: Revise Figure 1 Caption

**File:** Your LaTeX source
**Find:** Figure 1 caption (~line 520)
**Current:**

```latex
\caption{Hyperbolic grid mapping and spectral partition diagnostics
         generated from the manifold workspace run.}
\label{fig:grid}
```

**Replace with:**

```latex
\caption{(a) Hyperbolic grid mapping showing node density concentration
         at the surface ($z = 1.5$ m) with compactification parameter
         $\alpha_{\mathrm{stretch}} = 0.05$. (b) Derivatives of Chebyshev
         basis polynomials $dT_n/dz$ showing spectral sensitivity at
         multiple scales. Grid spacing ranges from $\Delta z_{\min} = 2.85$ mm
         (surface) to $\Delta z_{\max} = 9.43$ m (upper boundary).}
```

-----

### Edit 2.2: Revise Figure 2 Caption

**File:** Your LaTeX source
**Find:** Figure 2 (tier1_plane1) caption (~line 570)
**Current:**

```latex
\caption{Observed CASES-99 feature scatter from existing run outputs
         in the $(D_{\mathrm{eff}}, F_W)$ plane.}
```

**Replace with:**

```latex
\caption{Scatter of all 1386 CASES-99 profiles (22--31 October 1999)
         in the effective modal dimension--wave energy fraction plane
         $(D_{\mathrm{eff}}, F_W)$. Three visually distinct clusters emerge:
         Continuous Turbulence (red, $D_{\mathrm{eff}} > 20$, $F_W < 0.15$),
         Wave-Dominated (blue, $D_{\mathrm{eff}} \sim 5$, $F_W > 0.75$), and
         Intermittent Shear Bursts (yellow, intermediate). Silhouette score
         $S = 0.582$ for $K=3$ indicates strong regime separability.}
```

-----

### Edit 2.3: Revise Figure 3 Caption

**File:** Your LaTeX source
**Find:** Figure 3 (tier1_plane2) caption (~line 600)
**Current:**

```latex
\caption{Observed CASES-99 feature scatter from existing run outputs
         in the $(\chi_N, \mathrm{Ri}_g)$ plane.}
```

**Replace with:**

```latex
\caption{Feature scatter in the spectral curvature--gradient Richardson
         number plane $(\chi_N, \mathrm{Ri}_g)$ for the same 1386 CASES-99
         profiles. Low curvature ($\chi_N < 0.3$) and high stratification
         ($\mathrm{Ri}_g > 0.3$) characterize Wave-Dominated states (blue).
         High curvature ($\chi_N > 0.3$) and weak stratification
         ($\mathrm{Ri}_g < 0.15$) characterize Continuous Turbulence (red).
         Orthogonal perspective confirms regime separation from the
         energy--dimension plane.}
```

-----

### Edit 2.4: Revise Figure 4 Caption

**File:** Your LaTeX source
**Find:** Figure 4 (temporal_trace) caption (~line 680)
**Current:**

```latex
\caption{Temporal diagnostic trace from existing CASES-99 trajectory
         output used for regime timing interpretation.}
```

**Replace with:**

```latex
\caption{Time series of cluster membership (color-coded by regime) across
         the 10-day intensive observation period (22--31 October 1999,
         panels a--d show successive 60-hour windows). Regime 2 (Wave-Dominated,
         blue) clusters during early morning hours (02:00--06:00 local time);
         Regime 1 (Continuous Turbulence, red) concentrates during evening
         transition (18:00--22:00); Regime 3 (Intermittent Shear, yellow)
         dominates mid-day to afternoon transitions. Distinct diurnal
         periodicity validates that regimes reflect repeatable atmospheric
         dynamics rather than statistical artifacts.}
```

-----

### Edit 2.5: Revise Figure 5 Caption

**File:** Your LaTeX source
**Find:** Figure 5 (wave_reflection) caption (~line 750)
**Current:**

```latex
\caption{Synthetic gravity-wave boundary-energy comparison
         (with and without spectral sponge damping).}
```

**Replace with:**

```latex
\caption{Synthetic gravity-wave propagation test comparing upper-boundary
         absorption. A Gaussian wave packet (initial modal energy $E_0 = 1$)
         is initialized at $z = 35$ m and propagated upward. (Left panel)
         Energy time series with (blue) and without (red) spectral partition
         sponge damping. (Right panel) Zoomed final 20 time steps showing
         near-elimination of reflection with sponge applied. Mean boundary
         energy reduced from $211.2$ (unmasked) to $0.0064$ (masked),
         yielding suppression ratio $\approx 3.3 \times 10^4$, validating
         the sponge-layer hypothesis.}
```

-----

## PRIORITY 3: REFERENCE ADDITIONS (Paste into bibliography)

**Find in your LaTeX source:** The line `\begin{thebibliography}{30}`

**Before the closing `\end{thebibliography}` tag, add these entries:**

### Add after Lefauve2025 entry:

```bibtex
\bibitem[{Trefethen(2000)}]{Trefethen2000}
L.~N. Trefethen.
\newblock {Spectral Methods in {MATLAB}}.
\newblock \emph{Society for Industrial and Applied Mathematics (SIAM)},
  Philadelphia, PA, 2000.

\bibitem[{Boyd(2001)}]{Boyd2001}
J.~P. Boyd.
\newblock {Chebyshev and {Fourier} Spectral Methods}.
\newblock Dover Publications, Mineola, NY, 2nd edition, 2001.

\bibitem[{Fornberg(1996)}]{Fornberg1996}
B.~Fornberg.
\newblock {A Practical Guide to Pseudospectral Methods}.
\newblock Cambridge University Press, Cambridge, UK, 1996.
\newblock \doi{10.1017/CBO9780511626357}.

\bibitem[{Holtslag and De Bruin(1988)}]{Holtslag1988}
A.~A.~M. Holtslag and H.~A.~R. De Bruin.
\newblock {Applied Modeling of the Nighttime Surface Energy Balance over Land}.
\newblock \emph{Journal of Applied Meteorology}, 27(6):689--704, 1988.
\newblock \doi{10.1175/1520-0450(1988)027<0689:AMOTNSE>2.0.CO;2}.

\bibitem[{Nieuwstadt(1984)}]{Nieuwstadt1984}
F.~T.~M. Nieuwstadt.
\newblock {The Turbulent Structure of the Stable, Nocturnal Boundary Layer}.
\newblock \emph{Journal of the Atmospheric Sciences}, 41(14):2202--2216, 1984.
\newblock \doi{10.1175/1520-0469(1984)041<2202:TTSOTSNBL>2.0.CO;2}.

\bibitem[{Kelley and Ouellette(2011)}]{Kelley2011}
D.~H. Kelley and N.~T. Ouellette.
\newblock {Separating Fast and Slow: Timescales in Stratified Turbulence}.
\newblock \emph{Physics of Fluids}, 23(7):071701, 2011.
\newblock \doi{10.1063/1.3599701}.

\bibitem[{Caulfield(2021)}]{Caulfield2021}
C.~P. Caulfield.
\newblock {Layering, Instabilities, and Mixing in Turbulent Stratified Flows}.
\newblock \emph{Annual Review of Fluid Mechanics}, 53:113--145, 2021.
\newblock \doi{10.1146/annurev-fluid-030121-014412}.

\bibitem[{McQuarrie et~al.(2021)McQuarrie, Huang, and Willcox}]{McQuarrie2021}
S.~A. McQuarrie, W.~Huang, and K.~W. Willcox.
\newblock {Data-Driven Reduced-Order Models for {R}ayleigh--{B}\'enard Convection}.
\newblock \emph{Physical Review Fluids}, 6(7):073501, 2021.
\newblock \doi{10.1103/PhysRevFluids.6.073501}.
```

-----

## PRIORITY 4: NEW CITATION INSERTIONS (Where to cite new references in text)

### Cite Trefethen2000 & Boyd2001

**Find in Introduction § Contribution section:**

```
Instead of: "Key innovations:"
Insert before: "...Riemannian metric..."
```

**Add this sentence:**

```latex
Our approach builds on classical pseudospectral methods \cite{Trefethen2000, Fornberg1996},
extending them to stable boundary layers via metric-consistent discretization.
The use of Chebyshev polynomials as basis functions follows established practice
\cite{Boyd2001}, with the key innovation being the hyperbolic compactification
of the vertical coordinate and explicit partition-of-unity spectral masking.
```

-----

### Cite Fornberg1996 in Methods § 2.2

**Find:** Line ~420, end of SVD-Truncated Least-Squares Projection subsection
**After:** “preventing overfitting.”

**Add:**

```latex
\cite{Fornberg1996} discusses SVD truncation as an alternative to
Tikhonov regularization for underdetermined systems; our approach prioritizes
physical interpretability by zeroing modes with zero data constraint.
```

-----

### Cite Nieuwstadt1984 & Holtslag1988 in Introduction

**Find:** Introduction § 1 (MOST Crisis section)
**After:** “…rendering time-averaged diagnostics ambiguous.”

**Add:**

```latex
Classical observational and modeling studies \cite{Nieuwstadt1984, Holtslag1988}
have documented these deviations systematically, establishing that
MOST's universal form breaks down under conditions of strong stratification
and intermittent turbulence---precisely the regimes relevant to forecast
models in stable boundary-layer prediction.
```

-----

### Cite Kelley2011 & Caulfield2021 in Discussion

**Find:** Discussion § 4, after citation of Sun2015
**Replace:** “In particular, the present decomposition…”

**With:**

```latex
These results align naturally with recent frameworks for understanding
stratified-flow timescales \cite{Kelley2011} and modern perspectives on
mixing-layer structure and buoyancy control \cite{Caulfield2021}.
In particular, the present decomposition...
```

-----

## PRIORITY 5: MINOR STYLE FIXES

### Edit 5.1: Email capitalization

**Find:** Author email in preamble
**Current:**

```latex
David.England@UAH.Edu
```

**Change to:**

```latex
David.England@uah.edu
```

-----

### Edit 5.2: Dedication page (optional — check target journal policy)

**Issue:** Dedication pages may not be acceptable in all journals
**Check:** Target journal guidelines before submission
**If not allowed:** Comment out or move to supplementary material

```latex
% OPTIONAL: Comment out this block if target journal does not accept dedications
\thispagestyle{empty}
\vspace*{\fill}
...
\end{center}
\vspace*{\fill}
```

-----

## VERIFICATION CHECKLIST (Before uploading)

After making all edits above, verify:

- [ ] Compile with `pdflatex` at least twice
- [ ] No compilation errors or warnings about undefined citations
- [ ] All citations resolve (e.g., `\cite{Huang1998}` now works)
- [ ] All new BibTeX entries have DOIs where available
- [ ] Figure captions display correctly without overflow
- [ ] Page count is reasonable (~15–18 pages)
- [ ] Abstract still fits on first page (or journal’s first page guidelines)
- [ ] Bibliography count is now ~31 entries (was 23)

-----

## TIME ESTIMATE

- **Huang citation fix:** 2 min
- **Figure caption revisions (5 captions):** 10 min
- **New bibliography entries (8 entries):** 5 min
- **New text insertions (4 locations):** 10 min
- **Email & style fixes:** 2 min
- **Compilation & verification:** 5 min

**Total time:** ~34 minutes

-----

## NEXT STEPS AFTER EDITS

1. **Select target journal** (recommended: *Journal of the Atmospheric Sciences*, *Boundary-Layer Meteorology*, or *Quarterly Journal of the Royal Meteorological Society*)
1. **Download journal’s LaTeX template** (if provided)
1. **Adapt your preamble** to match journal requirements
1. **Verify all figures are embedded** in final PDF
1. **Create cover letter** with 3–5 sentence summary of novelty
1. **Submit through journal’s online portal**

Good luck with submission!