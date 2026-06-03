The **Cooperative Atmosphere-Surface Exchange Study–1999 (CASES-99)** was a comprehensive six-week field campaign conducted from late September to November 1999, centered near **Purcell, Oklahoma**, and **Leon, Kansas**. The project was designed as an intensive investigation of the **nighttime stable boundary layer (SBL)**, specifically focusing on the linkages between the atmosphere and the Earth's surface and the physical processes associated with evening and morning transition regimes.

### Objectives and Instrumentation
The campaign aimed to provide a detailed time history of **internal gravity waves (IGWs)**, Kelvin-Helmholtz shear instabilities, and turbulence events to evaluate their relative contributions to heat, moisture, and momentum fluxes. To achieve this, a massive array of **in situ and remote sensors** was deployed:
*   **Tower Complex:** A central **55–60 m main tower** was surrounded by two sets of three 10 m towers at radii of 100 m and 300 m. These towers were equipped with **CSAT3 sonic anemometers** and thermocouples at multiple levels to measure 3D wind and temperature fluctuations.
*   **Mobile and Remote Platforms:** The campaign utilized aircraft platforms, **microbarographs**, a **Tethered Lifting System (TLS)** involving kites and balloons, and various remote sensing tools including **lidars, radars, sodars, and scintillometers**.
*   **Data Resolution:** High-rate sensors like 3D sonic anemometers provided data resampled to **20 Hz**, while slower sensors for radiation and humidity were sampled at **1 Hz**.

### Significant Atmospheric Phenomena
Research independent of England’s manifold analysis identified several critical features of the nocturnal SBL during the campaign:
*   **Up-Gully Flow Surges:** Observations showed that when near-surface wind directions shifted through a specific "up-gully" orientation relative to a nearby terrain feature, it produced a **localized upward surge of vertical velocity** reaching at least 55 meters. This surge generated **outward-propagating gravity wave packets** with wavelengths between 100 and 250 meters.
*   **"Dirty" Waves:** During the night of October 20, 1999, researchers documented complex temperature and wind speed oscillations, referred to as **"dirty" waves**, which occurred very close to the surface (0.5 to 2.3 meters).
*   **MOST Departures:** One of the campaign's central findings was the confirmation that standard **Monin–Obukhov Similarity Theory (MOST)** fails under very stable conditions, as it cannot account for the rapid switches and hysteresis caused by intermittent turbulence and wave activity.
*   **Intermittent Turbulence:** Data from Oct 22–31 highlighted periods where sudden **turbulent bursts** triggered rapid changes in the Obukhov length, making traditional time-averaged diagnostics ambiguous.

### Site Characteristics and Environmental Context
The CASES-99 site was characterized by **gently rolling grassland** with an aerodynamic roughness ($z_0$) of approximately 0.01 meters. The **intensive observation period (IOP)** from October 22–31 was specifically marked by **strong radiative cooling**, high-amplitude stable inversions, and consistent dry conditions with a dew-point depression exceeding 15K. These conditions allowed fine-scale stable dynamics and coherent gravity wave activity to dominate the boundary layer, providing a "natural laboratory" for SBL research.

The summary provided captures the physical realities of the CASES-99 campaign that directly motivate the mathematical architecture of your pseudospectral pipeline. By contextualizing the site's environmental conditions, instruments, and atmospheric phenomena, we can trace exactly why your code requires things like hyperbolic stretching, multi-brand quality flag processing, and wave reflection suppression.

The physical details of the campaign map directly to the design choices in your Julia software pipeline:

### 1. Grid Dense Lower Resolution vs. The Tower Setup

* **The Physical Reality:** CASES-99 observed highly localized phenomena like "dirty" waves acting extremely close to the surface ($0.5$ to $2.3\text{ m}$) and a roughness length ($z_0$) of $0.01\text{ m}$. Standard Monin-Obukhov Similarity Theory (MOST) fails here because these thin, near-surface layers decouple from the air above.
* **The Code Implementation:** This is why your `UnifiedManifoldWorkspace` uses a fractional hyperbolic stretching parameter ($\alpha = 0.05$). Standard linear or cosine Chebyshev grids would space nodes too widely near the ground. By compressing the grid downward, your pipeline places sub-centimeter node resolution right between $1.5\text{ m}$ and $5\text{ m}$. This lets the pseudospectral matrices resolve the intense thermal gradients of those "dirty" waves without needing hundreds of computational modes.

### 2. Up-Gully Flow Surges and Spectral Partitioning

* **The Physical Reality:** The terrain-induced "up-gully" flow surges triggered upward vertical velocities reaching past the $55\text{ m}$ tower height, generating outward-propagating internal gravity wave (IGW) packets.
* **The Code Implementation:** In your `scripts/Report.jl`, you partitioned the Chebyshev modes into Meso ($\psi_M$), Wave ($\psi_W$), and Turbulent ($\psi_T$) windows. The Wave window (centered near mode 7) is engineered to isolate these exact non-turbulent IGW packets. Furthermore, because these waves propagate upward, your code treats the expanding grid spaces near the $55\text{ m}$ boundary as a natural **sponge layer** paired with the high-frequency filter ($\psi_T$) to absorb wave energy and prevent unphysical numerical wave reflections back down the tower.

### 3. Intermittent Turbulence, Hysteresis, and Data Filtering

* **The Physical Reality:** The Intensive Observation Period (IOP) from October 22–31 was dominated by strong radiative cooling and sudden turbulent bursts. These bursts create rapid switches and hysteresis, causing data drops or instrument ringing on the towers.
* **The Code Implementation:** This explains why generic NetCDF parsers fail and why you needed `clean_by_ncar_quality_flags` inside your pipeline.
* The Campbell CSAT3 sonics at the $1.5\text{ m}$, $5\text{ m}$, $30\text{ m}$, and $50\text{ m}$ levels required bitfield checks (`diag_Zm != 0`) to catch transducer blockage from dew or frost caused by that intense radiative cooling.
* The ATI sonics at $10\text{ m}$, $20\text{ m}$, $40\text{ m}$, and $55\text{ m}$ required checking the sub-sampling counters (`usamples_Zm < 10`) because rapid turbulent bursts could cause intermittent data dropouts within a single 20 Hz frame sequence.

### 4. Continuous Vectoring vs. Overwriting Reports

* **The Physical Reality:** CASES-99 ran continuously for weeks, meaning phenomena like the Richardson number ($Ri_f$) and effective dimension ($D_{\text{eff}}$) evolve over long, non-stationary historical timelines.
* **The Code Implementation:** Because your pipeline loops over every single day sequentially or in parallel (`make run-all-parallel`), storing data with a generic timestamp (like `TimeIdx` from 1 to 288) would cause separate days to overwrite each other. By adding the `FileDate` column to your master trajectory CSV and appending a `day_suffix` to your markdown summaries and Plots.jl figures, your software architecture treats each day of the campaign as an independent, non-stationary slice of a larger continuous manifold.
