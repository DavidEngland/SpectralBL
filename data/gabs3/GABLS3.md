GEWEX Atmospheric Boundary-Layer Study (GABLS): 3rd Case Study Briefing

Executive Summary

The GEWEX Atmospheric Boundary-Layer Study (GABLS) focuses on improving the representation of stable boundary layers in atmospheric models through global intercomparisons. While the first two GABLS cases utilized idealized conditions and prescribed surface temperatures, the 3rd case study utilizes a real-world dataset from the Cabauw meteorological site in the Netherlands. This case, spanning from July 1 to July 2, 2006, aims to evaluate Single Column Models (SCMs) against observations of decoupling, inertial oscillation, low-level jets, and the morning transition to convective conditions. Participants are required to run SCMs with full soil-vegetation and radiation modules, adhering to strictly prescribed initial conditions and forcing tendencies to ensure a valid intercomparison.

Project Overview and Objectives

GABLS provides a platform for the atmospheric boundary layer research community to benchmark Single Column Models, including research-grade models and those derived from operational weather and climate models.

* Previous Studies:
  * Case 1: An idealized study over snow with prescribed surface temperatures.
  * Case 2: Based on the CASES 99 experiment, also featuring prescribed surface temperatures.
* The 3rd Case Rationale: Earlier studies demonstrated that a lack of interaction with the surface and the simplicity of boundary conditions made it difficult to compare results with real-world observations. The Cabauw site was selected for its flat topography and high-quality observational data, making it ideal for studying the full diurnal cycle of the boundary layer.
* Modeling Constraints: Participants must use their full SCM suite. Operational institutes are requested to use their standard vertical resolution. Crucially, participants are instructed not to tune their models specifically to the Cabauw case.

Meteorological Case Context

The study is centered on a 24-hour period at the Cabauw site (51.9711° N, 4.9267° E, -0.7 m ASL).

Synoptic Situation

The simulation begins at 12:00 UTC on July 1, 2006, and concludes at 12:00 UTC on July 2, 2006.

* Pressure Systems: A stationary high-pressure system resides over Scandinavia, north of Cabauw.
* Wind Conditions: A relatively constant Easterly geostrophic wind of approximately 8 m/s exists. The high-pressure axis is tilted south, causing geostrophic wind to decrease with height.
* Atmospheric Disturbances: A dry air mass passes over the site at the start of the period. Between 00:00 and 03:00 UTC, a small synoptical disturbance is advected over the site, resulting in variations in temperature, humidity, and wind.

Boundary Layer Development

* Initial State: The case begins with a convective boundary layer height of approximately 1400 m.
* Nocturnal Transition: After sunset, a turbulent stable boundary layer develops. Decoupling and inertial oscillation lead to the formation of a low-level jet at a height of 200 m.
* Morning Transition: Sunrise initiates convective mixing, resulting in a fast-growing convective boundary layer that reaches approximately 1900 m by noon on the second day.

Prescribed Model Forcings and Parameters

Surface and Vegetation Parameters

To ensure consistency across models, specific surface parameters are prescribed:

* Albedo: 0.23 (consistent across NIR and VIS bands).
* Emissivity: 0.99.
* Vegetation: 100% grass with a Leaf Area Index (LAI) of 2.
* Soil Composition: Prescribed as clay (45% clay, 8% organic matter, 0% sand).
* Field Capacity: 0.47 m^3/m^3.

Roughness Lengths

The study distinguishes between momentum and heat roughness lengths based on regional obstacles (e.g., tree rows) versus local grassland:

* Momentum (z_{0m}): Regional roughness is set at 0.15 m (derived from wind gusts), which is significantly higher than the local grassland roughness of 0.03 m.
* Heat (z_{0h}): Set at 0.0015 m. This reflects a 10:1 ratio compared to local grassland roughness, as pressure drag from obstacles does not affect temperature transport in the same way it affects momentum.

Energy Budget and Evaporation

Evaporation from the grassland is the dominant factor in the surface energy budget.

* Sensible Heat Flux: Typically 100 W/m^2 during the day.
* Bowen Ratio: Modellers must adjust soil water content to achieve a Bowen ratio of 0.33 at the start of the simulation.

Initial Conditions and Profiles

Atmospheric Profiles (12:00 UTC, July 1, 2006)

Initial values are derived from De Bilt radiosoundings and Cabauw tower observations. To avoid model errors associated with the initial dry air mass advection, the humidity content was manually increased based on surrounding soundings to prevent an unrealistically dry boundary layer.

Height (m)	Temperature (°C)	Specific Humidity (kg/kg)
0	27.0	9.3 E-3
10	26.4	8.5 E-3
205	24.3	8.0 E-3
1800	9.0	7.5 E-3
2200	9.0	2.0 E-3
12000	-61.0	0.01 E-3

Soil Temperature Profile

The initial soil temperature profile is based on observed mean values:

Depth (m)	Temperature (°C)
0.00	23.4
0.04	21.4
0.20	17.9
1.00	12.2
2.00	10.0
inf	10.0

Dynamic Tendencies (Advection)

Because SCMs cannot naturally simulate horizontal disturbances, advective tendencies are prescribed based on 3D NWP model hindcasts and observations.

* Vertical Movement (Omega): Prescribed constant between 1500–5000 m (0.12 Pa/s until 17:00 UTC, then 0.00 Pa/s).
* Horizontal Tendencies: Prescribed constant between 200–1000 m. This includes specific schedules for temperature dynamic tendency (T_{adv}), specific humidity dynamic tendency (q_{adv}), and horizontal wind dynamic tendency (U_{adv}, V_{adv}).

Administrative Information

Key Project Personnel

* Fred Bosveld: KNMI (bosveld@knmi.nl)
* Cisco de Bruijn: KNMI (cisco.de.bruijn@knmi.nl)
* Bert Holtslag: Wageningen University and Research

Project Timeline

* February 1, 2008: Release of case set-up.
* May 1, 2008: Deadline for results to be included in preliminary analysis.
* November 15, 2008: Final deadline for sending results.
* June 26–27, 2009: GABLS Workshop in Boulder, Colorado, USA.
