# GABLS3 Flux Predictive Baseline

This report fits (1) a global weighted linear baseline and (2) a regime-conditioned sub-model (one fit per WindowTag). Lagged inertial-oscillation features (ΔU/Δt, ΔV/Δt, InertialMagDt) included. Sample weights are inversely proportional to window sample count.

## Window Counts

| WindowTag | n_samples |
|---|---:|
| OTHER | 132 |
| HOUR8 | 6 |
| HOUR9 | 6 |

## Target: MomUFluxTruth

### Global Weighted Linear Baseline

- wRMSE: 0.019601
- wMAE: 0.010796
- wR2: 0.585372
- wBias: -0.000000

#### By-Window Metrics (Global Model)

| WindowTag | wRMSE | wMAE | wR2 | wBias | n_samples |
|---|---:|---:|---:|---:|---:|
| OTHER | 0.032458 | 0.020029 | 0.494077 | -0.002101 | 131 |
| HOUR8 | 0.009273 | 0.008569 | 0.175063 | 0.004121 | 6 |
| HOUR9 | 0.004265 | 0.003860 | -19.392672 | -0.002036 | 6 |

#### Coefficients (Global Model)

| Feature | Coefficient |
|---|---:|
| Intercept | -18.37818978 |
| D_eff | 5.27614341 |
| F_W | 0.24957178 |
| chi_N | -36.56380969 |
| Ri_g | -0.01093707 |
| Ri_b | -0.00000735 |
| UShearMean | -0.08747052 |
| VShearMean | -0.46561010 |
| Deff_x_ShearMag | -33.95837928 |
| ShearMagMean | 124.74217711 |
| DeltaU_Dt | -0.00283516 |
| DeltaV_Dt | -0.00050222 |
| InertialMagDt | -0.00397990 |

### Regime-Conditioned Sub-Model

#### By-Window Metrics (Regime Model)

| WindowTag | wRMSE | wMAE | wR2 | wBias | n_samples |
|---|---:|---:|---:|---:|---:|
| OTHER | 0.031404 | 0.021043 | 0.526387 | -0.000000 | 131 |
| HOUR8 | 0.003463 | 0.002844 | 0.884974 | -0.000000 | 6 |
| HOUR9 | 0.000683 | 0.000624 | 0.477309 | -0.000000 | 6 |

#### Per-Tag Coefficients (Regime Model)

**WindowTag = HOUR8** *(restricted: N < 15, features: Intercept, D_eff, ShearMagMean)*

| Feature | Coefficient |
|---|---:|
| Intercept | 35.82180254 |
| D_eff | -9.82686628 |
| ShearMagMean | 2.77622544 |

**WindowTag = HOUR9** *(restricted: N < 15, features: Intercept, D_eff, ShearMagMean)*

| Feature | Coefficient |
|---|---:|
| Intercept | 4.06707012 |
| D_eff | -1.11658584 |
| ShearMagMean | 0.20286233 |

**WindowTag = OTHER**

| Feature | Coefficient |
|---|---:|
| Intercept | -6.47812950 |
| D_eff | 1.89854298 |
| F_W | -0.09500269 |
| chi_N | -19.61830615 |
| Ri_g | 0.00563733 |
| Ri_b | -0.00000940 |
| UShearMean | -1.00169261 |
| VShearMean | -0.81072929 |
| Deff_x_ShearMag | 6.62660415 |
| ShearMagMean | -22.45787336 |
| DeltaU_Dt | -0.00139666 |
| DeltaV_Dt | -0.00088775 |
| InertialMagDt | -0.00408844 |


## Target: HeatFluxTruth

### Global Weighted Linear Baseline

- wRMSE: 0.013257
- wMAE: 0.008412
- wR2: 0.816470
- wBias: 0.000000

#### By-Window Metrics (Global Model)

| WindowTag | wRMSE | wMAE | wR2 | wBias | n_samples |
|---|---:|---:|---:|---:|---:|
| OTHER | 0.021925 | 0.017068 | 0.757546 | 0.000308 | 131 |
| HOUR8 | 0.005885 | 0.004707 | -13.619964 | -0.002093 | 6 |
| HOUR9 | 0.003768 | 0.003527 | -141.963861 | 0.001787 | 6 |

#### Coefficients (Global Model)

| Feature | Coefficient |
|---|---:|
| Intercept | 31.49592880 |
| D_eff | -8.33986703 |
| F_W | -0.19472959 |
| chi_N | -30.21481955 |
| Ri_g | -0.01297616 |
| Ri_b | 0.00000154 |
| UShearMean | 3.59125035 |
| VShearMean | 2.07176798 |
| Deff_x_ShearMag | 50.24241572 |
| ShearMagMean | -187.21212807 |
| DeltaU_Dt | 0.00216426 |
| DeltaV_Dt | -0.00013480 |
| InertialMagDt | 0.00135560 |

### Regime-Conditioned Sub-Model

#### By-Window Metrics (Regime Model)

| WindowTag | wRMSE | wMAE | wR2 | wBias | n_samples |
|---|---:|---:|---:|---:|---:|
| OTHER | 0.021103 | 0.016525 | 0.775386 | 0.000000 | 131 |
| HOUR8 | 0.000939 | 0.000803 | 0.627787 | -0.000000 | 6 |
| HOUR9 | 0.000127 | 0.000114 | 0.836665 | -0.000000 | 6 |

#### Per-Tag Coefficients (Regime Model)

**WindowTag = HOUR8** *(restricted: N < 15, features: Intercept, D_eff, ShearMagMean)*

| Feature | Coefficient |
|---|---:|
| Intercept | 6.58923443 |
| D_eff | -1.80742892 |
| ShearMagMean | 0.44508145 |

**WindowTag = HOUR9** *(restricted: N < 15, features: Intercept, D_eff, ShearMagMean)*

| Feature | Coefficient |
|---|---:|
| Intercept | 2.01341700 |
| D_eff | -0.55374471 |
| ShearMagMean | 0.13777265 |

**WindowTag = OTHER**

| Feature | Coefficient |
|---|---:|
| Intercept | 22.33907207 |
| D_eff | -5.81418302 |
| F_W | -0.08916573 |
| chi_N | -33.94134622 |
| Ri_g | -0.04341073 |
| Ri_b | -0.00000102 |
| UShearMean | 6.23009938 |
| VShearMean | 2.96182111 |
| Deff_x_ShearMag | 27.03005121 |
| ShearMagMean | -105.04040157 |
| DeltaU_Dt | 0.00121332 |
| DeltaV_Dt | 0.00021651 |
| InertialMagDt | 0.00222133 |


## Target: MomVFluxTruth

### Global Weighted Linear Baseline

- wRMSE: 0.009354
- wMAE: 0.004113
- wR2: 0.104119
- wBias: -0.000000

#### By-Window Metrics (Global Model)

| WindowTag | wRMSE | wMAE | wR2 | wBias | n_samples |
|---|---:|---:|---:|---:|---:|
| OTHER | 0.016055 | 0.009288 | 0.121324 | 0.000325 | 131 |
| HOUR8 | 0.001552 | 0.001245 | -1.954584 | -0.000201 | 6 |
| HOUR9 | 0.001899 | 0.001845 | -17.834713 | -0.000121 | 6 |

#### Coefficients (Global Model)

| Feature | Coefficient |
|---|---:|
| Intercept | 3.30603232 |
| D_eff | -0.93294743 |
| F_W | -0.00624530 |
| chi_N | 4.99854187 |
| Ri_g | -0.00921372 |
| Ri_b | 0.00000076 |
| UShearMean | -0.38343986 |
| VShearMean | -0.43418636 |
| Deff_x_ShearMag | 14.06604832 |
| ShearMagMean | -51.46563736 |
| DeltaU_Dt | -0.00009479 |
| DeltaV_Dt | 0.00051630 |
| InertialMagDt | -0.00020871 |

### Regime-Conditioned Sub-Model

#### By-Window Metrics (Regime Model)

| WindowTag | wRMSE | wMAE | wR2 | wBias | n_samples |
|---|---:|---:|---:|---:|---:|
| OTHER | 0.015777 | 0.009979 | 0.151473 | -0.000000 | 131 |
| HOUR8 | 0.000347 | 0.000252 | 0.851827 | -0.000000 | 6 |
| HOUR9 | 0.000371 | 0.000319 | 0.282087 | -0.000000 | 6 |

#### Per-Tag Coefficients (Regime Model)

**WindowTag = HOUR8** *(restricted: N < 15, features: Intercept, D_eff, ShearMagMean)*

| Feature | Coefficient |
|---|---:|
| Intercept | 2.82919429 |
| D_eff | -0.77629766 |
| ShearMagMean | 0.22797396 |

**WindowTag = HOUR9** *(restricted: N < 15, features: Intercept, D_eff, ShearMagMean)*

| Feature | Coefficient |
|---|---:|
| Intercept | 0.95379960 |
| D_eff | -0.26040931 |
| ShearMagMean | 0.01183056 |

**WindowTag = OTHER**

| Feature | Coefficient |
|---|---:|
| Intercept | -0.92895289 |
| D_eff | 0.28361644 |
| F_W | 0.10851620 |
| chi_N | -2.58408508 |
| Ri_g | -0.02788360 |
| Ri_b | 0.00000023 |
| UShearMean | -0.63304430 |
| VShearMean | -0.64537941 |
| Deff_x_ShearMag | 2.95276211 |
| ShearMagMean | -10.70036039 |
| DeltaU_Dt | -0.00009173 |
| DeltaV_Dt | 0.00073040 |
| InertialMagDt | -0.00038021 |

