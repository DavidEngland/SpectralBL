# GABLS3 Flux Predictive Baseline

This report fits a weighted linear baseline with sample weights inversely proportional to window sample count.

## Window Counts

| WindowTag | n_samples |
|---|---:|
| OTHER | 132 |
| HOUR8 | 6 |
| HOUR9 | 6 |

## Target: MomUFluxTruth

### Weighted Global Metrics

- wRMSE: 0.021272
- wMAE: 0.012244
- wR2: 0.521267
- wBias: -0.000000

### By-Window Metrics

| WindowTag | wRMSE | wMAE | wR2 | wBias | n_samples |
|---|---:|---:|---:|---:|---:|
| OTHER | 0.034864 | 0.022089 | 0.422532 | -0.004371 | 132 |
| HOUR8 | 0.011218 | 0.010702 | -0.207249 | 0.005514 | 6 |
| HOUR9 | 0.004020 | 0.003941 | -17.122938 | -0.001143 | 6 |

### Coefficients

| Feature | Coefficient |
|---|---:|
| Intercept | -11.41361990 |
| D_eff | 3.35825832 |
| F_W | 0.23274886 |
| chi_N | -32.71087362 |
| Ri_g | -0.00700234 |
| Ri_b | -0.00000329 |
| UShearMean | -2.67592217 |
| VShearMean | -1.81683099 |
| Deff_x_ShearMag | 0.66292760 |

## Target: HeatFluxTruth

### Weighted Global Metrics

- wRMSE: 0.014616
- wMAE: 0.009685
- wR2: 0.781106
- wBias: 0.000000

### By-Window Metrics

| WindowTag | wRMSE | wMAE | wR2 | wBias | n_samples |
|---|---:|---:|---:|---:|---:|
| OTHER | 0.023278 | 0.017582 | 0.729186 | 0.002440 | 132 |
| HOUR8 | 0.008318 | 0.006966 | -28.215264 | -0.003667 | 6 |
| HOUR9 | 0.005465 | 0.004508 | -299.708126 | 0.001227 | 6 |

### Coefficients

| Feature | Coefficient |
|---|---:|
| Intercept | 16.01180852 |
| D_eff | -4.14711474 |
| F_W | -0.17556385 |
| chi_N | -28.42146430 |
| Ri_g | -0.02077156 |
| Ri_b | -0.00000489 |
| UShearMean | 6.61433158 |
| VShearMean | 3.30468204 |
| Deff_x_ShearMag | -1.52765954 |

## Target: MomVFluxTruth

### Weighted Global Metrics

- wRMSE: 0.009633
- wMAE: 0.004090
- wR2: 0.087223
- wBias: -0.000000

### By-Window Metrics

| WindowTag | wRMSE | wMAE | wR2 | wBias | n_samples |
|---|---:|---:|---:|---:|---:|
| OTHER | 0.016526 | 0.009477 | 0.100674 | 0.000641 | 132 |
| HOUR8 | 0.001607 | 0.001311 | -2.168227 | -0.000190 | 6 |
| HOUR9 | 0.001647 | 0.001482 | -13.178952 | -0.000452 | 6 |

### Coefficients

| Feature | Coefficient |
|---|---:|
| Intercept | -0.94544198 |
| D_eff | 0.21683818 |
| F_W | -0.00278455 |
| chi_N | 5.62176052 |
| Ri_g | -0.00945808 |
| Ri_b | -0.00000071 |
| UShearMean | 0.47328793 |
| VShearMean | -0.11244670 |
| Deff_x_ShearMag | -0.16832090 |
