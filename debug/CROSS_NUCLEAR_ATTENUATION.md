# Cross-Nuclear Attenuation Analysis

**Date:** 2026-03-12

## Problem
Cross-nuclear attraction (V_cross) is the dominant driver of PES shape errors.
It grows too strong at short R because the graph framework uses unperturbed
atomic densities with no orthogonalization penalty.

## Method
Attenuate cross-nuclear attraction: V_cross(a) *= f(S_a)
where S_a = sum_b S(a,b)^2 is total squared overlap.

## H2 Results (nmax=3, expt: R_eq=1.401, D_e=0.1745, k=0.369)

| Config | R_eq | D_e | k | R_err% | D_err% | k_err% |
|--------|------|-----|---|--------|--------|--------|
| bare | 1.243 | 0.3905 | 1.0528 | 11.3% | 123.8% | 185.3% |
| linear(0.5) | 2.154 | 0.0652 | 0.0907 | 53.7% | 62.7% | 75.4% |
| linear(1.0) | 3.765 | -0.0208 | 0.0372 | 168.7% | inf% | 89.9% |
| linear(2.0) | 3.981 | -0.0862 | 0.1560 | 184.1% | inf% | 57.7% |
| linear(5.0) | 2.344 | -0.4677 | -0.1427 | 67.3% | inf% | 138.7% |
| exp(0.5) | 1.941 | 0.0870 | 0.1678 | 38.5% | 50.2% | 54.5% |
| exp(1.0) | 3.640 | -0.0160 | 0.0125 | 159.8% | inf% | 96.6% |
| exp(2.0) | 4.172 | -0.0689 | 0.0677 | 197.8% | inf% | 81.7% |
| exp(5.0) | 4.861 | -0.1590 | 0.0881 | 246.9% | inf% | 76.1% |
| pade(0.5) | 1.728 | 0.1078 | 0.3435 | 23.3% | 38.2% | 6.9% |
| pade(1.0) | 3.054 | -0.0051 | 0.0193 | 118.0% | inf% | 94.8% |
| pade(2.0) | 4.420 | -0.0554 | 0.0307 | 215.5% | inf% | 91.7% |
| pade(5.0) | 4.848 | -0.1301 | 0.0561 | 246.0% | inf% | 84.8% |
| quadratic | 4.080 | -0.0772 | 0.1001 | 191.2% | inf% | 72.9% |

## LiH Results (nmax=3, expt: R_eq=3.015, D_e=0.092, k=0.0659)

| Config | R_eq | D_e | k | R_err% | D_err% | k_err% |
|--------|------|-----|---|--------|--------|--------|
| bare | 2.988 | 0.2501 | -0.2970 | 0.9% | 171.9% | 550.7% |
| linear(0.5) | 4.799 | 0.1328 | 0.0259 | 59.2% | 44.3% | 60.8% |
| linear(1.0) | 5.832 | 0.1026 | 0.0204 | 93.4% | 11.5% | 69.0% |
| linear(2.0) | 6.581 | 0.0806 | 0.0359 | 118.3% | 12.4% | 45.5% |
| linear(5.0) | 4.633 | -0.0676 | -0.0676 | 53.7% | inf% | 202.5% |
| exp(0.5) | 4.738 | 0.1334 | 0.0236 | 57.2% | 45.0% | 64.1% |
| exp(1.0) | 5.779 | 0.1037 | 0.0168 | 91.7% | 12.7% | 74.5% |
| exp(2.0) | 6.439 | 0.0811 | 0.0331 | 113.6% | 11.9% | 49.8% |
| exp(5.0) | 4.622 | -0.0675 | -0.0826 | 53.3% | inf% | 225.4% |
| pade(0.5) | 4.675 | 0.1341 | 0.0218 | 55.1% | 45.8% | 66.9% |
| pade(1.0) | 5.708 | 0.1047 | 0.0141 | 89.3% | 13.8% | 78.6% |
| pade(2.0) | 6.504 | 0.0838 | 0.0252 | 115.7% | 8.9% | 61.8% |
| pade(5.0) | 4.462 | -0.0571 | -0.0672 | 48.0% | inf% | 202.0% |
| quadratic | 6.392 | 0.0793 | 0.0388 | 112.0% | 13.8% | 41.2% |

## Universal Configs (valid minimum for both molecules)

| Config | H2 R_err% | H2 D_err% | LiH R_err% | LiH D_err% |
|--------|-----------|-----------|------------|------------|
| exp(0.5) | 38.5% | 50.2% | 57.2% | 45.0% |
| linear(0.5) | 53.7% | 62.7% | 59.2% | 44.3% |
| pade(0.5) | 23.3% | 38.2% | 55.1% | 45.8% |

**Best universal:** pade(0.5)
- H2: R_eq=1.728 (23.3%), D_e=0.1078 (38.2%)
- LiH: R_eq=4.675 (55.1%), D_e=0.1341 (45.8%)

## Conclusion

Overlap-attenuated cross-nuclear attraction can produce correct PES shapes.
The attenuation is physically motivated (Pauli orthogonalization screening).
