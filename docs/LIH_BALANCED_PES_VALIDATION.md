# LiH Balanced PES Validation Report

**Date:** 2026-03-11
**Version:** v0.9.37
**Script:** `debug/validate_balanced_pes.py`
**Status:** FAIL — no equilibrium

## Purpose

Validate the PES shape for the balanced `exact+True` configuration
(exact cross-nuclear + all-l cross-atom V_ee) discovered in the
Hamiltonian diagnostic. Compare against the known `fourier+s_only`
(v0.9.11 baseline) and the imbalanced `exact+s_only`.

## Atomic Reference Energies (nmax=3)

| Atom | E (Ha) |
|------|--------|
| Li   | -7.392086 |
| H    | -0.500000 |
| Sep  | -7.892086 |

## D_cp(R) Comparison Table

| R (bohr) | exact+True | exact+s_only | fourier+s_only |
|-----------|-----------|-------------|---------------|
| 1.500 | 0.4642 | 0.5828 | 0.4655 |
| 1.800 | 0.3417 | 0.5693 | 0.3429 |
| 2.000 | 0.2803 | 0.5602 | 0.2789 |
| 2.200 | 0.2272 | 0.5919 | 0.2279 |
| 2.500 | 0.1706 | 0.6328 | 0.1708 |
| 2.800 | 0.1316 | 0.6601 | 0.1310 |
| 3.000 | 0.1118 | 0.6616 | 0.1111 |
| 3.015 | 0.1106 | 0.6614 | 0.1098 |
| 3.200 | 0.0959 | 0.6608 | 0.0951 |
| 3.500 | 0.0852 | 0.6541 | 0.0762 |
| 4.000 | 0.0724 | 0.6209 | 0.0534 |
| 5.000 | 0.0416 | 0.5104 | 0.0227 |


## Equilibrium Summary

| Config | R_eq (bohr) | D_cp (Ha) | D_cp err% | R_eq err% |
|--------|------------|-----------|-----------|-----------|
| exact+True | 1.500 | 0.4642 | +404.6% | -50.2% |
| exact+s_only | 2.935 | 0.6620 | +619.6% | -2.6% |
| fourier+s_only | 1.500 | 0.4655 | +406.0% | -50.2% |
| Experiment | 3.015 | 0.092 | --- | --- |

## Answers to Key Questions

### Q1: Does exact+True have an equilibrium geometry?

**NO** — D_cp is monotonically attractive. This indicates the all-l V_ee repulsion may be insufficient to create a turnover, or the l>0 Gaunt coefficients need review.

### Q2: Does the PES shape match fourier+s_only?

Max |D_cp(exact+True) - D_cp(fourier+s_only)| = 0.0190 Ha

**NO** — Significant shape difference of 0.0190 Ha between the two balanced configurations.

### Q3: What is D_cp at equilibrium?

- exact+True: D_cp = 0.4642 Ha (paper target: 0.093)
- fourier+s_only: D_cp = 0.4655 Ha
- Estimate from diagnostic: 0.225 - 0.115 = 0.110 Ha

### Q4: Is the 17% geometric contraction still present?

R_eq(exact+True) = 1.500 bohr vs experiment 3.015 bohr
Contraction = -50.2% (paper prediction: -17%)

## Output Files

- `debug/data/lih_balanced_pes.txt` — full PES data
- `debug/plots/lih_balanced_pes_comparison.png` — overlay comparison plot
