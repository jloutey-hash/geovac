# Track AD: PK-Free Diagnostic and Ab Initio R_ref

**Date:** 2026-03-31
**Version:** v2.0.17

## Part 1: PK-Free Diagnostic

**Question:** Is the 25% residual l_max drift in composed LiH from PK or from the composed-geometry separation itself?

### Results (adiabatic, l_dependent channel mode)

| pk_mode | l_max | R_eq (bohr) | R_eq error % | D_e (Ha) |
|---------|-------|-------------|--------------|----------|
| none | 1 | 2.000* | -33.7% | 0.4455 |
| none | 2 | 2.000* | -33.7% | 0.4188 |
| ab_initio | 1 | 2.500 | -17.1% | 0.2024 |
| ab_initio | 2 | 3.278 | +8.7% | 0.1911 |
| ab_initio | 3 | 3.422 | +13.5% | 0.1397 |

*Grid-edge minimum — PES is monotonically attractive.

### Drift rates
- ab_initio PK (adiabatic): +0.461 bohr/l_max (l_max 1→3)

### Verdict
**The residual drift is NOT from PK.** Without PK, the PES is monotonically attractive at ALL l_max values — there is no equilibrium geometry. PK is essential for creating the equilibrium. The drift with increasing l_max is from the adiabatic approximation (differential angular correlation in higher-l channels adds binding preferentially at large R).

## Part 2: Ab Initio R_ref

### R_ref Candidates from Core Screening (Li, Z=3)

| Candidate | Value (bohr) | Source |
|-----------|-------------|--------|
| r_inv2 | 0.2712 | <1/r^2>-weighted core radius |
| r_pk_width | 0.3836 | 1/sqrt(B) from PK Gaussian |
| r_steepest | 0.3575 | max |dZ_eff/dr| |
| r_half_screen | 0.5006 | Z_eff = Z-1 |
| r_median | 0.5010 | N_core = 1 |
| r_avg | 0.5767 | mean core radius <r> |
| r_rms | 0.6774 | RMS core radius |
| r_peak | 0.8167 | peak of r^2*n(r) |

### R_ref Sweep at l_max=2

| R_ref (bohr) | R_eq (bohr) | R_eq error % | Effect |
|-------------|-------------|--------------|--------|
| None (baseline) | 3.278 | +8.7% | No R-dep PK |
| 0.5 (r_half_screen) | 3.278 | +8.7% | Zero effect |
| 1.0 | 3.278 | +8.7% | Zero effect |
| 2.0 | 3.278 | +8.7% | Zero effect |
| 3.0 | 2.844 | -5.7% | Pulls R_eq in |
| 5.0 | 2.500 | -17.1% | Too strong |
| 10.0 | 2.000 | -33.7% | Approaching PK-free |

### Analysis

All core-derived R_ref candidates (0.27 - 0.82 bohr) are far below the
minimum R in the PES scan (2.0 bohr). Since w_PK(R) = min(1, R/R_ref),
when R_ref < R_min_scan, w_PK = 1 everywhere and R-dependent PK has
ZERO EFFECT.

The R_ref value needed to improve R_eq at l_max=2 is ~2.5 bohr (between
the 2.0 and 3.0 sweep points). This is the EQUILIBRIUM DISTANCE ITSELF
— a molecular property, not derivable from atomic core data without
circular reasoning.

### Verdict

**No ab initio R_ref achieves <3%.** The R-dependent PK formula
w_PK(R) = delta_{l,0} * min(cap, R/R_ref) requires R_ref comparable
to R_eq (~3 bohr) to have any effect. All core-derived candidates are
~0.3-0.8 bohr, 4-10x too small. The formula is structurally unable to
connect atomic core properties to the molecular equilibrium scale.

The l_max drift is from the adiabatic approximation, not from PK. The
2D solver (Track A) is the correct fix path — it reduces drift by 4x by
bypassing the adiabatic bottleneck. Further improvement requires coupled-
channel Level 4 or non-adiabatic corrections, not PK modifications.

## Files Created

- `geovac/pk_rref.py` — Ab initio R_ref module (R_ref candidates + w_PK function)
- `tests/test_pk_rref.py` — 15 tests (all passing)
- `geovac/composed_diatomic.py` — Modified: added pk_rref parameter
- `debug/track_ad/` — Diagnostic scripts and results
