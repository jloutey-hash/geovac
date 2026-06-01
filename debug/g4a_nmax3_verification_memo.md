# G4a n_max=3 Verification Memo

**Date:** 2026-05-31
**Extends:** Sprint G4a (same date)

## Goal

Extend the Connes Standard Model axiom verification from n_max in {1, 2} to n_max=3. The full almost-commutative triple T = T_GV(n_max) x T_F^SM with A_F = C + H + M_3(C), dim H_F = 32.

## Dimensions

| n_max | dim_GV | dim_F | dim_H | Build time |
|:-----:|:------:|:-----:|:-----:|:----------:|
| 1 | 4 | 32 | 128 | 0.01s |
| 2 | 16 | 32 | 512 | 0.06s |
| 3 | 40 | 32 | 1280 | 0.61s |

## Axiom Results

### Bit-exact (zero residual at all n_max)

| Axiom | n_max=1 | n_max=2 | n_max=3 |
|:------|:-------:|:-------:|:-------:|
| J^2 = -I (KO-dim 1) | 0 | 0 | 0 |
| JD = +DJ | 0 | 0 | 0 |
| D Hermitian | 0 | 0 | 0 |
| gamma^2 = I | 0 | 0 | 0 |
| [a, JbJ^{-1}] = 0 (order-zero) | 0 | 0 | 0 |
| [[D,a], JbJ^{-1}] = 0 (order-one) | 0 | 0 | 0 |

All six load-bearing axioms pass at machine-zero (literal 0.0 in float64) at n_max=3. No finite-resolution degradation.

### Structural non-zero (documented, not violations)

| Check | n_max=1 | n_max=2 | n_max=3 | Normalized |
|:------|:-------:|:-------:|:-------:|:----------:|
| {gamma, D} Frobenius | 33.9 | 103.7 | 220.0 | /dim: 0.265, 0.203, 0.172 (decreasing) |
| [J, gamma] Frobenius | 22.6 | 45.3 | 71.6 | /sqrt(dim) = 2.000 exactly |

These are NOT axiom violations. The truthful Camporesi-Higuchi Dirac is chirality-diagonal, so {gamma, D} != 0 is structurally expected (KO-dim 1 is odd -- no chirality grading axiom). The [J, gamma] residual equals exactly 2*sqrt(dim_H) at every n_max, a structural constant reflecting the independent Z_2 grading (same finding as Sprint G3). Per-element {gamma, D} norm decreases with n_max.

## Gauge Group (n_max=3)

| Sector | Present? |
|:-------|:--------:|
| U(1) | Yes |
| SU(2) lepton | Yes |
| SU(2) quark | Yes |
| SU(3) color | Yes |
| Lepton-quark mixing | No (correct) |
| **Gauge group** | **U(1) x SU(2) x SU(3)** |

Full SM gauge content confirmed at n_max=3, identical to n_max=2.

## Falsifier

| Config | Higgs zero? | higgs/gauge ratio |
|:-------|:-----------:|:-----------------:|
| Yukawa imposed | No | 0.003 |
| Zero Yukawa | Yes | -- |

POSITIVE-THIN verdict extends to n_max=3: Higgs sector is non-trivial when Yukawa is imposed, trivially zero otherwise. GeoVac admits the SM but does not select the Yukawa.

## Verdict

All load-bearing Connes axioms are bit-exact zero at n_max=3 (dim_H=1280), confirming the Sprint G4a POSITIVE-THIN result. No axiom violation scales with cutoff; the construction is structurally exact at finite n_max. Full U(1) x SU(2) x SU(3) gauge group confirmed.

## Files

- Driver: `debug/g4a_nmax3_verification.py`
- Data: `debug/data/g4a_nmax3_verification.json`
