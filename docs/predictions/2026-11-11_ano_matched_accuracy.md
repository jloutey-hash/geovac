# Pre-registered prediction: the ANO-family matched-accuracy ladder for water (AM)

**Registered 2026-11-11**, before any scored water energy of the new rungs is computed. PI-directed ("Ok, yeah.
Let's do that stuff.") after record SZ found ANO-R1 23.6 mHa below cc-pVDZ at 24 functions.
**Commit status:** committed on its own, before any computation it scores.
**Rests on:** Paper 62 Tables dzp, matched and samesize (registry `p62_dzp_fc` -76.281739, `p62_ccpvdz_fc`
-76.241025, `p62_tzpol_fc` -76.266498, `p62_aug_fc` -76.273743, `p62_ccpvtz_fc` -76.332191, `p62_anor1_fc`
-76.264667, `p62_tzpol_pauli_fc` 505,467, `p62_aug_pauli_fc` 982,785); records MA and SZ.

## Why

Paper 62's matched-accuracy statement (at matched accuracy the Gaussian bases cost at least 9.2x the Pauli terms of
the Slater DZP) was measured on the cc-pV ladder only. ANO-R1 closes more than half of the equal-size gap, so the
ANO family may reach the Slater DZP sooner.

## The ladder (fixed before any scored run)

ANO-R1 is a column truncation of ANO-R2: same primitives, and R1's contraction columns equal R2's first columns
(checked 2026-11-11 on H and O, every shell). Every rung is ANO-R2 (Basis Set Exchange, revision 2021-09-25,
downloaded 2026-11-11 into `debug/lit_scan/dl_bse_samesize/`, excluded from git) kept to its first n contracted
functions per shell. The truncation is ours; it is how ANO sets are meant to be used. Driver
`debug/ano_matched_accuracy.py`; files `debug/data/basis/ano-r2_A{24,35,41,61}_H_O_Ne.json`.

| Rung | O | H | Water functions | Same shell composition as | Status |
|:--|:--|:--|:--:|:--|:--|
| A24 = ANO-R1 | [3s2p1d] | [2s1p] | 24 | cc-pVDZ | known: -76.264667 (record SZ) |
| A35 | [3s2p2d] | [2s2p] | 35 | cc-pVDZ+TZpol | scored here |
| A41 | [4s3p2d] | [3s2p] | 41 | aug-cc-pVDZ | scored here |
| A61 = ANO-R2 | [4s4p2d1f] | [3s2p1d] | 61 | none | neon control only; not scored |

Target: the Slater DZP's fc-CCSD(T), -76.281739 Ha. Scored energy: frozen-core CCSD(T), same geometry, solver and
Pauli counter as records MA and SZ.

## What I know before comparing (disclosed)

- **Truncation control (measured 2026-11-11):** water A24 gives fc-CCSD(T) -76.2646673, ANO-R1's registered
  -76.264667 to 3e-7 (inside its six-decimal rounding). So the truncation code reproduces the published ANO-R1.
- **Neon controls (non-scored, 2026-11-11; `debug/data/ano_matched_accuracy/matched_accuracy_controls.jsonl`),
  fc-CCSD(T) relative to Ne cc-pVDZ -128.678836, beside the cc-pV rungs of the same size (record MA):**

| Ne rung | Functions | RHF | fc-CCSD(T) | vs cc-pVDZ (mHa) | cc-pV rung of that size (mHa) |
|:--|:--:|:--|:--|:--:|:--|
| A24 | 14 | -128.521875 | -128.700170 | -21.3 | cc-pVDZ 0 |
| A35 | 19 | -128.521875 | -128.713428 | -34.6 | cc-pVDZ+TZpol -26.2 |
| A41 | 23 | -128.529145 | -128.775688 | -96.9 | aug-cc-pVDZ -30.5 |
| A61 | 33 | -128.529640 | -128.809554 | -130.7 | cc-pVTZ (30 fns) -123.6 |

  A41's fourth s and third p contracted functions are correlating functions on neon, which aug-cc-pVDZ's diffuse
  functions are not; hence the factor of three at equal size.
- **Neon-to-water transfer of the step from cc-pVDZ, measured on other bases:** cc-pVDZ+TZpol 0.97, aug-cc-pVDZ 1.07,
  ANO-R1 1.11, cc-pVTZ 0.74. Range used below: 0.74 to 1.11.
- The Slater DZP is 40.7 mHa below cc-pVDZ on water. Scaled from neon, A35 lands 25.6 to 38.4 mHa below cc-pVDZ (short
  of the target) and A41 lands 71.7 to 107.5 below (past it).

## Predictions

Notation: W(B) = water fc-CCSD(T) in basis B; D = -76.281739 (the Slater DZP).

- **AM1 (gated).** W(A35) > D: the 35-function rung does not reach the Slater DZP. Probability 0.75.
- **AM2 (gated).** W(A41) < D: the 41-function rung does. Probability 0.9.
- **AM3.** Pauli counts (fc, identity included) equal those of the cc-pV rung with the same shell composition:
  A35 505,467 on 68 qubits, A41 982,785 on 80 qubits.
- **AM4 (reported).** Each ANO rung is below the cc-pV rung of the same size: W(A35) < -76.266498 and
  W(A41) < -76.273743.

Gate: GO if AM1 and AM2 both hold.

## What each outcome means for Paper 62

- GO: the matched-accuracy crossing in the ANO family lies between 35 and 41 functions. The cheapest Gaussian rung
  measured to reach the Slater DZP then costs what aug-cc-pVDZ costs (9.2x the terms), but the bracket's lower side
  is A35 at 4.7x, so "at least nine times" no longer follows; Paper 62 states the bound as between 4.7x and 9.2x.
- AM1 fails (A35 reaches the target): the matched-accuracy cost is at most 4.7x the terms; Paper 62's
  matched-accuracy headline is restated at that figure and raised to the PI.
- AM2 fails: the ANO crossing lies beyond 41 functions and the cc-pV statement stands for both families.

Scored by: `python debug/ano_matched_accuracy.py water A35` and `... water A41`.

## Outcome (scored 2026-11-11, after registration commit 49c51a5e)

Rows: `debug/data/ano_matched_accuracy/` (same driver, solver and Pauli counter as records MA and SZ).

| Rung | Functions | RHF | fc-CCSD(T) | vs D (mHa) | vs cc-pVDZ (mHa) | Pauli (fc) | Qubits | 1-norm (fc) |
|:--|:--:|:--|:--|:--:|:--:|:--:|:--:|:--|
| A35 | 35 | -76.053557 | -76.287896 | -6.2 | -46.9 | 505,467 | 68 | 1632.93 |
| A41 | 41 | -76.057955 | -76.312401 | -30.7 | -71.4 | 982,785 | 80 | 2639.94 |

- **AM1 MISSED.** W(A35) = -76.287896 < D: the 35-function ANO rung already reaches the Slater DZP, by 6.2 mHa.
  The neon transfer read 25.6 to 38.4 mHa below cc-pVDZ; the measured step is 46.9 (transfer 1.36, outside the
  0.74-1.11 range taken from other bases).
- **AM2 HELD.** W(A41) = -76.312401, 30.7 mHa below D.
- **AM3 HELD.** 505,467 terms on 68 qubits and 982,785 on 80, equal to the cc-pV rungs of the same shell composition.
- **AM4 HELD.** A35 is 21.4 mHa below cc-pVDZ+TZpol; A41 is 38.7 mHa below aug-cc-pVDZ.

**Gate: not GO (AM1 failed).** Per the pre-registered reading: the matched-accuracy Gaussian cost is at most 4.7x the
Pauli terms of the Slater DZP (505,467 / 107,382), not at least 9.2x; the crossing in the ANO family lies between 24
and 35 functions and is not located more finely here. Paper 62's matched-accuracy headline is to be restated at that
figure; raised to the PI.
