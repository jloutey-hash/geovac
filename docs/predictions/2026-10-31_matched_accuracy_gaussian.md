# Pre-registered prediction: the matched-accuracy Gaussian comparison for water (MA)

**Registered 2026-10-31**, before any water energy of the two new bases is computed. PI delegated the
registration ("We can download that data. Seems fine to me." / "Register is your call tho."). Owed by Paper 62
(Sec. qubit, Sec. limits): the equal-size statement there is not a matched-accuracy one.
**Commit status:** committed on its own, before any computation it scores.
**Rests on:** Paper 62 Tables I and V (cc-pVDZ fc-CCSD(T) -76.241025, Slater DZP -76.281739, cc-pVTZ -76.332191,
all 12-Gaussian evaluator, own solver, R_OH 1.809 bohr, 104.5 deg; registry `p62_*`); CHANGELOG v5.20.11.

## Why

Paper 62 shows the Slater DZP 40.7 mHa below cc-pVDZ at the same 24 functions and the same Pauli count. That is an
equal-size comparison. The matched-accuracy question is the other way round: which Gaussian basis first reaches the
Slater DZP's energy, and what does it cost in qubits, Pauli terms and 1-norm? cc-pVTZ (58 functions) overshoots by
50 mHa. Two intermediate Gaussian rungs are added so the crossing is bracketed from both sides.

## The ladder (fixed before any scored run)

Same geometry, solver, frozen-core convention and Pauli counter as Paper 62 (`tests/test_water_ccsd.py`). Driver
`debug/matched_accuracy_gaussian.py`; basis files `debug/data/basis/{cc-pvdz,ccpvdz_tzpol,aug-cc-pvdz,cc-pvtz}_H_O_Ne.json`
(re-serialized from the Basis Set Exchange, ccRepo/Grant Hill revision 2018-10-03, downloaded 2026-10-31 into
`debug/lit_scan/dl_bse_augccpvdz/`, excluded from git).

| Basis | Water functions | Status |
|:--|:--:|:--|
| cc-pVDZ | 24 | known: fc-CCSD(T) -76.241025 (Paper 62) |
| cc-pVDZ+TZpol | 35 | scored here; cc-pVDZ's s and p shells with cc-pVTZ's polarization shells (O 2d, H 2p); our construction |
| aug-cc-pVDZ | 41 | scored here; the standard augmented set |
| cc-pVTZ | 58 | known: fc-CCSD(T) -76.332191 (Paper 62) |

Target: the Slater DZP's fc-CCSD(T), -76.281739 Ha (24 functions). Scored energy: frozen-core CCSD(T). RHF,
all-electron, Pauli count and 1-norm (fc and all-electron) reported for both new rungs.

## What I know before comparing (disclosed)

- No water energy of cc-pVDZ+TZpol or aug-cc-pVDZ has been computed with this solver.
- **Neon controls (non-scored, measured 2026-10-31; `debug/data/matched_accuracy_controls.jsonl`):**

| Ne basis | Functions | RHF | fc-CCSD(T) | RHF gain over cc-pVDZ (mHa) | fc gain (mHa) |
|:--|:--:|:--|:--|:--:|:--:|
| cc-pVDZ | 14 | -128.488776 | -128.678836 | | |
| cc-pVDZ+TZpol | 19 | -128.488776 | -128.705023 | 0.00 | 26.19 |
| aug-cc-pVDZ | 23 | -128.496350 | -128.709295 | 7.57 | 30.46 |
| cc-pVTZ | 30 | -128.531862 | -128.802454 | 43.09 | 123.62 |

  On the closed-shell atom the second polarization set is pure correlation (RHF unchanged to 1e-11) and worth 26
  mHa; augmentation is worth 30 mHa, a quarter of it RHF. Water needs 40.7 mHa over cc-pVDZ to reach the Slater DZP.
- The water augmentation gain will differ from neon's: water has diffuse functions on three centres and a bond
  region, and its RHF is further from its limit (41 mHa vs neon's cc-pVDZ, which is closer). No transfer ratio
  atom-to-water has been measured for augmentation, so the predictions below are directional.
- Timing: neon cc-pVTZ (30 functions) took 341 s of integrals; water at 41 functions is expected near 20-30 min.

## Predictions

Notation: W(B) = water fc-CCSD(T) in basis B. Target T = -76.281739.

**MA1 (gated).** W(aug-cc-pVDZ) > T: augmentation at double-zeta does not reach the Slater DZP.
- Falsifier: W(aug-cc-pVDZ) <= T.
- Reason: neon's augmentation gain is 30.5 mHa against the 40.7 needed; water would have to gain a third more than
  neon from the same step.
- Probability 0.65.

**MA2 (gated).** W(cc-pVDZ+TZpol) > T: a second polarization set does not reach the Slater DZP either.
- Falsifier: W(cc-pVDZ+TZpol) <= T.
- Reason: neon gains 26.2 mHa, all correlation; the Slater DZP's lead over cc-pVDZ is 29 mHa of SCF plus 12 of
  correlation, and polarization functions buy little SCF energy.
- Probability 0.8.

**MA3 (reported, not gated).** The RHF share of each new rung's fc gain: TZpol <= 0.3 (neon 0.00), augmentation
<= 0.5 (neon 0.25). Direction: both rungs are correlation-dominated steps, unlike the Slater DZP step, which is
SCF-dominated (0.71).
- Probability 0.7 each.

**MA4 (reported, not gated).** The two rungs are ordered as on neon: W(aug-cc-pVDZ) < W(cc-pVDZ+TZpol).
- Probability 0.6 (on neon the margin is 4 mHa; water's three centres could reverse it).

Also reported, not predicted: the Pauli counts and 1-norms at 35 and 41 functions (the count is fixed by orbital
count and symmetry and can be computed before the energy; it is not a prediction), and the interpolated
matched-accuracy cost, i.e. the Pauli count and 1-norm at the energy T, read linearly in energy between the two
rungs that bracket T. Under MA1 and MA2 the bracket is aug-cc-pVDZ to cc-pVTZ.

**Gate.** GO if MA1 and MA2 both hold. BORDERLINE if one misses. STOP if both miss (then a Gaussian basis of 35-41
functions matches the Slater DZP, and the matched-accuracy cost is at most that basis's).

## Caps

- One geometry, one method, one molecule, as in Paper 62. The intermediate rung cc-pVDZ+TZpol is our construction,
  not a standard basis; it is there to bracket, not to recommend.
- cc-pVDZ's polarization is the published one here (Paper 62 also reports it re-optimized, 1.2 mHa lower); the
  Slater target keeps its water-optimized polarization. The comparison is standard Gaussian bases against the
  Slater DZP as built.
- A linear interpolation of Pauli count in energy between two rungs is a reading, not a measurement of any basis.

---

(Appended after comparing. The predictions above are not edited.)

## Disposition (2026-10-31)

Runs: `python debug/matched_accuracy_gaussian.py water ccpvdz_tzpol`, `water augccpvdz`, after the record was
committed (6b6d85c2). Rows in `debug/data/matched_accuracy.jsonl`, log `debug/data/matched_accuracy.log`. Integrals
213 s and 293 s, solves 27 s and 51 s, under the p-function load.

| Basis | Functions | RHF | fc-CCSD(T) | ae-CCSD(T) | vs target T (mHa) | fc gain over cc-pVDZ (mHa) | RHF share |
|:--|:--:|:--|:--|:--|:--:|:--:|:--:|
| cc-pVDZ | 24 | -76.026795 | -76.241025 | -76.243143 | +40.71 | | |
| cc-pVDZ+TZpol | 35 | -76.031742 | -76.266498 | -76.269454 | +15.24 | 25.47 | 0.19 |
| aug-cc-pVDZ | 41 | -76.041422 | -76.273743 | -76.276008 | +8.00 | 32.72 | 0.45 |
| Slater DZP (target) | 24 | -76.055573 | -76.281739 | -76.285720 | 0 | 40.71 | 0.71 |
| cc-pVTZ | 58 | -76.057162 | -76.332191 | | -50.45 | 91.17 | 0.33 |

- **MA1: HELD.** W(aug-cc-pVDZ) = -76.273743, 8.00 mHa above T.
- **MA2: HELD.** W(cc-pVDZ+TZpol) = -76.266498, 15.24 mHa above T.
- **MA3 (reported): HELD** both. RHF share 0.19 (TZpol, bound 0.3) and 0.45 (augmentation, bound 0.5). Water's
  augmentation step is more SCF-heavy than neon's (0.45 vs 0.25), as the disclosed caveat allowed.
- **MA4 (reported): HELD.** aug-cc-pVDZ below cc-pVDZ+TZpol by 7.2 mHa (neon 4.3).
- **Gate: GO.**

Reported, not predicted (Jordan-Wigner, RHF canonical orbitals, identity in N_Pauli and out of the 1-norm):

| Basis | Q fc | N_Pauli fc | 1-norm fc (Ha) | Q ae | N_Pauli ae | 1-norm ae (Ha) |
|:--|:--:|:--:|:--:|:--:|:--:|:--:|
| Slater DZP | 46 | 107,382 | 706.3 | 48 | 128,793 | 754.6 |
| cc-pVDZ+TZpol | 68 | 505,467 | 1,814.9 | 70 | 570,822 | 1,871.5 |
| aug-cc-pVDZ | 80 | 982,785 | 1,901.1 | 82 | 1,091,076 | 1,958.3 |
| cc-pVTZ | 114 | 3,983,652 | 7,653.5 | 116 | 4,281,643 | 7,858.6 |

- The crossing lies between aug-cc-pVDZ and cc-pVTZ, 0.137 of the way in energy. Linear interpolation in energy
  puts the matched-accuracy Gaussian cost near 1.39 million terms, a 1-norm near 2,690 Ha and about 85 qubits:
  13 times the Slater DZP's terms and 3.8 times its 1-norm. The lower bound needs no interpolation: aug-cc-pVDZ,
  still 8 mHa short of the target, already costs 9.2 times the terms (982,785 vs 107,382), 2.7 times the 1-norm
  (1,901 vs 706) and 80 qubits against 46.
- External sanity (orientation only, not a cited reference): the aug-cc-pVDZ fc-CCSD(T) energy agrees with the
  commonly quoted value near -76.274 at the equilibrium geometry.

Scope, stated so it is not over-read: one geometry, one method, standard Gaussian bases as published against the
Slater DZP as built (water-optimized polarization). The interpolation is a reading. The result upgrades Paper 62's
qubit statement from equal-size to matched-accuracy: at matched or better accuracy the Gaussian basis costs at least
nine times the Pauli terms and 2.7 times the 1-norm.

## Addendum (2026-11-11): the "at least nine times" reading is superseded

The bracket above is for the cc-pV family only. A later registered run on the ANO family
(`docs/predictions/2026-11-11_ano_matched_accuracy.md`, record AM) found a 35-function truncation of ANO-R2
that is 6.2 mHa below the Slater DZP at 4.7 times its Pauli terms and 2.3 times its 1-norm. The matched-accuracy
cost is therefore at most 4.7 times the terms, not at least nine. The measurements above stand; only the closing
"at least nine times" sentence is withdrawn. Paper 62, Sec. "Matched accuracy", is the current statement.
