# Pre-registered prediction: same-size Gaussian control for Paper 62's water result (SZ)

**Registered 2026-11-11**, before any water energy in the three new bases is computed. PI-directed ("Ok, yeah,
let's do that.") after the PM raised that Paper 62 tests its 24-function Slater bases only against the cc-pV
family, while the benchmarking rule asks for the strongest available baseline.
**Commit status:** committed on its own, before any computation it scores.
**Rests on:** Paper 62 Tables I and V and the ladder (registry `p62_ccpvdz_fc` -76.241025, `p62_dzp_fc`
-76.281739, `p62_r5_fc` -76.244112, `p62_r5h_fc` -76.248015; R_OH 1.809 bohr, 104.5 deg; own solver,
`tests/test_water_ccsd.py`); record MA (`docs/predictions/2026-10-31_matched_accuracy_gaussian.md`).

## Why

cc-pVDZ is built so that DZ, TZ, QZ converge systematically, not to give the lowest energy at its size. If a
Gaussian basis of the same size and the same qubit count reaches the Slater energies, Paper 62's equal-size gain and
its matched-accuracy cost ratio (aug-cc-pVDZ, 9.2x the Pauli terms) were measured against a weak baseline.

## The bases (fixed before any scored run)

All 24 functions on water (counted, below). Downloaded from the Basis Set Exchange API on 2026-11-11 into
`debug/lit_scan/dl_bse_samesize/` (excluded from git) and re-serialized by `debug/samesize_gaussian_control.py prep`
into `debug/data/basis/{pcseg-1,def2-svp,ano-r1}_H_O_Ne.json`.

| Basis | Built for | Water functions |
|:--|:--|:--:|
| pcseg-1 (Jensen) | DFT/HF energies, segmented | 24 |
| def2-SVP (Weigend-Ahlrichs) | HF/DFT, economical | 24 |
| ANO-R1 (Zobel-Widmark-Veryazov) | correlated atomic natural orbitals | 24 |

ANO-R0 was downloaded first and is the minimal member (7 functions on water); it was replaced by ANO-R1 before any
energy was computed. Scored energy: frozen-core CCSD(T), same geometry, solver and Pauli counter as record MA.

## What I know before comparing (disclosed)

- **cc-pVDZ control (measured 2026-11-11):** water RHF -76.026795, fc-CCSD(T) -76.241025, 107,382 Pauli terms
  (fc, 46 qubits), 1-norm 674.29. All match Paper 62 and the registry to the printed digits.
- **Normalization control:** on Ne every basis has unit overlap diagonal to 3e-16; ANO-R1's overlap is the identity
  to 1e-10 (its contractions are orthogonal natural orbitals, so the general contraction is read correctly).
- **Neon controls (non-scored, measured 2026-11-11; `debug/data/samesize_control/matched_accuracy_controls.jsonl`),
  against Ne cc-pVDZ RHF -128.488776, fc-CCSD(T) -128.678836 (record MA):**

| Ne basis | RHF | fc-CCSD(T) | RHF vs cc-pVDZ (mHa) | fc vs cc-pVDZ (mHa) |
|:--|:--|:--|:--:|:--:|
| pcseg-1 | -128.445448 | -128.622494 | +43.3 (worse) | +56.3 (worse) |
| def2-SVP | -128.376407 | -128.564755 | +112.4 (worse) | +114.1 (worse) |
| ANO-R1 | -128.521875 | -128.700170 | -33.1 (better) | -21.3 (better) |

  This contradicts the PM's statement to the PI on 2026-11-11 that energy-built bases "routinely" beat cc-pVDZ on
  total energy: pcseg-1 and def2-SVP spend fewer primitives on the core and lose on the atom. Only ANO-R1 wins.
- Atom-to-water transfer ratios measured in Paper 62 for other basis steps: 0.50, 0.52, 0.64. Neon is not oxygen,
  and the hydrogen functions also change, which the neon control does not see, so predictions are directional.
- The Slater DZP's polarization exponents (O 3d, H 2p) were optimized on water's frozen-core CCSD (Paper 62); the
  rungs R5 and R5H set nothing on water. The fair comparison for "nothing tuned on water" is R5/R5H.

## Predictions

Notation: W(B) = water fc-CCSD(T) in basis B, mHa below cc-pVDZ (positive = better than cc-pVDZ).
Known: W(DZP) = 40.7, W(R5H) = 7.0, W(R5) = 3.1.

- **SZ1.** W(pcseg-1) < 0 and W(def2-SVP) < 0: both land above (worse than) cc-pVDZ, def2-SVP the worse of the two.
- **SZ2.** W(ANO-R1) lies in [5, 25] mHa (central guess about 12, from neon's 21.3 times a transfer of 0.5-0.65).
- **SZ3 (the decision).** The Slater DZP stays below the best same-size Gaussian by more than 10 mHa:
  W(DZP) - max over the three of W(B) > 10.
- **SZ4.** ANO-R1 beats the atom-set rung R5H: W(ANO-R1) > 7.0.
- **SZ5.** All three bases give exactly 107,382 Pauli terms on 46 frozen-core qubits (Paper 62's statement that the
  count is fixed by orbital count and the C2v pattern).

## What each outcome means for Paper 62

- SZ3 holds: the equal-size Slater DZP advantage survives the strongest same-size Gaussian tested; Paper 62 adds
  that baseline and the margin it measures.
- SZ3 fails: the matched-accuracy headline (aug-cc-pVDZ at 9.2x) was priced against the wrong baseline and must be
  re-priced against the best same-size Gaussian; raised to the PI.
- SZ4 holds: "set on the free atoms, below cc-pVDZ" stays true but is no longer below the best same-size Gaussian;
  Paper 62's atom-set claim must say so. SZ4 fails: the parameter-free result is stronger than expected.
- SZ5 fails: the Pauli-count statement in Paper 62 is wrong as stated; raised to the PI.

Scored by: `python debug/samesize_gaussian_control.py water {pcseg1,def2svp,anor1}`.

## Scored (2026-11-11, after the registration commit e00868ff)

Rows: `debug/data/samesize_control/matched_accuracy.jsonl`. W = mHa below cc-pVDZ's fc-CCSD(T) -76.241025.

| Basis | RHF | fc-CCSD(T) | W (mHa) | RHF part (mHa) | Pauli (fc, 46 qubits) | 1-norm (Ha) |
|:--|:--|:--|:--:|:--:|:--:|:--:|
| pcseg-1 | -76.003683 | -76.217263 | -23.8 | -23.1 | 107,382 | 674.10 |
| def2-SVP | -75.961010 | -76.174699 | -66.3 | -65.8 | 107,382 | 677.37 |
| ANO-R1 | -76.047820 | -76.264667 | +23.6 | +21.0 | 107,382 | 625.36 |
| (cc-pVDZ, control) | -76.026795 | -76.241025 | 0 | 0 | 107,382 | 674.29 |
| (Slater DZP, Paper 62) | -76.055573 | -76.281739 | +40.7 | +28.8 | 107,382 | 706.3 |
| (R5H, Paper 62) | | -76.248015 | +7.0 | | | |

- **SZ1 HELD.** pcseg-1 -23.8, def2-SVP -66.3; both worse than cc-pVDZ, def2-SVP the worse.
- **SZ2 HELD, at the edge.** ANO-R1 +23.6 in [5, 25]. The central guess (about 12) missed: the neon-to-water
  transfer is about 1.1, not 0.5-0.65.
- **SZ3 HELD.** DZP - ANO-R1 = 17.1 mHa > 10. The Slater DZP's margin over the best same-size Gaussian tested is
  17.1 mHa, not the 40.7 measured against cc-pVDZ.
- **SZ4 HELD.** ANO-R1 beats R5H by 16.6 mHa: the rung set on free atoms is below cc-pVDZ but above ANO-R1.
- **SZ5 HELD.** 107,382 Pauli terms in every basis.
- **Unregistered observation.** ANO-R1 is more accurate than cc-pVDZ with a lower 1-norm (625 vs 674 Ha), so Paper
  62's "the 1-norm rises with accuracy" holds across its own kits but not across Gaussian families.

**Correction (2026-11-11, same day, after scoring).** SZ4 was registered against R5H (7.0 mHa below cc-pVDZ), but
Paper 62's best atom-set rung is R5H with the 3d also set on the atom (record
`docs/predictions/2026-11-01_scale_ladder_3d.md`): -76.260011, 19.0 mHa below cc-pVDZ. Against that rung ANO-R1 is
ahead by 4.7 mHa, not 16.6. SZ4 holds as registered; the PM's first report to the PI quoted the 16.6 figure and was
corrected in the same session. Paper 62 carries the 4.7 mHa figure.
