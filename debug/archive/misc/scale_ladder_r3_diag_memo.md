# R3 non-monotone step of the scale ladder: diagnostic (2026-10-30)

PI: "Let's do the R3 non-monotone diagnostic." Diagnostic only (feedback_diagnostic_before_engineering).
Nothing is tuned on water; every substituted exponent is a value the ladder already holds.
Record: `docs/predictions/2026-10-30_scale_ladder.md` (SL, GO). Driver: `debug/scale_ladder.py`.
This memo's hypotheses were written before any new number below was computed.

## The fact to explain

R2 -> R3 splits the O p exponent from the O s-valence exponent. The free O atom (radial HF, the ladder's
objective) gains 3.4 mHa. Water loses 2.85 mHa at fc-CCSD(T). From the stored rows (`debug/data/scale_ladder.jsonl`):

| Piece (water, R2 -> R3) | Change (mHa) |
|:--|--:|
| RHF | +0.13 |
| fc-CCSD correlation | +2.57 |
| (T) | +0.15 |
| fc-CCSD(T) | +2.85 |

Exponents: R2 s-valence = p = 3p = 3d = 2.28124. R3: s-valence 2.35807, p = 3p = 3d = 2.21747, 1s core 7.81.

## Hypotheses (stated before computing)

- **H-A (virtual drift, expected):** the ladder's objective is HF on the atom, which sees the occupied 2p and
  only weakly the 3p, and does not see the 3d at all (`radial_atom` has no d). R3 lets the p exponent go diffuse
  for the occupied shell; the 3p and 3d, which carry most of the CCSD correlation at 24 functions, ride along and
  get worse for correlation, which wants compact correlating functions. Prediction: putting the 3d (and 3p) back to
  2.28124 with everything else at R3 recovers most of the 2.85 mHa on water; and a one-centre closed-shell probe
  in the same Gaussian-fitted basis shows the same correlation loss from R2 to R3 exponents, so the effect is not
  molecular transfer at all.
- **H-B (RHF transfer):** ruled out by the table; RHF moves 0.13 mHa.
- **H-C (Gaussian fit artifact):** ruled out by construction. `tests/test_water_ccsd.sto_shape` fits one shape per
  (l, n, n_gauss) and `sto_shell` scales it by zeta^2, so the fit error is the same at every exponent.
- **H-D (s-valence shift):** the s-valence exponent tightening 2.281 -> 2.358 costs water. Tested by the same
  substitution in the other direction (s back to 2.28124, p family at 2.21747).

## Runs

Script `debug/scale_ladder_r3_diag.py`. (a) `oion`: O^2- (closed shell, 10 electrons, one centre, the O
primitive set only, 15 functions) fc-CCSD(T) through `test_water_ccsd.run` at the R2 and R3 exponent sets and the
substitutions. (b) `water <name>`: R3dB (3d at 2.28124), R3vB (3p and 3d at 2.28124), R3sB (s-valence at
2.28124). Rows: `debug/data/scale_ladder_r3_diag.jsonl`.

## Results

### One-centre probes (done before the water substitution runs finished)

Neutral O atom, triplet M_S = 1 determinant, UHF reference, fc-CCSD(T) through the same spin-orbital CCSD
(`debug/scale_ladder_r3_diag_atom.py`, `debug/data/scale_ladder_r3_diag_atom.json`). The ladder's own atom
objective is the radial average-of-configuration HF, so the UHF step differs from its 3.4 mHa.

| Set | UHF | fc-CCSD(T) corr | fc-CCSD(T) total | corr vs R2 (mHa) |
|:--|--:|--:|--:|--:|
| R2 | -74.555126 | -0.105160 | -74.660286 | |
| R3 | -74.557968 | -0.102436 | -74.660404 | +2.72 |
| R3dB (3d back at 2.28124) | -74.558026 | -0.104591 | -74.662618 | +0.57 |
| R3vB (3p and 3d back) | -74.559040 | -0.104834 | -74.663874 | +0.33 |
| R3sB (s-valence back) | -74.555983 | -0.102706 | -74.658690 | +2.45 |

O^2- closed shell (`scale_ladder_r3_diag.py oion`, `debug/data/scale_ladder_r3_diag_oion.json`): correlation
-162.41 (R2), -161.92 (R3), -163.29 (R3dB), -161.78 (R3vB), -162.33 (R3sB) mHa. Its RHF moves 34 mHa between the
sets because the anion wants other exponents, so it is a poor stand-in; the neutral atom is the probe that counts.

**Reading.** On the atom the R3 step gains 2.84 mHa at UHF and loses 2.72 mHa of correlation: a wash (-0.12 mHa) at
fc-CCSD(T). Water lost 2.72 mHa of correlation (CCSD 2.57 + (T) 0.15) and 0.13 mHa at RHF. So the correlation loss
transfers one to one and the HF gain does not transfer at all. The step is non-monotone on water only because the
ladder's objective is HF on the atom, and the two halves of the step transfer differently. About 80 percent of the
correlation loss is the 3d exponent (a passenger of the 3p in the ladder design, invisible to `radial_atom`), 10
percent the 3p, 10 percent the s-valence shift.

**Prediction for the water substitution runs (written before reading them):** if the correlation loss is the same
one-centre effect, water R3dB recovers about 2.1 mHa of the 2.85 (fc-CCSD(T) near -76.1106), R3vB about 2.4 mHa
(near -76.1109), R3sB about 0.3 mHa (near -76.1088), with RHF within 0.2 mHa of R3's -75.890159 in all three.

### Water substitution runs (`debug/data/scale_ladder_r3_diag.jsonl`, about 1150 s each under load)

| Run | RHF | fc-CCSD | fc-CCSD(T) | vs R3 (mHa) | RHF part | corr part |
|:--|--:|--:|--:|--:|--:|--:|
| R3 (reference, SL row) | -75.890159 | -76.103503 | -76.108484 | | | |
| R3dB (3d back at 2.28124) | -75.890568 | -76.105643 | -76.110616 | -2.13 | -0.41 | -1.72 |
| R3vB (3p and 3d back) | -75.896834 | -76.113248 | -76.118395 | -9.91 | -6.68 | -3.24 |
| R3sB (s-valence back) | -75.885064 | -76.097878 | -76.102842 | +5.64 | +5.10 | +0.55 |

**Scoring the prediction.** R3dB: held (predicted -76.1106, got -76.110616; RHF moved 0.41, just outside the 0.2
clause). R3vB: the correlation clause held (recovered 3.24, predicted about 2.4 or more), the RHF clause failed by
6.7 mHa. R3sB: the correlation clause missed by 0.85 mHa in the wrong direction, the RHF clause failed by 5.1 mHa.

**Exponent sensitivities at the HF level, atom (UHF) against water (RHF), in mHa:**

| Move | O atom UHF | water RHF |
|:--|--:|--:|
| 3d 2.217 -> 2.281 (at R3) | -0.06 | -0.41 |
| 3p 2.217 -> 2.281 (at R3dB) | -1.01 | -6.27 |
| s-valence 2.281 -> 2.358 (p family at 2.217) | -1.99 | -5.10 |
| p family 2.281 -> 2.217 (s at 2.281) | -0.86 | +5.23 |

## Verdict

The R3 step is non-monotone on water for two separate reasons, one per half of the step.

1. **Correlation half, transfers one to one.** The ladder's objective is HF on the atom. It does not see the 3d at
   all and sees the 3p only through the occupied 2p, so at R3 the correlating functions drift diffuse with the 2p.
   The atom loses 2.72 mHa of fc-CCSD(T) correlation; water loses 2.72 mHa. About 80 percent of it is the 3d
   (a passenger of the 3p by the ladder's design), the rest the 3p and the s shift. Putting the 3d back on water
   recovers 2.13 of the 2.85 mHa, as the atom predicted.
2. **HF half, does not transfer.** The atom gains at HF from both moves of the step (p family diffuser, s valence
   tighter). Water's RHF loses 5.2 mHa from the first and gains 5.1 mHa from the second, so the atom's HF gain
   cancels to 0.13 mHa on water. Water's RHF is two to six times more sensitive than the atom's to the 3p and
   s-valence exponents, and of the opposite sign for the p family: the free atom's 2p^4 wants a more diffuse p than
   the O in water's bonds. A step that moves exponents in two directions has no single HF transfer ratio.

Together: on the atom the step is a wash at fc-CCSD(T) (-0.12 mHa), and water sees only the half that transfers,
the loss. This is a diagnosis of the ladder's objective, not of the primitive family. Not done here, and a PI call
because it changes the ladder's rule: setting the 3d exponent by an atomic correlated criterion rather than tying it
to the 3p. R4 and R5, whose gains are dominated by HF-driven radial splits (the p split carries 0.90 on both), are not
touched by this; their transfer ratios 0.735 and 0.515 remain the record's.

H-A confirmed for the correlation half. H-B is replaced by the cancellation mechanism above (RHF transfer is not
absent; two large opposite pieces cancel). H-C ruled out by construction. H-D: the s shift is worth 5.1 mHa to water
at RHF, in the atom's direction, and is what cancels the p-family loss.
