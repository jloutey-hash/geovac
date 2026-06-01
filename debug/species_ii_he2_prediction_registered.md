# PRE-REGISTERED PREDICTION — He2 responsiveness test (frozen before compute)

**Registered:** 2026-05-30, before any He2 computation. Confinement-reframe charter §6
bar: predict past the known. This file fixes the falsifier so it cannot be softened
post-hoc (W3 selection-bias discipline).

## The question the 5-system fit CANNOT answer
Responsiveness law (bit-exact, 5 systems): **1 active electron pair → FROZEN**
(entanglement R-independent over the PES); **≥2 pairs → RESPONSIVE** (entanglement
sweeps as the bond forms/breaks). H2/NaH/KH frozen (1 pair); LiH/MgH2 responsive (2 pair).

Confound: every ≥2-pair system in the fit had **INTERACTING** pairs (LiH core+bond
overlap; MgH2 two bonds share the metal). So the data cannot distinguish:
  (A) the law is about pair **COUNT** (≥2 → responsive), vs
  (B) the law is about pair **INTERACTION** (≥2 *interacting* → responsive).

## Why He2 splits them
He2 = two closed He 1s² pairs, van der Waals bound (~1e-7 Ha) = "very unstable" at the
literal extreme. It is the first ≥2-pair system whose two pairs are **NON-interacting**
(two independent closed shells, negligible overlap in the binding region). Count says
responsive; interaction says frozen.

## FROZEN PREDICTIONS (committed)

- **Metric:** S_mode = von Neumann entanglement entropy of the He|He mode-reduced ground
  density matrix (degeneracy-robust ensemble RDM), mirroring `ensemble_mode_entropy` used
  for the 5 fit systems. **Responsiveness = RANGE of S_mode over the R-sweep.**
- **R-sweep:** He–He distance over a wide window spanning the vdW/large-R region down to
  compressed/small R (target window ~ R ∈ [1.5, 8.0] bohr; exact grid recorded at run).
- **Classification thresholds (same spirit as the 5-system discriminator):**
  - FROZEN  if S_mode range over the LARGE-R (vdW) sub-window (R ≳ 4 bohr) < 0.02.
  - RESPONSIVE if S_mode range over that same large-R sub-window > 0.05.
  - (0.02–0.05 = ambiguous band; report as such, no spin.)

- **MY BET (the mechanistic refinement, B):** He2 is **FROZEN at large R** (S_mode range
  over R ≳ 4 bohr < 0.02), and any response appears only on COMPRESSION (small R) and
  TRACKS inter-shell overlap, not pair count. Reading if confirmed: responsiveness law
  refines to "≥2 **interacting** pairs"; the chemistry-wall cure must induce inter-pair
  correlation, not merely add pairs.

- **FALSIFIER of my bet (i.e., pair-COUNT law A wins):** He2 is RESPONSIVE across the
  large-R region (S_mode range over R ≳ 4 bohr > 0.05) even where the shells barely
  overlap. Reading if this happens: pair count alone drives responsiveness; my
  interaction-mechanism story is wrong.

## Honest caveats fixed in advance
- Energy will be wrong (framework cannot bind He2; W1e). We read the ENTANGLEMENT response,
  not binding energy — per the H2 precedent (contaminated energy, clean FROZEN signal).
- If the GS is degenerate, use the ensemble (degeneracy-robust) RDM, never an arbitrary
  eigenvector (the fabrication-incident trap).
- Report S_mode(R) AND inter-shell overlap(R) so "tracks overlap" is checkable, not asserted.
- One run, read back from JSON, no result written before the numbers exist.
