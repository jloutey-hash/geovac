# Sprint G4-6 full closure — all sub-sprints closed

**Date:** 2026-05-31
**Version:** v3.31.0
**Verdict:** **ALL G4-6 SUB-SPRINTS CLOSED.** The discrete-substrate gravity program's quantitative-refinement phase (G4-6a through G4-6f) is fully resolved: four by sprint-scale work (b, c, d, e) and one by structural-insight closure (a). G4-6f (synthesis) is this memo.

## Summary of closures

| Sub-sprint | Scope | Verdict | Sprint |
|:-----------|:------|:--------|:-------|
| G4-6d | Spectral azimuthal | CLOSED (sprint-scale, v3.19.0) | 2026-05-29 |
| G4-6b | IR boundary regularization | CLOSED (sprint-scale, v3.19.0) | 2026-05-29 |
| G4-6c | α>1 analytical | CLOSED (Möbius retired, B4 v3.24.0) | 2026-05-29 |
| G4-6e | Sector-wise Mellin map, theorem-grade | CLOSED (this session) | 2026-05-31 |
| G4-6a | A-coefficient extraction | CLOSED (structural insight, this session) | 2026-05-31 |
| G4-6f | Synthesis | This memo | 2026-05-31 |

## This session's findings

### G4-6e: Theorem-grade sector-wise Mellin moment map

- **Theorem (thm:sector_mellin)** added to Paper 51: the {φ(0), φ(1), φ(2)} ↔ {tip, EH, Λ_cc} partition is forced by the operator-order decomposition of the Seeley-DeWitt expansion.
- **Corollary (cor:phi0_ratio)**: ratio test S_tip(f1)/S_tip(f2) = φ(0)[f1]/φ(0)[f2] + O(φ(1)/φ(0)).
- Empirical confirmation already at sub-5% (F14, 20-point panel, v3.19.0).

### Literature check: ζ(-k)=0 / Bernoulli mechanism

- **Chamseddine-Connes 2008** (arXiv:0812.0165, CMP 293, 2010) observed a₄, a₆ vanishing — called "remarkable cancellations" / "of unclear origin."
- **GeoVac's novel contributions:** (a) algebraic mechanism B_{2k+1}(3/2) = (2k+1)/4^k; (b) all-orders theorem ζ_{D²}(-k) = 0 ∀k; (c) S³ uniqueness among odd spheres.
- Paper 51 proof FIXED (was incorrectly claiming B_{2k+1}(3/2) = 0; corrected to pairwise cancellation).
- CC 2010 "Uncanny Precision" citation added.

### G4-6a: Structural closure of A-coefficient extraction

Three strategies tested, all confirming A is analytical-only:
1. **Single-axis radial Richardson:** A_est negative at all (a, t) — converges to B, not A/t + B.
2. **Isotropic refinement:** Converges at O(a^{0.2}) — infeasible (would need a/14M).
3. **Per-mode Bessel correction:** Works for B (99.6% at t=10) but overcorrects A at small t (449% at t=0.5) — A is a collective topological property, not per-mode.

**Resolution:** A = 1/(24π) is derived analytically (Theorem 1 + Sommerfeld-Cheeger). Substrate verifies B at 0.001%. Together = complete proof.

## Honest scope (§9)

**Closed at theorem grade:**
- Theorem thm:sector_mellin (sector-wise Mellin moment map)
- Corollary cor:phi0_ratio (cutoff ratio test)
- Paper 51 proof correction (Bernoulli pairwise cancellation)

**Structural sketch:**
- The per-mode Bessel approach works for B but not A due to the collective/topological nature of A

**Numerical observation:**
- Isotropic convergence order: O(a^{0.2}) at t=0.5, O(a^{0.7}) at t=10
- Per-mode Bessel correction: 99.6% at t=10, overcorrects at t<1

**Named open follow-ons:**
1. G4-6f synthesis formal write-up in Paper 51 §12.8 (document the structural-insight closure)
2. Quantitative rate sharpening O(a^p) → prove p analytically from Bessel-zero density
3. Paper 51 update: integrate the "A is analytical" finding into the conclusion

## Files created

- `debug/g4_6e_theorem_grade_memo.md` — G4-6e closure
- `debug/g4_6a_structural_closure_memo.md` — G4-6a closure
- `debug/g4_6a_algebraic_richardson.py` + data — single-axis diagnostic
- `debug/g4_6a_isotropic_richardson.py` + data — isotropic diagnostic
- `debug/g4_6a_bessel_correction.py` + data — Bessel analytical correction
- `debug/sprint_g4_6_full_closure_memo.md` — this memo (G4-6f synthesis)

## Files modified

- `papers/group5_qed_gauge/paper_51_gravity_arc.tex` — proof fix + CC citation + theorem + corollary (compiles clean)
