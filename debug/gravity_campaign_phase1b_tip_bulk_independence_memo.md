# Gravity Campaign — Task 1 (tip/bulk independence) verdict

**Date:** 2026-05-29
**Type:** Diagnostic synthesis (no new computation; combines G4-6b + G4-5d existing data).
**Verdict:** **SEPARABLE, structurally.** The conical-tip / replica entropy coefficient extracts independently of the non-converging bulk Weyl A-coefficient. Direction F (full-spectral-action north star, multi-month bulk A-convergence) is **AXED**. Campaign proceeds on Direction R (entropy/tip north star).

## The gate

Does the entropy coefficient (B_substrate, the thing that becomes +1/6 in S = A/4) extract independently of the bulk Weyl A-coefficient (1/24π), which does NOT converge under single-axis substrate refinement (G4-6a thread-8 clean negative)?

## The two load-bearing facts (from existing data)

**(1) G4-6b — B is clean, measured where A is negligible.** B_substrate is read at large t (t=10), where the bulk contribution A/t ≈ 0.0013 is <1% of B ≈ 0.17. Result: B = 0.163 vs continuum 1/6, −2.3%, R-independent at R≥10; Richardson → B_∞ = 0.173 (+3.96%, finite-a discretization). The only entanglement ever seen was a *joint-fit artifact* (fitting A and B together on the small-t panel is ill-conditioned; the UV undershoot in A poisons both). Resolution: don't joint-fit — measure B at large t, extract A from small-t residuals. The contamination is B→A (A is the casualty), not A→B.

**(2) G4-5d — they are different Mellin moments (the structural reason).** Sector-by-Mellin-moment map:
- entropy / topological tip ↔ **φ(0)** (zeroth, unweighted moment)
- Einstein–Hilbert (the A-coefficient region) ↔ φ(1)
- cosmological constant ↔ φ(2)

The tip term tip(t) = (dK/dα)|_{α=1} − K_disk(t) is *bounded* (~0.15, peaks +0.163 at t≈10), formed as the difference of two large ~530 values at small t. **The divergent bulk A-coefficient cancels in that difference by construction**: at α=1 the wedge derivative's bulk piece equals the disk's bulk piece, so they subtract before the entropy is formed. The entropy (φ(0)) and the bulk A (φ(1)/φ(2)) are orthogonal Mellin moments.

## Why this is structural, not a measurement convenience

The entropy is not merely *measurable* in a regime where A is small — it is the φ(0) moment of a quantity (the replica difference) from which the bulk A has *already cancelled*. The bulk A-coefficient's non-convergence is a non-convergence of the absolute K_disk(t) at small t; that absolute value never enters the entropy, because the entropy is the difference dK/dα − K_disk in which the bulk leading-order cancels.

## Honest caveats

- The entropy retains a *mild cutoff dependence* (φ(0) is log-divergent in continuum, regulated by the substrate t-grid bounds). This is Class-1 calibration (expected, per G8 / Paper 51 §10), and is a separate matter from the bulk A *non-convergence*. Cutoff dependence ≠ inheriting the A problem.
- B measured at α=1 only (the entropy is the replica derivative at α=1, so this is the right point). α>1 behavior of B is the separate Möbius/conical thread (G4-6c), non-gating for the entropy.

## Consequence for the campaign

- **Direction F axed.** The only thing that would have forced the full-spectral-action / multi-month bulk-A-convergence route was the entropy inheriting the bulk non-convergence. It does not.
- **R1 (lock the entropy coefficient) is mostly pre-completed** by G4-6b + G4-5d: B = 0.163 (2.3%), R-independent, structurally the φ(0) moment with bulk cancelled. R1 collapses to consolidation.
- **R2 (propinquity convergence proof on the tip/φ(0) term) is the core deliverable** — does discrete S_tip → A/4 with a rate, via the Paper 38/40 machinery.
- **R3** bounds the bulk A-coefficient as a separate, named scope boundary (cosmological-constant-order; bounded, not solved). The α>1 scope-limited closure folds in here.

## Cross-references

- `debug/g4_6b_ir_boundary_first_move_memo.md` — B clean at large t, joint-fit artifact resolution
- `debug/g4_5d_cutoff_dependence_memo.md` — sector-by-Mellin-moment map (φ(0)/φ(1)/φ(2))
- `debug/gravity_campaign_phase1_moebius_sign_memo.md` — Phase 1 (Möbius sign false alarm)
