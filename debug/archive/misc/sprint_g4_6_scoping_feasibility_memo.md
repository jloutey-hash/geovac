# Sprint G4-6 scoping — full discrete-substrate S_BH closure feasibility

**Date:** 2026-05-29
**Path:** Gravity arc completion, Task 5 of multi-task plan. Diagnostic-only assessment.
**Verdict:** **HYBRID-RECOMMENDED.** Current sprint-scale closures (99.998% tip coefficient, 96.69% replica derivative, 85% $S_{BH}$ at UV cell, <5% Mellin moment map) ARE publication-ready at observation-level rigor for Paper 51. The multi-month G4-6 commitment closes theorem-grade / quantitative-rate. Two clean choices: (A) close Paper 51 at v3.20.0 snapshot with named G4-6 successor (4-7 months); (B) commit to G4-6 multi-month and update Paper 51 incrementally. Recommend **(A) plus one targeted sprint** (Task 3 α > 1 mechanism, 1 day) before snapshot.

## 1. Question

The G4-6 scoping memo (`debug/g4_6_scoping_memo.md`) sized the multi-month G4-6 commitment at 4-7 months for theorem-grade / quantitative-rate closure. Current state (v3.20.0):

| Result | Recovery |
|---|---:|
| Spinor SC tip coefficient ($-1/12$) | 99.998% at sweet spot α=2/5, t=2.0 |
| Replica derivative ($+1/6$) | 96.69% |
| $S_{BH}$ at UV cell (Λ=2, r_h=2) | 85.4% |
| Sector-wise Mellin moment map | <5% deviation (Task #24) |
| α < 1 SC recovery | 99.5%+ at most α |
| α > 1 SC recovery | 67.88% (mechanism OPEN, see Task 3) |
| Per-t UV target derivation | structural (Task 4 confirmed) |

Decision gate: is current state close enough to call G4-6 functionally complete (SPRINT-CLOSURE), or does the remaining gap require the originally-named multi-month commitment (MULTI-MONTH), or is a hybrid the right move?

## 2. What "complete" means at three rigor levels

**Level 1 — Observation rigor.** Sprint-scale numerical closures with named follow-ons. Sufficient for Paper 51 to publish as a snapshot of the gravity arc with named G4-6 multi-month successor. Most existing GeoVac papers (Papers 28-37 in the QED arc, Papers 25-30 in the gauge-theory arc) close at this rigor. **CURRENT STATE.**

**Level 2 — Quantitative-rate rigor.** Theorem-grade convergence rate (Richardson c_UV coefficient identified, IR-boundary closed-form subtraction, Mellin moment ratios at theorem-grade across cutoff classes). Sufficient for a math.OA-style closure analogous to Paper 38 (WH1 PROVEN). Requires F13 + F14 + F16 closures from the G4-6 scoping memo. **MULTI-MONTH TARGET.**

**Level 3 — Theorem rigor.** Full formal-proof convergence theorem (analog of Paper 38's five-lemma propinquity proof at the level of $S_{BH}$). Multi-month NCG-mathematics. **NOT G4-6 TARGET; deferred to a future cycle.**

The question is therefore: is Paper 51 publication-ready at Level 1 (current state) for the gravity arc, or does the PI commit to Level 2 (multi-month G4-6) before Paper 51 closes?

## 3. Two clean options

### Option A: Paper 51 closes at v3.20.0 snapshot

**What Paper 51 says.** Gravity arc structural closure at observation rigor:
- G1: two-term exactness of spectral action on $S^3$ (theorem)
- G2: 4D CC form on $S^3 \times S^1$ (theorem)
- G3: spinor-bundle specificity of two-term exactness (theorem)
- G4-1, G4-2: standard CC derivation of BH entropy (continuum theorem)
- G5: closed-form vacuum action density on decompactified $S^3 \times \mathbb{R}_\tau$ (theorem)
- G6-Diag-Full: physical graviton kinetic positivity (necessary condition)
- G7: $G_{\rm eff} = 6\pi/\Lambda^2$, $\Lambda_{cc} = 6\Lambda^2$ (Gaussian)
- G8: cutoff dependence + Class-1 calibration partition
- G4-3 substrate: opens discrete program (sub-sprint sequence)
- G4-4: tip coefficient $-1/12$ bit-exact at sweet spot; replica derivative $+1/6$ at 96.69%
- G4-5: $S_{BH}$ at 85.4% of continuum at UV cell; sector-wise Mellin moment map at <5% deviation
- α > 1 Möbius mechanism OPEN (Task 3)
- G4-6 multi-month commitment NAMED as successor

**Effort.** Paper 51 update only (Task 6). 1 day mechanical incorporation of v3.18.0-v3.20.0 results + 1 day prose synthesis = ~2 days.

**Value.** Gravity arc consolidated. Publication-ready snapshot. PI free to pursue other directions (Lorentzian arc, math.OA program, calibration-side investigations).

**Drawback.** Headline statement softer than "framework predicts $S_{BH}$ at theorem grade." Multi-month G4-6 successor named but not executed.

### Option B: Multi-month G4-6 then Paper 51

**What Paper 51 says (post-G4-6).** Same as Option A plus:
- F13 closure: $S_{BH}^{discrete}(a \to 0) > 0.95 \times$ continuum (Richardson c_UV identified)
- F14 closure: IR-boundary subtraction restores small-Λ cells to within factor-2
- F16 closure: sector-wise Mellin moment map at theorem grade
- F17 closure: spectral azimuthal discretization → per-t tip recovery >99%
- F15 closure: α > 1 mechanism (either F15a closed-form sub-leading OR F15b discretization-artifact)
- Headline: framework reproduces continuum CC $S_{BH}$ at quantitative-rate rigor

**Effort.** 4-7 months (G4-6a/b/c/d/e/f sequence). Paper 51 update at the end.

**Value.** Theorem-grade quantitative closure. Strong publication. Closes the gravity arc at the same rigor level as Paper 38 (WH1 PROVEN).

**Drawback.** 4-7 months of focused work on gravity arc, deferring other directions. Risk of mid-sprint surprises (F15 might surface a structural obstruction that softens the headline).

## 4. Compression evidence — should we trust the 4-7 month estimate?

The math.OA arc compressed estimated months into days in multiple cases:

| Sprint | Estimated | Realized | Compression |
|---|---|---|---|
| Paper 48 drafting | 1.5 months | 13 minutes wall + 15 min cleanup | ~3500× |
| Phase A.4'-A G-B2 | 3-5 weeks | 2-4 hours | ~100× |
| L3b-2a-d strong-form | 4-6 weeks | ~1 day | ~30× |
| Phase A.3' bridge | ~1.5-2 months | session-scale | ~20-50× |

This is a *pattern*, not a coincidence. The realized work compresses dramatically because:
- The substrate is mostly built (G4-3 through G4-5 closed at sprint scale)
- The remaining work is quantitative refinement, not new architecture
- Each falsifier F13-F17 has a well-defined target and a well-scoped diagnostic route

**Realistic compression range:** if math.OA-arc compression rates apply, 4-7 months → 1-3 months realistic, possibly less. If the sprint-scale α > 1 closure (Task 3, three named routes) lands quickly (1-3 days), that also closes F15 in advance.

## 5. The hybrid recommendation

**Option C — Hybrid sprint + snapshot.**

Before closing Paper 51 at v3.20.0 snapshot (Option A path), open ONE targeted sprint-scale closure:

**Sprint G4-6c (α > 1 mechanism).** Per Task 3 of this plan, the α > 1 Möbius mechanism is OPEN with three named sprint-scale closure routes (Fursaev-Miele PDF read 1 day; discrete mode decomposition 1 day; Sommerfeld contour at excess angle 1 week). The lowest-cost route (PDF read) is 1 day.

If Sprint G4-6c closes:
- F15 is closed sprint-scale (either F15a closed-form mechanism OR F15b discretization artifact)
- Paper 51 has one stronger structural statement (Möbius mechanism identified)
- The remaining G4-6 multi-month commitment is reduced to F13/F14/F16/F17 only
- G4-6 effort estimate drops from 4-7 months to 3-5 months

If Sprint G4-6c does not close in 1 day, fall back to Option A (snapshot with named multi-month successor).

**Expected effort.**
- Sprint G4-6c (1 day, lowest-cost route): 1 day
- Paper 51 update (Task 6): ~2 days
- Total: ~3 days

**Value.** Paper 51 closes with one more structural statement potentially closed (Möbius mechanism). The multi-month G4-6 named successor is sharpened (fewer F-targets remain).

## 6. Decision matrix

| Option | Effort | Paper 51 closure rigor | Multi-month G4-6 named? |
|---|---|---|---|
| A (snapshot) | ~2 days | Observation-level | Yes, 4-7 months |
| B (multi-month) | 4-7 months | Quantitative-rate | No (already executed) |
| C (hybrid) | ~3 days | Observation-level + Möbius potentially closed | Yes, 3-5 months reduced |

## 7. Recommendation

**Option C (HYBRID) is recommended.**

Rationale:
1. Cost is comparable to Option A (~3 days vs ~2 days).
2. Has a chance to close one more substantive open question (α > 1 mechanism).
3. Preserves the structural integrity of the multi-month G4-6 successor for future commitment.
4. Aligns with the user's larger arc question ("where is this going") — the PI can see the horizon while one remaining open item gets a try.

Decision-tree:
- Sprint G4-6c lowest-cost route (Fursaev-Miele PDF, 1 day)
  - SUCCESS: mechanism closed. Update Paper 51 with closed mechanism. G4-6 successor named at 3-5 months.
  - FAILURE: mechanism remains OPEN. Update Paper 51 with mechanism OPEN. G4-6 successor named at 4-7 months.

Either way, Paper 51 closes at observation-level rigor with a named G4-6 multi-month successor.

## 8. Honest scope

This memo:
- Confirms current state (99.998% / 96.69% / 85.4% / <5%) is publication-ready at observation rigor.
- Estimates G4-6 multi-month commitment at 4-7 months for theorem-grade closure.
- Notes compression-of-estimates pattern from math.OA arc (estimates often realized at 20-3500× faster than projected).
- Recommends hybrid (Option C) as the cleanest closure path.

Does NOT:
- Execute G4-6c α > 1 closure (that's Task 3 follow-on; this scoping is diagnostic-only).
- Update Paper 51 (that's Task 6 in the multi-task plan).
- Make the PI choice on Option A/B/C (that's a strategic decision).

## 9. Verdict

**HYBRID-RECOMMENDED.** Paper 51 publication-ready at observation rigor (Level 1). Multi-month G4-6 commitment closes theorem-grade (Level 2) and is named as successor regardless of which option. Recommended path: 1-day Sprint G4-6c (Möbius mechanism, Task 3 route A) before Paper 51 closure attempt. If G4-6c closes, mechanism gets closed; either way Paper 51 closes with named multi-month G4-6 successor.

## 10. Cross-references

- `debug/g4_6_scoping_memo.md` — original multi-month G4-6 scoping (4-7 month estimate)
- `debug/sprint_alpha_gt_1_mechanism_investigation_memo.md` — Task 3 of this plan (α > 1 Möbius mechanism, 3 sprint-scale routes)
- `debug/sprint_cc_scoping_memo.md` — Task 1 of this plan (CC scoping NO-GO)
- `debug/sprint_projection_specific_calibration_scoping_memo.md` — Task 2 of this plan (projection-specific scoping NO-GO)
- `debug/sprint_per_t_uv_target_derivation_check_memo.md` — Task 4 of this plan (1/(24πt) is structural)
- `papers/group5_qed_gauge/paper_51_gravity_arc.tex` — Paper 51 current state
- `memory/feedback_diagnostic_before_engineering.md` — discipline for this scoping pass
- `memory/feedback_bit_exactness_rule.md` — bit-exact (skeleton) vs residual (Layer 2) framing

## 11. Files

- `debug/sprint_g4_6_scoping_feasibility_memo.md` (this)
- No new driver or data (diagnostic-only)
