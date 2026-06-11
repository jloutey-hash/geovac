# Sprint v3.19.0 follow-on thread — Umbrella synthesis (tasks #24-28)

**Date:** 2026-05-29
**Path:** Gravity arc, sprint-scale follow-on thread on the v3.19.0 second-push (G4-5 Möbius + Mellin map + UV bracketing + IR cure + G4-6 scoping). Five tasks executed single-thread in ~70 minutes (vs ~1.5 week estimate at task creation).
**Verdict:** **MIXED — three POSITIVE, one PARTIAL, one MIXED — with two substantive structural corrections to v3.19.0.**

## 1. Thread summary

| # | Task | Effort | Verdict | Substantive content |
|---|------|--------|---------|---------------------|
| 24 | F14 20-pt Mellin panel | 5 min | **POSITIVE strict-5%** (was PARTIAL 47.9%) | Sector-wise Mellin moment map closed at sub-5% across 3 cutoffs × 3 Λ |
| 25 | α > 1 validation at α ∈ {4, 5, 10} | 5 min | **POSITIVE-EMPIRICAL-LOCK** | Möbius extrapolates at mean 1.71%; α=10 lands at -0.032% (essentially exact) |
| 26 | Fursaev-Solodukhin literature grounding | 30 min | **MIXED — v3.19.0 mechanism FALSIFIED** | hep-th/9512134 is Preitschopf's "Octonions"; Möbius has no published-literature support; mechanism OPEN |
| 27 | Geometric-mean azimuthal | 30 min | **PARTIAL-OVERSHOOT with substantive structural finding** | F6 bit-exact, edge ratio 2/π verified, GM is genuine bracket interior at every t but UV bracket too wide for sub-percent closure |
| 28 | Wedge-spectral-density UV target | 5 min | **POSITIVE-CLOSED-FORM-IDENTIFIED** | Per-t UV target is **1/(24πt)** from Dowker 1977 + Cheeger 1983; v3.19.0 "spectral overshoot" was normalization artifact |

Memo paths:
- `debug/g4_5d_F14_20pt_panel_memo.md`
- `debug/alpha_gt_1_moebius_validation_4_5_10_memo.md`
- `debug/fursaev_solodukhin_1995_grounding_memo.md`
- `debug/g4_5a_geomean_azimuthal_memo.md`
- `debug/wedge_spectral_density_per_t_uv_target_memo.md`

## 2. The two substantive structural corrections to v3.19.0

### 2.1 Task #26: v3.19.0 Track 5 mechanism attribution FALSIFIED

The v3.19.0 sub-agent claimed the Möbius factor F(α) = α/(2α-1) is the "Fursaev-Solodukhin spinor double-cover correction" with claimed reference arXiv:hep-th/9512134. Task #26 literature grounding found:

- **hep-th/9512134 is Preitschopf's "Octonions and Supersymmetry"** — not Fursaev-Solodukhin. Fabricated arXiv ID.
- **The actual spinor-on-cone paper is Fursaev-Miele 1996** (hep-th/9605153, Nucl. Phys. B 484, 697), not Fursaev-Solodukhin 1995.
- **Fursaev-Miele 1996 states spin 1/2 "resembles the scalar case"** — antisymmetric Cheeger-like with no Möbius modification at α > 1.
- The Möbius α/(2α-1) match (empirically robust, locked by task #25 at sub-2%) has **no published-literature mechanism**.

**Action taken:** Paper 51 §12.7.7 mechanism paragraph revised in-place (mechanism OPEN, retracting the Fursaev-Solodukhin attribution). Bibliography extended with Fursaev-Miele 1996 and Dowker 1994. Memory file `memory/alpha_gt_1_moebius_closed_form.md` revised.

**The curve-fit-audit memory rule (`feedback_audit_numerical_claims.md`) worked exactly as designed:** caught a confabulated mechanism claim that the v3.19.0 sub-agent produced without verification. The empirical 6-point Möbius match is preserved as POSITIVE-EMPIRICAL-LOCK; only the mechanism attribution was confabulated.

### 2.2 Task #28: v3.19.0 Track 2 "spectral overshoot" was a normalization artifact

The v3.19.0 Track 2 framing positioned FD (1.31%) and spectral (813%) as bracketing the true per-t UV target with FD undershoot and spectral overshoot. Task #28 derived the actual per-t UV target from Dowker 1977 + Cheeger 1983:

**The continuum prediction is tip^cont(t) = 1/(24πt) + constant + O(t).**

The leading $1/(24\pi)$ coefficient was derived by differentiating Dowker's spinor cone formula $-(1/12)(1/\alpha-\alpha)$ at α=1. The "constant" approaches +1/6 (Lichnerowicz/Seeley-DeWitt) at intermediate t before exponential decay at very large t.

**Empirical verification (re-analyzing v3.19.0 Track 2 + task #27 data):**

| Scheme | tip(t=a²) | vs +1/6 (v3.19.0 framing) | vs true UV target 5.305 |
|---|---:|---:|---:|
| FD | 0.0022 | 1.31% | **0.04%** |
| GM | 0.5932 | 355.9% | **11.2%** |
| Spectral | 1.3556 | 813% | **25.6%** |

**All three discretizations UNDERSHOOT the true UV target.** Spectral is closest at 25.6%, GM intermediate at 11.2%, FD essentially fails at 0.04%. The v3.19.0 "spectral overshoot" was relative to the IR Lichnerowicz baseline +1/6 = 0.167, off by a factor of 31.83 from the true UV target 5.305 at t = a².

This reframes the G4-6c azimuthal-refinement sub-sprint target: **recover the 1/(24πt) UV divergence**, not land at +1/6. Substrate UV refinement (a → a/2 with proportional N_ρ increase) at multi-week compute cost is the structural fix.

## 3. The three POSITIVE closures

### 3.1 Task #24: F14 20-point Mellin panel POSITIVE-strict-5%

8-point panel (G4-5d-refined v3.19.0 Track 1) gave PARTIAL with max deviation 47.9% on the sharp-cutoff channel at Λ ∈ {1, 2}. Two cures applied:
1. Denser t-grid: 20 log-spaced points + 3 sharp-edge anchors at t = 1/Λ² for Λ ∈ {0.5, 1, 2}.
2. Explicit analytical edge insertion at t = 1/Λ² for sharp cutoff (turning log-trapezoidal-across-discontinuity into clean analytical truncation).

Result: max deviation drops 47.9% → **4.82%**, mean 18.2% → **2.67%**. Sector-wise Mellin moment map (tip ↔ φ(0), EH ↔ φ(1), Λ_cc ↔ φ(2)) closed POSITIVE at strict 20% gate, in fact at sub-5%.

Sanity vs G4-5d at t = 1: bit-exact (diff = 0.0).

### 3.2 Task #25: α > 1 Möbius validation POSITIVE-EMPIRICAL-LOCK

Tested v3.19.0 Track 5 Möbius slope = -(1/12)·α/(2α-1) at α ∈ {4, 5, 10} (not in original fit set {1.5, 2, 3}):

| α | rel err |
|---|---:|
| 4 | +2.78% |
| 5 | +2.32% |
| **10** | **-0.032%** |

Mean 1.71% (better than original 2.3% on fit set). The α = 10 result at -0.032% is essentially exact at the asymptote where F(α) → 1/2. The rel err **decreases** with α (3.3% at α=1.5 → 0.032% at α=10), confirming the Möbius asymptote is structurally locked and not a 3-point fit coincidence.

### 3.3 Task #28: Per-t UV target POSITIVE-CLOSED-FORM-IDENTIFIED

Derived 1/(24πt) closed form from standard continuum theory (Dowker 1977 + Cheeger 1983). Empirical fit of tip(t) = A/t + B on small-t data gives A_spec = 0.003017 (22.7% of continuum A_cont = 0.013263) — confirming spectral has the right 1/t shape at intermediate t but undershoots at substrate UV cells due to discrete cutoff effects.

## 4. The PARTIAL finding (task #27)

Geometric-mean azimuthal discretization gives **F6 bit-exact at α = 1** (max_diff = 0.0), edge ratio 0.631 matches predicted 2/π = 0.637 (0.8% slack from finite N_φ). GM is **genuinely a bracket interior at every t** (FD ≤ GM ≤ Spec throughout the trace). However at t = a² the GM recovery is 355.9% relative to +1/6, outside the [50%, 200%] gate window — PARTIAL-OVERSHOOT.

In light of task #28's correction, the [50%, 200%] gate itself was based on the IR Lichnerowicz baseline +1/6. Against the true UV target 5.305, GM lands at 11.2% — undershoot. The PARTIAL verdict against the gate is structurally correct, but the deeper structural finding (GM is a genuine bracket interior at every t, F6 bit-exact at α=1) survives intact.

**Cheap geometric-mean cure for sub-percent UV closure ruled out.** Substrate refinement (a → a/2) is the real fix.

## 5. Implications for G4-6

The v3.19.0 G4-6 scoping memo named four target structural closures (i)-(iv). This thread updates the status:

| Target | v3.19.0 status | Post-thread status |
|---|---|---|
| (i) Subleading O(a²/r_h²) UV | Open, multi-substrate Richardson extrapolation | **Reframed:** target is 1/(24πt) UV divergence recovery (task #28 closed-form). Multi-substrate refinement is the right structural direction. |
| (ii) Subleading O(r_h²/R²) IR | Open, boundary regularization | Unchanged. |
| (iii) α > 1 structural asymmetry | Closed by v3.19.0 Möbius | **Reinforced** by task #25 lock at α ∈ {4, 5, 10}. Mechanism reframed as OPEN (task #26). |
| (iv) Spectral azimuthal discretization | Open, DST/Fourier replacement | **Partial closure** by task #27 (GM as cheap interim) + task #28 (correct target identified). True closure needs substrate refinement, not just spectral azimuthal. |

**G4-6 estimate after this thread: 3-6 months** (unchanged from v3.19.0). G4-6a (multi-substrate UV foundation, the sequential gate) is now the load-bearing next sprint.

## 6. Honest scope

- **Closed at quantitative-rate level (NOT theorem-grade):** F14 Mellin map (#24), Möbius α > 1 extrapolation (#25), per-t UV target formula (#28).
- **Closed at "literature grounding" level:** task #26 — fabricated citation falsified, mechanism for Möbius REMAINS OPEN.
- **Partial verdict with substantive structural finding:** task #27 — GM is genuine bracket interior, F6 bit-exact, but UV bracket too wide for sub-percent closure.
- **Structural sketch:** none in this thread. All findings are either numerical observation or closed-form derivation from standard continuum results.
- **Named open follow-ons:**
  - Mechanism identification for Möbius α/(2α-1) (G4-6 sub-sprint)
  - Multi-substrate Richardson extrapolation a → a/2 to recover 1/(24πt) UV target (G4-6a, the sequential gate)
  - Wedge-spectral-density derivation at theorem-grade rigor (currently sketch-level; G4-6e)
  - Fursaev-Miele 1996 §III spinor calculation PDF-level verification (currently relying on search-result abstract quote)

## 7. What was modified

- **Production code:** none.
- **Tests:** none.
- **Papers:** Paper 51 §12.7.7 mechanism paragraph revised in-place (task #26 — retract Fursaev-Solodukhin attribution, frame mechanism as OPEN). Bibliography extended with two bibitems (`fursaev_miele1996`, `dowker1994`). Page count unchanged (25 → 25, three-pass clean).
- **Memory:** `memory/alpha_gt_1_moebius_closed_form.md` revised (task #26 — flag mechanism as OPEN, log v3.19.0 retraction; task #25 — log empirical validation closure).
- **Hard prohibitions:** none violated. No Paper 2 changes, no fitted parameters, no geometry-hierarchy modifications.

## 8. Substantive new content (durable)

Four findings worth keeping beyond this sprint:

1. **Sector-wise Mellin moment map at sub-5%** (task #24) — empirical closure of the v3.18.0 / v3.19.0 structural reframing of G8. The moment map (tip ↔ φ(0), EH ↔ φ(1), Λ_cc ↔ φ(2)) is now empirically locked across 9 channels (3 cutoffs × 3 Λ).

2. **Möbius α > 1 closed form empirically locked** (task #25) — six-point match (3 fit + 3 validation) at sub-2%, with α=10 essentially exact at -0.032%. Mechanism remains OPEN (task #26) but the empirical identification is structural.

3. **Per-t UV target = 1/(24πt) closed form** (task #28) — derived from standard continuum theory (Dowker 1977 + Cheeger 1983). Replaces the v3.19.0 "+1/6 baseline" reading at small t. Reframes the G4-6c sub-sprint target.

4. **v3.19.0 confabulation caught and corrected** (task #26) — the curve-fit-audit memory rule (`feedback_audit_numerical_claims.md`) successfully caught a fabricated arXiv citation and unsupported mechanism attribution. The discipline worked as designed.

## 9. Files

Drivers + data + memos for all five tasks:
- `debug/g4_5d_F14_20pt_panel.{py, _memo.md}` + `debug/data/g4_5d_F14_20pt_panel.json`
- `debug/alpha_gt_1_moebius_validation_4_5_10.{py, _memo.md}` + `debug/data/alpha_gt_1_moebius_validation_4_5_10.json`
- `debug/fursaev_solodukhin_1995_grounding_memo.md` (no driver/data, literature work only)
- `debug/g4_5a_geomean_azimuthal.{py, _memo.md}` + `debug/data/g4_5a_geomean_azimuthal.json`
- `debug/wedge_spectral_density_per_t_uv_target.{py, _memo.md}` + `debug/data/wedge_spectral_density_per_t_uv_target.json`
- `debug/sprint_v3_19_0_followon_thread_synthesis_memo.md` (this)

## 10. Cross-references

- v3.19.0 sprint (the predecessor): CHANGELOG.md v3.19.0, `debug/sprint_g4_5_parallel_push_2_synthesis_memo.md`
- v3.18.0 sprint (G4-5 first push): CHANGELOG.md v3.18.0, `debug/sprint_g4_5_parallel_push_2_synthesis_memo.md` references
- G4-6 scoping (the multi-month forward plan): `debug/g4_6_scoping_memo.md`
- Paper 51 §12.7.7 (the v3.19.0 paper extension, now revised per task #26): `papers/group5_qed_gauge/paper_51_gravity_arc.tex`
- Curve-fit audit rule: `memory/feedback_audit_numerical_claims.md`
- Diagnostic-before-engineering rule: `memory/feedback_diagnostic_before_engineering.md`
- α > 1 Möbius memory: `memory/alpha_gt_1_moebius_closed_form.md`
- Sector-wise Mellin moment map memory: `memory/sprint_g4_5_sector_mellin_map.md`
