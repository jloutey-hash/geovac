# Sprint B3 Phase 3, Sprint 3 — admissibility settled, band-stability exact, Bures refuted (2026-06-10)

**Goal:** the three named Sprint-3 items from the Sprint-2 memo: (1) cone-graded
admissibility on the wedge (multiplier system, not matrix-unit kicks); (2) the
convergence statement groundwork; (3) Bures transverse-Hessian positivity —
prove or refute.

**Verdict: ALL THREE SETTLED — admissibility is the flow-commutation grading
(the folding REORGANIZES the causal classes); the band-limited penalty
structure is bit-exactly cutoff-INDEPENDENT (convergence is exact, not
asymptotic, for fixed bands); Bures positivity REFUTED. One attribution
correction to Sprint 2 (period-π).**
Driver `debug/wh7_b3_phase3_sprint3.py` (deterministic seeds, JSON in
`debug/data/wh7_b3_phase3_sprint3.json`), falsifier
`tests/test_wh7_b3_phase3_sprint3.py` (14/14; Sprint-1/2 falsifiers 21/21
regression clean), Paper 45 Q1 extended + two in-place corrections,
25 pp GATE: PASS.

## 1. Admissibility ≠ cone grading (T1): the folding reorganizes the classes

Window wedge: plain-swap reflection R: m′ → −m′ (HemisphericWedge convention),
isometry V onto the +1 eigenspace (dim 9, unfolded weights 2|m′| =
{0,0,0,0,1,1,2,2,2}), K_W = V†|K|V. Folding the seven Phase-2 multiplier
classes produces two structural surprises:

- **The (2,1) spacelike class is annihilated** (fold norm < 10⁻¹², the
  Hermitized generator is R-odd). An entire causal class has no wedge shadow.
- **The (2,2) timelike top class becomes admissible.** Upstairs it has maximal
  transfer ‖[K,G]‖ = 4 and the largest linear penalty (median |D⁺| = 0.257);
  on the wedge [K_W, G_W] < 10⁻¹², evenness 2.5×10⁻¹⁵, median |D⁺| = 3×10⁻⁴.
  Mechanism: G connects the weight pair (−4, +4) (in 2m_j units ±2m′ = ±4);
  the m′ → −m′ identification maps both ends to the single wedge weight 4 —
  **maximal weight transfer folds to zero absolute-weight transfer.**

So "admissible = timelike/null class" is **refuted**: the even-penalty
(admissible) sector on the wedge is ker ad_{K_W} — the flow-commutation
grading — which contains the extreme-spacelike (m′=0, q_F = −1) classes AND
the folded timelike top class, while intermediate classes ((0.5,0.5),
(1.5,1.5), null) stay linear. The cone grading (Phase 2) and the admissibility
grading (Sprints 2–3) are genuinely different structures; Phase 2 already
said the cone is "a grading, not a signature," and the state level confirms
it is also not the penalty-parity grading.

**Magnitude organization (12-seed ensembles, full window):** median |D⁺| per
class ordered by Spearman rank against three predictors over the five
transfer>0 classes: transfer r = +0.97 (p = 0.005), q_F r = +0.90 (p = 0.037),
b r = +0.56. Weight transfer organizes at least as well as the causal ratio
(the two predictors are themselves correlated; the one discriminating pair at
fixed transfer 2 — null 0.108 vs (2,1) spacelike 0.063 — hints the cone adds
secondary structure; 12 seeds, anecdotal, not pinned).

## 2. Band-stability is EXACT (T2): the convergence statement collapses

Geovac wedge sweep n_max = 2..5 (`for_bisognano_wichmann`, dims 16→8, 40→20,
80→40, 140→70) with **fixed band-limited data**: reference perturbation and
both kicks supported on the n_fock ≤ 2 wedge labels (identical label set at
every cutoff, asserted).

- **Bit-exact cutoff independence:** E(0.2; Δw=0), E(0.2; Δw=2), and the
  linear coefficient D⁺(Δw=2) agree across all four cutoffs to 10⁻¹⁵–10⁻¹¹
  (the 10⁻¹¹ is the 1/ε₀ amplification of machine rounding). Why exact: in
  the D_max ratio σ^{−1/2}ρσ^{−1/2} the Gibbs normalization 1/Z cancels
  algebraically, and band-limited conjugations act as the identity off the
  band, so the ratio operator is 𝟙 on the bulk ⊕ a band block that depends
  only on band weights and band unitaries — strictly n_max-independent.
- **The type-III shadow:** the naive trace-norm orbit scale decays as
  exactly 0.180522/Z (constant to 10⁻⁶ × Z across all four cutoffs) — the
  wedge KMS state flattens and has no trace-norm limit. Only band-relative
  (D_max-based) objects survive; the D_max-matching interval recovery is
  cutoff-stable at 1.7×10⁻⁴ (grid-limited) at every n_max.
- **Evenness is cutoff-uniform** (10⁻¹⁵ at every n_max) — as it must be: the
  mechanism is an exact symmetry at each finite cutoff, not an asymptotic one.

**Consequence for the Phase-3 prize:** for FIXED band-limited data there is
nothing to converge — the interval-penalty structure is exact at every
cutoff. The genuinely open limit is **band exhaustion** (increasing bands,
n_fock ≤ B with B → ∞), which is where the Mondino–Sämann-shaped statement
lives. That is the remaining Phase-3 target, now sharply posed.

## 3. Bures positivity REFUTED (T3)

Adversarial panel: 2400 cells = 25 seeds × θ ∈ {0.1, 0.3, 0.6, 1.0} ×
t_total ∈ {0.3, 1.0, 2.5} × 8 random flow-commuting kicks (random Hermitian
blocks within K-weight eigenspaces; max ‖[K,G]‖ = 0 exactly) × sampled
ε ∈ {0.05, 0.2, 0.5}.

- Bures excess negative on **574/2400** cells (worst −2.3×10⁻²).
- The Sprint-2 m′=0 **class generators themselves** go negative on
  **134/600** continuity cells once θ and t_total widen beyond the Sprint-2
  panel (θ ≤ 0.6, t = 1.0, ε = 0.2 only).
- D_max contrast on the same panel: 499/2400 negative.

The Sprint-2 "Bures positive on low-transfer classes" observation was a
narrow-panel artifact — downgraded in Paper 45 in place. Standing verdict:
**no tested functional has a sign-definite even-sector penalty.** The
evenness (parity) of the admissible penalty layer is a theorem; its sign is
not. The chain inequality stays the only universal direction.

## 4. Correction to Sprint 2: period-π attribution

Sprint 2 attributed the wedge orbit period π to "positive-weight folding."
Wrong: the **full** Dirac space already has U_π = −𝟙 at every cutoff
(verified 10⁻¹⁵–10⁻¹⁶ at n_max = 2..5), because the all-odd 2m_j spectrum
makes every weight difference even. The period-π is the spinor
spin-statistics grading σ_π(F) = (−1)^{2b}F (Phase 1) restricted to a pure
half-integer sector — upstairs and on the wedge alike, no folding involved.
Paper 45 corrected in place; Sprint-2 falsifier comment fixed (assertion
unchanged and still true).

## Phase-3 state after Sprint 3

The three-layer architecture is now fully verified at finite cutoff with the
penalty layer mechanism-complete: **flow = time order** (intervals =
operational flow parameter, additive, cutoff-stable); **translation seminorm
= metric** (B1); **gradings**: cone = causal type of multiplier classes
(Phase 2), flow-commutation = admissibility/penalty parity (Sprints 2–3) —
two different gradings, related by folding in a now-computed way.

## Sprint 4 (named — the remaining prize)
- **Band exhaustion:** the interval-penalty structure for increasing bands
  B → ∞ at n_max = ∞-side scaling; identify the limit object and rate
  (Mondino–Sämann-shaped convergence). This is the one genuinely open item
  of the Phase-3 charter.
- Secondary: closed-form proof of the fold-transfer rule (which (b, m′)
  classes fold to zero / to commuting; the (2,1) annihilation and (2,2)
  conversion suggest a clean parity formula in (b, m′, m)).
