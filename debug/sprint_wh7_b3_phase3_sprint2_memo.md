# Sprint B3 Phase 3, Sprint 2 — evenness mechanism, cost-functional verdict, wedge substrate (2026-06-10)

**Goal:** the four named Sprint-2 items from the Sprint-1 memo: (1) ensemble
statistics for the excess-sign caveat; (2) cost-functional comparison (D_max vs
fidelity/Bures vs symmetrized divergences); (3) mechanism of the bimodal p;
(4) wedge restriction as the proper substrate + the operational interval
functional.

**Verdict: MECHANISM-CLOSED — one evenness symmetry explains both the Sprint-1
bimodal scaling AND the sign caveat; no cost functional is sign-stable; the
whole structure transfers bit-level to the HemisphericWedge substrate.**
Driver `debug/wh7_b3_phase3_sprint2.py` (deterministic seeds, JSON in
`debug/data/wh7_b3_phase3_sprint2.json`), falsifier
`tests/test_wh7_b3_phase3_sprint2.py` (14/14; Sprint-1 falsifier 7/7 regression
clean), Paper 45 Q1 extended, 25 pp GATE: PASS.

## 1. The mechanism (T3) — evenness from flow-commutation

The Sprint-1 "bimodal kick-cost scaling" and the frozen sign caveat are one
fact. Let ω_s = U_s ω₀ U_s† be the modular orbit, s2(ε) = e^{iεG} ω_{t/2}
e^{−iεG} the kicked midpoint, c any **unitarily invariant** cost functional.

**If [G, K] = 0** (the kick conserves K-weight), then e^{iεG} commutes with
U_s, and conjugating the two legs gives the exact identity

  c12(ε) = c(ω̃₀(−ε), ω_{t/2}) = c23(−ε),  with ω̃₀(ε) = e^{iεG} ω₀ e^{−iεG},

so the kick excess E(ε) = c12(ε) + c23(ε) − 2·c(ω₀, ω_{t/2}) is an **even
function of ε** — quadratic leading order, for *every* unitarily invariant
cost, regardless of smoothness. Verified bit-level on the m′=0 classes:
[K,G] ~ 10⁻¹⁵, leg-evenness residual 1.2–1.4×10⁻¹⁴ across the full two-sided
ladder ε ∈ ±[10⁻⁴, 0.4], D± < 10⁻³. Crucially the per-leg first-order terms do
**not** vanish (ρ-slot slope −0.380 and +0.056 on the two m′=0 classes;
analytic perturbation formula = numerical quotient to 3×10⁻⁴) — the selection
is **cancellation between legs**, not per-leg vanishing. The Sprint-1 suspicion
("D_max max-eigenvalue non-smoothness + top-eigenspace weight selection") is
retired: the real mechanism is smoothness-independent and cost-universal.

**If [G, K] ≠ 0**, the identity breaks at O(ε) and a smooth **signed** linear
term appears (D⁺ ≈ D⁻ ≠ 0 on all five m′≠0 classes; evenness residuals
3×10⁻² – 1×10⁻¹). Two corollaries:

- **The Sprint-1 caveat is explained.** The excess sign at small ε is the sign
  of a generic linear coefficient — reference-state dependent by construction.
  Nothing was wrong with the states; the observable is genuinely signed.
- **The bimodal p was a window artifact.** The ε^{1.2} production fits were
  log-log fits of signed-linear + quadratic mixtures over ε ∈ [0.05, 0.4]
  conditioned on positive excess. (0.5, 0.5) sat in the "quadratic" group only
  because its linear coefficient is small (|D| ≈ 0.027 vs 0.09–0.22); the true
  dichotomy is m′ = 0 vs m′ ≠ 0, witnessed exactly: ‖[K, G_class]‖ = 2|m′|
  identically (1, 2, 2, 3, 4 on the five broken classes, 10⁻¹⁵ on m′=0).

Flow-translation invariance alone (no commutation needed) gives
c12(0) = c23(0) — observed at 6×10⁻¹⁵, including for the asymmetric D_max.

## 2. Cost-functional comparison (T2) — chain status

96-cell panel (Sprint-1 design: 7 classes × {raw, tangent-projected} × 4 ε +
40 random kicks) at the primary configuration, five functionals:

| functional | chain violations | baseline deficit | status |
|:-----------|:----------------:|:----------------:|:-------|
| D_max (Datta) | 0/96 | +2.185 | theorem |
| Bures angle | 0/96 | +0.194 | metric, theorem |
| trace distance | 0/96 | +0.207 | metric, theorem |
| Umegaki | **90/96** | **−0.098** | fails generically |
| Jeffreys | **88/96** | **−0.172** | fails generically |

Sprint-1's 0/96 Umegaki violations were **reference-state luck**: at this
sprint's primary seed the Umegaki baseline deficit is already negative on the
*on-orbit* triple. The Paper 49 Datta-replacement lesson
(`m1_datta_max_divergence_replacement`) reproduces generically on
modular-orbit triples, not only on TICI cocycle triples.

## 3. Ensemble sign statistics (T1) — the distributional closure

504 cells: 24 seeds × 3 perturbation strengths θ ∈ {0.1, 0.3, 0.6} × 7
classes, excess at ε = 0.2, all five functionals. Fraction of cells with
positive excess:

| functional | overall | m′=0 + (0.5,0.5) classes | |m′| ≥ 1 classes |
|:-----------|:-------:|:------------------------:|:----------------:|
| D_max | 0.69 | 0.71–0.86 | 0.60–0.62 |
| Umegaki | 0.90 | 0.96–0.99 | 0.68–0.97 |
| Bures | 0.93 | **1.00 / 1.00 / 1.00** | 0.75–0.97 |
| trace | 0.89 | 0.82–0.99 | 0.69–0.97 |
| Jeffreys | 0.90 | 0.96–0.99 | 0.71–0.97 |

**No functional achieves state-independent positive excess** — and after §1
this is structural, not accidental: on weight-transferring kicks the leading
term is linear with generic sign, for every unitarily invariant cost. The
answer to Sprint-1's "identify the state classes where positivity holds" is
that positivity is **not a state-class property; it is a kick-class
property** (evenness ⇒ leading-quadratic; sign then set by the Hessian).
Two structured observations, frozen as observations (audit rule — panel-wide,
not proven): Bures excess is positive on *every* ensemble cell for the three
low-transfer classes; D_max is the *least* sign-stable functional overall
(0.69) even though its chain inequality is the theorem anchor — chain
soundness and sign stability are independent virtues, and no tested
functional has both.

## 4. Wedge substrate (T4) — the structure transfers, period halves

`for_bisognano_wichmann(n_max=3)` → dim 40 full, dim 20 wedge; unfolded
K_W = diag(two_m_j > 0), **odd positive integers**, so the wedge KMS state at
β = 1 is a genuine Boltzmann state (unlike the full-space delta-like
structure — `modular_hamiltonian.py` docstring). Anchors: KMS 7×10⁻¹⁶,
flow-invariance 3×10⁻¹⁸. **The orbit period is π, not 2π**: weight
*differences* are even and U_π = −1 exactly (6×10⁻¹⁶) — the positive-weight
folding halves the thermal-time period; intervals on the wedge are defined
mod π.

- D_max chain: 0/76 violations (graded kicks + 40 random).
- **Dichotomy transfers bit-level.** Kicks graded by K-weight transfer
  Δw ∈ {0, 2, 4} (only even transfers exist on the wedge):
  Δw = 0 → [K_W, G] = 0, evenness 4×10⁻¹⁴, p = 1.99/2.01/1.99;
  Δw ≥ 2 → ‖[K_W,G]‖ = Δw, evenness broken at O(ε), signed linear terms
  (D⁺ from −0.115 to +0.117 — both signs at one reference state, as the
  mechanism predicts).
- **Operational interval functional.** ℓ̂(ω_a, ω_b) = argmin_τ
  d_tr(σ_τ(ω_a), ω_b): recovers the true flow separation at 3.7×10⁻⁴
  (600-point grid + parabolic refinement) and is additive on ordered triples
  at 1.1×10⁻³. The interval is a **state-pair functional**, not a
  parametrization choice — orbit injectivity (min trace-distance 0.026 on
  (0, π)) is what makes it well-defined.

## 5. Architecture fixed for the rest of Phase 3

Interval = flow parameter (operational, additive, mod π on the wedge). Costs
= penalty layer with a now-proven grading: weight-conserving (admissible)
deformations pay an even/quadratic penalty universally; weight-transferring
deformations produce signed first-order cost changes that no unitarily
invariant functional stabilizes. The chain inequality (D_max / Bures / trace)
remains the universal "detour never costs less" direction.

## Sprint 3 (named)
- Cone-graded admissibility on the wedge: connect the Δw grading to the
  Phase-2 causal classes (needs the multiplier system on the wedge, not just
  matrix-unit kicks); decide whether "admissible = timelike/null class"
  reproduces the even-penalty layer.
- The convergence statement (the Phase-3 prize): n_max → ∞ behavior of the
  interval functional + penalty layer, Mondino–Sämann-shaped.
- Bures transverse-Hessian positivity on low-transfer classes: prove or
  refute the panel-wide observation (§3).
