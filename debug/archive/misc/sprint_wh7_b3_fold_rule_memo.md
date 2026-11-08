# Sprint B3 Phase 3, Sprint 3b — the fold-transfer rule in closed form (2026-06-10)

**Goal:** the named Sprint-4 secondary from the Sprint-3 memo: closed-form proof
of which (b, m′) multiplier classes fold to zero / to flow-commuting generators
under the wedge reflection. Exact CG arithmetic throughout (discrete-for-skeleton
rule: Layer-1 claim ⇒ exact arithmetic, no float sweeps as evidence).

**Verdict: CLOSED-FORM-PROVEN — one CG-symmetry lemma + a parity rule; and the
two Sprint-3 headline fold facts are exactly characterized as j_max = 1
WINDOW-EDGE effects. The invariant content survives; the dramatic version
((2,1) annihilated, (2,2) timelike → admissible) does not generalize.**
Driver `debug/wh7_b3_phase3_sprint3b_fold_rule.py` (sympy CG, exact; JSON in
`debug/data/`), falsifier `tests/test_wh7_b3_fold_rule.py` (9/9, 4 s), Paper 45
admissibility passage reworked in place, 25 pp GATE: PASS.

## 1. The lemma (block reflection identity, exact)

For the plain-swap wedge reflection R: m′ → −m′ (HemisphericWedge convention)
and the window compression C^b_{μ′μ} with matrix elements
N·⟨b μ′; j₂ m₂′ | j₁ m₁′⟩⟨b μ; j₂ m₂ | j₁ m₁⟩, N = √((2j₂+1)(2b+1)/(2j₁+1)):

  [R C^b_{μ′μ} R]_(j₁,j₂) = (−1)^{b+j₂−j₁} [C^b_{−μ′,μ}]_(j₁,j₂)

— the CG symmetry ⟨j₁,−m₁′|b,μ′;j₂,−m₂′⟩ = (−1)^{b+j₂−j₁}⟨j₁ m₁′|b,−μ′;j₂ m₂′⟩
read blockwise. Verified exactly (sympy, 50 entries over all seven class
generators); conventions validated against the numeric b1 substrate at 2×10⁻¹⁵.
The exponent b+j₂−j₁ is an integer by the CG integer-ladder rule, so the phase
is always well-defined.

## 2. The parity rule (μ = 0 column)

For Hermitized G = C^b_{μ′,0} + C† the lemma gives each (j₁,j₂) block a
definite R-conjugation parity

  ε(j₁,j₂) = (−1)^{b+j₂−j₁+μ′},

and **the wedge folding annihilates a block iff ε = −1** (an ε = +1 block
folds to its P-restriction — nonzero, Frobenius-partial; the complement lives
in the antisymmetric wedge sector). Verified blockwise-exact for all five
μ = 0 classes; exact rational Frobenius² fold ratios:

| class | blocks (ε) | fold F² ratio | K_W-commutes |
|:------|:-----------|:---:|:---:|
| (1, μ′=0) | (0,1)+, (1,0)+, (1,1)−, (½,½)− | 6/19 | yes |
| (2, μ′=0) | (1,1)+ | 5/6 | yes |
| (2, μ′=1) | (1,1)− | **0** (annihilated) | — |
| (1, μ′=1) null | (0,1)−, (1,0)−, (1,1)+, (½,½)+ | 11/19 | no |
| (2, μ′=2) | (1,1)+ | 1/2 | yes (j≤1 only) |

Mixed classes (μ′, μ both ≠ 0): R maps the (μ′,μ) component to (−μ′,μ),
disjoint from the span of G — no parity eigenstructure, strictly partial folds:
(½,½) → 3/8, (3/2,3/2) → 1/4 (exact). All ratios are rationals — π-free
Layer-1 content, as the skeleton requires.

## 3. The window-edge theorem (corrects the Sprint-3 reading)

Both Sprint-3 headline fold facts are **j_max = 1 edge effects**, proven
exactly at j_max = 3/2:

- **(2,1) revives.** At j ≤ 1 the b = 2 system has only the (1,1) block
  (ε = −1 ⇒ total annihilation). At j_max = 3/2 the half-integer blocks
  (½,3/2) and (3/2,½) enter with ε = +1: the folded generator is nonzero
  with exact fold ratio **8/41**, supported on exactly those blocks.
- **(2,2) stops commuting.** At j ≤ 1 the range forces the single mirror
  transition m′ = −1 → +1 (the only m′ with m′+2 in range), so the folded
  generator is purely mirror and flow-commuting. At j_max = 3/2 non-mirror
  transitions (m′ = −3/2 → ½, −½ → 3/2, …) enter: 8 non-mirror entries,
  [K_W, G_W] ≠ 0, mirror Frobenius fraction exactly **9/25**.

**Invariants (the stable statements):** (i) μ′ = 0 generators fold to
weight-diagonal (flow-commuting, admissible) operators at every tested
window; (ii) the mirror component (transitions m₁′ = −m₂′) is always
admissible; (iii) the per-block parity rule itself is window-independent —
what changes with the window is only which blocks exist. Sprint 3's
"folding reorganizes the cone classes" survives in this per-block form;
the specific annihilation/conversion statements were window artifacts and
Paper 45 has been qualified in place.

## 4. Status of the B3 Phase-3 arc after 3b

Sprints 1–3b complete. The single remaining charter item is **band
exhaustion** (the Mondino–Sämann-shaped limit over increasing bands, on top
of the Sprint-3 bit-exact fixed-band stability) — multi-week, unchanged by
this sprint. No new follow-ons opened.
