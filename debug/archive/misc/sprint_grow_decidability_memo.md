# Sprint: the (AA|BB) accidental-zero sweep + the cross-class decidability wall

**Date:** 2026-08-15. **Branch:** work/sparsity-boundary.
**Drivers:** `debug/bet2_aabb_accidental_zero_sweep.py`,
`debug/bet2_nmax3_sample.py`. **Paper:** 58 (§census outlook, corrected here).

Two bets against Paper 58's `g` row, one empirical and one structural. Together
they settle how far the (AA|BB) "counted = true density" result generalizes, and
they correct one over-optimistic sentence in the certified paper.

---

## Bet 2 — the (AA|BB) "no accidental zeros" is robust, not census-specific

Paper 58 DECIDED the (AA|BB) block on ONE config (Z_A=3, Z_B=1, n_max=2, R=3):
195/195 nonzero, 0 accidental, 0 missed-symmetry. This bet asked: general, or a
property of that config?

Two facts, read off the machinery, focus the test: the *permitted* set (195) and
the *Gaunt symmetry-zero* count are both pure angular data, hence charge- and
R-independent. So the **only** quantity that can move is the count of *accidental
radial zeros*. Sharpest hiding places: equal charges (degenerate cross-center
rates — and N2/F2 in the paper's own swap table are homonuclear), rate-coincidence
charges (4:2 puts a 2s_A exponent = 1s_B exponent), and radial nodes (n_max=3).

**Result — CLEAN everywhere tested.** Same Lindemann decision (group by exp rate,
require every rational coefficient vanish), positive control (an M-violating
quartet is decided *zero*) passing in every config:

| axis | swept | accidental zeros |
|:--|:--|:--:|
| charge | 3:1, 1:1, 2:2, 2:1, 4:2 | 0 (all 195 each) |
| R | 1, 2, 3, 4, 5, 5/2, 7/3 | 0 |
| basis | n_max=2 (all 195) + n_max=3 (100 node-bearing sample, 2 configs) | 0 |

~590 independent decisions, aimed squarely at cancellation (equal charge,
rate-coincidence, nodes, special R), zero hits. So the "counted = true density on
(AA|BB)" upgrades from a single data point to a class property.

*Edges:* n_max=3 is a sample (100 of 7425 node-bearing permitted), so corroborating
not exhaustive; still (AA|BB) only; empirical over a grid, not a theorem. But see
Bet 3 — the theorem underneath is functional independence.

---

## Bet 3 — the cross-class DECISION is transcendence-blocked, not "ordinary work"

Went to supply the independence argument Paper 58 says would decide the 80% cross
bulk, and which it calls "ordinary work rather than obstructions in principle."
**It is not ordinary work.** Read the actual seed structure off the closed forms:

- exp terms e^{-lambda R} (lambda rational) — Lindemann handles these;
- **log** content is `log(2), log(5), log(7)` — transcendental *constants*, not
  functions of R — Baker's theorem decides their independence. **Not a wall.**
- **E_1(k R)** functions (k rational) — present in hybrid AND exchange;
- **standalone gamma** — verified `EulerGamma * e^{-p1-p2}/(2 p1 p2)` as a free
  term in the exchange building block, and present at **every** tau
  (`ordered_xi_general`, tau=0..3), so it survives assembly, does not cancel.

The census decides vanishing at a *fixed algebraic R*, which needs the seed
*values* independent. exp (Lindemann) and log-constants (Baker) are fine, but:

1. **E_1-value independence at algebraic points is an open problem** of
   transcendence theory (Schanuel-adjacent) — blocks every cross class.
2. **gamma is not known even to be irrational** — blocks the exchange class
   specifically, on top of the E_1 wall.

**Crux: weight-one =/= decidable.** (AA|BB) was decidable because it is the *one*
class where pi cancels AND no E_1/gamma appears — pure {exp}, Lindemann. The
instant a single E_1 or gamma enters, you leave Lindemann and hit open
transcendence. That is a **categorical** line — (AA|BB) | {hybrid, exchange} — not
the "harder but not categorically different" gradient the paper described. The
paper conflated "weight-one filtration" (true, the seeds don't escalate in weight)
with "decidable" (false — weight-one already contains E_1 and gamma).

**The achievable substitute (unconditional).** The seed *functions*
{1, e^{-lambda R}, E_1(mu R), ln R} are linearly independent over Q(R) — a
*theorem* (Liouville-Rosenlicht: E_1 is the standard non-elementary integral,
transcendental over the exp field), corroborated by a full-rank 7x7 Wronskian
probe. So each cross-class entry is nonzero *as a function of R* unless its grouped
coefficients vanish → any zero is isolated at special R. With the paper's own
Gaussian corroboration (permitted entries generically nonzero) plus Bet 2's clean
grid, the honest cross-class position is: **functionally nonvanishing
(unconditional), generically nonzero (measured), pointwise-decidable only
*conditional* on a Schanuel-type hypothesis** for the E_1 and gamma values.

---

## The strategic question (PI): advance E_1/gamma transcendence, or postulate?

**Do not advance it.** gamma-irrationality has been open since Euler; E_1-value
independence is at the frontier. Neither is a structure-search-with-verifier
problem (the class where this workflow has an edge, e.g. Bet 1's closed form) —
they are deep-Diophantine, no verifier, field-defining-hard. Wrong bet.

**Postulate it — and this is the framework's own discipline, not a cop-out.**
Conditional-on-Schanuel results are standard practice across arithmetic geometry.
More to the point, Paper 18's taxonomy rule is *pin the irreducible transcendental
to its projection and stop chasing it*; E_1 and gamma are already tagged
embedding-tier (E_1 *is* the Neumann kernel Q_0; gamma is E_1's boundary constant).
Chasing their transcendence would violate the corpus's own rule. So: pin them,
state cross-class decidability as *conditional* on a named value-independence
hypothesis, keep the functional-nonvanishing theorem + Bet 2 empirics as the
unconditional floor. Caveat: the E_1 wall is cleanly Schanuel-adjacent, but
standalone gamma is thornier (gamma-irrationality is not a clean Schanuel
corollary), so the exchange class's conditional decision rests on a heavier
hypothesis than the hybrid's — state that, don't gloss it.

---

## Bet 1 — the 3-centre one-electron integral closes at weight 1, γ-free

Companion result from the same session (the "does weight-1 survive a third centre?"
question of build plan §10.4). The 3-centre one-electron integral
⟨χ_Y|−Z_X/|r−X||χ_Z⟩ (all three centres distinct) is the two-electron *exchange*
assembly with electron 2's density replaced by a point charge at X: the η-integral
becomes P_τ(η_X), the ξ-integration collapses to a *fixed* ordered split at ξ_X, and
the prefactor loses one electron's a³·2π. Neumann kernel + per-τ weight are taken
verbatim from the validated `exchange_value`, so no constant is re-derived.

Validated three independent ways: an X-centred spherical-grid 3-D reference (the
1/|r−X| singularity cancelled by the s² Jacobian) — 1s×1s = −0.341962, 2p₀×1s =
−0.186875, 2p₊₁×2p₀ = −0.027107; the numerical Neumann assembly (matches ~1e-7,
τ-sum converged by τ_max=10); and the symbolic closed form (matches, built through
the engine's validated weight-1 moments `finite_power_exp`/`upper_integral`/
`log_shift_moment`).

**Function content {exp, E₁, ln} — weight one, no dilogarithm, γ-FREE.** Sharper
than the two-centre exchange: that class's Euler γ came from its ξ=1 endpoint, and a
source pinned at ξ_X > 1 never reaches it, so γ drops out — a *proper subset* of the
exchange seed set {exp, E₁, ln, γ}. So weight-1 (the arc's central structural
property) survives a third centre.

The object closed is water's **one-body three-centre V_ne block** — the LARGER of
water's two missing capabilities (−6.92 Ha vs the two-body −2.35 Ha, §10.4). Honest
scope: (a) σ=0 built + weight-inspected; σ≠0 predicted identical by the increment-3d
pole-absorption argument, σ=−1 reference in hand. (b) Does NOT solve water — the
two-body 3-centre ERI (T2) remains the genuine wall, and closed-form value is
all-or-nothing at the tensor level, so a Hamiltonian with only the one-body block
closed is neither accurate nor decidable. (c) Decidability still E₁-blocked (Bet 3
wall) but γ-free — one transcendence-wall better than the exchange class.

Drivers: `debug/bet1_three_center_1e_{reference,assembly,symbolic}.py`. Backing:
`tests/test_two_center_eri_aabb.py::test_three_center_1e_closes_weight_one_gamma_free`.
Build plan §10.4 RESULT.

## Disposition

- **Paper 58 corrected** (certified; PI standing authorization to fix on finding a
  problem): the (AA|BB) decided claim now carries the grid robustness; the "ordinary
  work" paragraph is replaced by the categorical-line + functional-substitute +
  conditional-decision picture.
- Backing: `tests/test_two_center_eri_aabb.py` gains a non-census grid-robustness
  check and a gamma-survives-every-tau pin.
- No `g`-row *tier* change (AA|BB stays decided, cross classes stay counted). What
  changed is the *characterization* of why the cross classes are counted: open
  transcendence, not missing labour.
