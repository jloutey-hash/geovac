# Pre-registered prediction: does the artanh term drop from the Slater slice? (SA1, SA2)

**Registered:** 2026-10-20, before any computation that scores it (PI-directed: "Let's do that",
on the proposal to move the period-ladder work from values to whole functions).
**Prompted by:** Paper 59 `sec:density_ring`, the Slater-slice paragraph: "whether that term drops,
as it does for the chain at the fingerprint fibre, was not tested."
**Extends:** `docs/predictions/2026-10-19_p3e2_widened.md` (P3e-3, HELD) and the chain's order-four
result in Paper 59 `sec:density_ring` (guard `tests/test_paper59_density_ring.py`, the
`test_chain_partial_section_and_order_four_factor` and `test_chain_lies_on_the_order_four_factor`
tests).
**Commit status:** committed on its own, before any computation that scores it.

## The question

On the slice, G(t) = int_0^inf x K1(tx) J(x)^2 dx is annihilated by W = R3 R1 (theta + 2), order
five. Sol(W) contains t^-2 and a second solution H of R1 (theta + 2) that carries artanh t and has a
pole of order 3 at t = 1. Is G built without H? In operator terms: is G annihilated by an operator
of order four, so that its minimal annihilator has order four, not five?

## What I know before computing (disclosed)

- Only invariant line of Sol(R1 (theta + 2)) is t^-2, because H picks up a multiple of t^-2 around
  t = 1. R3 is irreducible (equivalent to F_3's operator). So any invariant subspace of Sol(W) that
  contains G has dimension 4 or 5, and dimension 4 needs a 4-dim invariant subspace containing t^-2.
  (Reasoning by hand.)
- Such a subspace exists if and only if there is a rational partial section S of order at most two
  with R1 S = 1 modulo R3 (right division). Then P = (theta + 2) - S R1 (theta + 2) has order four,
  kills t^-2, and W = Q1 P. This is the chain's construction (`debug/nosplit_order4.py`).
- If P exists and every local solution of P at t = 1 is free of logarithms, then ker P is exactly the
  set of solutions of W that are single-valued around t = 1. G is analytic on Re t > 0, so P G = 0
  would follow exactly. This is the chain's argument (`debug/nosplit_membership.py`).
- At the chain's fibre the answer was yes: S1 exists and the chain lies on P.
- G is analytic at t = 1 (the integral converges for all t > 0). Hand estimates: G ~ 9 pi^3/(8 t^2)
  as t -> 0 and G ~ 6 pi / t^4 as t -> infinity. Not used in scoring.
- Nothing about S for the Slater pieces has been computed.

## Model

- Pieces: R1 = (1 - t^2) theta - (1 + 7 t^2), theta + 2, and R3 as hard-coded in
  `tests/test_paper59_chain_slater_slice.py`. All coefficients are even in t, so work in x = t^2,
  delta = x d/dx, theta = 2 delta, as the chain guard does.
- **Section search:** the chain guard's exact machinery (`_split_bounds`, `_split_search`), with
  Top = R3 and Bot = R1, over Q, inside the pole and degree bound set by the exact local exponents of
  R3 R1 and R3. Nothing is computed modulo a prime.
- **Membership:** P G at t = 1/2 and 13/10 by the guard's quadrature of G and its theta-derivatives,
  at 25 digits, divided by the size of the largest term.
- **Exact membership:** local exponents and log structure of P at x = 1.

## Predictions

**SA1 (the section exists).** R1 S = 1 mod R3 has an exact rational solution S within the bound.
Then W = Q1 P exactly, with P of order four killing t^-2.
- *Falsifier:* the exact system is inconsistent within the bound. Then G's minimal annihilator has
  order five and H is present, provided the bound is valid.

**SA2 (G lies on P).** P G / (largest term) < 1e-18 at 25 digits at both t, and P has no logarithm
in its local solutions at t = 1, so P G = 0 holds exactly.
- *Falsifier:* P G is clearly nonzero (> 1e-6 relative), or P has a logarithm at t = 1 and the
  numerical test fails. SA2 is void if SA1 fails.

Confidence: moderate. The reason is the chain analogy. The analyticity argument decides SA2 only once
SA1 holds; nothing known forces SA1.

## Novelty status

- **SA1, SA2: PROBABLY KNOWN** by analogy. Paper 59 shows the same outcome for the Delta^-3 chain at
  the fingerprint fibre. The Slater weights change R1 and R3, so it is not a corollary, but the
  expectation comes from our own result. Searched: Paper 59, the 2026-10-17 and 2026-10-19 records,
  the BFKNS banana review cited in the 2026-10-19 record. Nothing addresses this object.
  A hit is a **CHECK**. A failure is a **MISS**, and it would be a new structural difference between
  the physical and the Delta^-3 weights, recorded as such.

## Controls

- **Positive, exact:** the same search code on the chain pieces must return the chain's recorded S1
  (`S1_DELTA` in the chain guard), and on the chain's planted split it must find a solution.
- **Negative, exact:** S with one coefficient changed by 1% must leave a nonzero remainder.
- **Negative, numerical:** P built from that perturbed S must give P G > 1e-6 relative at both t.
- **Bound check:** the exponent table is printed in full and compared with the hand values at x = 0
  and x = 1 before the search runs.

## Honest caps

- SA1 rests on the bound. A bug in the exponents could shrink it. The table is printed and checked.
- This is one slice, (t;1,1). The generic fibre is not covered.
- A held SA2 says only that H is absent. It does not identify G among the solutions of P.

## Result

(To be appended after scoring; the predictions above are never edited.)
