# Pre-registered prediction: a proof that the slice and the chain lie on P (MP1, MP2)

**Registered:** 2026-10-20, before any computation that scores it (PI-directed: "Love proofs. Let's do
it", on the proposal to find the condition that puts G on P).
**Prompted by:** `docs/predictions/2026-10-20_slice_artanh.md`, disposition. P G = 0 was measured to
1e-26, but the analyticity argument fails at t = 1 (P has a log there) and at t = 2 (P keeps a
half-integer exponent). The same holds for the chain.
**Extends:** that record (SA1, SA2) and Paper 59 `sec:density_ring`.
**Commit status:** committed on its own, before any computation that scores it.

## The question

Is there a local condition at one point that cuts Sol(W) down to exactly ker P, and that G satisfies
for a reason we can prove? Then P G = 0 is a theorem, not a measurement.

## What I know before computing (disclosed)

- Derived by hand, nothing computed:
  - J' = 2x K1 + x^2 K0 > 0 and J(0) = 0, J(inf) = 3 pi/2. So 0 <= J <= 3 pi/2. With
    int_0^inf y K1(y) dy = pi/2 and K1 > 0, this gives 0 < G(t) <= 9 pi^3 / (8 t^2) for all t > 0.
  - Chain: chi(x) = int_0^x s K1(s) ds rises from 0 to pi/2, and int_0^inf K0 = pi/2. So
    0 < C(t) = (8 t^3/pi^3) int chi^2 K0(tx) dx <= t^2 for all t > 0.
- From the earlier bound tables (computed, `debug/data/slice_artanh_A.json`): at x = 0, R3 R1 has
  exponents -5/2, -1/2, 1/2, 3/2 for the slice and 1/2, 3/2, 3/2, 5/2 for the chain. Adding theta + 2
  (solution x^-1) or theta - 2 (solution x) gives W's exponents at x = 0, by hand: slice
  {-5/2, -1, -1/2, 1/2, 3/2}, chain {1/2, 1, 3/2, 3/2, 5/2}. Not yet computed for W itself.
- P's exponents at x = 0 have not been computed for either object.

## The argument to be tested

At a regular singular point, the solutions of W that are O(|x|^r |log x|^K) for some K form a subspace
whose dimension is the number of exponents with real part at least r. Take r = -1 for the slice and
r = 1 for the chain. If W has exactly one exponent below r there, that subspace has dimension 4. If
every exponent of P at x = 0 is at least r, then ker P lies inside it, and the two are equal. The
bound puts G (or C) in the subspace, so P G = 0 (or P C = 0) exactly.

## Predictions

**MP1 (slice).** W's exponents at x = 0 are {-5/2, -1, -1/2, 1/2, 3/2}, with -5/2 simple. P's exponents
at x = 0 are {-1, -1/2, 1/2, 3/2}, so -5/2 is the one dropped. Then P G = 0 is proved.
- *Falsifier:* P keeps -5/2 at x = 0. Then the argument fails at x = 0, and the fallback below is used.

**MP2 (chain).** W's exponents at x = 0 are {1/2, 1, 3/2, 3/2, 5/2}, with 1/2 simple. P's exponents are
{1, 3/2, 3/2, 5/2}. Then P C = 0 is proved.
- *Falsifier:* P keeps 1/2 at x = 0.

**Fallback, registered now.** If either fails at x = 0, try x = infinity with the large-t facts
G = O(t^-4) and C = O(t^-?) only if they can be proved by a bound like the ones above. No other point
is tried.

Confidence: moderate. P drops one of W's five solutions, and the dropped one must be large at some
point; there is no reason yet that this point is x = 0.

## Novelty status

- **MP1, MP2: NOT FOUND** in the literature for these objects. But the outcome they prove was measured
  first (SA2 and Paper 59's chain numerics), and the method is the standard boundary-condition argument
  for Fuchsian equations. Decided now: a hit is scored **CHECK**, as a tier upgrade from measured to
  proved, not as a new outcome. A failure is a **MISS**.

## Controls

- **Exponent code:** the exponents of R3 R1 at x = 0 from the same code must reproduce the bound tables
  above, and t^-2 (slice) and t^2 (chain) must be exponents of P.
- **Negative:** P built from the perturbed section (s0 times 1.01) is not a right factor of W, so the
  argument does not apply to it; it is used only as the numerical control already on record.
- **Bounds:** the two inequalities are checked numerically at several t as a sanity check, not as proof.

## Honest caps

- The dimension statement for solutions of bounded growth is standard Frobenius theory. It is used, not
  re-proved.
- The proof covers the one slice (t;1,1) and the chain at the fingerprint fibre only.

## Result

(To be appended after scoring; the predictions above are never edited.)

## Disposition (scored 2026-10-20)

Script `debug/slice_membership_D_x0.py`, data `debug/data/slice_membership_D.json`, log
`debug/slice_membership_D.log`.

**Controls.** The exponents of R3 R1 at x = 0 reproduce the bound tables (slice -5/2, -1/2, 1/2, 3/2;
chain 1/2, 3/2, 3/2, 5/2). P keeps x^-1 (slice) and x (chain). The bounds hold at t = 0.05, 0.3, 1, 3:
G t^2 / (9 pi^3/8) = 0.923, 0.616, 0.242, 0.050, and C / t^2 = 0.805, 0.435, 0.156, 0.033.

**MP1: FAILED.** W at x = 0 is {-5/2, -1, -1/2, 1/2, 3/2} as predicted, but P is {-5/2, -1, -1/2, 1/2}.
P drops the smallest solution, not the largest, so no growth condition at x = 0 gives ker P.
Class: **MISS**.

**MP2: FAILED.** W at x = 0 is {1/2, 1, 3/2, 3/2, 5/2} as predicted, but P is {1/2, 1, 3/2, 3/2}. Same
pattern. Class: **MISS**.

**Fallback: FAILED, by inspection, no computation needed.** P kills t^-2 (slice) and t^2 (chain). These
are the fastest-growing solutions at x = infinity. G = O(t^-4) excludes t^-2, so ker P is not contained
in the set of solutions that obey the bound, and no growth condition at infinity gives ker P either.

**Next.** A parity argument at infinity, found while scoring, is registered separately in
`docs/predictions/2026-10-20_slice_parity_proof.md` before it is computed.
