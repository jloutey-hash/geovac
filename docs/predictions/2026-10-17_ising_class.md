# Pre-registered prediction: the Ising-class moments on the period ladder

**Registered:** 2026-10-17, before any computation in this file (PI-directed: "let's tear the
ising-class", after the ring extensions of `2026-10-17_ring_extensions.md`).
**Extends:** Paper 59 `sec:four_density` (the tree ladder F_n and its operators); seam row J10.
**Commit status:** committed on its own, before any computation.

## The bridge

Paper 59's tree ladder is F_n(t) = int_0^inf K0(x)^(n-1) K0(t x) dx. At t = 1 it is the
Ising-class moment c_{n,0} = int_0^inf K0(x)^n dx of Bailey, Borwein, Broadhurst and Glasser
(J. Phys. A 41, 205203 (2008)). So the ladder says something about these values, and the values
are a test of the ladder that was not built for it.

## What I know before computing (disclosed, from memory, not yet checked against sources)

c_{1,0} = pi/2; c_{2,0} = pi^2/4; c_{3,0} is a rational multiple of Gamma(1/3)^6 / (pi 2^(2/3)),
a complex-multiplication elliptic value; c_{4,1} = int x K0^4 = 7 zeta(3)/8; the n = 5, 6 moments
are tied by Broadhurst to L-values of modular forms of levels 15 and 6 (cited in Paper 59). The
predictions below for n <= 4 are therefore informed. The predictions for n >= 5 and the parity
rule are not.

## Conventions

Rungs as in `2026-10-16_linking_ceiling.md` and `2026-10-17_density_loop.md`: Tate (rational
combinations of products of pi and odd zeta values) < elliptic < K3 < the four-density rung.
"Singular point" means a singular point of F_n's operator from the Mellin route of Paper 59.

## IS1 (main, structural): the parity rule

The finite nonzero singular values of t^2 for F_n's operator are exactly the thresholds of the
banana with n-1 unit masses, t^2 = (n-1-2j)^2 for 0 <= j < (n-1)/2. Hence t = 1 is a singular
point if and only if n is even. Recorded values for n = 2..5: {1}, {4}, {1, 9}, {4, 16}.
Prediction for n = 6, 7, 8: {1, 9, 25}, {4, 16, 36}, {1, 9, 25, 49}. Confidence: high for the
pattern, moderate that no extra true singular point appears at larger n.

Test: build F_n's operator for n = 6, 7, 8 by the Mellin route, check it on F_n by quadrature at
>= 3 values of t to >= 25 digits, and read off the true (non-apparent) singular points and the
local exponents at t = 1.

## IS2 (values, odd n): odd moments stay on their rung

At odd n the point t = 1 is ordinary, so c_{n,0} is a value of the rung's period at a generic
point and does not drop. Predictions: c_{3,0} is elliptic and not Tate (informed); c_{5,0} and
c_{7,0} have no integer relation with the Tate basis of weight <= n (pi^a zeta(b) products) nor
with the n = 3 elliptic constant times Tate numbers, at >= 200 digits.

## IS3 (values, even n, low confidence): even moments drop

At even n the point t = 1 is a threshold, and the value is a period at a singular point, which is
where values degenerate. Predictions: c_{2,0} is Tate (informed). c_{4,0} drops below the
four-density rung: it is an L-value of a weight-4 modular form attached to the conifold at
x = 1/4, and it is not Tate. c_{6,0} drops likewise, to something below its rung. Confidence: low;
"drops" for n = 6 has no sharp operational form beyond the PSLQ bases tried, which must be
stated.

## Tests and decision rules

1. Moments c_{n,0} for n = 1..8 to >= 250 digits by two independent routes (direct quadrature,
   and the Mellin-route recurrence or a series).
2. PSLQ against each stated basis at the stated precision. A "no relation" needs a positive
   control: the same call recovers c_{2,0} = pi^2/4, c_{4,1} = 7 zeta(3)/8, and a planted relation
   of the same height.
3. Literature values used as checks must be verified against the primary source, not memory.
4. Bounded searches state their bounds.

## Falsifiers

1. IS1: a true singular point at t = 1 for odd n, or none for even n, or a singular value off the
   threshold list, with the operator checked on F_n.
2. IS2: a relation between c_{5,0} or c_{7,0} and the stated lower bases, with controls passing.
3. IS3: c_{4,0} found Tate; or no weight-4 L-value relation for c_{4,0} found where the conifold
   reading says one should be (reported as not confirmed, not as falsified, if the L-value
   search is incomplete).

## Honest caps, stated before the result

- Single-scale values only (t = 1). Two-scale moments are Paper 59's main object and are untouched.
- PSLQ negatives are evidence at a precision and basis, not proofs.
- The parity rule is about the operators. The value-level readings IS2 and IS3 are heuristics on
  top of it, and the known n <= 4 values informed them.
