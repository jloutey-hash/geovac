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

## Result (2026-10-17, same session, after commit 169a6ce8)

### IS1: the operator part holds for n = 2..8 [SYMBOLIC]

With `_t_ode(n)` and `_local_exponents` of `tests/test_paper59_four_density_rung.py` (the operators
for n = 2..5 match `F_OPS` exactly; all are Fuchsian):

| n | finite nonzero singular t^2 | exponents at each | t = 1 singular |
|:-:|:--|:--|:-:|
| 2 | 1 | 0, 0 | yes |
| 3 | 4 | 0, 1/2, 1 | no |
| 4 | 1, 9 | 0, 1, 1, 2 | yes |
| 5 | 4, 16 | 0, 1, 3/2, 2, 3 | no |
| 6 | 1, 9, 25 | 0, 1, 2, 2, 3, 4 | yes |
| 7 | 4, 16, 36 | 0, 1, 2, 5/2, 3, 4, 5 | no |
| 8 | 1, 9, 25, 49 | 0, 1, 2, 3, 3, 4, 5, 6 | yes |

No other finite singular points. The singular values are exactly the thresholds (n-1-2j)^2 and
t = 1 is singular iff n is even, as registered. The pattern of exponents was not registered: at
even n one exponent, (n-2)/2, is doubled (one logarithm); at odd n one is the half-integer
(n-2)/2 (a square-root branch). Owed by the registration: the operators for n = 6, 7, 8 checked
on F_n by quadrature (dispatched).

Observation, not registered: F_n(t) is real-analytic on t > 0, since every t-derivative of the
integrand stays integrable. So at each threshold F_n is the solution without the logarithm or the
branch. At even n, t = 1 is such a point, and c_{n,0} is the value there of the solution analytic
at the singular point.

### IS1 quadrature, IS2, IS3: confirmed for the registered value cases (the conifold attachment untested); parity finds no support at n = 8

Memo `debug/sprint_ising_class_memo.md`; scripts `debug/ising_class_*.py`; data
`debug/data/ising_class_*`. Moments c_{n,0}, n = 2..8 (at n = 1 route B reaches 93-101 digits), by two independent routes (mpmath
tanh-sinh with an own K0; exact log-series on [0,1] plus arb on [1,R]) agree to >= 286 digits.
PSLQ controls: every run recovers c_{2,0} = pi^2/4 and c_{4,1} = 7 zeta(3)/8 and a planted
relation of height bound/100; every relation found re-checked at full precision.

- **IS1: CONFIRMED.** The operators for n = 6, 7, 8 annihilate F_n at t = 0.5, 0.7, 1.3, 1.7 to
  1.1e-42 (dps 45); adding 1 to one coefficient gives 7e-7 to 1e-1. No source term.
- **IS2: CONFIRMED within bounds.** c_{3,0} = 3 Gamma(1/3)^6/(32 pi 2^(2/3)), not Tate (bound
  1e32); it also equals (pi^2/2) L(48.3.e.a, 1). c_{5,0} and c_{7,0}: no relation with Tate,
  Tate times c_{3,0}, or homogeneous bases with five modular L-values (bounds 1e14 to 1e21 for
  n = 5, 1e7 to 1e15 for n = 7, over all bases; 1e14 to 1e17 and 1e7 to 1e9 over the lower
  bases that Paper 59 quotes).
- **IS3: CONFIRMED for the registered value cases.** c_{4,0} is not Tate (bound 1e24). PSLQ finds
  c_{4,0} = 4 pi^2 L(8.4.a.a, 2), f = eta(2 tau)^4 eta(4 tau)^4 of weight 4 and level 8, and
  c_{6,0} = 108 pi^2 L(24.6.a.a, 4), weight 6 and level 24; both to 5e-261. (The disclosed memory above placed n = 6 at level 6; the form is of level 24.) The conifold
  attachment was not tested.
- **Not registered, and the rule did not anticipate it: n = 8.** c_{8,0} has no Tate relation and
  no single L-value relation found (bounds 1e7 to 1e12). So "even moments drop" finds no support as a
  general rule (a bounded search, evidence and not a disproof). The organizing quantity is the size of the period matrix of the Kloosterman
  motive (Lim, Tu, Yu, arXiv:2609.13576): c_{n,0} is a single L-value where that matrix is 1x1,
  n = 1, 2, 3, 4, 6, and sits inside a 2x2 determinant for n = 5, 7, 8. Parity matched only the
  cases it named.

Checked by me, independently of the agent's code: c_{4,0} by direct quadrature and by BBBG's
(pi^4/4) 4F3(1/2,1/2,1/2,1/2; 1,1,1; 1) both equal 4 pi^2 L(f, 2) to 2e-38, with L(f, 2) built
from the q-expansion of eta(2 tau)^4 eta(4 tau)^4 and the Fricke split, sign +1; the same code
at the wrong level misses by 33%. Lim-Tu-Yu exists (title and authors confirmed by search; its
theorem statements were read by the agent, not by me).

Prior art. The values are not ours: c_{3,0} and c_{4,0} (as a 4F3) are Bailey-Borwein-
Broadhurst-Glasser's (verified in their Table 1); the L-value forms of c_{3,0} and c_{4,0} are
collected by Lim-Tu-Yu, who credit them to earlier work, and the numerical c_{6,0} form is theirs. What the ladder adds is the reading: t = 1 is a threshold
of F_n's operator exactly at even n, with one logarithm there, and a square-root branch at the
thresholds of odd n; and the value-level parity heuristic built on it finds no support at n = 8.

Guarded (2026-10-18, its own step): the threshold table and exponents for n <= 8, the n = 6, 7, 8 operators on F_n, and the c_{3,0}, c_{4,0}, c_{6,0} identities, in `tests/test_paper59_four_density_rung.py`, fire-tested 14/14. The 24.6.a.a coefficients the guard uses are LMFDB's; the Kloosterman-route a_p of this run match them for every prime 5 <= p <= 97.

Unidentified: c_{5,0}, c_{7,0}, c_{8,0}. Stopped for feasibility: the weight-3 level-525 form and
the symmetric-square L-function the literature attaches to the n = 5, 7 determinants.
