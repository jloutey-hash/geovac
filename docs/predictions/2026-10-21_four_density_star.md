# Pre-registered prediction: the real three-dimensional four-density star (FD1-FD3)

**Registered:** 2026-10-21, before any computation in this file (PI-directed: "let's try to work it
down", the punchlist of the 2026-10-21 high-level review, item "the real three-dimensional four-density
reduction").
**Prompted by:** Paper 59 `sec:four_density` states as a cap that "the real three-dimensional
four-density integral is assumed, not shown, to reduce into the one-dimensional Bessel-moment family of
Eq. (eq:tree); the analogous step is checked only at n = 3". The OBSERVATION that correlated integrals
linking four densities climb past the elliptic rung rests on that step.
**Extends:** `2026-10-16_linking_ceiling.md`, `2026-10-17_density_loop.md`, `2026-10-19_p3e2_widened.md`.
**Commit status:** committed on its own, before any computation in this file.

## The question

Take the fingerprint fibre, densities rho~_v = Delta_v^(-3), Delta_v = sqrt(c_v k^2 + 1), and a star: one
central density joined to three leaves by three Coulomb kernels 1/k^2 in three dimensions. Is the
minimal annihilator of its one-scale slice built only from the operators of the tree ladder, F_2, F_3 and
F_4 (Paper 59 eq:four_density_op, AESZ #58), and pieces of order one, with F_4's piece present?

## What I know before computing (disclosed)

- n = 3, the chain (a star with two leaves): C(t) = (8 t^3/pi^3) int_0^oo chi(x)^2 K0(t x) dx with
  chi(x) = int_0^x s K1(s) ds, scales (1,1,t), t on the central density (`debug/density_loop_B_ode.py`,
  Paper 59 `sec:density_ring`). chi obeys theta chi = x^2 K1 and ((theta-2)^2 - 1 - x^2) theta chi = 0.
  The Mellin operator of C has order seven and divides as N2 R3 R1 (theta-2), with R3 equivalent to F_3's
  operator; W = R3 R1 (theta-2) annihilates C (derived, `2026-10-21_wg_closure.md`).
- By hand, not computed: in position space the central density is proportional to t^3 K0(t r) and each
  leaf's potential to chi(r)/r. So the three-leaf star on the slice (t;1,1,1) is
  S(t) = N t^3 int_0^oo chi(x)^3 K0(t x) x^(-1) dx, with N fixed by the same derivation that gives
  8/pi^3 at n = 3.
- By hand, not computed: chi's solution space is spanned by 1, int x K1 and int x I1. So chi^3 lies in the
  symmetric cube, of dimension 10 = 1 + 2 + 3 + 4. The 4-dimensional block is the symmetric cube of the
  Bessel pair, the block that gives F_4 for int K0^3 K0(t x) dx. I expect the Mellin operator to have order
  at most 11.
- F_4's operator, its irreducibility (Fuchs exponent argument) and its identification with AESZ #58 are
  in Paper 59 `sec:four_density` [SYMBOLIC].
- Nothing in this file has been computed.

## Model

S(t) as above. The annihilator of chi^3 in theta_x is the symmetric cube of chi's operator
(`density_loop_B_ode.sym_power` with N = 3, m = 3). The Mellin route is `mellin_t_ode`, with the weight
x^(-1) absorbed into the shift alpha and E = M[K0], E(sigma + 2) = sigma^2 E(sigma). Factoring, equivalence
maps and checks follow the slice scripts (`debug/p3e2w_slice_*.py`, `debug/slice_artanh_*.py`,
`debug/density_loop_D*.py`). Exact arithmetic over Q for operators and maps; mod-p only where stated.
Numerical checks by quadrature of theta_t^j applied under the integral, at 32 digits or more.

## Predictions

**FD1 (the F_4 piece).** The Mellin operator L of S has order at most 11. It divides exactly into pieces,
and one piece of order four is equivalent, by an exact map rational in t, to F_4's operator
(eq:four_density_op, in t or in x = 1/(4 t^2)).
- Falsifier: no piece equivalent to F_4's operator.

**FD2 (S needs it).** The right factor of L below the F_4 piece does not annihilate S: applied to S it
is nonzero at two values of t, by more than 1e-5 of the largest term.
- Falsifier: it annihilates S to 1e-20 (then S sits lower, on the K3 rung or below).

**FD3 (nothing above).** The minimal annihilator of S, found as the shortest right factor of L that
annihilates S (residual at most 1e-25 of the largest term at 32 digits at four values of t, with a
perturbation of one part in 200 in one coefficient giving above 1e-5), has composition factors each
equivalent to F_4's, F_3's or F_2's operator or of order one.
- Falsifier: the minimal annihilator needs a factor of order at least two that is equivalent to none of
  F_2, F_3, F_4.

## Novelty status

- FD1, FD3: NOT FOUND. Paper 59's targeted search found no source placing a molecular integral on these
  operators. A hit on both is a **FINDING**: it removes the "assumed, not shown" cap for the star. The
  search is repeated at scoring.
- FD2: PROBABLY KNOWN in kind (the 1D tree identity, Paper 59 eq:tree, puts the holomorphic period there);
  a hit is a CHECK.

## Controls

- Positive: the same code with m = 2 and alpha for n = 3 reproduces the chain's registered order-seven
  operator; with K0^3 in place of chi^3 it reproduces F_4's operator.
- Positive: the derivation of N reproduces 8/pi^3 at n = 3, and S at one value of t agrees with an
  independent momentum-space evaluation of the star to 1e-8 or better.
- Negative: the perturbation control of FD3.

## Caps stated in advance

The star only; the four-density path (two internal densities) is not tested, and in three dimensions the
1D argument that the shape does not matter does not carry over. One slice (t;1,1,1). The fingerprint
weights, not Slater ones. The annihilation is measured, not derived.

## Feasibility rule

If the symmetric-cube Mellin operator or its factorization cannot be completed within the session's
memory, the computation stops and the record says STOP for feasibility, unscored, with the stage reached.

## Disposition (appended 2026-10-21, after scoring)

Computed by a sub-agent from this record; scripts `debug/four_density_star_*.py`, data and logs
`debug/data/four_density_star_*`. PM checks: theta - 3 divides every coefficient of L; E2 kills
t^3 (1 + theta) applied to the Legendre period (series to t^80); the final scoring script
`four_density_star_S_w8score.py`, rerun, reproduces the exact divisions (L = N o W8, W9 = Z o W8) and
all residuals.

**Controls.** m = 2 reproduces the chain's registered order-seven operator; K0^3 reproduces F_4's
operator and eq:four_density_op; the normalization is (2/pi)^(leaves+1), 8/pi^3 for the chain and
N = 16/pi^4 for the star; S agrees with an independent momentum-space evaluation to 1.7e-21 (t = 1.6)
and fully at 20 digits (t = 0.7), and the same code with two leaves reproduces the chain to 18 digits.

**FD1: PARTLY.** L has order 19 (t-degree 18), above the registered 11; no source term (L S at most
6e-45 of its largest term at 45 digits, at t = 0.5, 1.3, 2.5, 3.5). L divides exactly as
N11 o R4pp o X1 o E2 o (theta - 3): theta - 3 (solution t^3); E2 = (theta-3)^2 - t^2 theta (theta-2),
equivalent to the Legendre operator of F_2 by f -> t^3 (1 + theta) f; X1 of order one; R4pp of order four,
equivalent to F_4's operator by an exact map of order three (numerators of degree at most 17,
denominator (t^2-1)^2 (t^2-9)^3), found mod p, rebuilt by CRT, verified exactly. The falsifier does not
fire. Novelty search repeated at scoring (2026-10-21): no source places this integral on AESZ #58. A
related prior instance exists: J. Matuzas, "A Calabi-Yau period in a four-electron correlation integral:
modular physical fibre and a hypergeometric endpoint", Zenodo 10.5281/zenodo.23218203 (2026-10-07): a ring
of Coulomb links plus one diagonal with Slater exponents, on AESZ-330 (landing page read; full text not
read). Class: **MISS** for the order count; **FINDING** for the F_4 piece, narrowed by that related work.

**FD2: HELD.** Below the F_4 piece, R4 = X1 o E2 o (theta - 3) leaves 4.0e-3 to 2.19 of the largest term
on S; in the second exact series, B5 leaves 6.3e-3 to 0.24. Class: **CHECK**.

**FD3: HELD (minimality bounded).** W8 = R4pp o X1 o E2 o (theta - 3), order 8, annihilates S: residual
at most 2.4e-45 at 45 digits at four t; a perturbation of 1/200 in the theta^2 coefficient of R4pp gives
1.5e-4 to 0.69. Composition factors: F_4, F_2 (Legendre) and two of order one; no F_3 (K3) factor. Lower
right factors fail on S. Cap: minimality bounded mod p in the box (t^2-1)^4 (t^2-4)^3 (t^2-9)^4
(5t^2+4)^3 (5t^4+2t^2+2)^3, numerators t^-10 to t^74; other order-two subs not enumerated. Class:
**FINDING**. The star skips the K3 rung, and its one-leaf moment is elliptic where the chain's is
algebraic.
