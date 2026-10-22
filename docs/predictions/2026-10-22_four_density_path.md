# Pre-registered prediction: the real three-dimensional path of four densities (FP1-FP3)

**Registered:** 2026-10-22, before any computation in this file (PI-directed: "Let's do those next steps",
item two of the 2026-10-21 overview).
**Prompted by:** Paper 59 `sec:four_density` states that the real three-dimensional integral reduces into
the family of eq:tree "only for the star ... for the path of four densities the step is assumed". FD1-FD3
(`2026-10-21_four_density_star.md`) settled the star. This record tests the path.
**Extends:** `2026-10-21_four_density_star.md`, `2026-10-16_linking_ceiling.md`,
`2026-10-17_density_loop.md`.
**Commit status:** committed on its own, before any computation in this file.

## The question

Take the fingerprint fibre, densities rho~_v = Delta_v^(-3), Delta_v = sqrt(c_v k^2 + 1), and a path
a - b - c - d: four densities joined by three Coulomb kernels in three dimensions, two of them internal.
Is the minimal annihilator of its one-scale slice built only from the operators of the tree ladder,
F_2, F_3 and F_4 (Paper 59 eq:four_density_op, AESZ #58), and pieces of order one, with F_4's piece
present?

## What I know before computing (disclosed)

- The star (FD1-FD3): S = (16/pi^4) t^3 int chi^3 K0(t x) dx/x, chi = int_0^x s K1(s) ds; Mellin
  operator of order 19; minimal annihilator W8 with composition factors F_4, Legendre and two of order
  one, no K3 piece. My order counts ran low twice (GS1: 17 against <= 10; FD1: 19 against <= 11), so
  no order is predicted here.
- By hand, not computed. In position space a density of scale a is proportional to a^3 K0(a r), and the
  potential of a leaf to chi(a r)/r. Put t on the internal density b, scales (a, b, c, d) = (1, t, 1, 1).
  The angular average of 1/|r - r'| is 1/max(r, r'), and r r'/max(r, r') = min(r, r'), so
  P(t) = N t^3 int_0^oo int_0^oo min(r, r') K0(t r) chi(r) K0(r') chi(r') dr dr'
       = N t^3 int_0^oo K0(t r) chi(r) H(r) dr,  H(r) = int_0^oo min(r, r') K0(r') chi(r') dr',
  with H'' = -K0 chi and H independent of t. N is fixed by the derivation that gives 8/pi^3 for the
  chain and 16/pi^4 for the star. The same steps give the chain, t^3 int K0(t r) chi(r)^2 dr, and the
  star, as checks.
- So P is a single Mellin integral against K0(t r), with g = chi H holonomic: theta(theta - 1) H =
  -r^2 K0 chi. The naive dimension of g's solution space is 3 x 8 = 24, against 10 for the star.
- In one dimension the tree identity (eq:tree) makes the holomorphic period shape-independent; in three
  dimensions the shape enters through H, an extra integration.
- Nothing in this file has been computed.

## Model

P(t) as above. The annihilator of g in theta_r from the exact relations of chi, K0 and H (tensor and
composition, as in `debug/generic_slice_A.py` and `tests/test_paper59_generic_slice.py`). The Mellin route
is `mellin_t_ode` with alpha, beta for the weight, as for the star. Exact arithmetic over Q for operators
and maps; mod p only where stated, with multi-prime CRT and an exact check before any piece is used.
Numerical checks by quadrature of theta_t^j applied under the integral, at 32 digits or more; H by its
two one-sided integrals.

## Predictions

**FP1 (the F_4 piece).** The Mellin operator L of P divides exactly into pieces, and one piece of order
four is equivalent, by an exact map rational in t, to F_4's operator (eq:four_density_op).
- Falsifier: no piece of L equivalent to F_4's operator.

**FP2 (P needs it).** The right factor of L below the F_4 piece does not annihilate P: applied to P it is
nonzero at two values of t, by more than 1e-5 of the largest term.
- Falsifier: it annihilates P to 1e-20.

**FP3 (nothing above).** The minimal annihilator of P, found as the shortest right factor of L that
annihilates P (residual at most 1e-25 of the largest term at 32 digits at four values of t, with a
perturbation of one part in 200 in one coefficient giving above 1e-5), has composition factors each of
order one or equivalent to F_2's, F_3's or F_4's operator.
- Falsifier: the minimal annihilator needs a factor of order at least two equivalent to none of
  F_2, F_3, F_4.

## Novelty status

- FP1, FP3: NOT FOUND. Paper 59's targeted search and the FD1 search found no source placing a molecular
  integral on these operators; the nearest related work is Matuzas 2026 (a different four-electron
  integral, on AESZ-330). A hit on both is a **FINDING** that removes the "assumed for the path" cap.
- FP2: PROBABLY KNOWN in kind (eq:tree); a hit is a CHECK.

## Controls

- Positive: the same reduction with one leaf removed gives the chain's registered order-seven operator;
  the code that builds g, applied to chi^3/r in place of chi H, gives the star's order-19 operator.
- Positive: P at one value of t agrees with an independent momentum-space evaluation of the path to 1e-8
  or better, which also fixes N.
- Negative: the perturbation control of FP3.

## Caps stated in advance

One slice, (1, t, 1, 1), t on an internal density. The fingerprint weights, not Slater ones. The
annihilation is measured, not derived. Minimality is bounded, not proved, unless stated otherwise at
scoring.

## Feasibility rule

If the annihilator of g, its Mellin operator or the factorization cannot be completed within the
session's memory, the computation stops and the record says STOP for feasibility, unscored, with the
stage reached. A switch to the leaf slice (t, 1, 1, 1), where the kernel is chi(t r) and g = K0 H, is
allowed only as a new registration.
