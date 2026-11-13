# Pre-registered prediction: the five-density star reaches rung five (SA)

**Registered 2026-11-13**, before any quantity below is computed. PI-directed (item #3, option A after v6.7.21:
"attainment at n = 5").
**Commit status:** committed on its own, before any computation it scores.
**Rests on:** Paper 63 Prop. tree_reduction and Prop. tree_rung (v6.7.20, v6.7.21); prop:tree_bound (F_n's operator
irreducible for n <= 8).

## Why

Prop. tree_rung bounds the three-dimensional Coulomb tree at the fingerprint weights: every composition factor of
the fibre's operator has order at most n, and the pieces of order n are F_n's operator (even shifts) or an odd
moment's. Whether F_n's piece is actually present (attainment) was computed only for n = 3, 4. This record tests an
argument that would prove it for stars with n odd and n <= 7, at no machine-search cost.

## The argument (to be scored, not assumed)

Star of five on (t;1,1,1,1): S5(t) = (2/pi)^5 t^3 int_0^oo K0(t x) chi(x)^4 x^-2 dx, chi = int_0^x s K1. Write
chi = pi/2 - R, R(x) = int_x^oo s K1 ~ sqrt(pi/2) x^(1/2) e^-x. On [1, oo) the integrand splits into pieces decaying
as e^(-k x), k = 0..4; only R^4 decays as e^(-4x), and [0, 1] gives a function analytic near t = -4. So near t = -4
the fibre's only singular part is t^3 int_1^oo K0(t x) R^4 x^-2 dx, whose leading behaviour is a square-root branch
(t + 4)^(-1/2), monodromy eigenvalue -1. The pieces of order <= 4 in Prop. tree_rung are smooth at t = -4 and the
order-one pieces have constant coefficients, so their local monodromy there is trivial; an extension of such
pieces has only eigenvalue 1. Hence some composition factor has eigenvalue -1 at t = -4, it must be an order-5
piece, and for the star every order-5 piece has an even shift (four leaf factors r^2 K0 and the root's r^2), so by
Prop. tree_rung it is F_5's operator. The same count gives the exponent (4 - n)/2 for a star of n: half-integer
for odd n, an integer (no eigenvalue -1) for even n.

## Predictions (scored by `debug/star5_attainment.py`, exact where marked)

- **SA1 (exact).** phi(L_4) (F_5's operator) has leading coefficient proportional to (t^2 - 4)(t^2 - 16), and its
  local exponents at t = -4 include a half-integer.
- **SA2 (exact).** For c <= 3 and every shift lambda in -6..12, the leading coefficient of phi(L_c(theta - lambda))
  is nonzero at t = +-4.
- **SA3 (measured, the star of five).** On the boundary value t = -(4 - eps) + i0 the branch is the imaginary part
  -pi t^3 int I0((4-eps) x) R^4 x^-2 dx. With I5(eps) = int_1^oo I0((4-eps) x) R(x)^4 x^-2 dx, eps^(1/2) I5(eps)
  tends to C5 = pi^2/(8 sqrt 2) = 0.8724 (from the leading asymptotics), within 5% at eps = 1e-3, approaching
  monotonically as eps decreases through 1e-1, 3e-2, 1e-2, 3e-3, 1e-3.
- **SA4 (exact, n = 7).** phi(L_6) has local exponents at t = -6 including a half-integer.

**Control (measured, the star of four, the even case).** I4(eps) = int_1^oo I0((3-eps) x) R(x)^3 x^-1 dx grows as a
pole: eps I4(eps) tends to C4 = (pi/2)^(3/2) / sqrt(6 pi) = 0.4535 within 5% at eps = 1e-3, and eps^(1/2) I4(eps)
grows without bound. This checks that the measurement distinguishes a square-root branch from an integer exponent;
it says nothing against attainment at n = 4, which Paper 63 computed directly.

## Verdict rules

- **GO:** SA1-SA3 hold and the control behaves as stated. Paper 63 gains attainment for the star of five, and of
  seven if SA4 holds, as [SYMBOLIC PROOF + MEASURED].
- **BORDERLINE:** SA3's limit is nonzero but off C5 by more than 5%, or not monotone; the branch stands, the
  amplitude reading does not.
- **STOP:** SA1 fails (F_5 has no half-integer exponent at -4, so the argument is wrong somewhere), SA2 fails, or
  eps^(1/2) I5 tends to zero.

## Result (scored 2026-11-13, after the record was committed)

Driver `debug/star5_attainment.py`, output `debug/data/star5_attainment.json`.

- **SA1: HELD.** phi(L_4) has leading coefficient -(t-4)(t-2)(t+2)(t+4); its local exponents at t = -4 are
  0, 1, 3/2, 2, 3 (half-integer 3/2).
- **SA2: HELD.** No zero of the leading coefficient at t = +-4 for c <= 3, lambda in -6..12.
- **SA3: MISSED as registered (11.5% at eps = 1e-3, registered 5%); monotone and nonzero.** sqrt(eps) I5 =
  1.39973, 1.28302, 1.14957, 1.03844, 0.97256 at eps = 1e-1, 3e-2, 1e-2, 3e-3, 1e-3, against C5 = 0.87236. The
  approach is slower than predicted (an O(eps^(1/2)) relative correction).
- **SA4: HELD.** phi(L_6) has local exponents 0, 1, 2, 5/2, 3, 4, 5 at t = -6.
- **Control: HELD.** eps I4 = 0.66981, 0.56026, 0.50219, 0.47243, 0.46110 (C4 = 0.45345, 1.7% at 1e-3);
  sqrt(eps) I4 grows from 2.12 to 14.58.

**Verdict: BORDERLINE** under the rules: the branch stands (nonzero, monotone, and the even-n control shows a
pole instead), the amplitude reading at the registered precision does not.

**Post-hoc, not scored.** At eps = 1e-3, 3e-4, 1e-4, sqrt(eps) I5 = 0.97256, 0.92872, 0.90534; a fit linear in
sqrt(eps) extrapolates to 0.87457, +0.25% from C5. For n = 7 the leading coefficients of phi(L_c), c <= 5, are
prod (t^2 - (c - 2j)^2) and have no zero at +-6 for any shift lambda in -8..14 (a shift never changes the top
weight, so this holds for every lambda).

**Post-hoc, not scored (self-check, 2026-11-13).** The reviewer noted that the argument does not need the
eigenvalue -1: pieces nonsingular at t = -(n-1) form connections there, so any non-holomorphic term forces
a c = n-1 piece. The R^(n-1) term gives eps^-3/2, eps^-1, eps^-1/2, -log eps, eps^1/2, eps log eps for
n = 3..8, so Paper 63 states attainment for every star with 3 <= n <= 8 (F_n irreducible there); the
control's pole is the n = 4 case.
