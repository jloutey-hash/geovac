# Pre-registered prediction: the order-two piece C2b of the correlated triangle at gamma = 3/2 (EC1-EC3)

**Registered:** 2026-10-22, before any computation in this file (PI-directed: "Ok, let's do that next
registration", after the v5.21.5 overview named this as the cheapest next registration).
**Prompted by:** `2026-10-22_correlated_triangle_g32.md`, disposition: "C2b is an irreducible order-two
factor ... It is not identified." Paper 59 `sec:density_ring` carries the same open item.
**Extends:** `2026-10-22_correlated_triangle_g32.md` (GT1-GT3).
**Commit status:** committed on its own, before any computation in this file.

## The question

C2b is the order-two composition factor of L9, the annihilator of Y at gamma = 3/2, stored as
`factor_C2b` in `debug/data/correlated_triangle_g32_series_final.json` (x = t^2, theta = t d/dt,
x-degree 30). Is C2b the Picard-Fuchs operator of a family of elliptic curves?

## What I know before computing (disclosed)

From the GT run (`debug/data/correlated_triangle_g32_7_refine.log`), local data of C2b in x:
- x = 0: exponents 3/2, 7/2, with a log.
- x = 1/4 and x = 25/4: exponents -5, 0, with a log.
- x = infinity: exponents -11, -10, with a log.
- x = 4 and x = 9: exponents 0, 1, no log.
- The other roots of the leading coefficient have not been examined.
- In t the four points are 0, +-1/2, +-5/2, infinity, that is t = |1 - gamma| and 1 + gamma.

Derived before computing, from these data alone:
1. The earlier lead is wrong as stated. The PM note in the GT record called the pattern "Beauville-type".
   The cross ratio of 0, 1/4, 25/4, infinity has orbit {25, 1/25, -24, -1/24, 25/24, 24/25}. Beauville's
   six families have rational cross-ratio orbits only for -1 and -8. And at x = 0 the monodromy is minus a
   unipotent, not a unipotent. So C2b is not a Beauville family.
2. If C2b is elliptic, its j-map J has degree exactly 6. Suppose every singular point other than the four
   above is apparent. Then J has poles only at the four points, every point over j = 0 is ramified to an
   order divisible by 3, and every point over j = 1728 to an order divisible by 2. Riemann-Hurwitz then
   gives deg J = 12 - 6E, where E counts the extra ramification. One fibre is starred, since x = 0 carries
   minus a unipotent, and any quadratic twist keeps the number of minus signs odd. The Euler number
   deg J + 6 (number of stars) is a multiple of 12, so deg J is 6 mod 12. Hence deg J = 6 and E = 1.
3. The family would then be a rational elliptic surface with fibres I_n* at x = 0 and I_a, I_b, I_c at
   1/4, 25/4, infinity, with n + a + b + c = 6, Mordell-Weil rank 1, and a one-parameter moduli. That fits
   the singular points moving with gamma.
4. Pulled back to the t-line, the twist by t is absorbed. The result would be a semistable rational
   elliptic surface with six fibres at t = 0, +-1/2, +-5/2, infinity.

Nothing in this file has been computed.

## Model

- **Stage A:** exact local data of C2b at every root of its leading coefficient and at infinity.
- **Stage B:** monodromy of C2b. Numerical analytic continuation at 60 digits or more. Loops around
  x = 0, 1/4, 25/4 from one base point, with step radii kept below the distance to every root of the
  leading coefficient. Traces of the Fricke words g_i, g_i g_j, g_i g_j g_k and g_1 g_3 g_2.
- **Stage C:** every rational J of degree 6 with poles at the four points. Write J = P^3 / D and
  J - 1728 = c Q^2 / D, with deg P = 2, deg Q = 3 and D = x^n (4x - 1)^a (4x - 25)^b. Solve exactly for
  all ten splittings of 6 into four positive widths.
- **Stage D:** for each candidate J, the Picard-Fuchs operator of the family y^2 = 4X^3 - gX - g,
  g = 27J/(J - 1728), pulled back to t and twisted by t. Its monodromy traces are compared with Stage B.
  For a match, an exact map A of order at most one over Q-bar(t) is sought with C2b o A = 0 modulo that
  operator, by the method of the GT2 map search.

## Predictions

**EC1 (local structure, CHECK).** Every singular point of C2b other than x = 0, 1/4, 25/4 and infinity is
apparent: integer exponents and no log. Each of the four points carries one Jordan block. The monodromy
has determinant 1.
- Falsifier: a further point with a log or with non-integer exponents.

**EC2 (integral monodromy).** The traces of all Fricke words in Stage B are rational integers, to 1e-30.
- Falsifier: one trace more than 1e-20 from an integer at 60-digit precision, with the run at 80 digits
  agreeing. Then C2b is not the period of any family of elliptic curves, over any field.

**EC3 (the elliptic family).** At least one candidate J of Stage C has monodromy traces equal to those of
C2b, and for it an exact map A exists.
- Pass: J and A exhibited; the identity is exact; a 1/1000 change in one coefficient of A breaks it.
- Falsifier: EC1 and EC2 pass, Stage C is complete, and no candidate's traces match C2b's. By point 2 above
  the search is then exhaustive, so C2b is not elliptic.
- UNCLEAR: traces match but the map is not found in the box (t-powers in [-30, 30], numerator x-degree
  at most 120).

Not scored, reported: the widths (n, a, b, c); whether J lies in Q(x) or needs a number field; the
Mordell-Weil generator; how the family behaves as gamma tends to 1.

## Novelty status

- EC2 and EC3: to be searched at scoring. A hit on both is a **FINDING** of moderate weight. It would
  place every piece at gamma = 3/2 as algebraic, elliptic or K3, with no new period type.
- EC1: a CHECK.

## Controls

- **Positive, J-solve:** Legendre's family. Its period operator is gauged by a random order-one operator
  over Q(x). The pipeline must find integral traces, recover j = 256 (l^2 - l + 1)^3 / (l^2 (l - 1)^2)
  from the three-pole version of Stage C, and find the exact map.
- **Positive, Stage C:** an equation P^3 - Q^2 = lambda D planted from a known sextic must be recovered
  with its widths.
- **Negative:** a Heun operator with singular points 0, 1/4, 25/4, infinity, exponent difference zero at
  each, and the accessory parameter set to 7/3. It must fail EC2. By point 1 above no such Heun operator is
  elliptic, whatever its accessory parameter.

## Caps stated in advance

- One gamma value, the fingerprint fibre, the exponential geminal, apex placement, as for GT.
- EC3's exhaustiveness rests on point 2, which rests on EC1.
- Integrality at 1e-30 on finitely many traces is numerical evidence for EC2, not a proof. The exact map
  in EC3 is the proof.

## Feasibility rule

If Stage C or the map search cannot be completed within the session's memory, the computation stops. The
record then says STOP for feasibility, unscored, with the stage reached.

## Disposition

**Computed:** 2026-10-22, after the registration commit 751a89fa. Scripts `debug/c2b_ident_*.py`, data and logs
`debug/data/c2b_ident_*`. The feasibility rule did not fire.

**EC1: PASS (CHECK).** The leading coefficient of C2b in d/dx form is x^2 (4x - 1)(4x - 25) times a squared
quartic, a sextic and a degree-14 factor. At the quartic's roots the exponents are 1, 2; at the sextic's and the
degree-14 factor's roots they are 0, 2. There is no log at any of them, at six primes each. So all three are
apparent. The four points carry one Jordan block each, minus a unipotent at x = 0. Every exponent sum is an
integer, so the monodromy has determinant 1.

**EC2: PASS.** Traces of the Fricke words at 60 and 80 digits:
- g0 = -2, g1 = 2, g2 = 2 (loops about x = 0, 1/4, 25/4);
- g0g1 = g0g2 = 0, g1g2 = -2, g0g1^-1 = g0g2^-1 = -4, g1g2^-1 = 6;
- g0g1g2 = 10, g0g2g1 = 2.
The largest distance to an integer is 9.8e-51 at 60 digits and 4.3e-71 at 80. The two runs agree to 1.9e-51.

**EC3: PASS (FINDING, with EC2).** Stage C solved all ten splittings exactly, by resultants, with numerical roots
of the exact univariate polynomials at 80 digits. It found 73 solutions; 40 are real and 3 are rational:
- widths (1, 2, 2, 1): J1 = (16x^2 + 152x + 25)^3 / (x (4x - 1)^2 (4x - 25)^2);
- widths (2, 1, 1, 2): J2 = (16x^2 - 88x + 25)^3 / (x^2 (4x - 1)(4x - 25));
- widths (2, 1, 1, 2): J3 = 16 (16x^2 - 68x + 25)^3 / (81 x^2 (4x - 1)(4x - 25)).

J1 - 1728 = (4x - 5)^2 (16x^2 - 616x + 25)^2 / (x (4x - 1)^2 (4x - 25)^2), a square over D, as the form requires.
Its one extra critical point is x = -5/4, a simple one, as point 2 of the registration required.

For J1 the pulled-back j-line operator has tr^2/det equal to C2b's tr^2 on every Fricke word, to 4.9e-78 at
50 digits. J2 matches too, to 2.0e-68. Both J's have the widths the traces imply. Over all 73 solutions, at 30
digits, exactly these two match. Every other solution, J3 included, differs by at least 20 in some tr^2
(`debug/data/c2b_ident_D_compare_30.json`).

The exact map. The source operator is the j-line hypergeometric operator pulled back by J1 and multiplied by
x^(1/2) P^(-1/4), P = 16x^2 + 152x + 25. It is (256x^5 - 1344x^4 - 1680x^3 + 500x^2) y'' +
(256x^4 + 640x^3 - 2480x^2) y' + (-352x^2 + 380x + 125) y, with one apparent point, x = -5/4.
A = (N0 + N1 d/dx) / ((4x - 1)^5 (4x - 25)^5 (4x + 5)), with deg N0 = 22 and deg N1 = 23. C2b o A = 0 modulo
the source, exactly; the nullspace has dimension 1; a 1/1000 change in N0 breaks it. The other twist,
S = (4x - 1)(4x - 25), has no map in the same box. A is injective, since the source is irreducible (J is not
constant) and A is a nonzero operator of order one.

The family. With R = (4x - 5)(16x^2 - 616x + 25), P^3 - R^2 = 1728 x (4x - 1)^2 (4x - 25)^2 exactly, and
E_x: Y^2 = X^3 - 3x^2 P X - 2x^3 R has j = J1. Its discriminant is 2985984 x^7 (4x - 1)^2 (4x - 25)^2, with
orders (0, 0, 1) at infinity. The fibres are I1* at x = 0, I2 at 1/4 and 25/4, I1 at infinity: a rational
elliptic surface, with Mordell-Weil rank 8 - 5 - 1 - 1 = 1 by Shioda-Tate. Over t, the e = 1 model has six fibres I2, at t = 0, +-1/2, +-5/2,
infinity. The source solution is a period: x omega_2(E_x) / source = 1.81379936423421785059407825764 at
x = 0.004, 0.01, 0.02, 0.04, which is pi/sqrt 3 to 30 digits. The pi here is the normalization of the
period integral, Paper 34's period projection, not a physical constant.

So C2b's solutions are, through the exact map A, periods of an explicit family of elliptic curves over Q(x).
Every piece of L9 at gamma = 3/2 is algebraic, elliptic or K3.

**Controls.**
- Positive, J-solve: PASS. Legendre's j = 256 (x^2 - x + 1)^3 / (x^2 (x - 1)^2) was recovered exactly at
  widths (2, 2, 2), with its isogenous forms at (1, 1, 4), (1, 4, 1) and (4, 1, 1).
- Positive, Stage E: PASS. The gauged Legendre operator gives integer traces to 1.1e-37, and the map search
  finds the planted order-one map, nullity 1.
- Positive, Stage C: PASS. A planted P was recovered exactly, the only solution.
- Negative: PASS. The Heun operator with accessory parameter 7/3 has non-integer traces, such as
  g0g1 = -8.2001..., and some complex ones. It fails EC2.

**Not scored.**
- Widths (1, 2, 2, 1) for J1. J2 has the swapped widths (2, 1, 1, 2), as a 2-isogeny does. By the trace match it
  has the same projective monodromy, so it is isogenous up to a twist.
- J1 lies in Q(x).
- The Mordell-Weil generator was not computed.
- As gamma tends to 1, the points 1/4 and 25/4 move to 0 and 4 and the surface degenerates; not computed.

**Novelty.** NOT FOUND, in two searches: one for a three-body exponentially correlated integral with an
elliptic Picard-Fuchs piece, one for the j-invariant or discriminant. A FINDING of moderate weight.

**Deviations.**
- Stage C used resultants and numerical roots in place of a Groebner basis. The basis stalled; the elimination
  is exact up to the root-finding, and every rational root was checked exactly.
- The map's denominator includes 4x + 5, the source's apparent point. The registered box named only
  t-powers and the numerator degree. The first search, without that factor, had nullity 0. A Wronskian bound
  shows A can have a simple pole there.
- The registered map box allowed t-powers in [-30, 30]. The map found has none beyond the twist x^(1/2) built
  into the source.

**Caps.** As registered: one gamma, the fingerprint fibre, apex placement. Integer traces are shown on
finitely many words; the exact map is the proof. L9's minimality is bounded, as in GT, and does not affect
this result.
