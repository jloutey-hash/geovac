# Pre-registered prediction: a loop of linked densities does not climb the period ladder

**Registered:** 2026-10-17, before any computation in this file (PI-directed: "Yes please",
answering "Want me to write the pre-registration for the loop test?").
**Extends:** `2026-10-16_linking_ceiling.md` (PL1-PL4) and `2026-10-13_three_electron_k3.md`
(P3a-P3e); Paper 59 `sec:three_electron`, `sec:four_density`; atlas doors 17-20.
**Commit status:** to be committed on its own, before any computation, if the PI authorizes it.
If computation starts uncommitted, the Result section must say so.

## The question

Every rung on the ladder so far is a tree: densities joined by kernels with no closed ring.
The tree identity (PL1) works because on a tree every momentum integral collapses. A ring of
kernels leaves one momentum that does not collapse, the loop momentum. Feynman integrals climb
faster with loops than with legs. Does a ring of linked densities climb faster than a chain?

The case that matters for the code is the triangle: three electrons, each pair joined by a
kernel. The LiH explicitly correlated engine carries one such class, the three-body triangle
<f_ik f_jk / r_ij>, the one class that does not separate (CHANGELOG v5.15.13).

## Setup

The Coulomb triangle T = <rho_a(1) rho_b(2) rho_c(3) / (r12 r13 r23)>: three densities,
three Coulomb kernels. In momentum space, with the Fourier convention of Paper 59 eq:momentum,
each kernel carries a momentum k_e, each density sees L_v = the signed sum of the momenta on its
kernels, and L_a + L_b + L_c = 0. Three kernel momenta map onto two independent L's, so one
three-dimensional momentum is left over: the loop.

Fingerprint fibre (the convention of the PL1 3-D check): densities spherical at zero
separation, with transform rho~_v(k) = Delta_v(k)^(-3), Delta_v = sqrt(c_v k^2 + 1),
a_v = c_v^(-1/2). The slower root Delta_v^(-1) is not used: by power counting the triangle then
diverges logarithmically at large momentum (dimension 6 - 3 - 3 = 0), while the chain converges.

## Predictions

**PT1 (the loop integrates out exactly, symbolic).** For Coulomb kernels the loop integral is
the massless one-loop triangle in three dimensions. It is "unique" (three propagator indices 1,
summing to D = 3), so by the star-triangle relation

    int d^3 l / (|l - x1|^2 |l - x2|^2 |l - x3|^2) = pi^3 / (|x1 - x2| |x2 - x3| |x3 - x1|),

with the normalization from pi^(D/2) prod_i Gamma(D/2 - alpha_i) / Gamma(alpha_i), alpha_i = 1.
So

    T = (2 pi^2)^(-3) pi^3 int d^3L_a d^3L_b rho~_a rho~_b rho~_c / (|L_a| |L_b| |L_c|),

over the same two momenta and the same three densities as the three-electron chain. For
spherical densities the magnitude measure 8 pi^2 p q w dp dq dw cancels the weight exactly:

    T = (1/pi) int_{|p-q| <= w <= p+q} rho~_a(p) rho~_b(q) rho~_c(w) dp dq dw.

At the fingerprint the w integral is elementary, int dw (c w^2 + 1)^(-3/2) = w / Delta(w), so

    T = (1/pi) int_0^inf int_0^inf dp dq Delta_a(p)^(-3) Delta_b(q)^(-3)
          [ (p+q)/Delta_c(p+q) - |p-q|/Delta_c(p-q) ].

This is an integral of an algebraic function on the chain's K3 surfaces (P3b: the double cover
of P^2 branched on six lines; the two terms are related by q -> -q), with boundary lines p = 0,
q = 0 and p = q.

**PT2 (main prediction: a loop does not climb the geometry).** Let G(t) be T at scales
(a_a, a_b, a_c) = (1, 1, t), and L3 the chain operator of Paper 59 (4 theta^3 - t^2 (theta+1)^3,
theta = t d/dt, the equation of F_3(t) = int K0(x)^2 K0(tx) dx). Three outcomes are possible.

- (A) G lies in the chain's D-module: G = (R0 + R1 D + R2 D^2) F_3 with rational R_i.
  Same rung, nothing added.
- (B) G = (R0 + R1 D + R2 D^2) F_3 + E, with E a nonzero elementary function (algebraic
  functions times polynomials in logarithms of algebraic functions). Same rung, plus a
  boundary or weight term.
- (C) The minimal operator of G has a right factor of order >= 3 that is not gauge-equivalent
  to L3, for example an order-4 factor with no invariant quadric, as F_4 has. The loop climbs.

**Predicted: (A) or (B), not (C).** Reason: PT1 puts the loop on the tree's base. Hunch, low
confidence: (B), because the triangle region adds boundary lines the chain's 1-D fibre lacks.
Also predicted: the singular points of G's operator lie in {0, +2, -2, infinity}, the chain's.

**PT3 (one-centre limit: a loop is not rational).** Take the one-centre 1s Slater densities
and the exponent triples of P3c / F5, zeta = (1,1,1), (1,2,3), (2,1,3), where the chain is
rational (23/54, 133/144, 289/225). Then the triangle is not rational at (1,2,3) and (2,1,3):
it carries at least one logarithm or inverse tangent of exponent ratios. Reason: with rational
transforms the inner w integral of PT1 leaves an arctan, where the chain's leaves a rational
function. No prediction at (1,1,1), which may collapse. This is a real prediction, not only
prior art. Remiddi's closed form has dilogarithms, but it is for densities carrying extra 1/r_i
factors, and differentiating in the exponents to reach plain 1s densities lowers the weight.

**PT4 (consequence, not tested here).** If PT2 holds, the ladder is counted by the number of
linked densities whatever the ring structure: a ring of three sits with the chain of three, and
it is trees that climb. The LiH triangle uses geminal kernels, whose loop is a massive triangle.
That is not "unique", so PT1's one-line step changes. The 3-D massive triangle is believed
elementary (Nickel 1978, not verified), which would give the same reading.

## Tests and decision rules

1. **PT1.** The 3-D loop integral at five random external-momentum triples against
   pi^3/(|x12||x23||x31|), relative agreement < 1e-8. Then the reduced 2-D formula for T against
   an independent route at one fingerprint scale set and one one-centre Slater triple, < 1e-8.
   The independent route is a direct position-space evaluation, or Remiddi's or Fromm-Hill's
   closed form if the primary source is retrieved and its convention matched.
2. **PT2.** Compute the minimal annihilator of G (series route, order <= 8, polynomial degree
   <= 30, exact or >= 60-digit coefficients). Decide (A)/(B)/(C) by D-module membership against
   F_3 (linear algebra for R0, R1, R2 of degree <= 12), local exponents, singular points, and the
   quadric counter of `debug/linking_ceiling_C_quadrics.py` on any new factor.
3. **PT3.** The one-centre values to >= 60 digits. PSLQ against Q with height bound 10^12 must
   find no relation at (1,2,3) and (2,1,3). Then report, as evidence only, the smallest basis
   among {1, log 2, log 3, pi, arctan values, dilogarithms} that fits.

## Controls (the test is void if one fails)

- The 3-D chain C(t), the Coulomb chain at the same fingerprint, goes through the same pipeline
  and must come out (A) or (B) with respect to F_3. This checks that the pipeline puts a known
  K3-rung object on the K3 rung, and fixes any gauge or normalization freedom before G is read.
- F_4(t) goes through the same pipeline and must come out (C): order 4, no quadric. This checks
  that the pipeline can see a climb when there is one.
- With one kernel dropped (the two-edge case), the momentum bookkeeping of PT1 must reproduce the
  one-centre chain values above, as a check of the conventions.

## Falsifiers

1. PT1: the loop integral or the reduced formula disagrees with its check.
2. PT2: outcome (C) for G, with both controls passing.
3. PT3: the triangle is rational at (1,2,3) or (2,1,3).

## Honest caps, stated before the result

- Fibre statements at zero separation with the fast-decay root, like every rung in Paper 59.
  Whether the Feynman integrations lower the class of the integrated value is PI item 1.
- Coulomb kernels only. The LiH geminal triangle is a consequence (PT4), not a test.
- Only the triangle. Larger rings, and rings with trees hanging off them, are not covered.
- Three-dimensional momenta here, while PL1 used the one-dimensional convention. That the 3-D
  four-density integral reduces into the 1-D family is still assumed (PL caps); PT2's chain
  control is the 3-D comparison at n = 3.

## Prior art, checked before registering

- Remiddi, Phys. Rev. A 44, 5492 (1991): the one-centre three-electron correlation integral with
  Slater functions in closed form, in rational functions, logarithms and dilogarithms (abstract;
  misprints noted by Sims and Hagstrom, 2003). Fromm and Hill (1987) treat the general case, and
  Harris, Phys. Rev. A 55, 1820 (1997) proves the two forms equivalent. Full text not read: the
  form and convention must be checked before any use in test 1.
- The star-triangle (uniqueness) relation in three dimensions: D'Eramo, Peliti and Parisi;
  Symanzik for general D. The pi^3 normalization above is my evaluation of the general formula,
  not a quoted source; test 1 checks it.
