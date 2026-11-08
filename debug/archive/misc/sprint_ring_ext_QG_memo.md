# Test QG: the exponential-geminal triangle on the period ladder (memo)

**Record:** `docs/predictions/2026-10-17_ring_extensions.md`, test QG (committed 8ba60ae5 before
any computation). Scripts `debug/ring_ext_QG_*.py`, data `debug/data/ring_ext_QG_*`. Nothing outside
`debug/` was touched. Nothing was committed.

**Outcome: G-chain at gamma = 1 (exact K3 factor). At gamma = 3/2 the K3 piece is also there, by
its local data and a rigidity argument, not by an explicit map. QG2, the main prediction, is
falsified.** The mass undoes the drop. Y does not sit on the elliptic rung with the Coulomb
triangle. It sits with the three-electron chain. All registered controls pass.

## Object and conventions

Y(t) = <rho_i rho_j rho_k exp(-gamma r_ik) exp(-gamma r_jk) / r_ij>, scales (1, 1, t), t on the
apex k. Fingerprint densities rho~ = (c k^2 + 1)^(-3/2), and in position space rho is proportional
to K0(a r). The scale t is an inverse length (rho_k ~ K0(t r)), as in G_mom of the density-loop
record. Fourier convention F~(k) = int e^{-ik.x} F d^3x, with f~ = 8 pi gamma/(k^2+gamma^2)^2.

## QG1: the massive loop is elementary [SYMBOLIC + MEASURED]

The closed form used is the solid-angle identity. With M_ij = (m_i^2 + m_j^2 + s_ij)/2 and
D = m1 m2 m3 + M12 m3 + M13 m2 + M23 m1:

    int d^3l / prod (|l-x_i|^2 + m_i^2) = pi^2 arctan(sqrt(det M)/D) / sqrt(det M).

It follows from Feynman parameters, since int_simplex (a^T M a)^(-3/2) = Omega/sqrt(det M) and
Omega is the solid angle of the Gram vectors. The (2,2,1) geminal loop is d/dm1^2 d/dm2^2 of this
form, taken exactly with a bi-dual number. Direct 3-D integration was done in spherical coordinates
with Gauss rules, in float64, independently of the closed form. At five random triples, for each of
gamma = 1 and 3/2, the agreement is:

- loop (1,1,1) with masses (gamma, gamma, 0);
- loop (2,2,1) with masses (gamma, gamma, 0);
- loop (1,1,1) with masses (gamma, gamma, gamma), used by the control.

The worst relative difference is **2.5e-14** (gate 1e-8: PASS). The massless limit of the closed
form gives pi^3/(|x12||x13||x23|) to 1.3e-12. Transcendentals: pi^2 is the loop measure
pi^(3/2) Gamma(3/2) (an M1-type volume factor); the arctan is the solid angle (weight one,
elementary).

**Citations.** Nickel, J. Math. Phys. 19, 542 (1978): **NOT VERIFIED**. A web search did not find
the paper, and the AIP URL I tried led to a different paper. The solid-angle formula
(Van Oosterom–Strackee) was used from memory and is **NOT VERIFIED** against its source. The
closed form is backed by the numerical check, not by either citation.

## Two independent routes and the normalisation

- **Route A** is the momentum triangle with the closed-form loop (a 2-D integral for Phi~).
- **Route B** is a "double bubble". The apex kite Phi~(w) is a two-loop self-energy with two
  bubbles, {rho_i, rho_j} and {f, f}, joined by the Coulomb line. Feynman parameters on each bubble
  and the position-space Coulomb kernel give an algebraic double integral:

      Phi~ = 2 gamma^2 int_0^1 int_0^1 du dv sqrt(u(1-u)) v(1-v) K / (M^3 N^5),

  with M = sqrt(1 + u(1-u) z), N = sqrt(gamma^2 + v(1-v) z), z = w^2, and K rational in
  (M, N, (u-v)^2 z).
- The two routes agree to **2e-16 to 4e-16** at three w and both gamma (dps 15).
- The Coulomb member of the same pipeline (gamma = 0, all three kernels Coulomb) reproduces the
  record's G(t) to **4.8e-15** (t = 0.7) and **2.6e-14** (t = 1.6).

Then Y(t) = (1/(4 pi^2)) t^3 int_0^inf z^(1/2) (z + t^2)^(-3/2) phi(z) dz, with phi(z) = Phi~(sqrt z).

**gamma -> 0 control (registered): PASS.** Put Yukawa kernels on all three edges. The gamma = 0
member equals G at both t (numbers above). At gamma = 0.1, 0.03 and 0.01 the deviation Y/G - 1 is
-0.456, -0.177, -0.064 (t = 0.7) and -0.378, -0.139, -0.050 (t = 1.6). It shrinks roughly like
gamma log(1/gamma) [MEASURED].

## The annihilator: route

phi(z) is analytic at z = 0. Its Taylor coefficients phi_n came from Cauchy's formula on |z| = 2;
the nearest singularity is at |z| = 4. Each point is a 2-D Gauss sum.

- gamma = 1: 600 coefficients at dps 330.
- gamma = 3/2: 600 coefficients at dps 330 and again at dps 520. The two agree to the expected
  precision.
- Exact check: phi_0 = 3 pi/16 at gamma = 1 (from the closed-form K at z = 0) matches.
- phi_n lies in pi Q.

Operators were guessed from the series (U - 1 rows plus a random normalisation, then residuals on
at least N/3 check rows). The rational recognition was then re-checked on all 600 coefficients.

The function Y has a log t part at t = 0: Y_log = t^3 sum binom(-3/2, n) phi_n t^(2n) (derived from
the kernel). The monodromy around t = 0 sends Y to Y + 2 pi i Y_log, so **any annihilator of Y
annihilates Y_log**. Y_log's minimal operator is therefore a right factor of Y's.

## gamma = 1: everything exact except the guessed phi-operator

1. **phi** has a 2nd-order operator, t-degree 4 in z [MEASURED: residual 2.9e-155 on all 600
   coefficients; order 1 is impossible from the exponents]:

       15360 th(th+1) + 64 z(128th^2+660th+505) + 16 z^2(72th^2+676th+1487) + 4 z^3(44th+151)
       - z^4(2th+7)(2th+9).

   - True singular points: z = 0, -4, infinity. The points z = 20 and z = -12 are apparent.
   - It is **equivalent to Legendre** (F_2's family) at argument t^2 = -z/4. The exact rational
     map has denominator (t^2-1)^4 (t^2-3)^3 (t^2+5)^3 and numerator degree 20. The operator
     identity holds in exact arithmetic [SYMBOLIC].
   - So the geminal kite is elliptic.
2. **Two routes to Y's operator.** (a) The Mellin route from phi's operator gives an order-6 upper
   bound. (b) Y_log's series gives an order-6 operator independently. Their exact greatest common
   right divisor is

       G4 (order 4, t^2-degree 4):
       -200960(th-5)(th-3)^2(th-1) + 64t^2(th-3)(1648th^3+6992th^2-9043th-1010)
       - 16t^4(932th^4+12704th^3+54000th^2+71790th-45) + 12t^6 th(20th^3+376th^2+1637th+1875)
       + 9t^8 th(th+2)(th+4)(th+6).

   G4 annihilates Y_log to 3e-155 on 600 coefficients. No order-3 annihilator exists up to
   t^2-degree 99.
3. **G4 on Y itself** (the full function, from a 32-digit grid). G4 Y relative to the largest term
   is 2.3e-26, 5.1e-27, 1.0e-27 and 9.3e-28 at t = 0.5, 0.8, 1.3 and 3.0. The ~3e-23 absolute
   residual is explained quantitatively by a 4e-30 small-w truncation constant: Q0(0) c0 matches.
   **No source term.** The quadrature check (h against 2h) is 3e-17 to 1e-15, so about 30 digits.
   The minimal annihilator of Y is G4.
4. **Exact peeling [SYMBOLIC].** G4 = R3 o R1, verified as an exact operator identity.
   - R1 = theta - r has the algebraic solution t^3 (33t^6 - 612t^4 + 5232t^2 - 49088)/(t^2-4)^(9/2).
   - **R3 is equivalent to F_3's operator 4th^3 - t^2(th+1)^3.** There is an exact gauge map A
     (order 2, denominator (9t^4+312t^2-12560)^2) with G4* o A = 0 modulo L3, in exact
     arithmetic; and R3 = R*, with Sol(R) = A Sol(L3). F_3's operator is self-adjoint up to sign.
5. **No split (D12-style, bounded).** f3 is not a right factor of G4: denominators up to
   (t^2-4)^4 (apparent)^2 (33t^6...), numerator degree up to 50. Positive control: the reversed
   composition R1 o R3 has its split found by the same search at (t^2-4)^2 (apparent)(33t^6...),
   degree 30. So the K3 piece is reached through an integral with an algebraic kernel, as in the
   chain.
6. **Composition factors: algebraic (order 1) and K3 (F_3), and nothing higher.** Y needs both:
   the minimal order is 4, Y is not algebraic, and there is no order-3 annihilator.
7. **Quadric counter.** R3 gives 0 and G4 gives 0 (F_3 control: 1; Legendre: 0). This is not a
   contradiction. The constant-coefficient quadric is not invariant under a gauge map with
   derivatives, and R3 reaches F_3 only through A. The exact map is the rung statement.
8. **Singular points:** t = 0, +-2, infinity. At gamma = 1 the three thresholds 2, 2 gamma and
   1 + gamma coincide at t = 2. The apparent points are the roots of 9t^4 + 312t^2 - 12560.

## gamma = 3/2

1. **phi**: an exact integer operator of order 6 and z-degree 23 (coefficients up to 10^30). Its
   residual is 2.6e-345 on 600 coefficients at dps 520 [MEASURED]. Orders 4 and 5 have none up to
   degree 79/65.
   - Singular points: z = -4, -9 = -4 gamma^2, -25/4 = -(1+gamma)^2, -1/4 = -(1-gamma)^2, and 0.
   - Apparent points: a degree-17 factor.
2. **Y_log**: an exact integer operator **L9** of order 9 and t^2-degree 25 (coefficients up to
   10^45). Its residual is 9e-345 on 600 coefficients. Order 8 has none up to degree 43.
   - Leading factor: (t^2-9)^2 (t^2-4)^2 (4t^2-25)(4t^2-1) times an apparent polynomial of degree
     19 in t^2.
   - Local exponents:
     - t = 2: {0,...,6, -1/2, -3/2}.
     - t = 3: {0,...,6, -3/2, -5/2}.
     - t = 5/2 and t = 1/2: {-1, 0,...,7}.
     - t = 0: {-2, 0, 1, 3, 3, 5, 5, 7, 9}.
     - infinity: {0, 2, 3, 4, 4, 6, 6, 6, 8}.
3. **L9 on Y itself**: 6.9e-28, 8.0e-29 and 2.3e-28 at t = 0.8, 1.3 and 1.7. With one coefficient
   perturbed by 1e-6 relative, the control gives 1.3e-10, 6.7e-12 and 5.5e-14. **No source.** At
   t = 4 both give 1e-19, the quadrature floor there, so t = 4 is not counted.
4. **Pieces found.**
   - (a) An exact order-2 right factor R2 [SYMBOLIC]: L9 = M7 o R2 exactly. R2 is reducible, with
     solutions t^-2(88t^2+939) (rational) and t^-2(11t^2-39)(t^2-4)^(-3/2) (algebraic).
     Quadric count 0.
   - (b) An order-2 right factor through the t = 3 = 2 gamma branch (exponent -3/2). Found mod p;
     not identified. It has no order-1 factor through that solution.
   - (c) An order-3 quotient N3, the minimal right factor of L9* through its exponent -1 solution,
     mod p, t^2-degree 34. Its local data match F_3's twisted by t:
     - singular points only t = 0, +-2 and infinity, plus an apparent factor of degree 66;
     - maximally unipotent at t = 0 (exponents -9, -5, -1, two log obstructions);
     - integer exponents at infinity;
     - a pseudo-reflection at t = 2 (exponents 0, 1, 13/2).

     The gamma = 1 control has the same data in the same mod-p pipeline (exponents -5, -3, -1;
     0, 1, 11/2), and there the exact F_3 map exists. A rank-3 Fuchsian system on P^1 minus
     {0, 4, infinity} in x = t^2 with a pseudo-reflection is hypergeometric, and is fixed by its
     eigenvalues at 0 and infinity (Beukers–Heckman rigidity; citation from memory, **NOT
     VERIFIED**). So **N3 is equivalent to F_3's K3 operator** (twist by t, which is rational in t).

     Tier: [MEASURED mod p + rigidity argument]. An explicit map was **not found** within
     denominators (apparent)^(1,2) (t^2-4)^(0,1) and numerator degree up to 146. The same blind
     pipeline finds the map for gamma = 1 (positive control). At gamma = 3/2 the apparent factor has
     degree 66 instead of 10, so the bound is the likely reason.
5. **Not established at gamma = 3/2.** The full composition-factor list is incomplete. Up to four
   orders between R2 and N3 are only partly identified: piece (b), and the rest of the order-6 part
   below N3. I did not prove that no factor above K3 occurs. The local data give no sign of one:
   every exponent is rational, and the only maximally unipotent block found is N3's. **G-climb is
   not excluded at gamma = 3/2.** The quadric count was not run on N3 (it exists only mod p).

## Against the registration

- **QG1: holds.**
- **QG2: falsified.** There is a K3 factor at gamma = 1 (exact) and very likely at gamma = 3/2.
  Outcome **G-chain** at gamma = 1. At gamma = 3/2: G-chain or higher, with G-climb not excluded.
- **Why the reasoning failed.** "An inverse tangent adds weight, not geometry" was the wrong frame.
  The loop is elementary, but the kite it sits in is elliptic: an exact Legendre map at gamma = 1.
  Attaching the apex density is a Stieltjes transform with a Beta-type Mellin multiplier. That turns
  a 2F1(1/2,1/2;1) family into F_3's 3F2(1/2,1/2,1/2;1,1;t^2/4), a Clausen/Hadamard mechanism. For
  the Coulomb triangle the star-triangle step removes the kite structure and gives G a lower rung.
  The massive loop has no star-triangle step.
- **Side prediction "new singular points that move with gamma": confirmed.** At gamma = 3/2 they are
  t = |1-gamma|, 2, 1+gamma and 2 gamma, which coincide at t = 2 for gamma = 1.
- **"Non-split attachment": confirmed at gamma = 1** (bounded, with a positive control).
- **Where the decision rule did not anticipate the outcome.** At gamma = 3/2 the K3 rung is fixed by
  local data plus rigidity, not by the rule's "rational gauge, exact verification of any map found".
  By the rule's letter, gamma = 3/2 is "K3-equivalent, map not exhibited".
- **The chain control (D9-D11 K3 piece):** the scan on the chain's W finds f3 as a quotient (pass).
  The other scan controls also pass: F_3 finds itself, Legendre finds itself, and F_4 finds neither.

## Defects found and fixed during the work

1. Python's max over flint balls returned the wrong element. This made the guesser reject a genuine
   operator (G4) and accept degenerate ones. Fixed by comparing midpoints as mpf.
2. A ball-arithmetic solve gave over-wide radii and false rejections. Switched to approximate solve
   plus residual checks.
3. Spurious "annihilators" appeared at the precision threshold when the fit used nearly all rows.
   Fixed by requiring at least N/3 check rows; genuine operators sit 15+ orders below the spurious
   level.
4. The rational recognition demanded 320-digit agreement when the null vector carries about 250.
5. A |z| = 3.5 Cauchy run for gamma = 3/2 was inaccurate (residual 1e-103 against 1e-155). It is
   kept as `..._rho3.5_UNRELIABLE.json` and not used.
6. The first gamma = 3/2 rung scan was underdetermined (520-term series, 375 columns). Rerun at 1300
   terms.

None of these changes a reported number.

## Secondary placement (t on i): STOP-feasibility

The kite bubble form for t on i agrees with the triangle route to 4e-16 at three p and both gamma.
For phi_1 at gamma = 1, 600 coefficients were computed; the singularity is at |z| = gamma^2 = 1. Its
minimal operator is order 7, z-degree 21 (residual 3e-154). It was not rationalised, and Y_i's
operator was not built. Reason: time; the operators grow fast. Not a result.

## Caps

- Fibre statements only. One geminal shape. Two values of gamma. Apex placement only.
- phi's operators are guessed and verified on 600 coefficients, not derived [MEASURED]. Peeling,
  maps and divisions are exact where stated.
- The gamma = 3/2 K3 identification rests on mod-p data plus rigidity. Its middle factors are open.
- Bounded searches state their bounds above.
