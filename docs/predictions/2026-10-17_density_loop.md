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

## Result (2026-10-17, same session, after commit ade983ab; scripts `debug/density_loop_{A_pt1_pt3,A2_pt3_exact,B_ode,B2_inhom,C_factor,D_gauge,D2_shared,D4_factors,D5_rational,D6_peel,D7_m5,D8_confirm}.py`, data and logs in `debug/data/density_loop_*`)

**PT1 HOLDS. PT3 HOLDS, stronger than predicted. PT2: the main claim holds (no climb), but the
registered three-way rule did not anticipate the outcome. The triangle lands below the chain's
rung, a fourth outcome. Read the PT2 row for how it falls under the rule as written.**

| Check | Outcome |
|:--|:--|
| PT1a, loop identity | the Feynman-parameter route equals pi^3/(\|x12\|\|x23\|\|x31\|) at five random triples, relative 6e-13 to 1.2e-12 (the quadrature's own level). |
| PT1b, reduced formula | the star-triangle reduction equals Remiddi's closed form (PPR 2004 eq. 7, read from the extracted text and confirmed by this agreement) to 8e-22 and 1.7e-21 at the three one-centre triples, and the position-space form below equals the 2-D momentum formula to 1e-21 at the fingerprint. |
| PT3, one centre | not rational: no relation with Q at height 10^12 at 100 digits, at (1,2,3) and (2,1,3). Stronger than predicted: genuine dilogarithms survive. Exact at (1,1,1): T = 11/3 - pi^2/8 - (3/2) Li2(-1/2) + elementary logs and log products (data file). At (1,2,3), PSLQ at 100 digits: 360 T = 3344 - 180 pi^2 + 1080 [Li2(1/5) + Li2(1/6) + Li2(1/3)] - 3360 ln2 - 525 ln5 + 540 ln^2 2 + 1080 ln2 ln3 + 1080 ln^2 3 - 2160 ln2 ln5 + 540 ln^2 5. The triangle is symmetric in its densities, so (1,2,3) = (2,1,3). Weight two, where the chain at the same triples is rational. |
| Position form | at the fingerprint the triangle is a three-Struve moment, G(t) = (pi/2) t^2 int x^2 psi(x)^2 psi(tx) dx with psi = L_{-1} - I_1, and the chain is C(t) = (8 t^3/pi^3) int chi(x)^2 K0(tx) dx with chi = int_0^x s K1(s) ds. |
| Operators (Mellin route) | triangle: L_G = -4 th (th-3)^2 + t^2 (5 th^3 - 7 th^2 - 5 th + 3) - t^4 th (th+2)(th+3), with L_G G = -3 t^2 exactly (18 digits at three t): an elementary remainder the regularized Mellin transform cannot see. Chain: an order-7 operator, residual 2e-19 and 4e-21. |
| **PT2, the triangle** | L_G is reducible. Its regular solution at t = 0 is algebraic, (1 - t^2/4)^(-3/2). The quotient is gauge-equivalent to Legendre's equation (F_2, the two-density elliptic rung): a relation q y = p0 f2 + p1 theta f2 is found, and the adjoint test finds f2 as a quotient. No relation to F_3's K3 family on either side. Quadric count 0. Composition factors: elementary t^2, algebraic (1 - t^2/4)^(-3/2), elliptic F_2. **No K3 piece and nothing new: the triangle sits on the elliptic rung, below the chain's K3.** Against the registered rule: not (A) and not (B), since G is not in F_3's module; by the letter it is (C), since L_G is an order-3 right factor not gauge-equivalent to L3; but L_G is reducible into lower pieces, so it is not (C) in the sense the record states ("the loop climbs"). Main claim, "not (C)" read as no climb: supported. |
| PT2 side prediction | singular points in {0, +-2, infinity}: **FAILS**. t = +-1 are genuine (exponents 0, 1, 0 at t = 1, with a log); they belong to the F_2 piece. |
| PT4 reading | "a ring of three sits with the chain of three": **not supported**; the ring sits lower. |
| Control, F_4 | no relation to f2 or f3 on either side (with rational coefficients); irreducible order 4 from PL2. Passes: the method sees a new rung as new. |
| Control, 3-D chain | its operator peels exactly as M5 . R1 . (theta - 2): bottom pieces elementary (t^2, from the leaf-free moment t^3 int K0(tx) dx = (pi/2) t^2; R1 has solution t^3/(1-t^2)^2, from the one-leaf moment, whose Mellin transform is a Gamma ratio). **F_3's K3 family is a sub-piece of M5** (rational-gauge relation, confirmed at 900 series terms). M5's top order-2 piece is not F_2 at argument t^2 or t^2/4: unidentified. So the method places the 3-D chain on the K3 rung, as the control requires in substance. The letter "(A) or (B)" is not established, because that order-2 top piece is not shown to be elementary. **Settled later the same day (steps D9-D12, below):** the top piece is not part of the chain, and the letter fails within the searched bounds. |
| Control, bookkeeping | the 2-D chain formula uses the convention of `test_one_centre_limit_is_rational` (passes: 23/54, 133/144, 289/225), and the chain's position and momentum forms agree to 1e-21. |

**Reading.** At the 3-D fingerprint, closing three densities into a ring does not raise the
rung; it lowers it. The loop integrates out exactly (star-triangle). Each density then enters
through psi = L_{-1} - I_1, a Bessel function with a constant source, and the product
structure that gives the chain its K3 piece does not survive: what is left is one elliptic
family, the same one an ordinary two-density repulsion integral carries, plus algebraic and
elementary pieces. At one centre the same ring raises the weight from rationals to dilogarithms
(Remiddi), with no new geometry. So in both limits a loop trades geometry for weight, or keeps
both low; in neither does it climb.

**Defects found and fixed during the computation, recorded so they are not repeated.**
(1) My first Mellin-to-ODE converter divided a common factor out of the theta-polynomials;
theta does not commute with t^2, so the factor must come out of the recurrence in s. The F_3
self-check passed only because F_3 has no common factor. (2) My first F_3 reference series used
argument t^2; the recurrence of 4 th^3 - t^2 (th+1)^3 gives t^2/4. The same slip was in Paper 59
(line 928) and in `2026-10-16_linking_ceiling.md`; both are corrected. (3) A first factor test
allowed only polynomial coefficients; the change of basis needs rational ones (it missed the
triangle's known F_2 piece until fixed). (4) Explicit zero entries in a sparse series tripped a
sanity check. None changes a reported number.

**Caps.** Fibre statements at zero separation with the fast-decay root, like every rung in
Paper 59; the integrated-value question (PI item 1) is untouched. Coulomb kernels only; the LiH
geminal triangle (PT4) is untested here (tested the same day: K3 at gamma = 1, record `2026-10-17_ring_extensions.md`). Only the triangle. The comparison is between the 3-D triangle
and the 3-D chain at the same fingerprint; the 1-D convention of PL1 was not used for the ring,
where it does not apply (the loop momentum is a flat direction of the kernel-free form). The
chain control's no-split result (D12) is bounded in the degrees and denominators of the map.
Annihilator and relation searches are bounded in order and degree, as stated in each script.

**Capture (same day).** Home: Paper 59 `sec:density_ring`; seam row J11 (F); atlas door 21.
Guard: `tests/test_paper59_density_ring.py`, written as its own step and fire-tested 9/9. While
writing it, two facts were re-derived by routes independent of the mod-p searches: L_G G = -3 t^2
on the momentum form with the t-derivatives taken analytically (exact at t = 0.5, 1.3, 2.5), and
in exact rational arithmetic L_G = M o (theta - r), r = 3t^2/(4 - t^2), with the map
f -> t^3 (1 + theta) f / (4 - t^2) sending Legendre solutions to solutions of M. The orders (1 and 2)
leave no room for an order-3 K3 factor in the triangle. The data file `density_loop_B2.json` was
stale (written with the pre-fix operator) and was regenerated: L_G G = -3 t^2 at all six t, and
L_C C below 1e-18. The chain's F_3 sub-piece was then still mod-p evidence; see the next paragraph.

**The chain control, settled (same day, steps D9-D12; scripts `debug/density_loop_D{9_exact_k3,10_k3_exact_check,11_chain_minimal,12_split}.py`, data `debug/data/density_loop_D{9,10,11_dps20,11_dps32,12}.*`).**
(1) The "mismatch" between my notes' M5 and the quotient computed with a monic first-order
factor is exactly a right multiplication: the monic route gives M5 o (1 - t^2), because
D7's R1 = (theta-3) - t^2 (theta+1) is (1 - t^2) (theta - r1). Checked exactly. (2) Exact
rational linear algebra on the F_3 series (with a control: F_3's own operator gives nullity 1;
a first try missed a cube in the series ratio and failed that control, so the control earned its
place) gives the gauge A = [t^3 (t^4 - 18t^2 - 4), t^3 (2t^4 - 43t^2 - 40), t (t^6 - 27t^4 - 96t^2 + 32)] / (t^2 - 4)^2
on (f, theta f, theta^2 f), and M5 o A = 0 modulo F_3's operator in exact arithmetic, with no
series. So M5 = N2 o R3 exactly, R3 of order 3 equivalent to F_3's operator (apparent
singularity at 4t^2 + 11 = 0). (3) N2 is regular at 0 (exponents 5, 7) and at t = +-2
(exponents -1, 0), and irregular at infinity: it is not an elliptic factor. And it is not
needed: W = R3 o R1 o (theta - 2) annihilates C, residual at most 4e-19 at 20 digits and 2e-31 at
32 digits at t = 0.5, 1.3, 2.5, 3.5 (falling with the working precision), while the same W with
one coefficient of R3 changed (198 -> 199) gives 1.9 to 604. The Mellin route's operator carried a
spurious top factor. (4) Does the K3 piece split off the elementary pieces? Outcomes (A) and (B)
both need it: if C = B F_3 + E with B rational and E elementary, then, because F_3's operator is
irreducible, W o B vanishes modulo it, which is a split. Exact search over B with denominators
t^a (t^2-1)^b (t^2-4)^c (4t^2+11)^e (a <= 8, b, c <= 3, e <= 2, numerators of degree <= 16)
finds none; the positive control, the reversed composition (elementary) o R3, finds its split.
So the chain is K3 reached through integrals with elementary kernels, not a rational
combination of F_3 plus an elementary term. **Verdict on the control:** in substance it passes
(the chain's only non-elementary piece is F_3's K3 family, now exact); by the letter "(A) or
(B)" it fails within the searched bounds. The pipeline's check that it puts a known K3 integral
on the K3 rung holds; the registered letter was too narrow for a chain whose leaves are iterated
integrals (chi = int s K1 ds). Paper 59 `sec:density_ring` carries this. Guard: three slow tests added to `tests/test_paper59_density_ring.py`, written after the capture as their own step and fire-tested 8/8 (a changed map coefficient, F_3's operator at argument t^2, a wrong R1, W without theta - 2, a wrong chi, a wrong t^3 shift, a split operator in place of W, a search blinded to degree 2).

**Follow-up (same day, three agents; memos `debug/sprint_{nosplit,minimality_audit,banana_dictionary}_memo.md`).**
(1) The no-split is now proved: a bound on the poles and degree of any splitting map, from the
local exponents, searched exhaustively in exact arithmetic (27 unknowns, rank 27 against 28),
with planted splits found by the same code. Finer: R3 o R1 splits, so the chain satisfies an
order-4 factor P of W (measured, P C / h at most 1e-33 at 30 digits, 9e-45 at 40; control 2e-2),
and the artanh solution is absent: C = c t^2 + t^2 int s^-3 (B' f) ds. Registered (A) and (B)
are refuted exactly, (B) for any elementary term. (2) The minimality audit finds no defect in any
Paper 59 rung claim. The triangle's source term is L_G applied to a constant, so L_G (G + 1) = 0,
and G needs its elliptic piece ((theta - r)(G + 1) is nonzero). (3) For n = 2 to 5 the tree ladder is the banana
ladder under a Hadamard product with the Legendre family (prior art: Almkvist-Zudilin), so the
triangle's elliptic piece is the n = 2 rung itself. Paper 59 `sec:density_ring` and
`sec:four_density` carry all three.

**Prior-art check, Avery (2026-10-17, PI's copy of the LaTeX source of J. Avery and J. S. Avery,
"Hyperspherical Harmonics and their Physical Applications").** Searched for three-electron,
Hylleraas, correlated, explicitly, products r12 r13 r23, star-triangle, uniqueness, Remiddi,
Fromm, Pachucki: no matches. Its interelectron repulsion integrals join two densities only
(ch. 9, Fourier and Sturmian; app. B, Slater orbitals); its "triangles" are mass-spring examples
(ch. 11-12); its "method of trees" (app. C, Aquilanti) builds hyperspherical harmonics and is
unrelated to trees of linked densities. The triangle is not in this book.
