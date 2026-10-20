# Pre-registered prediction: the correlated triangle, and a ring of four

**Registered:** 2026-10-17, before any computation in this file (PI-directed: "let's use some sub
agents to run as many of these next steps as is reasonable", after the ring result of
`2026-10-17_density_loop.md`).
**Extends:** `2026-10-17_density_loop.md` (PT1-PT4, Result, and "The chain control, settled");
Paper 59 `sec:density_ring`, `sec:four_density`, `sec:three_electron`; seam rows J9-J11.
**Commit status:** to be committed on its own, before any computation, if the PI authorizes it.
If computation starts uncommitted, the Result section must say so.

## Why these two

The Coulomb triangle sits on the elliptic rung, one below the chain of three, because a closed
loop of Coulomb kernels integrates out exactly in three dimensions (the star-triangle relation).
Two questions follow. (1) The LiH explicitly correlated engine meets the triangle with geminal
kernels, <f_ik f_jk / r_ij>, whose loop is massive and not "unique", so the exact step is gone.
Does that triangle still sit low? (2) Is "a ring sits one rung below" a rule, or a fact about
three? The smallest test is a ring of four Coulomb kernels.

## Conventions (as in PT1-PT3)

Fingerprint fibre: densities spherical at zero separation, rho~_v(k) = Delta_v(k)^(-3),
Delta_v = sqrt(c_v k^2 + 1), a_v = c_v^(-1/2); in position space rho_v is proportional to
K0(r / a_v). theta = t d/dt. Rungs: elementary or algebraic < elliptic (Legendre, F_2) <
K3 (F_3's operator 4 theta^3 - t^2 (theta+1)^3) < the four-density rung (AESZ #58).
"Sits on rung X" means: the minimal annihilator of the function, checked on the function
itself, has composition factors no higher than X, and at least one factor equivalent (rational
gauge, rational pullback allowed and stated) to X's operator. Lesson of 2026-10-17: an operator
from a Mellin recurrence is an upper bound; its minimality and any source term must be checked
on the function (the triangle's -3 t^2; the chain's spurious N2).

## QG: the correlated triangle

Object: Y(t) = <rho_i rho_j rho_k f(r_ik) f(r_jk) / r_ij> with the exponential geminal
f(r) = exp(-gamma r), at gamma = 1 and gamma = 3/2, scales (a_i, a_j, a_k) = (1, 1, t), t on the
apex k (the density carrying both geminals). Secondary: t on i. Route allowed: f(r) =
-d/dgamma [exp(-gamma r)/r], so Y is a mixed second gamma-derivative of the Yukawa-Yukawa-Coulomb
triangle with two masses set equal; the 3-D one-loop triangle with two massive and one massless
propagator is then the loop integral.

**QG1 (the massive loop is elementary, symbolic + measured).** The 3-D one-loop triangle with
masses (gamma, gamma, 0) is an elementary function of the three external momenta (inverse
tangents and square roots), in closed form; check it against direct 3-D integration at five
random external triples, relative < 1e-8. Prior art to verify, not assumed: Nickel, J. Math.
Phys. 19, 542 (1978).

**QG2 (main: the geminal does not raise the rung).** Y(t) sits on the elliptic rung or lower,
at both gamma: no factor equivalent to F_3's K3 operator and nothing higher. Reason: the loop
is still elementary, so after it the integrand lives on the same two momenta as the Coulomb
triangle, with an inverse tangent where there was a rational kernel; an inverse tangent adds an
integration, which can add weight, not geometry. Confidence: moderate. Also predicted: new
singular points that move with gamma (for example where c_k gamma^2 hits a fixed value), and
non-split attachment of the elliptic piece (an inverse tangent is an integral).

Outcomes: (G-low) elliptic or lower, as predicted. (G-chain) a K3 factor: the geminal triangle
sits with the chain, the mass undoes the drop. (G-climb) a factor above K3.

## QR: a ring of four

Object: S(t) = <rho_a rho_b rho_c rho_d / (r_ab r_bc r_cd r_da)>, four densities in a square of
Coulomb kernels, scales (1, 1, 1, t). Position-space route allowed: with K_v(x, y) =
int rho_v(z) / (|x - z| |z - y|) d^3z, S = int rho_a(x) rho_c(y) K_b(x, y) K_d(x, y), a
three-dimensional integral over |x|, |y| and the angle. The 3-D one-loop box with massless
propagators is not unique (indices sum to 4, not 3), so the triangle's exact step does not apply.

**QR1 (main, low confidence: a ring of n sits at the rung of a tree of n - 1).** S(t) has a
composition factor equivalent to F_3's K3 operator and nothing higher: the square sits on the
K3 rung, one below a tree of four. Reason: the only data point, the triangle, sits one below the
tree of its size; the guess is that closing a ring costs exactly one rung. This is a guess
from one point, registered so that it can fail.

Outcomes: (R-n-1) K3, as predicted. (R-ell) elliptic or lower: rings collapse to the bottom
whatever their size. (R-tree) the four-density rung: a ring of four sits with the tree of four.
(R-other) a factor not on the ladder.

**QR2 (one-centre limit, evidence only).** At one centre with 1s Slater densities the square is
not rational, by analogy with the triangle's dilogarithms. Report the smallest fitting basis.

## Tests and decision rules

1. QG1 as stated. For QR, the K_v kernel against direct 3-D quadrature at five (x, y) pairs,
   relative < 1e-8.
2. For each of Y and S: obtain an annihilating operator (Mellin route or series route), then
   check L f on the function at >= 3 values of t at >= 25 digits, including any source term;
   peel it exactly into factors (rational arithmetic); for each factor of order >= 2 test
   equivalence to Legendre and to F_3's operator (rational gauge, exact verification of any
   map found) and apply the quadric counter of `debug/linking_ceiling_C_quadrics.py`. Test
   whether the function needs every factor (the D11 check).
3. Bounded searches state their bounds. A "none" from a search needs a positive control that
   the same search finds a planted answer.

## Controls (the test is void if one fails)

- gamma -> 0 control for QG: with Yukawa kernels exp(-gamma r)/r on all three edges, the
  gamma -> 0 limit reproduces the Coulomb triangle G(t) of PT2 (numerically at two t).
- Square control for QR: with the fourth kernel removed (an open chain of four), the same
  pipeline must reproduce the tree-of-four rung of Paper 59 `sec:four_density`.
- The chain of three at the same fibre must still come out with its K3 piece (D9-D11).

## Falsifiers

1. QG2: a K3 factor or higher in Y at either gamma, with the controls passing.
2. QR1: the square sits on any rung other than K3, with the controls passing.

## Honest caps, stated before the result

- Fibre statements only; the integrated value is PI item 1, untouched.
- One geminal shape (exponential); the cusp-correct r exp(-gamma r) and the linear-exponential
  geminal of the LiH engine are not covered. Two values of gamma.
- One larger ring. Rings with trees attached are not covered.
- The ring-of-four computation may be too heavy for the precision a rung decision needs; a STOP
  for feasibility is an allowed outcome and must be reported as such, not as a result.

## Result (2026-10-17, same session, after commit 8ba60ae5)

### QR: STOP for feasibility. No rung decision.

Memo `debug/sprint_ring_ext_QR_memo.md`; scripts `debug/ring_ext_QR_{box,kernel,pw,mom,data,search,omega}.py`,
data `debug/data/ring_ext_QR_*`. QR1 is neither confirmed nor falsified; none of R-n-1, R-ell,
R-tree, R-other can be read, and the open-chain-of-four control was not carried to a rung.

What passed. Test 1: the kernel K_v has a one-dimensional form (a Feynman-parameter integral of
-z + (pi z/2)(H_1 - Y_1), Struve and Bessel), matching direct 3-D quadrature at five (x, y) and
two scales to 2e-16 .. 9e-14. The 3-D massless box reduces exactly to its four triangles
(at N = D + 1 the dimension-shift term vanishes): Box = pi^3 [sum over vertices of 1/(product of
the three edges there)] / (Ptolemy sum of opposite-edge products), matching Feynman parameters
at four configurations to 7e-16 (prior art not searched; likely known). The ring itself by two
independent routes, position-space partial waves and momentum space, agreeing to 2e-8:
S(0.7) = 0.022023792114497, S(1.6) = 0.04913982868458; S / t^2 tends to the chain of three at
(1,1,1) as t -> 0; route (a) reproduces the triangle and the chain of three, whose recorded
operators annihilate its data to 1e-14.

Why it blocked. After the loop the square keeps a weight that couples all its variables (the
triangle's was a constant), so no operator came out exactly, and creative-telescoping tools are
not installed. A float64 operator search failed its own positive controls (it did not recover the
chain's known K3 operator), so by test rule 3 its "none" for S is void. Unblock: 60+ digit data
(the evaluator is cancellation-free and ports to arb or double-double), or creative telescoping
on route (b)'s explicit four-dimensional algebraic integrand, or a closed form for the weight.

Correction to the conventions paragraph above: the fibre density with transform
(c k^2 + 1)^(-3/2), a = c^(-1/2), is proportional to K0(a r), not K0(r / a). The computation
used the transform, so nothing depends on the slip.

### QG: G-chain at gamma = 1. QG2, the main prediction, is falsified.

Computation began after commit 8ba60ae5. Memo `debug/sprint_ring_ext_QG_memo.md`; scripts
`debug/ring_ext_QG_*.py`, data `debug/data/ring_ext_QG_*`. All registered controls pass.

QG1 holds [SYMBOLIC + MEASURED]. The 3-D loop with masses (gamma, gamma, 0) is a solid angle,
pi^2 arctan(sqrt(det M)/D)/sqrt(det M), and its (2,2,1) mass-derivative is the geminal loop. Both
match direct 3-D integration at five triples for each gamma, worst relative difference 2.5e-14.
The Nickel 1978 reference was not found and is NOT VERIFIED; the closed form rests on the check.

Controls. The gamma -> 0 member (Yukawa on all three edges) reproduces G(t) to 5e-15 at t = 0.7
and 3e-14 at t = 1.6, and the gap closes steadily at gamma = 0.1, 0.03, 0.01. Two routes to the
kite (momentum triangle; a double-bubble algebraic double integral) agree to 4e-16. The rung scan
finds F_3 in F_3, Legendre in Legendre, neither in F_4, and the K3 piece in the chain's W.

gamma = 1 [SYMBOLIC, except where marked].
- The kite phi(z) has an order-2 operator [MEASURED on 600 coefficients, residual 3e-155], mapped
  exactly to Legendre at t^2 = -z/4. The kite is elliptic.
- The minimal annihilator of Y is an operator G4 of order 4: the exact common right divisor of the
  Mellin-route operator and of the operator of Y's log t part. G4 Y is 1e-27 to 2e-26 of the
  largest term at t = 0.5, 0.8, 1.3, 3.0 [MEASURED, 30-digit grid]; the small remaining absolute
  residual is a constant offset in Y of about 4e-30, matched by Q_0(0) at three t. No source term.
  No order-3 annihilator up to t^2-degree 99.
- G4 = R3 R1. R1 has the algebraic solution t^3(33t^6-612t^4+5232t^2-49088)/(t^2-4)^(9/2).
  R3 is the adjoint of an operator R with Sol(R) = A Sol(L3), A an exact order-2 map with
  denominator (9t^4+312t^2-12560)^2; L3 is anti-self-adjoint, so R3 is equivalent to F_3's
  operator. The K3 piece does not split off (bounded search; the same search finds the planted
  split in R1 R3). Composition factors: algebraic and K3, nothing higher; Y needs both.
- Singular points t = 0, +-2, infinity; at gamma = 1 the thresholds 2, 2 gamma and 1 + gamma
  coincide at t = 2.

Checked by me with the helpers of `tests/test_paper59_density_ring.py`, on the agent's printed
operators: the memo's G4 equals the stored one; G4*(A f) = 0 mod L3 exactly; A is injective;
G4* = Q1 R exactly with Q1 of order one; L3* = -L3; R1's solution is exact; G4 is right-divisible
by R1; changing one coefficient of A by 0.1% breaks the identity. Independently, phi(0) = 3 pi/16
at gamma = 1 by a position-space route the agent did not use (closed-form convolutions).

gamma = 3/2 [MEASURED mod p + rigidity]. Y_log has an exact integer operator L9 of order 9; L9 Y
is below 7e-28 at t = 0.8, 1.3, 1.7 while a 1e-6 perturbation gives 1e-10 to 5e-14. Singular
points t = 1/2, 2, 5/2, 3, that is |1-gamma|, 2, 1+gamma, 2 gamma: the side prediction that they
move with gamma holds. An order-3 quotient N3, found mod p, has F_3's local data twisted by t
(maximally unipotent at 0, a pseudo-reflection at t = 2, no other true singular points), so by
the rigidity of hypergeometric systems (Beukers-Heckman, used from memory, NOT VERIFIED) it is
equivalent to F_3. The explicit map was not found within the searched bounds; the same blind
pipeline finds it at gamma = 1. The middle factors are not all identified, so a factor above K3
is not excluded at gamma = 3/2.

Verdict against the registration. QG1 holds. QG2 is falsified: outcome G-chain at gamma = 1; at
gamma = 3/2, K3 by local data, G-climb not excluded. By the decision rule's letter, gamma = 3/2 is
"K3-equivalent, map not exhibited". Non-split attachment holds at gamma = 1.

Why the registered reasoning failed. It said an inverse tangent adds weight, not geometry. The
loop is elementary, but the kite around it is elliptic, and attaching the apex density, a
Stieltjes transform, turns the elliptic family into F_3's K3 family (a Clausen-type step). The
Coulomb triangle stays low only because its star-triangle step removes the kite. So the drop of
the ring is a fact about uniqueness, not about rings.

Not a result: the secondary placement (t on i) is STOP for feasibility.

Caps. Fibre only; one geminal shape (exponential, not the cusp-correct or linear-exponential
geminals of the LiH engine); two values of gamma; apex placement. The kite operators are guessed
from series and checked on 600 coefficients, not derived. gamma = 1 is the special value where
the thresholds coincide.


## Novelty classification (appended 2026-10-19)

Classes per `docs/predictions/README.md`, set after the fact by a literature audit (`debug/prediction_novelty_audit_memo.md`; this record predates the novelty field). FINDING = held, not found in the literature; CHECK = held, already known; MISS = failed.

- QG1: CHECK (Davydychev-Delbourgo 1998 Sec. 5.1, after Nickel 1978; now credited in Paper 59)
- QG2: MISS
- QR1: UNSCORED
- QR2: UNSCORED
