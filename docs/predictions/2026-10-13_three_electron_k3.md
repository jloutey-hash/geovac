# Pre-registered prediction: the three-electron chain integral is a K3, not a genus-2 curve

**Registered:** 2026-10-13, after the four-centre check (`2026-10-13_four_centre_genus.md`,
P1 confirmed) and BEFORE any three-electron computation (PI-directed: "go straight to P3").
**Supersedes:** P3 of the four-centre registration ("genus 2 first appears in three-electron
integrals"), which is re-stated here in a form that can fail.
**Rule being extended:** Paper 59 `sec:momentum`, `sec:genus` (fingerprint = the D = 0 fibre,
whose period identifies the curve).

## The integral

The simplest genuinely three-electron Coulomb integral, the *chain*
I = <rho_a(1) rho_b(2) rho_c(3) / (r12 r13)> = int rho_a(r) V_b(r) V_c(r) d^3r,
with rho_a, rho_b, rho_c two-centre 1s Slater densities. In momentum space
I ~ int d^3p d^3q  rho_b(p) rho_c(q) rho_a(p+q) / (p^2 q^2)  (times phases):
two momentum integrations, three densities.

## Why "genus" was the wrong question

Each two-centre density carries one dispersion root (Paper 59; confirmed for four centres):
Delta_b = sqrt(c_b p^2 + 1), Delta_c = sqrt(c_c q^2 + 1), Delta_a = sqrt(c_a |p+q|^2 + 1).
With two integration variables the natural object is a surface, not a curve.

## Predictions

**P3a (reduction).** At the fingerprint fibre (D = 0, all phases trivial, W = 0), the cos-angle
integral between p and q is elementary and keeps Delta_a algebraic (no log, no new root): the
fibre becomes a 2-D integral over scalar p, q >= 0 whose ONLY square roots are
sqrt(c_b p^2+1), sqrt(c_c q^2+1), sqrt(c_a (p +/- q)^2 + 1).

**P3b (the rung is a K3 surface).** Homogenized, the three quadrics factor into six lines in P^2.
The double cover of P^2 branched along six lines is a K3 surface (a Calabi-Yau 2-fold), with at
worst rational double points for every c_a, c_b, c_c > 0 (no point where four lines meet).
Consequences: the coordinate slices (p fixed) are genus-1 curves whose modulus moves with p
(non-isotrivial), and generic line slices are genus-2 curves. So the old hunch "genus 2" is right
only on generic slices; the invariant new rung is a K3.

**P3c (degeneration).** If every density is one-centre (rational transform, no root), the chain
integral is elementary plus logarithms: no K3, no elliptic.

**Pattern (conjecture, not tested here).** Two-centre ERI (one momentum, two densities): double
cover of P^1 on 4 points = elliptic curve (CY 1-fold), matching Paper 59. Chain with two momenta:
P^2 on 6 lines = K3 (CY 2-fold). Same ladder as the banana Feynman integrals (sunrise elliptic,
three-loop banana K3). Whether longer trees give CY n-folds is open; the angular reduction for
more than one angle is not elementary.

## Falsifiers (any one breaks P3a/P3b)

1. The angle integral produces a logarithm or a fourth independent root at D = 0.
2. The reduced 2-D fibre disagrees with a direct (p, q, cos theta) evaluation.
3. For some generic c's four of the six lines are concurrent (or two coincide), so the cover is
   not K3.
4. The p-fixed slice is isotrivial (modulus independent of p), i.e. the surface is a product or
   rational, not a genuine K3.
5. The full chain integral in momentum form disagrees with an independent position-space route.

## Result (2026-10-13, same day; `debug/three_electron_k3_check.py`, log + json in `debug/data/`)

All five falsifiers checked; none triggered.

| Check | Outcome |
|:--|:--|
| F1 angle integral | Elementary, and stronger than predicted: for every D, not just D = 0. Yukawa kernel gives -e^{-D Delta}/(D c p q); the D = 0 density 1/Delta^5 gives -1/(3 c p q Delta^3). No log. sympy residual 0. |
| F2 reduced vs direct fibre | Agree to 1e-15 to 4e-14 at three (c_a, c_b, c_c) triples. |
| F3 six lines | No four concurrent for any positive c. Triple points only on the loci sqrt(c_a)(sqrt(c_b) +/- sqrt(c_c)) = sqrt(c_b c_c) etc. (rational double points; still K3). |
| F4 slices | The p-fixed slice is genus 1 and its j-invariant moves with p (1.2e7 at p = 0 down to 1.75e3 at p = 5): non-isotrivial. Five random line slices have six distinct branch points: genus 2. |
| F5 one-centre degeneration | Momentum route equals the exact position-space value to 1e-15: 23/54, 133/144, 289/225 for zeta = (1,1,1), (1,2,3), (2,1,3). Stronger than predicted: rational, not even a log. |

**Addendum P3d (registered after F1-F5, before any period test).** The line configuration is
special: conjugate pairs whose three meeting points are collinear (all on w = 0), and only two
moduli survive (the ratios of c_a, c_b, c_c). A two-parameter K3 family generically has Picard
rank >= 18, and rank-18 K3s typically have periods that are *products of two elliptic periods*
(Shioda-Inose / Kummer of E x E'). Prediction: the D = 0 fibre is a product-of-elliptic-periods
object (plus lower-rung boundary terms), not a "new" K3 transcendental. Test: one-parameter
Picard-Fuchs order, or PSLQ of a high-precision fibre value against K K', K E', E E' at the
moduli of the two coordinate-slice curves. If P3d holds, the three-electron rung is
"elliptic squared", not a new class.

## P3d result (2026-10-13, after the registration commit 5542244e; `debug/three_electron_p3d_check.py`)

**P3d HOLDS for the holomorphic period, and the closed form is classical.** A 1-D Fourier
transform (Paper 59 eq:K0 at D = 0) sends the K3 period
P = int_{R^2} dp dq / sqrt((cb p^2+1)(cc q^2+1)(ca(p+q)^2+1)) to a three-scale Bessel moment,
P = 8/(pi sqrt(ca cb cc)) V3(1/sqrt ca, 1/sqrt cb, 1/sqrt cc), V3(a,b,c) = int_0^inf K0(at)K0(bt)K0(ct) dt.
Bailey (1936), restated as Bailey-Borwein-Broadhurst-Glasser, J. Phys. A 41, 205203 (2008),
arXiv:0801.0891, eq. (33): (2c/pi) V3 = K(k-)K(k+') + K(k+)K(k-'), with
k+- = (sqrt((c+a)^2-b^2) +- sqrt((c-a)^2-b^2))/(2c), a sum of products of two complete elliptic
integrals with moduli algebraic in the scales.

| Check | Outcome |
|:--|:--|
| Equal-mass anchor V3(1,1,1) | quad = Gamma form = (pi/2)K3K3' = eq. (33): all 30 digits agree (6.94882278107962978943643644547) |
| V3 quad vs eq. (33), four (ca,cb,cc) | agree to 1e-34, including three cases where k+- are complex conjugates |
| P direct 2-D quadrature vs Bessel reduction | agree to 2e-10 to 7e-10 (relative ~1e-12; scipy subdivision warning on the slow tails) |

**Reading.** The geometry climbs (curve to K3) but the period does not open a new class: the
three-electron chain's period is "elliptic squared", as Duhr (MPIM talk, 2024-08-19) reports for
three-loop banana maximal cuts with three equal masses (products of complete elliptic integrals
with algebraic arguments). Our case has three independent scales and still factorizes, by Bailey.
GeoVac's contribution here is the reduction (three-electron Slater chain to Bailey's V3), not
the closed form.

**Not settled (the next frontier, the analogue of Paper 59's MMV question).** The actual D = 0
fibre (F2's integrand) is not the holomorphic period: the two Coulomb kernels leave 1/(pq)
factors after the angle integral, i.e. iterated (antiderivative) structure, so the integral is
a relative period that may carry regulator-type content beyond K K', as the full banana
integral does (elliptic trilogarithms, Bloch-Kerr-Vanhove) beyond its maximal cut.

## Prior art to check

**Corpus first (searched 2026-10-13, PI prompt).** Paper 59 `sec:context` (around l.680) already
places the two-electron ERI next to the two-loop *sunrise* (elliptic) and cites
Bloch-Kerr-Vanhove (three-loop banana, K3) and Bonisch-Fischbach-Klemm-Nega-Safari (all-loop
banana). `debug/sprint_routeC_irregular_resurgent_memo.md` l.31 set the K3 banana objects aside:
"a DIFFERENT geometry (K3 twofold), NOT our genus-1 two-scale sunrise". That verdict stands for
two electrons; the three-electron chain is the first GeoVac object on the K3 side. Correct
dictionary: n momentum integrations here <-> (n+1)-loop banana (CY n-fold), i.e.
ERI <-> sunrise, chain <-> three-loop banana. `debug/lit_scan/elliptic_eri_feynman_memo.md`
items 12-15 found no chemistry use of these periods (as of that scan). Paper 59 also records
Broadhurst's "two-scale moment -> product of K's" mechanism, the P3d mechanism one level down.
Known, to verify against primary sources: the equal-mass three-loop banana K3 has Picard rank 19
and periods that are symmetric squares of elliptic periods, which supports P3d.

**Literature.** Banana / sunrise Calabi-Yau ladder (Bloch-Vanhove, Bloch-Kerr-Vanhove, Bonisch-Fischbach-Klemm-
Nega-Safari); double covers of P^2 on six lines (Matsumoto-Sasaki-Yoshida, E(3,6) hypergeometric);
three-electron Slater integrals (Fromm-Hill 1987 for one-centre r12 r23 r31, dilogarithms;
Harris; Remiddi). Not checked before registering.
