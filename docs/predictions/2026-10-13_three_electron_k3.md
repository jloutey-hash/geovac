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

## P3e registration: the true fibre (registered 2026-10-13 after 2d540a32, before computing it)

**Setup (derivation, to be checked numerically).** With rho = (c k^2+1)^{-5/2} and
g = (c_a u^2+1)^{-3/2}/(3 c_a), the true D = 0 fibre (F2's integrand) folds to
F = -(1/2) PV int_{R^2} rho_b(p) rho_c(q) g(p+q) / (p q), and a 1-D Fourier transform gives
F = 4/(27 pi c_a^2) int_0^inf x K1(a x) J(b x) J(c x) dx,  a,b,c = c^{-1/2},
J(X) = int_0^X z^2 K2(z) dz = 3 L(X) - 3 X K0(X) - X^2 K1(X),  L(X) = int_0^X K0.
L is the iterated piece (the Coulomb kernels' 1/(pq)). Writing Ki(bx) = bx int_1^inf K0(bxu) du
turns F into a two-fold mass integral of three-Bessel moments, i.e. of Bailey K K' products:
the Coulomb kernels act as integrations over the K3 family, as Feynman parameters do for
Paper 59's ERI. The residues of F's 2-form along p = 0 and q = 0 are the Paper 59 ERI curves
y^2 = (c_c q^2+1)(c_a q^2+1) and (c_b p^2+1)(c_a p^2+1).

**P3e-1 (equal scales collapse).** At c_a = c_b = c_c, integration by parts swaps every Ki
against a Bessel of the same argument, so F lies in the Q-span of {c30, c32, Tate numbers}
(c3k = int x^k K0^3): no iterated content. Falsifier: PSLQ at 60 digits finds no such relation.

**P3e-2 (generic scales do not collapse).** At unequal rational masses F is NOT in the
Q-span of {the K3 periods (int K0K0K0 and its three first mass-derivative moments), pi x the
two residue-curve periods (int K0K0, int x K1 K0 for each), Tate numbers built from the
mass logs}. Prediction: PSLQ at 100 digits finds no relation with coefficients below 10^6.
Falsifier: a relation that survives a jump to 140 digits. Honest cap: a null PSLQ is evidence,
not proof; a positive one would be decisive the other way.

## P3e result (`debug/three_electron_p3e_check.py`; values cached in `debug/data/`)

| Check | Outcome |
|:--|:--|
| 1-D formula vs F2's 2-D fibre | agree at all three (c_a,c_b,c_c) triples to 3e-15 to 4e-14 (the 2-D accuracy) |
| P3e-1, equal scales | **holds.** 216 G(1,1,1) = 80 c30 - 6600 c32 + 81 pi^3, identical at 40 and 70 digits. Only K3 periods plus pi^3. |
| K3 block at masses (1,2,4) | two exact relations among the six moments int x^k K K K (one at 40 digits; V = 13 D1 + 4 D2 - 26 S12 - 44 S13, identical at 110 and 150), leaving four independent periods {V, D1, D2, S12}: rank 4 = 22 - 18, consistent with Picard rank 18. |
| P3e-2, generic, test as registered | **INCONCLUSIVE: the registered test cannot see a known relation.** At 150 digits, all three PSLQ calls (basis alone, control, G) return 21-term relations with coefficients 1e5-1e6. That is the chance-relation scale for 21 numbers at 150 digits (about 10^(150/21) = 10^7), so all three are noise. The positive control shows the registered cap of 1e6 was set too low. Every piece of Gctrl reduces exactly to {V, D1, D2, S12} at masses (1,2,4) (x^3 K1K0K0, x^4 K1K1K0, x^4 K1K0K1, x^5 K1K1K1; PSLQ at 60 digits). The assembled reduction is Gctrl = (11050792 V - 27272112 D1 - 52735824 D2 + 87003024 S12)/363825, with a residual of 1.6e-149 at 150 digits (`debug/data/three_electron_p3e_control_reduction.json`). Genuine relations at these masses need coefficients around 1e8, so a null at 1e6 is not evidence either way. Two basis facts were also checked: P_uv = pi K((v-u)/(v+u))/(u+v) exactly, and Q_uv = -dP/du is the quasi-period (not elementary), so the basis is the right kind. |
| P3e-2, generic, powered rerun (amended after the control, before G was examined) | **HOLDS (evidence, not proof).** Values were recomputed by an independent route: one shared double-exponential grid for every integral (`debug/three_electron_p3e_fast.py`). It agrees with the adaptive 150-digit values to 3e-149, and the step-h/step-2h difference is 4e-264. That took 48 minutes, where the adaptive route needed more than 9 hours just for the basis. In the 25-element search at 320 digits (cap 1e11), the control's exact relation above is recovered. G's returned relations at 260 and 320 digits are different dense 25-term vectors with coefficients 1e9-1e10, i.e. chance relations. |
| P3e-2, smaller bases (more PSLQ headroom; `debug/three_electron_p3e_subspans.py`) | In every sub-basis the control is recovered with its exact coefficients at both 260 and 320 digits. G has no relation that reappears at both precisions. K3 + {pi^3, pi^2, pi, 1}: nothing below about 1e24. Adding the two residue curves' periods and quasi-periods, with and without a factor pi: none below 1e12 at 260 digits. Adding log2: none below 1e10 (260) or 1e13 (320). Adding log2 and log3: none below 1e9 (260). Every relation returned was a dense vector at the cap that did not reappear at the other precision. **Reading:** at equal scales the fibre collapses onto K3 periods plus pi^3; at scales (1,2,4) it does not, so the generic fibre carries content beyond the K3 periods, the residue curves' elliptic periods and simple Tate numbers. This is the expected iterated, mixed-motive extension, the three-electron analogue of what Paper 59 found one level down. **Caps:** a null PSLQ is not a proof. One mass triple was tested. The basis omits candidates such as dilogarithms of mass ratios, the (b,c) residue pair, and products of elliptic periods (which the iterated structure could produce). |

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
