# QR: the Coulomb ring of four at the fingerprint fibre (2026-10-17)

Registration: `docs/predictions/2026-10-17_ring_extensions.md` (committed 8ba60ae5 before any
computation). Object: S(t) = <rho_a rho_b rho_c rho_d / (r_ab r_bc r_cd r_da)>, scales (1,1,1,t).
Prediction QR1 (low confidence): the K3 rung, one below the tree of four.

## Verdict: STOP-feasibility

No annihilator of S(t) was obtained, so no rung is decided: not R-n-1, R-ell, R-tree or
R-other. Built: the object (two independent evaluators that agree), the K kernel (step 1,
passed), and a closed form for the loop. Blocked: the annihilator. The loop does not collapse
as the triangle's did, and the float64 search used here cannot find an operator even as large
as the chain of three's (calibrated below). What would unblock it is listed at the end.

Convention note. The registration writes "rho_v proportional to K0(r / a_v)" with
a_v = c_v^(-1/2). The Fourier transform of (c k^2 + 1)^(-3/2) is a^3 K0(a r)/(2 pi^2), so the
density is K0(a_v r), with mass a_v, and "scale t" means rho_d = t^3 K0(t r)/(2 pi^2). This is
the convention of the density-loop code (G_mom uses c = 1/t^2), and it is the one used here.
The registration's K0(r/a_v) is a slip in wording, not a different object.

## 1. The K_v kernel (registration test 1): PASS

K(x, y) = int rho(z) / (|x - z| |z - y|) d^3z. A 1-D form, from a Fourier transform, a Feynman
parameter u, and the q-integral in closed form:

    K(x,y) = (2/pi) int_0^{pi/2} dth  a Im F1(a (b - i X)) / X,   u = sin^2 th,
    X = |(1-u) x + u y|,  b = |y - x| sqrt(u(1-u)),
    F1(z) = int_0^inf e^{-zq} (1+q^2)^{-3/2} dq = -z + (pi z/2) (H_1(z) - Y_1(z))

(H_1 Struve, Y_1 Bessel; F1 checked against quadrature at complex z; the limit rho -> delta
gives 1/(|x||y|) exactly). Against direct 3-D quadrature in prolate spheroidal coordinates with
foci x, y (scipy, float64), at five random (x, y) and two scales (a = 1, 0.6): relative
agreement 2e-16 to 9e-14 (`ring_ext_QR_kernel.py`, `data/ring_ext_QR_kernel.json`).
Requirement < 1e-8: met.
[Transcendentals: the pi's are the Fourier/angular measure, the Paper 18 calibration class.
Struve and Bessel functions come from the density's transform. No new constants.]

## 2. The loop in closed form: the 3-D massless box

In density-momentum space, S = (2 pi)^(-12) (4 pi)^4 int d^3L_a d^3L_b d^3L_c prod rho~_v Box.
Here Box is the one-loop massless box in three dimensions. Its "dual points" P_0..P_3 are
spaced so that consecutive distances are the density momenta. Since N = D + 1, the
dimension-shift coefficient (N - D - 1) is zero. The box therefore reduces exactly to its four
triangles: Box = sum_i b_i Tri_i, b = S^{-1} 1, S_ij = |P_i - P_j|^2. This is Melrose's
mechanism, re-derived here by integration by parts in Feynman parameters. Each triangle is
unique (star-triangle), Tri = pi^3/(d d d). sympy collapses the sum to

    Box = pi^3 [ sum_i 1 / prod_{j != i} d_ij ] / (d_01 d_23 + d_02 d_13 + d_03 d_12).

The numerator is a sum over vertices of 1/(product of the three edges at that vertex). The
denominator is the Ptolemy sum of opposite-edge products, which never vanishes. The Ptolemy
determinant in the b_i (zero on planar cyclic quadrilaterals) cancels out. Checked against a
direct Feynman-parameter evaluation (primary sectors, x = y^2) at four random configurations:
relative 0 to 7e-16 (`ring_ext_QR_box.py`, `data/ring_ext_QR_box.json`).
Tier: SYMBOLIC (sympy) + MEASURED. Prior art: probably known, since 3-D one-loop integrals are
standard in ABJM-type amplitude work. I did not search. NOT VERIFIED; do not present it as new
without a search.

Consequence. Use the rim lengths l_a..l_d and the diagonals e (P0P2), f (P1P3). As in the
triangle, the rim lengths cancel the box's edge factors:

    S = (1/(2 pi^3)) int d^4 l  prod_v rho~_v(l_v)  Omega(l),
    Omega = int de int_0^{2pi} dphi [ (l_b l_c + l_a l_d)/e + (l_a l_b + l_c l_d)/f ] / (l_a l_c + l_b l_d + e f),

with phi the dihedral angle of the folded quadrilateral. Omega is homogeneous of degree 0;
numerically Omega(2l) = Omega(l) to 15 digits. Its small-l_d limit is 4 pi l_d l_b/(l_a l_c),
also checked numerically. For the triangle the analogous shell function is a constant on the
triangle-inequality region, which is why the triangle collapsed to a 2-D formula. Here Omega is
not constant. It has log singularities on {l_a = l_d, l_b = l_c} and {l_a = l_b, l_c = l_d},
where a diagonal can vanish. Its inner phi-integrals are elliptic integrals of the third kind on
a curve that depends on (e, l). Whether the e-integration leaves Omega elementary was not
decided. So S is an integral over at least four dimensions with a weight that does not look
elementary. The triangle was an integral over two dimensions with a constant weight. This is
the structural reason the triangle's exact step does not carry over.

## 3. Two independent evaluators of S(t): they agree

(a) Position-space partial waves (`ring_ext_QR_pw.py`). Expand 1/|x - y| = sum_l g_l P_l,
g_l = r<^l / r>^(l+1). Around a ring of n the angular integrals give (4 pi)^n / (2l+1)^(n-1),
so Ring_n = (4 pi)^n sum_l (2l+1)^(1-n) R_l. Because the kernels are semi-separable, R_l is a
sum over the n! radial orderings of nested 1-D integrals. Each nesting is an exponential
smoothing in log r, done panel by panel with exact exponential weights (Gauss-Legendre,
float64). The l-tail is fitted in (2l+1)^(-k), k >= 2n-2, and summed with the Hurwitz zeta
function. Resummed, the l-sum is a Legendre chi function of the ratio r_min/r_max: chi_2 for the
triangle, chi_3 (a weight-3 polylogarithm) for the square.
- Validation on the triangle, against the record's 25-digit values from independent Struve and
  momentum routes: G(0.7) relative 1.2e-15, G(1.6) relative 8e-16. On the chain of three
  (l = 0, open), C(0.7) and C(1.6) agree to float64.
- The record's operators applied to this evaluator's theta-derivatives at 64 points,
  t in [0.25, 4]: (theta - 2) L_G G = 0 to a maximum normalised residual of 2.3e-15, and the
  chain's W C = 0 to 9.2e-15. This checks the theta-derivative weights and the evaluator
  together. It does not reach the bimodal orderings (two minima and two maxima around the
  ring), which first appear at n = 4; route (b) covers those.
- S(0.7) = 0.022023792114497034, S(1.6) = 0.0491398286845813. Both are stable to 1e-16 under
  l_max 80 -> 120, three tail fits, and a refined grid. The l = 0 term is 99.5% of S.

(b) Momentum space, kernel momenta on the edges (`ring_ext_QR_mom.py`). The partial waves of the
fingerprint transform are geometric:
(1 + c(k^2 + k'^2) - 2ckk'x)^(-3/2) = (h/(ckk'))^(3/2) (1 - h^2)^(-1) sum (2l+1) h^l P_l(x),
with h = z - sqrt(z^2 - 1) and z = (c(k^2+k'^2)+1)/(2ckk'). This is the classical generating
function (1 - h^2)(1 - 2xh + h^2)^(-3/2). For any ring the angular weight is (4 pi)^n (2l+1),
so the l-sum resums to (1 + H)/(1 - H)^2, H = prod h_v:

    Ring_n = (4pi)^(2n) (2pi)^(-3n) int_{R+^n} d^n k  prod_v w_v (1+H)/(1-H)^2,
    w_v = [1/Delta_v(k - k') - 1/Delta_v(k + k')] / (2 c_v k k'),
    h_v = (Delta_v(k + k') - Delta_v(k - k'))^2 / (4 c_v k k').

A 4-D tensor Gauss rule converges as N^-3.5. At N = 96 the raw values are 7e-7 and 1.0e-6
from route (a). A three-point power extrapolation agrees with (a) to 1.7e-8 (t = 0.7) and
2.1e-8 (t = 1.6). The triangle through the same code converges toward the record's values.

Sanity check: S(t)/t^2 -> 0.1557888..., the chain of three at scales (1,1,1), as t -> 0
(0.1151, 0.1305, 0.1407 at t = 0.1, 0.05, 0.025, with a t log t correction). This is the
expected leading term: at large separation the two kernels at d give 1/R^2 and leave the open
chain a-b-c.

## 4. The annihilator: why it is blocked

The triangle's operator came from the Mellin route on a 1-D Bessel moment. That moment exists
because the loop collapses into one factor per density. For the square, every route set up here
leaves the loop coupled: the l-sum (or chi_3) couples all four radii; the box's Ptolemy
denominator couples all six edges of the tetrahedron; (1 + H)/(1 - H)^2 couples all four kernel
momenta. On the Mellin route, replacing the t-density by a power makes its vertex elementary.
The vertex then becomes a modified edge of a ring of three,
E_l = (2l+1)/(sigma-1) g_l [max^(sigma-1)/(2l+2-sigma) - min^(sigma-1)/(sigma+2l)].
But the l-sum still carries l-dependent rational weights, so a recurrence in sigma needs
creative telescoping over l and the radii. That is not feasible by hand. No creative-telescoping
software is installed (sympy and mpmath only; no Sage, flint, Maple or Mathematica).

The numerical route was calibrated, and it fails the controls (`ring_ext_QR_data.py`,
`ring_ext_QR_search.py`; `data/ring_ext_QR_data.npz`, `data/ring_ext_QR_search.json`).
I computed theta^j f for j <= 10 at 64 points, t in [0.25, 4], for S, the open chain of four
C4, and the controls G and C3. The search covered every (order r <= 10, t^2-degree d) with at
most 44 unknowns. It trained on t <= 1.6 and tested by extrapolation to t >= 2.4.
- G (15 unknowns): its known operator is the best candidate (extrapolation residual 1.8e-10),
  but only marginally ahead of spurious candidates (1e-9 to 1e-8).
- C3 (24 unknowns): its known operator W is NOT found. At 24 unknowns the training matrix
  already has several singular values near 1e-16, because smooth columns are nearly dependent.
  The true null vector is not isolated in float64.
- C4 and S: nothing in range beats the spurious level (about 1e-8).
So the float64 search reaches about 15 to 20 unknowns, below the chain of three's operator,
which carries the K3 rung. A "none" for S is uninformative: the same search says "none" for a
known K3 object. The positive control fails, so by registration rule 3 no negative may be read.
The open-chain-of-four control was evaluated (C4 and all its theta-derivatives) but could not
be run to a rung, for the same reason. Its exact operator by D-finite closure has order
<= 18 (y1 = (2/pi) chi, y2'' = -(4/pi^2) K0 chi, y3'' = -(2/pi) K0 y2,
C4 = (2t^3/pi) int r K0(tr) y3); it was not built.

## 5. What would unblock the rung test

1. Precision. The smallest natural singular values of the smooth-column matrices fall roughly
   like 10^(-0.8 n) with n unknowns (G, n = 15: 1e-14; C3, n = 24: 1e-16). If S's operator is
   about as large as the chain's Mellin operator (order about 7, t^2-degree about 6, n about
   56), the data need about 60 digits, more if it is larger. The evaluator is all-positive,
   nested and cancellation-free, so it ports directly to multiprecision. In mpmath, the float64
   run (535 s for 704 values) would take days. A vectorised arb/flint or double-double port
   would make it practical.
2. Exact route. Creative telescoping (ore_algebra in Sage, Koutschan's HolonomicFunctions,
   Lairez's periods) or Griffiths-Dwork on the 4-D algebraic integrand of route (b). That
   integrand is explicit and algebraic in (k, t), with t only in w_d and h_d.
3. A closed form for the shell function Omega(l). If it is elementary, S is a 4-D integral of
   algebraic-times-log functions.
4. The open chain of four's exact operator (D-finite closure, order <= 18). It is the control
   for route 2, and it would test Paper 59's cap that the 3-D four-density integral reduces
   into the 1-D family.

## QR2 (one centre, 1s Slater)

Not run, because it is not cheap. R_l is rational for each l, but the 60+ digits an
integer-relation search needs require an exact or high-order tail of the l-series
(terms ~ l^-6). No evidence offered.

## Tiers and caps

- K kernel 1-D form: SYMBOLIC derivation + MEASURED (5 points, < 1e-13).
- 3-D box closed form: SYMBOLIC (sympy) + MEASURED (4 configurations, < 1e-15). Prior art NOT
  VERIFIED.
- Ring representations (partial waves; kernel momenta with (1+H)/(1-H)^2; the Omega form):
  derived here, MEASURED by two-route agreement at two t (2e-8) and on the triangle (1e-15).
  The Omega form was evaluated numerically, but S was not computed through it.
- Annihilator and rung: NOT OBTAINED. No reading among R-n-1/R-ell/R-tree/R-other is
  supported. QR1 is neither confirmed nor falsified.
- Fibre statements only; float64 throughout. Computation started after commit 8ba60ae5.

## Files

`debug/ring_ext_QR_{box,kernel,pw,mom,data,search,omega}.py`;
`debug/data/ring_ext_QR_{box,kernel,mom,search}.json`, `debug/data/ring_ext_QR_data.npz`.
Nothing outside debug/ was modified; nothing was committed.
