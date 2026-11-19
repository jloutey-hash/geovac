# Cold review of `debug/tree_attain_l_memo.md` (l > 0 part of Paper 63 conj:tree_attain)

Date 2026-11-19. Reviewer did not write the memo. Own code only (none of the memo's model, quadrature or Gaunt
routines): `debug/tree_attain_l_review_symbol.py` (exact operator for B1, sympy linear algebra over Q(x)),
`debug/tree_attain_l_review_numeric.py` (relocation identity, three root types, two cutoffs),
`debug/tree_attain_l_review_tailfast.py` (top constants from exterior tails). Data
`debug/data/tree_attain_l_review_{symbol,numeric,tailfast}.json`. No paper, test or tracked file edited.

## Verdicts

**Main claim (homogeneous P_v, degree <= 2, 3 <= n <= 8, Pi != 0 => F_n in the fibre's D-module): SOUND WITH
CAVEATS.** I broke no step. Caveats 1-4 below govern the wording of a paper sentence; none blocks the statement.

**Text correction (2), "exponent shifts by the degrees of the derivatives" is wrong for centre derivatives:
CORRECT.** Exact identity, checked numerically to 1e-9.

**Text correction (1), principal symbol is the wrong condition for non-homogeneous P_v: CORRECT AS A STATEMENT
ABOUT THE LEADING COEFFICIENT; NOT a counterexample to the conjecture, and "Pi is the right quantity" overstates.**
Pi != 0 is a sufficient condition for this proof, not a necessary one; in the memo's own example 1 the conclusion
almost certainly holds although Pi = 0 (see attack 5). The memo itself says "neither case is claimed"; the paper
edit must keep that.

## Setup, made exact before judging

Densities are 3-D functions rho_v = P_v(d) a^3 K0(a|r|)/(2 pi^2): K0 is the fingerprint *radial profile* of
each density, not a Yukawa kernel; the kernels are Coulomb. The root density is P_0(d)[t^3 K0(t|r|)], and the
"T" kernel K0(tx) is that same profile restricted to the radius after angular integration
(`prop:tree_reduction`: fibre = 4 pi int r^2 <rho_root Psi> dr). So the relocation is a genuine 3-D integration
by parts of P_0(d) off the root density onto Psi = product of the root's children's potentials, followed by the
angular average. The memo has this right.

## Attacks

| # | attack | result |
|:-:|:--|:--|
| 1 | root move: boundary terms at 0 and oo | tried, held |
| 2 | O(log^k) at 0 for degree <= 2 sources; what the transfer needs | tried, held |
| 3 | step (i): twists, L-dependent Poisson, Gaunt; parity; z^-1 | tried, held (parity is a Z/2 grading, provable) |
| 4 | step (ii): exponent and constant C | tried, held (independent asymptotics, two trees) |
| 5 | non-homogeneous examples | integrals verified; NOT counterexamples to the conjecture |
| 6 | independent least-order a for B1, exact | tried, held (identical order, degree, weight, symbol) |
| 7 | m != 0 densities | accept as "covered by the argument, not computed" |

**1. Relocation.** At most one derivative ever sits on K0 in a surface term: for degree 2, the first
integration by parts leaves eps^2 * (d K0 ~ 1/eps) * Psi, the second eps^2 * (K0 ~ log eps) * (d Psi). The
feared eps^2 * eps^-2 * Psi(0) = O(1) term would need d^2 K0 on the sphere, which no integration by parts
produces. Distributionally d_i d_j K0 has no delta part either: K0 ~ -log r and its gradient is O(1/r), so the
delta coefficient lim eps^2 * oint n_i d_j K0 is O(eps). With Psi = O(log^k) and d Psi = O(r^-1 log^k) the two
surface terms are O(eps log^k eps) -> 0; at infinity K0(tr) decays. My numerics (own radial potentials, log grid
4e5 points, Simpson):

| root (n = 3) | Psi(0) | rel. difference direct vs relocated, t = 0.7, 1.5, 2.6 | r0 dependence |
|:--|:-:|:--|:--|
| p_z root - s - p_z | 0 | 2e-15 .. 1e-14 | none |
| h_2 root - p_z - p_z (B1) | 0 | 4e-13 .. 3e-11 (grid accuracy) | none |
| lap root, two s leaves | 1 | abs. diff 4.31e-8, 4.24e-7, 2.21e-6 at r0 = 1e-8; 1000x smaller at r0 = 1e-11 | linear in r0 |

The lap difference equals 4 pi t^3 r0 Psi(0) to three digits at every t: it is exactly the cutoff surface term
oint Psi d_n K0 over the excised ball, O(eps), and vanishes in the limit. Small wording point: the memo calls
its lap offset "t-independent"; in my normalisation the absolute offset scales as t^3 (the t^3 prefactor). Not a
defect of the proof.

**2. Regularity at 0.** `prop:tree_rung` needs boundary terms of T(theta g) = -(1+theta_t)T(g) to vanish, i.e.
x K0(tx) g(x) -> 0 for every g in A G; T(z g) needs nothing. With sources O(r^-2 log^k) (the worst at degree 2:
h_l(d)K0 = (-1)^l h_l(r^) K_l, lap K0 = K0 - K1/r), channel potentials are O(log^{k+1}) for L = 0 and O(log^k)
for L >= 1 (both halves of the Green's function integral are O(log^k) for L >= 1; I checked by hand). G = x^2 *
<P_0(-d) Psi> is then O(log^k): second derivatives of Psi are O(r^-2 log^k) (product rule, and lap Phi = -W).
theta and z do not worsen a log-power expansion. Condition met. This uses that every function has a
sum c r^alpha log^j r expansion at 0, which `prop:tree_rung_l` already asserts; the memo inherits it.
Degree 3 would give r^-3 sources and break both this and `prop:tree_rung_l`, as the memo says.

**3. Step (i).**
- *Symbol unaffected by the l > 0 additions.* In the memo's filtration F_w the kept potentials carry weight 0, so
  theta: Phi -> theta Phi and theta(theta Phi) = -theta Phi + L(L+1) Phi - r^2 W act at weight 0 on the kept
  factor (the r^2 W part drops to the next filtration piece). The integer twists only add the weight-0 scalar p
  (theta r^p = r^p (theta + p)). So gr carries theta = Sym^k[[0,Z],[1,0]] (x) 1 and the characteristic
  polynomial is chi_k^(2^#kept): roots at Theta = (k-2j)X, |k-2j| <= n-2. Held. No z^-1 is ever needed inside
  a: A acts on r^p B for any integer p; negative even twists in the pure part are cleared by z^J outside a, and
  phi(z)^J only adds a power of t to the leading coefficient.
- *Parity.* The memo argues it; it is in fact a grading and can be proved in three lines: give a model element
  r^p * prod(Phi_{y,L_y} or theta Phi_{y,L_y}) sitting in channel L the degree p + sum L_y - L mod 2. theta and
  z preserve it; Gaunt keeps L = l_w + sum L_c mod 2; a density h_l(d) lap^j has twist -l-2j in channel l
  (degree 0); d_z moves twist and channel by one each; opening Phi_{y,L} replaces it by r^2 W_{y,L}, of the same
  degree. G is homogeneous of degree 0 (channel 0, start p = 0), hence so is a G, and on the pure top (no kept
  potentials, channel 0) degree 0 means p even. Why it is needed: q must lie in C[z]<theta> for
  T(q K0^{n-1}) = phi(q) F_n; an odd twist would give an odd moment (the banana pieces of `prop:odd_banana`).
  This is the same bookkeeping as `prop:odd_banana` (i), applied after relocation. Held; I recommend the paper
  state it as a grading, not as an observation.
- Least-order a is computed mod one prime for n = 4, 5 and exactly only for n = 3; the argument, not the
  computations, carries n = 6..8. Same status as the s-type proposition.

**4. Step (ii).** h_l(d) K0 = (-1)^l h_l(r^) K_l(r) exactly, and K_l ~ K0 ~ (pi/2x)^(1/2) e^-x, so a centre
derivative changes the angular factor, never the power of x; lap K0 = K0 - K1/r likewise. Independent
asymptotics (own exterior tails E = Phi - moment term, Gauss-Laguerre on [x, oo), Richardson x = 50..400):

| tree | Pi | C predicted | C computed | rel | local exponent of G^top e^{2x}, x 400 -> 800 |
|:--|:-:|:--|:--|:-:|:-:|
| B1: h_2 root - p_z - p_z | 8/15 | pi/15 = 0.2094395102 | 0.2094395098 | 2e-9 | 0.9995 |
| D: s root, two h_2 leaves | 1/5 | pi/10 = 0.3141592654 | 0.3141592431 | 7e-8 | 0.990 |

The exponent is (5-n)/2 = 1 with O(1/x) drift, not 1 + (degrees). B1 also confirms the factor
(-(n-1))^{d_0} in Pi: dropping it (principal-symbol integral alone) predicts pi/60, off by 4.

**5. Non-homogeneous examples.** Verified: <mu^2 - 1/3> = 0 while int omega_z^2 = 4 pi/3; <1 - omega_z> = 1
while int omega_z = 0. Both are also obvious once one notes that next to s densities only the leaf's L = 0
channel enters the fibre:
- leaf d_z^2 - 1/3: L = 0 channel = (1/3) lap K0 - K0/3 = -K1/(3r) = (1/3) r^-2 theta K0, nonzero, tail one
  power of x lower. The fibre equals that of an s-type tree with a twisted radial leaf. Step (i) is untouched,
  and step (ii) runs with exponent (5-n)/2 - 1 and a nonzero constant; the epsilon-expansion of
  int x^{(2-n)/2} e^{-eps x} is non-analytic for every n (eps^-1/2, -log eps, eps^1/2, eps log eps, ...). So the
  conjecture's CONCLUSION very plausibly holds here (by the same proof one order down), although Pi = 0.
  This is no counterexample to the conjecture as written (hypothesis true, conclusion true).
- leaf 1 + d_z: fibre = the s-type fibre (L = 0 channel is K0), so it reaches the rung by `prop:star_attain`,
  though the conjecture's hypothesis fails. A failed hypothesis does not falsify an "if" statement.
I also tried to make the conjecture false (principal integral != 0 but fibre identically 0) with degree <= 2
leaves on an s root: a single non-s vertex kills the fibre only if its L = 0 part c0 K0 + c2 lap K0 vanishes,
forcing c0 = c2 = 0 and hence a zero principal integral; two leaves (d_z + a, d_z + b) give
a b phi0^2 + phi1^2/3, never identically zero. No counterexample found. Verdict: the principal symbol is not the
leading asymptotic coefficient for non-homogeneous P_v (lower-degree parts sit at the same exponential order),
so the conjecture's stated hypothesis is the wrong *mechanism*; whether the conjecture is true for
non-homogeneous P_v is open, and the examples suggest it is (both reach the rung).

**6. Independent operator for B1.** My own model, built by hand from the Poisson equations and the d_z ladder
(<h_2(d)(f P_2)> = (1/5) r^-2 (theta+1)(theta+3) f, G = (1/5)(theta+1)u, u = (theta+3)Phi_2), eight coordinates
over Q(x): least theta-order 5, z-degree 4, weight 9, coefficients even in x, symbol term t^w * lead(phi(a)) =
-4t(t-1)^2(t+1)^2, a G pure and nonzero with even x-powers only (lowest x^2, so even J = 0 works). Identical to the
memo's row B1. One top constant (B1) recomputed by my own quadrature: attack 4.

**7. m != 0.** The radial operators depend only on L (gradient of f Y_LM has radial factors (theta - L) f/r and
(theta + L + 1) f/r for every M), Gaunt coefficients with M obey the same parity rule, and Pi is a sphere
average that does not care about axial symmetry. Accept as "the proof covers it; computations are axial only".
The paper must not say "checked" for m != 0.

## Required caveats (for the paper sentence)

1. Scope: homogeneous P_v (each P_v a sum of h_l(d) lap^j with l + 2j = d_v fixed), d_v <= 2, centre derivatives
   only, 3 <= n <= 8, t on the root, other scales 1, condition Pi != 0 (equivalently, for homogeneous P_v, the
   conjecture's angular integral != 0). Tier: [SYMBOLIC PROOF] for the argument, [MEASURED] for the checks.
2. State that the relocation (step 0) is the one new step and that its boundary terms vanish because a degree-2
   derivative puts at most one derivative on K0 in any surface term.
3. Computations: exact (over Q) only at n = 3 (two trees, now three with this review's B1 recomputation); mod one
   prime at n = 4, 5; none at n >= 6; axial densities only.
4. Parity: state it as a Z/2 grading (proof above), or cite `prop:odd_banana` (i) and note it applies after
   relocation.
5. Non-homogeneous P_v: say the conjecture's principal-symbol hypothesis is not the leading coefficient there and
   that Pi != 0 suffices for the proof; do not say "the conjecture fails" or "Pi is the right condition".
   Pi = 0 cases are not covered.

## Proposed scope wording

"For homogeneous centre-derivative densities (each P_v a combination of h_l(d) nabla^{2j} with l + 2j = d_v <= 2),
3 <= n <= 8, t on the root and every other scale one, the conjecture holds: if the angular integral is nonzero,
F_n lies in the D-module of the fibre. The proof moves the root's derivatives onto the rest of the integrand by
integration by parts in three dimensions, after which the three steps of Proposition star_attain apply with the
exponent (5-n)/2 unchanged. For non-homogeneous P_v the principal symbol does not give the leading coefficient;
the proof needs instead Pi = <P_0((n-1)omega) prod_{v != root} P_v(-omega)> != 0, which is sufficient but not
necessary. The operator of step (i) is computed exactly for n = 3 and modulo a prime for n = 4, 5, axial densities
only."
And replace "for l > 0 the exponent of the non-holomorphic term shifts by the degrees of the derivatives" by
"centre derivatives leave the exponent unchanged (h_l(d)K0 = (-1)^l h_l(r^)K_l, with the tail of K0); only scale
derivatives raise it, by one power of x each."

## Wrong answers a guard should reject

- Relocation: a direct-minus-relocated difference that does not shrink with the cutoff r0 (an O(1) boundary term),
  for a lap root with Psi(0) != 0; the right difference is 4 pi t^3 r0 Psi(0) -> 0.
- B1 symbol with a root at t = +-2 (= n-1), or with root multiplicities other than 0^1, (+-1)^2; theta-order != 5.
- B1 top constant pi/60 (principal-symbol integral without the (-(n-1))^{d_0} factor) instead of pi/15; D top
  constant other than pi/10.
- Top exponent 1 + d (shifted by derivative degrees: 3 for B1, 5 for D) instead of (5-n)/2 = 1.
- A pure coordinate of a G with odd r-power.
- For leaf d_z^2 - 1/3 next to s densities: Pi = 0 and principal integral 4 pi/3 (a guard that takes either
  equal); any assertion that this leaf falsifies the conjecture.
