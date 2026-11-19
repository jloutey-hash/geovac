# tree_attain_l: the l > 0 part of Paper 63 conj:tree_attain

Date 2026-11-19. Scripts `debug/tree_attain_l_model.py` (module model), `debug/tree_attain_l_algebra.py` (step (i),
mod p), `debug/tree_attain_l_numeric.py` (relocation, transfer), `debug/tree_attain_l_tail.py` (step (ii));
data `debug/data/tree_attain_l_{algebra,numeric,tail}.json`. No paper, test or tracked file edited.

## Verdict

**GO, for homogeneous densities; one new step, read by nobody else yet.** Goal 1 (leaf-only) and goal 2 (root
with l > 0) both close, by the same proof. The obstacle the paper names for goal 2 ("a root with l > 0 makes the
fibre a sum of several transforms") disappears: move the root's derivatives onto the rest of the integrand by
three-dimensional integration by parts. The fibre is then one transform t^3 T(G), and the three steps of
`prop:star_attain` (as written in the paper since v6.7.39: (i) operator with good symbol, (ii) split on [1, oo),
(iii) simplicity of D F_n) run with G in place of H.

## Statement proved

Take a Coulomb tree of n densities, 3 <= n <= 8, at the fingerprint weights, t on one vertex (the root), every
other scale 1, each density rho_v = P_v(d) a_v^3 K0(a_v r)/(2 pi^2) with P_v a polynomial in the center
derivative d **homogeneous of degree d_v <= 2** (any combination of harmonics h_l(d) and lap with l + 2j = d_v).
If int_{S^2} prod_v sigma_v(-omega) d omega != 0 (the conjecture's hypothesis), then F_n lies in the D-module of
the fibre, so the minimal operator of the fibre has a composition factor equivalent to the operator of F_n. For
n >= 9 the two modules meet in a nonzero submodule.

More generally, without homogeneity, the same proof needs

  Pi := < P_0((n-1) omega) prod_{v != root} P_v(-omega) >_{S^2} != 0,        (full polynomials, lap -> |omega|^2)

which for homogeneous P_v equals (-(n-1))^{d_0} <prod_v sigma_v(-omega)>. For non-homogeneous P_v the
conjecture's hypothesis (principal symbol) and Pi != 0 are different conditions in both directions (e.g. a leaf
d_z^2 - 1/3 next to s densities: principal-symbol integral 4 pi/3, Pi = 0; a leaf 1 + d_z next to s
densities: principal-symbol integral 0, Pi = 1). Neither case is claimed. If the paper keeps "principal symbol",
the conjecture's statement should either restrict to homogeneous P_v or use Pi.

**Vanishing angular integral: no claim.** When Pi = 0 the leading coefficient of the top term cancels and step
(ii) as written gives nothing. Step (i) does not care (the cancelling root shape B2 has a good operator, below),
so the proof would go through whenever some later coefficient of the top term's asymptotic series is nonzero; the
records' eps^{-1/2} for B2 suggests that, but nothing is claimed here for Pi = 0.

## Proof

Hat normalisation: non-root rho_v = P_v(d) K0(r); channel-L potential Phi_{w,L} of the subtree of w solves
(theta - L)(theta + L + 1) Phi_{w,L} = -r^2 W_{w,L}, theta = r d/dr, W_{w,L} the channel-L part of
rho_w prod_{c child} Phi_c. Constants (2 pi^2, 4 pi/(2L+1)) are dropped; they are nonzero and do not move
modules.

**(0) Relocation.** The fibre is int rho_0 Psi d^3r, Psi = prod_{w child of root} Phi_w (Fubini, as in
`prop:tree_reduction`). Integrate by parts in three dimensions:

  int [P_0(d) t^3 K0(t|r|)] Psi d^3r = t^3 int K0(t|r|) [P_0(-d) Psi] d^3r = 4 pi t^3 T(G),
  G(x) = x^2 < P_0(-d) Psi >_{S^2}(x),  T(g) = int_0^oo K0(tx) g dx.

Boundary terms: on the sphere |r| = e, K0 ~ log e, d K0 ~ 1/e, Psi = O(log^k e) and d Psi = O(e^-1 log^k e)
(potentials of sources O(r^-2 log^k), the worst a degree-2 density gives, are O(log^{k+1}) in every channel,
as in the proof of `prop:tree_rung_l`); each term is e^2 x O(e^-1 log^k e) -> 0. At infinity K0(t r) decays for
Re t > 0. Channel by channel, d_z [f P_L] = (L+1)/(2L+1) r^-1 (theta - L) f P_{L+1} + L/(2L+1) r^-1
(theta + L + 1) f P_{L-1}, and lap [f P_L] = r^-2 (theta(theta+1) - L(L+1)) f P_L, so G is a finite sum of
r-twisted polynomials in theta applied to channel parts of Psi: an element of the same function module, with
G = O(log^k) at 0 (degree <= 2 gives at most r^-2, cancelled by x^2). So every function in A G is O(log^k) at 0
and the transfer T(b g) = phi(b) T(g), b in A = C[z]<theta>, holds with vanishing boundary terms, as in
`prop:tree_rung`. The root's own derivatives are gone: what is left is an s-type root over a module with l > 0
content. (For the conjecture's class there is no d_a; a root d_t^i would be handled by the remark that D(d^i y)
and D y have the same factors of order > 1.)

**(i) Operator with good symbol.** The module model is the memo's (`debug/tree_attain_memo.md`, Lemma 1) with
three changes, none of which touches the symbol. Coordinates: (k, i, kept, p) with r^p K0^{k-i} D^i prod_kept
(Phi_{y,L} or theta Phi_{y,L}), D = theta K0, k the number of opened vertices, p in Z bounded below.
- Radial parts of densities lie in r^{-d_v} B (Hobson: h_l(d) f = h_l(r)(r^-2 theta)^l f; lap = r^-2
  theta(theta+1)), so products of k of them lie in r^lambda Sym^k B.
- Potential pieces: theta (theta Phi_{y,L}) = -theta Phi + L(L+1) Phi - r^2 W_{y,L}, with W_{y,L} (Gaunt
  re-expansion, constant coefficients) in the next filtration step.
- Twists r^p: theta r^p = r^p (theta + p).
Filter by k. On the piece N^(k)/N^(k+1) give z weight 2, theta weight 1, F_w = span{r^{e+2m} K0^{k-i} D^i pi :
2m + i <= w}: finite-dimensional, bounded below. The scalars p, L(L+1), -1 and the potential map
Phi -> theta Phi have weight 0, so on gr theta acts as Sym^k [[0, Z],[1, 0]] tensor 1 and gr is killed by
chi_k^{2^{#kept}}, chi_k = prod_j (Theta - (k-2j) X). The good-filtration argument of Lemma 1 (checked by the
review, attack 1) then gives a in A, a product of operators of symbol chi_k^m, k <= n-2, with a G in the pure
top N^(n-1), the C[z]-span of r^p K0^{n-1-i} D^i.

Parity: give Phi_{y,L} parity L and r^p parity p. Densities have twist parity = l (d_v = l + 2j), Gaunt keeps
L = l_w + sum L_c mod 2, Poisson keeps it, d_z changes twist and channel by one each, and the root takes channel 0
of something whose model elements all have twist parity = channel parity. So every coordinate of G, hence of
a G, has p even, and z^J a G = q K0^{n-1}, q in A, for J large (K0^{n-1} is cyclic in Sym^{n-1} B over C[z],
integer pivots). By the symbol rule (*), phi(z^J a) has leading coefficient a power of t times prod (t + k - 2j)
over k <= n-2: nonzero at t = +-(n-1). (Odd p never needs an odd element of the ring: r^p B with p odd is an
A-module like any other.)

**(ii) The split on [1, oo).** For x > 0, Phi_{y,L} = (M_{y,L}/(2L+1)) x^{-L-1} + E_{y,L}, M = int_0^oo s^{L+2} W,
E the exterior tail, which obeys the same equation. Expanding recursively, G = sum_S q_S G^S, constants q_S.
For |S| <= n-2, G^S lies in a module of the kind in (i) with pieces k <= |S| (the moment terms are twisted
solid harmonics, and P_0(-d) of them stays so), so int_1^oo K0(tx) G^S dx obeys an equation with singular points
at 0 and integers |t| <= n-2 and a right side (boundary terms at x = 1) holomorphic on the log cover: its
continuation through the upper half plane is holomorphic at p = -(n-1). The term S = all non-root vertices is,
by the local rule (leading tail -s_0/k^2 in every channel, so pointwise in omega with finitely many channels,
uniformly), and d -> -(n-1) omega on e^{-(n-1)x} for the root's -d,

  G^top = C x^{(5-n)/2} e^{-(n-1)x} (1 + O(1/x)),
  C = (-1)^{n-1} (pi/2)^{(n-1)/2} Pi / prod_{u != root} k_u^2.

The exponent is the s-type one: center derivatives do not shift it (h_l(d) K0 = h_l(-omega) K_l at leading
order, and K_l, K0 have the same leading tail). The paper's sentence "the exponent of the non-holomorphic term
shifts by the degrees of the derivatives" is true only of scale derivatives d_a (each adds one power of x); it
is not true of the center derivatives the conjecture is about. With C != 0 (Pi != 0) the rest of step (ii) is
verbatim: the imaginary part at t = -(n-1-eps)+i0 has integrand C (2 pi (n-1))^{-1/2} x^{(4-n)/2} e^{-eps x}(1 +
O(1/x)), of one sign, so J (n <= 6) or J' (n = 7, 8) is unbounded and T(G) is not holomorphic at p.

**(iii)** As in the paper: if phi(q) F_n were zero, T(G) would solve phi(z^J a) y = 0, regular at p, and be
holomorphic there. So phi(q) F_n is a nonzero element of D T(G) = D f and of D F_n, which is simple for n <= 8.

## Checks

**Step (i), least-order a with a G pure, mod 2147483629** (`tree_attain_l_algebra.py`). In every row: the
symbol at (1, -t) splits into integer roots of absolute value <= n-2, is nonzero at t = -(n-1), a G has no
non-pure coordinate left and is nonzero, and every pure coordinate has even r-parity.

| tree (root = 0) | n | theta-order | z-degree | weight | roots of the symbol (multiplicity) |
|:--|:-:|:-:|:-:|:-:|:--|
| control: s star of 3 | 3 | 7 | 5 | 13 | 0^3, +-1^2 (memo: -4t^3(t-1)^2(t+1)^2) |
| control: s path of 3, end | 3 | 6 | 3 | 8 | 0^2, +-1^2 (memo: weight 8) |
| A: p_z - s(root) - p_z | 3 | 7 | 5 | 13 | 0^3, +-1^2 |
| B1: d(root) - p_z - p_z | 3 | 5 | 4 | 9 | 0, +-1^2 |
| B2: cancelling root - p_z - p_z (Pi = 0) | 3 | 5 | 3 | 7 | 0, +-1^2 |
| C: p_z(root) between p_z and s | 3 | 11 | 14 | 31 | 0^3, +-1^4 |
| D: d - s(root) - d | 3 | 7 | 5 | 13 | 0^3, +-1^2 |
| E: s(root) - p_z - p_z | 3 | 6 | 3 | 8 | 0^2, +-1^2 |
| F: star 4, s root, leaves p_z, p_z, s | 4 | 32 | 122 | 254 | 0^10, +-1^7, +-2^4 |
| G: star 4, p_z root, leaves p_z, s, s | 4 | 31 | 105 | 219 | 0^9, +-1^7, +-2^4 |
| H: path 4, s root at end, p_z, p_z inside | 4 | 12 | 13 | 30 | 0^4, +-1^2, +-2^2 |
| I: star 4, s root, leaves d, d, s | 4 | 32 | 133 | 276 | 0^10, +-1^7, +-2^4 |
| J: path 4, d root at end, p_z, p_z, s | 4 | 11 | 16 | 35 | 0^3, +-1^2, +-2^2 |
| K: lap root, s and d_z^2 nested, s leaf | 4 | 21 | 53 | 111 | 0^5, +-1^4, +-2^4 |
| P: path 5, d root at end, p_z, s, p_z, s | 5 | 19 | 48 | 99 | 0^3, +-1^4, +-2^2, +-3^2 |

(The stars of four with two different leaf types need order 32 against 16 for the s star of four: identical
leaves make the product symmetric and the module smaller.) Exact over Q for A and B1
(`tree_attain_l_numeric.py`, N3): symbols -4 t^3 (t-1)^2 (t+1)^2 and -4 t (t-1)^2 (t+1)^2, matching the mod-p
roots.

**Model validation** (float64 log grid, 4e5 points; independent closed forms for the densities, channels of the
potentials by numerical mu-quadrature, not the model's Gaunt code):
- N1: model density channels of d_z, d_z^2, h_2, lap and the cancelling root applied to K0, against
  hand-differentiated closed forms (d_z K0 = -mu K1, d_z^2 K0 = mu^2 K0 + (2mu^2-1)K1/r, lap K0 = K0 - K1/r):
  worst relative difference 1.5e-15.
- N2, the relocation (direct: root derivative on K0(tr), no model; relocated: int K0(tr) G dr, G the model
  element evaluated), at t = 0.7, 1.5, 2.6: A, C, G (p_z roots, s root) agree to 2e-15 to 4e-15; B1, J (d roots,
  nested) to 2e-14 to 3e-14; B2 and K (roots containing lap) to 1e-9 to 7e-9, a t-independent offset that is the
  grid cutoff r_0 = 1e-9 (offset 8.88e-11 at r_0 = 1e-9, 8.88e-13 at r_0 = 1e-11: the monopole boundary term
  ~ r_0 Psi(0), which vanishes as r_0 -> 0, as the proof says).
- N3, the transfer phi(a) T(G) = T(a G) with the exact a, phi applied under the integral: B1 at t = 1.5:
  -19.95348076605 both sides (5e-15); t = 2.5: -2.501091081271 both (2e-15). A at t = 1.5: -3771.468070448 vs
  ...446 (6e-13); t = 2.5: 40.32393308153 vs ...142 (3e-12, largest term 2.9e3).

**Step (ii): the top constant** (`tree_attain_l_tail.py`, N4; nested Gauss-Laguerre in exponentially scaled
variables, 48 nodes, channels by 16-point mu-quadrature, root derivatives via theta E and the channel ODE;
Richardson on x = 100 .. 1600; no model and no local rule):

| tree | n | Pi = <P_0((n-1)w) prod P_v(-w)> | prod k_u^2 | C predicted | C computed | rel |
|:--|:-:|:-:|:-:|:--|:--|:-:|
| A | 3 | 1/3 | 1 | 0.523598775598 | 0.523598775598 | 5e-14 |
| B1 (d root) | 3 | 8/15 | 4 | 0.209439510239 | 0.209439510239 | 2e-13 |
| F | 4 | 1/3 | 1 | -0.656233747738 | -0.656233747739 | 5e-13 |
| G (p_z root) | 4 | -1 | 1 | 1.96870124322 | 1.96870124322 | 8e-13 |
| J (d root, nested) | 4 | 6/5 | 36 | -0.0656233747738 | -0.0656233747742 | 5e-12 |
| K (lap root, d_z^2 nested) | 4 | 3 | 4 | -1.47652593241 | -1.47652593241 | 4e-13 |
| L: fork 5, d root, p_z inside | 5 | 32/15 | 64 | 0.0822467033424 | 0.0822467033427 | 4e-12 |
| M: star 5, p_z root, leaves p_z, d, s, s | 5 | -8/15 | 1 | -1.31594725348 | -1.31594725348 | 3e-12 |
| N: spider 6, d leaf at the branch, p_z inside | 6 | 2/15 | 225 | -0.00183255032972 | -0.00183255032972 | 2e-12 |

Consistency with the records: A/(s chain) = 1/3 (the p_z chain's "one third"); B1's Pi = 8/15 = (n-1)^2 x 2/15,
the record's 2/15 being stated with the root's t^l factor removed.

**Step (ii): the branch** (N5; J(eps) on a log grid to x = 2e4, eps = 1e-1 .. 1e-3, least squares with the next
terms of the expansion): F (n = 4, pole) predicted singular coefficient -0.15114995, fitted -0.15112068 (2e-4);
M (n = 5, eps^{-1/2}, p_z root) predicted -0.46525761, fitted -0.46528585 (6e-5). The fit precision is limited
by the float grid; the constants themselves are checked to 1e-12 above.

## What remains

- **Non-homogeneous P_v.** Proved under Pi != 0; the conjecture's principal-symbol hypothesis is a different
  condition there (examples above). Not claimed in either direction.
- **Pi = 0.** Not claimed. Step (i) holds (B2 computed); step (ii) would need the first nonzero coefficient of the
  top term's series, and a proof that one is nonzero.
- **Scale derivatives d_a** (in `prop:tree_rung_l` but not in the conjecture): the same proof should work with
  the exponent raised by the top d_a-degree I of the non-root vertices (each d_a gives -x at leading order), and
  a root with a single power d_t^i; a root mixing powers of d_t gives a sum of t^{3-i} T(b_i H_i) with different
  powers of t, which this memo does not handle. Not checked.
- **n <= 8** is inherited from the simplicity of D F_n (Sec. four_density; n = 6, 7, 8 by the exponent argument).
- **No exact check** of step (i) at n >= 6 with l > 0, and only mod one prime at n = 4, 5; exact over Q only for
  A and B1 (n = 3). Non-axial densities (m != 0, e.g. p_x) are covered by the proof but not computed.

## Where this could be wrong

1. Step (0) is new. Its boundary terms are argued from the O(log^k) behaviour of potentials at the origin, which
   rests on degree <= 2; the numerical relocation check covers p_z, d, lap and cancelling roots on nested trees
   (agreement to 2e-15 .. 3e-14, the lap cases to the cutoff, scaling as the cutoff). A degree-3 root would
   break it (r^-3 sources are not integrable), as it breaks `prop:tree_rung_l`.
2. Step (i) reuses Lemma 1's good-filtration argument, reviewed once. The l > 0 changes (twists, channel
   potentials, the extra weight-0 scalars) are argued not to touch gr; fifteen trees, mod one prime, agree.
   A wrong reduction mod p would need the least order to drop modulo that prime.
3. The parity bookkeeping (pure top has even twist) is argued and observed (parities [0] in every row); if it
   failed, z^J a G would not be a C[z]-combination and q would not exist in A.
4. The local rule's uniformity in omega uses that each potential has finitely many channels (degree <= 2 at
   each of n vertices gives L <= 2n); fine for fixed trees, not a uniformity statement in n.
5. The model is a free-generator model; relations among the true functions only make the module a quotient, and
   N3 confirms the model numerically for two n = 3 trees (one nested, d root).
6. The claim that the paper's "exponent shifts by the degrees of the derivatives" is wrong for center
   derivatives rests on the leading tail h_l(-omega) K_l; the table's exponent (5-n)/2 for p_z and d roots and
   leaves, fitted to 1e-12 in C with that exponent fixed, supports it.
7. The non-homogeneous examples are mine and computed only at the level of the angular integral.

## Files

- `debug/tree_attain_l_model.py`, `debug/tree_attain_l_algebra.py`, `debug/tree_attain_l_numeric.py`,
  `debug/tree_attain_l_tail.py`
- `debug/data/tree_attain_l_algebra.json`, `debug/data/tree_attain_l_numeric.json`,
  `debug/data/tree_attain_l_tail.json`
- `debug/tree_attain_l_memo.md`
