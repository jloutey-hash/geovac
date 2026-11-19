# tree_attain: the s-type part of Paper 63 conj:tree_attain

Date 2026-11-18. Scripts `debug/tree_attain_algebra.py`, `debug/tree_attain_tail.py`; data
`debug/data/tree_attain_algebra.json`, `debug/data/tree_attain_tail.json`. No paper, test or tracked file edited.

## Verdict

**GO, with one stated condition.** For every Coulomb tree of n s-type densities at the fingerprint weights, t on
any one vertex, 3 <= n <= 8, the minimal operator of the fibre has a composition factor equivalent to the operator
of F_n. In fact F_n itself lies in the module of the fibre: F_n = B f for an operator B with rational coefficients.

The condition: the proof below is **not** "the star proof run as written". Two things had to be added.

1. The step the paper names as not carried out (the split on [1, oo)) is carried out here (Lemma 2). It works as the
   paper expected.
2. The last inference of the star proof ("quotients and extensions of such modules are connections there ... so a
   piece with c = n-1 must survive") does not hold as written for even n, for stars as well as trees. Composition
   factors over C(t) do not control singular points: d/dt composed with (t-p) d/dt has two factors, both equivalent to
   d/dt, and the solution log(t-p). A top piece of the filtration can be nonzero as a module over operators regular at
   p and still vanish over C(t). For odd n the inference can be rescued (factors with meromorphic local solutions give
   integer exponents, and the branch is a half-integer power); for even n (a pole at n = 4, logarithms at n = 6, 8)
   it cannot, as written. n = 4 is covered by the paper's computed factorisations. Lemma 1 below replaces that step
   for every n, stars included.

Both lemmas are short and are proved in this memo; they are new text, read by nobody else yet. That is the
remaining risk, and it is why this is not reported as an unqualified GO.

## What the star proof uses

(a) fibre = (2/pi)^n t^3 T(H), T(g) = int_0^oo K0(tx) g dx, H = x^(3-n) chi^(n-1) (Prop. tree_reduction);
(b) chi = pi/2 - R, so on [1, oo) H is a sum of terms R^k x^(3-n), k <= n-1;
(c) terms with k <= n-2 carry at most n-2 Bessel factors and their transforms on [1, oo) are analytic at t = -(n-1);
(d) the term k = n-1 gives, at t = -(n-1-eps)+i0, -i pi int_1^oo I0((n-1-eps)x) R^(n-1) x^(3-n) dx, integrand
    c_n x^((4-n)/2) e^(-eps x)(1+O(1/x)), c_n != 0: not holomorphic at eps = 0;
(e) "so a piece with n-1 Bessel factors survives" (the inference criticised above);
(f) that piece has even shift (Prop. odd_banana (i)) so is equivalent to F_n's operator, irreducible for n <= 8
    (Sec. four_density; n = 6, 7, 8 "by the exponent argument", guard tests/test_paper63_four_density_rung.py).

## Setting

Hat normalisation: rho = K0, (Phi W)(x) = (1/x) int_0^x W s^2 ds + int_x^oo W s ds, so theta(theta+1) Phi W = -x^2 W
(theta = x d/dx). Root the tree at the vertex carrying t. For a non-root vertex u: W_u = K0 prod_{w child of u} P_w,
P_u = Phi W_u, k_u = number of vertices in the subtree of u. H = x^2 prod_{w child of root} P_w, and
fibre f = (2/pi)^n t^3 T(H). A = C[z]<theta>, z = x^2; phi(z) = t^-2 theta_t^2, phi(theta) = -1 - theta_t;
T(bH) = phi(b) T(H) for b in A, Re t > 0 (Prop. tree_rung, "the transfer"; every theta^j H is O(x^2 log^m x) at 0).
Write D = theta K0 = -x K1, theta D = z K0. p = -(n-1).

Weights: theta has weight 1, z weight 2. The top-weight part of b in A is its symbol sigma(b)(X, Theta) (z -> X^2).
Symbols multiply. phi(z^k theta^j) has order 2k+j in theta_t with leading coefficient (-1)^j t^(-2k), so

  (*)  phi(b) has order w = weight(b) and leading coefficient t^(-w) sigma(b)(1, -t).

chi_c := prod_{j=0..c} (Theta - (c-2j) X). [EXACT, script] the top-weight part of L_c (the annihilator of K0^c) is
chi_c for c = 1..7; this is the paper's leading coefficient prod (t^2 - (c-2j)^2), read through (*).

## Lemma 1 (an operator with good symbol that sends H to n-1 Bessel factors)

**Statement.** There is a in A, a product of operators with symbols chi_c^m, 0 <= c <= n-2, such that
a H = sum_i g_i(z) K0^(n-1-i) D^i with polynomial g_i. Hence, by (*), phi(a) has leading coefficient a power of t times
a product of factors (t + c - 2j), 0 <= c <= n-2: it does not vanish at t = +-(n-1).

**Proof.** *The function module.* For a set S of non-root vertices with S + root connected, let dS be the children of
S + root not in S, and N_S the C[z]-span of K0^(|S|-i) D^i prod_{y in dS} pi_y, pi_y in {P_y, theta P_y}. H is in
N_empty. theta acts as a derivation, theta K0 = D, theta D = z K0, and
theta(theta P_y) = -theta P_y - z K0 prod_{c child of y} P_c, whose last term replaces y by its children and adds
one K0: it lies in N_{S+y}. So N^(k) = sum_{|S| >= k} N_S is an A-module, N^(k) contains N^(k+1), and
N^(n-1) = N_{all} is the pure span of K0^(n-1-i) D^i.

*The pieces.* N^(k)/N^(k+1) is the image of the sum over |S| = k of Sym^k B tensor (C[z] p_y + C[z] p'_y, theta p = p',
theta p' = -p') over y in dS. Each such tensor product has a finite filtration with pieces Sym^k B tensor (rank one,
theta = mu), mu an integer: free over C[z] on e_i = K0^(k-i) D^i, theta e_i = mu e_i + (k-i) e_(i+1) + i z e_(i-1).

*One piece.* Let P be such a piece, F_w P the span of z^m e_i with 2m + i <= w (finite-dimensional, zero for w < 0).
theta F_w is in F_(w+1), z F_w in F_(w+2). On gr P the class of theta acts C[Z]-linearly as Sym^k of [[0, Z], [1, 0]],
with characteristic polynomial chi_k; by Cayley-Hamilton chi_k kills gr P. So the lift X of chi_k (weight k+1) has
X F_w in F_(w+k). Take g in F_(w0) P and G_w = A g meet F_w P. gr^G(A g) is a submodule of gr P, which is finitely
generated over the Noetherian C[Z, Theta]; so it is generated in degrees <= w2, and by induction on w,
G_w = sum_{v <= w2} F_(w-v)A G_v. G_(w2) is finite-dimensional, so it lies in F_(w1)A g, and G_w is in F_(w+w1)A g
for all w. Then X^m g is in G_(w0+mk), so X^m g = b g with weight(b) <= w0 + w1 + mk < m(k+1) once m > w0 + w1:
a_g = X^m - b kills g and has symbol chi_k^m.

*Assembling.* Go down the filtration: lift the class of H in the top piece, kill it, the result lies one step lower,
repeat. All pieces above N^(n-1) have k <= n-2. The product a has the stated symbol and a H is in N^(n-1). QED.

**[EXACT, script]** For each tree below the least-order a (over C(z), cleared of denominators) was computed by exact
linear algebra in this module model; its symbol at (1, -t) factors over Q, and all roots are integers of absolute
value <= n-2, as the lemma says (the least-order operator's symbol divides X^(2m) times the lemma's).

| tree, root | n | theta-order of a | weight | t^w x leading coefficient of phi(a) | a H pure and nonzero |
|:--|:-:|:-:|:-:|:--|:-:|
| star 3, centre | 3 | 7 | 13 | -4 t^3 (t-1)^2 (t+1)^2 | yes |
| path 3, end | 3 | 6 | 8 | 4 t^2 (t-1)^2 (t+1)^2 | yes |
| star 4, centre | 4 | 16 | 56 | c t^6 (t-2)^2 (t-1)^3 (t+1)^3 (t+2)^2 | yes |
| path 4, end | 4 | 12 | 32 | c t^4 (t-2)^2 (t-1)^2 (t+1)^2 (t+2)^2 | yes |
| path 4, inner | 4 | 23 | 105 | c t^7 (t-2)^4 (t-1)^4 (t+1)^4 (t+2)^4 | yes |
| path 5, end | 5 | 20 | 100 | c t^4 (t-3)^2 (t-2)^2 (t-1)^4 (t+1)^4 (t+2)^2 (t+3)^2 | yes |

(c a nonzero integer.) The exact operator for the fork of five and the path of five rooted at its centre did not
finish inside the 20-minute cap (the path of four rooted inside took 645 s); for those trees Lemma 1 rests on its
proof, and only Lemma 2 is checked numerically.

**[NUMERICAL]** The identity phi(a) T(H) = T(aH), with phi applied under the integral sign and both sides nonzero:
star 3 at t = 1.5, 2.5: -9.3648159755e4 vs -9.3648159753e4 (1.3e-11), -1.3409060459e3 vs ...58e3 (2.2e-11);
path 3 rooted at an end: 33.601343683 both sides (6e-13), 2.1830686306 (2e-12). For n >= 4 the coefficients of a are
too large for float64 (star 4: terms 1e66 against a result 1e54), so the identity is not checked numerically there;
the transfer itself is Prop. tree_rung's and is guarded.

## Lemma 2 (the split on [1, oo): the step the paper left open)

For x > 0, P_u = Q_u/x + E_u with Q_u = int_0^oo W_u s^2 ds and
E_u(x) = -(1/x) int_x^oo W_u(s) s (s-x) ds, and theta(theta+1) E_u = -x^2 W_u (the same Poisson equation, since
theta(theta+1) kills 1/x). Expanding every potential, recursively,

  H = sum_S q_S H^S,   q_S = prod_{y in dS} Q_y,   H^S = x^(2-m_S) prod_{w in S, child of root} E_w^(S),

where S runs over the same sets as in Lemma 1, E_w^(S) is the exterior tail of W_w^(S) = K0 x^(-m) prod E_c^(S)
(children c in S; m counts children not in S). H^S is O(e^(-|S| x) x^M).

**(i) |S| <= n-2.** H^S lies in a module of exactly the kind of Lemma 1 built on the subtree S (tails in place of
potentials, an overall integer power of x, which shifts theta and not the symbol). All its pieces, the bottom one
included, are Sym^k B twists with k <= |S| <= n-2. So there is a_S with a_S H^S = 0 and symbol a product of chi_k^m,
k <= n-2. Integrating by parts on [1, oo), phi(a_S) int_1^oo K0(tx) H^S dx = boundary terms at x = 1, which are
Laurent polynomials in t times theta_t^j K0(t). The equation has leading coefficient nonzero away from t = 0 and the
integers |t| <= n-2, and a right side analytic on the cover of C minus 0. So int_1^oo K0(tx) H^S dx continues along
any path avoiding those points and is holomorphic at p.

**(ii) S = all non-root vertices (one term, q_S = 1).** By the local rule (a source s_0 x^a e^(-kx)(1+O(1/x)) has
exterior tail -(s_0/k^2) x^a e^(-kx)(1+O(1/x)); one line of Watson's lemma applied to the formula for E),

  H^top = C x^((5-n)/2) e^(-(n-1)x) (1+O(1/x)),   C = (-1)^(n-1) (pi/2)^((n-1)/2) / prod_{u != root} k_u^2  != 0.

int_1^oo K0(tx) H^top dx is holomorphic for Re t > -(n-1) on the sheet reached through the upper half plane, and at
t = -(n-1-eps)+i0 equals a function analytic at eps = 0 minus i pi J(eps),
J = int_1^oo I0((n-1-eps)x) H^top dx, integrand C (2 pi (n-1))^(-1/2) x^((4-n)/2) e^(-eps x)(1+O(1/x)), of one sign
for large x. So J (n <= 6) or its first eps-derivative (n = 7, 8) is unbounded as eps -> 0+.

**(iii)** int_0^1 K0(tx) H dx is analytic on the cover of C minus 0.

**Conclusion (S).** Continued from Re t > 0 through the upper half plane, T(H) is not holomorphic at t = -(n-1):
every other term is, and the top term or its derivative is unbounded. This holds for every s-type tree and every
n >= 2 (take more derivatives for n >= 9).

**[NUMERICAL, tail script]** (a) P_u - Q_u/x equals the tail formula on 1 <= x <= 30 to 9e-14 at every vertex.
(b) e^((n-1)x) H^top against the 7-term local-rule series (the Tail class of tests/test_paper63_local_rule.py,
whose lead is asserted equal to C exactly): agreement 3e-10 to 4e-9 at x = 55, falling as x^-7; against the
leading term alone the difference is 0.5% to 6%, falling as 1/x. (c) the leading singular coefficient of J, fitted
on eps = 1e-1 .. 1e-4 with the next terms of the expansion in the basis:

| tree, root | n | prod k_u^2 | singular term | predicted C Gamma(a+1)/sqrt(2 pi (n-1)) | fitted |
|:--|:-:|:-:|:--|:--|:--|
| path 3, end | 3 | 4 | eps^-3/2 | 0.0981748 (pi/32) | 0.0981748 |
| star 4 | 4 | 1 | eps^-1 | -0.45345 | -0.45345 |
| path 4, end | 4 | 36 | eps^-1 | -0.0125958 | -0.0125958 |
| path 4, inner | 4 | 4 | eps^-1 | -0.113362 | -0.113362 |
| star 5 | 5 | 1 | eps^-1/2 | 0.872358 (pi^2/(8 sqrt 2)) | 0.872356 |
| path 5, end | 5 | 576 | eps^-1/2 | 0.00151451 | 0.00151451 |
| path 5, centre | 5 | 16 | eps^-1/2 | 0.0545224 | 0.0545225 |
| fork 5 (0-1, 1-2, 1-3, 3-4), root 0 | 5 | 64 | eps^-1/2 | 0.0136306 | 0.0136306 |
| same fork, root at the branch vertex | 5 | 4 | eps^-1/2 | 0.21809 | 0.21809 |
| spider 6 (0-1, 1-2, 1-3, 3-4, 3-5), root 0 | 6 | 64 | log eps | 0.00245212 | 0.00245346 |

Constants on record: the star of five's pi^2/(8 sqrt 2) and the star of four's (pi/2)^(3/2)/sqrt(6 pi) are the
k_u = 1 cases. The s chain's 1/pi^2 (test_T3, chain leaf - root - leaf) is the star of three: c_3 Gamma(3/2) = pi/8
in the hat normalisation, times (2/pi)^2 (4 pi/(2 pi^2)) = 8/pi^3. The new content is the factor 1/prod k_u^2: each
vertex's potential is taken at the exponential rate of its whole subtree.

## The theorem

**Claim.** s-type Coulomb tree, n densities, fingerprint weights, slice (t; 1, ..., 1) with t on any vertex. Then
(1) for every n >= 2, the module D_t f of the fibre meets D_t F_n in a nonzero submodule;
(2) if the operator of F_n is irreducible (n <= 8), F_n is in D_t f, so the minimal operator of f has a composition
factor equivalent to the operator of F_n.

**Proof.** Take a from Lemma 1 and write a H = q K0^(n-1), q in A (K0^(n-1) is a cyclic vector of Sym^(n-1) B over
C[z], Prop. tree_rung "the pieces"). The transfer holds on both modules, so
phi(a) T(H) = T(aH) = phi(q) F_n. If this were zero, T(H) would solve phi(a) y = 0, whose leading coefficient does
not vanish at p, and would be holomorphic at p along every path: contradiction with (S). So phi(q) F_n is a nonzero
element of D_t T(H) = D_t f (t^3 is a unit) and of D_t F_n. For n <= 8, D_t F_n is simple, so it is generated by
phi(q) F_n and lies in D_t f. QED.

Every piece with n-1 Bessel factors has even shift here by construction (a is in A, which contains only z = x^2),
so Prop. odd_banana (i) is not needed for the s-type case.

**Induction statement, written out.** For a rooted subtree with root u and every admissible S inside it: (I1) the
function module generated by W_u (or P_u) is filtered by the number of opened vertices, top-count piece = pure
products of k_u Bessel functions, every other piece a quotient of Sym^k B twists with k < k_u; (I2) P_u = Q_u/x +
E_u with E_u = sum_S q_S E_u^(S), E_u^(S) of rate |S|, and E_u^full ~ -(A_u/k_u^2) x^(-k_u/2) e^(-k_u x), where
A_u = sqrt(pi/2) prod_{w child} (-A_w/k_w^2). Both pass from children to parent by multiplying by K0 and applying
the Poisson equation.

## Range of n and why

- 3 <= n <= 8: full statement. The only input that stops at 8 is the irreducibility of F_n's operator, taken from
  the paper (n = 2..5 by the Fuchs-relation count, n = 6, 7, 8 by the exponent argument; guard
  tests/test_paper63_four_density_rung.py). Not re-derived here.
- n >= 9: part (1) only (some composition factor of F_n's operator is present).
- Any root; the conjecture's slice is covered in full for s-type.

## What is not proved

- l > 0. Needed: (i) Lemma 1 with multipole potentials (order-one pieces r^L, r^(-L-1): symbol Theta, so the symbol
  argument should carry over; the pure top then has Laurent coefficients and an even shift by Prop. odd_banana (i));
  (ii) (S) with the leading coefficient times the angular integral, nonzero by hypothesis; (iii) for a root with
  l > 0 the fibre is a sum of several transforms E_i T(r^mu g_i), so one operator a no longer suffices and a common
  left multiple with controlled symbol is needed. Not attempted.
- Minimality or the full factor list of any tree's operator. Only the presence of the F_n factor.
- Nothing about the integrated value, other kernels (Yukawa, geminal), other weights, or unequal non-root scales
  (the rates then differ and the top rate is sum of scales; the argument looks the same but was not checked).

## Where this could be wrong

1. Lemma 1's symbol step is my own write-up of the standard good-filtration argument, for the filtered ring
   C[z]<theta> with weights (2, 1). The exact checks confirm its conclusion on six rooted trees (n = 3, 4 and
   the path of five from an end); no branched non-star tree and nothing at n = 6, 7, 8 was checked exactly. A reviewer should check the claim "gr^G(A g) finitely generated gives
   G_w in F_(w+w1) A g".
2. The function module is modelled abstractly (K0, D, P_y, theta P_y as free generators). Relations among these
   functions would only make the true module a quotient, which is harmless for Lemma 1, and the n = 3 numerical
   identities confirm the model for two trees; for n >= 4 the model is checked only through its rules.
3. (S) uses "holomorphic at p" for the lower terms along the specific path through the upper half plane; the
   equations for the lower terms have singular points at 0 and at integers |t| <= n-2, all avoided. If a lower term
   had rate n-1 the argument would fail; rates are |S| by construction.
4. The O(1/x) in the top integrand must be uniform in small eps; it is (I0's expansion at argument >= n-2).
5. Irreducibility for n = 6, 7, 8 is inherited, as is the transfer's boundary behaviour at 0.
6. The numerical fits in the table confirm the leading constant to 5 to 6 digits (0.05% for the n = 6 log); they
   are float64 on a grid, with the 7-term series beyond x = 50.
7. My reading that the star proof's last step is incomplete for even n is a reading of two sentences; if the authors
   meant the argument of Lemma 1 it was not on the page.

## Files

- `debug/tree_attain_algebra.py`, `debug/data/tree_attain_algebra.json` (+ `.log`)
- `debug/tree_attain_tail.py`, `debug/data/tree_attain_tail.json` (+ `.log`)
- `debug/tree_attain_memo.md`
