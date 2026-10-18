# Minimality audit of Paper 59's operator-based rung claims (2026-10-17)

**Task.** Check every claim in Paper 59 of the form "f satisfies L" or "f's rung is X because L
is/contains Y", and seam rows J4, J7, J9, J10, J11, for the two defects the ring sprint exposed:
(a) a source term the regularized Mellin transform cannot see (L f = pure powers, not 0);
(b) a rung read off a non-minimal (reducible) operator.

**Verdict: GO. No DEFECT.** Two claims that rested on bounded searches are now settled without a
bound (irreducibility of the F_4 and F_5 operators). Two checks the ring and chain readings needed,
and did not carry, were run and pass. One small sharpening for the ring is reported (Section 4).

Scripts (all new, in `debug/`): `minimality_audit_fuchs.py`, `minimality_audit_ring.py`,
`minimality_audit_ring_num.py`, `minimality_audit_chain.py`, `minimality_audit_pf.py`.
Data: `debug/data/minimality_audit_{fuchs,ring,ring_num_dps15,chain,pf}.json`.
Nothing outside `debug/` was modified. Nothing was committed.

## 1. The claims, and how each operator was obtained

| # | Claim | Owner | How L was obtained | L f checked on f? | Minimality |
|:--|:--|:--|:--|:--|:--|
| 1 | Fourth-order eq:pf annihilates L(D) = sqrt(c1) N(D) | `sec:reduction` | Integration by parts on the integrand's own first-order relation Q f' + Q' f/2 = 0, boundary zero at the branch point (symbolic) | Yes: residual 1e-24 (paper); re-run here, Section 3 | Irreducible: proved (Fourier-Laplace of a simple D-module) and exhibited (integer monodromy) |
| 2 | Legendre modulus operator kills K(rho), K(1-rho) | `sec:obstruction` | Classical | Yes, 1e-25 | Order 2, classical |
| 3 | In-module source identity eq:inmodule | `sec:obstruction` | Fit with held-out residual | Measured on S itself | Not an annihilator claim; outside both defects |
| 4 | F_n(t) = int K0^(n-1) K0(tx) dx satisfies the order-n Mellin-route operator, n = 2..5 | `sec:four_density` | **Mellin route** (defect (a) applies) | Yes: driver `linking_ceiling_B_ode.py` takes theta^j under the integral analytically and evaluates at t = 0.3, 0.7, 1.6, 30 digits; relative residual 1e-28 to 1e-32 for every n. Guard repeats n = 4 at t = 0.7 (1e-15) | Bounded search before (degree 40 driver, 12 guard); **now proved irreducible for n = 2..5, Section 2** |
| 5 | eq:four_density_op is the n = 4 operator in x = 1/(4t^2) (AESZ #58); 0 quadrics; Sp4/SL4 | `sec:four_density` | Symbolic identity with #4 | via #4 | Proved irreducible (Section 2), so the quadric count is a statement about F_4's own D-module |
| 6 | n = 5: one quadric, orthogonal rank five | `sec:four_density` | Mellin route | via #4 | Proved irreducible (Section 2); "new relative to n = 4" stays open as stated |
| 7 | F_3's operator is the symmetric square (3F2 at t^2/4) | `sec:four_density`, `sec:density_ring` | Mellin route | via #4 | Proved irreducible (Section 2) |
| 8 | Triangle: L_G G = -3t^2, (theta-2) L_G G = 0 | `sec:density_ring` | Mellin route, source found after | Yes: 6 t at 20 digits (B2) and analytic theta on the momentum form (guard, 2 t) | Reducible by design; minimal content checked in Section 4 |
| 9 | Triangle pieces: algebraic, elliptic (Legendre via exact map), elementary; "sits on the elliptic rung" | `sec:density_ring` | Exact division + exact map | The elliptic piece's presence was **not** checked; run here (Section 4) | M irreducible (equivalent to Legendre). Presence of M in G's minimal operator: **checked here, holds** |
| 10 | Chain control: W = R3 R1 (theta-2) annihilates C, N2 spurious, only non-elementary piece is F_3's K3 family | `sec:density_ring` | Mellin route, then exact peeling | Yes: 4 t at 20 and 32 digits, residual <= 2e-31, perturbed control 1.9 to 604 | R3 equivalent to F_3's operator (irreducible). That C actually reaches R3 was **not** checked; run here (Section 5), holds |
| 11 | Tree identity eq:tree | `sec:four_density` | Fourier / delta functions | n/a (not an operator) | n/a |
| 12 | Three-electron chain is a K3; period = V3 bar = Bailey's K K' sum; P3e collapse | `sec:three_electron` | Geometry + value identities + PSLQ | n/a | n/a |
| 13 | Four centres add no genus | `sec:four_centre` | Momentum reduction, symbolic radicand | n/a | n/a |
| 14 | One-centre ring values (Remiddi), (1,2,3) relation | `sec:density_ring` | Value identities + PSLQ | n/a | n/a |

The 2026-10-13 records (four-centre genus, three-electron K3) contain no operator-route claims:
their rungs come from the branch locus and from value identities (Bailey's formula, PSLQ). The
2026-10-16 record's operators are #4 to #7.

## 2. Irreducibility of the F_n operators without a degree bound

Method (`minimality_audit_fuchs.py`). For a Fuchsian operator L of order n with singular set S,
a right factor of order k has, at each p in S, exponents forming a sub-multiset of L's exponents
(the valuation multiset of a subspace of local solutions), and each extra (apparent) singularity
adds at least 1 to Fuchs' relation. Hence a factor of order k needs
sum over p in S of (smallest sum of k exponents at p) <= (|S| - 2) k(k-1)/2.
If this fails for every k, L is irreducible, with no bound on degrees or apparent singularities.
The script computes the indicial polynomials exactly (sympy), checks Fuchs' relation for L itself,
and, as a positive control, must flag the known reducible L_G as "factor possible".

| Operator | Singular points | Exponents | Fuchs sum | Smallest sum vs bound, k = 1, 2, 3, 4 | Verdict |
|:--|:--|:--|:--|:--|:--|
| F_2 (Legendre) | 0, +-1, oo | 0,0 / 0,0 / 1,1 | 2 = 2 | 1 > 0 | irreducible |
| F_3 | 0, +-2, oo | 0,0,0 / 0,1/2,1 / 1,1,1 | 6 = 6 | 1 > 0; 3 > 2 | irreducible |
| F_4 (AESZ #58, in t) | 0, +-1, +-3, oo | 0,0,2,2 / 0,1,1,2 (each conifold) / 1,1,1,1 | 24 = 24 | 1 > 0; 6 > 4; 13 > 12 | irreducible |
| F_5 | 0, +-2, +-4, oo | 0,0,1,2,2 / 0,1,3/2,2,3 / 1,1,1,1,1 | 40 = 40 | 1 > 0; 6 > 4; 14 > 12; 25 > 24 | irreducible |
| L_G (control) | 0, +-1, +-2, oo | 0,3,3 / 0,0,1 / -3/2,0,1 / 0,2,3 | 12 = 12 | -3 <= 0; 2 <= 4 | factor possible (correct: it factors) |

Working in t, not x = 1/(4t^2), matters: a degree-two pullback could in principle split an
irreducible operator; the t-operator is the one F_n satisfies, and it is irreducible. Each F_n is
nonzero and satisfies its operator on the actual function (#4), so each operator is F_n's minimal
operator. Defect (b) cannot apply to the F_n rungs, and the quadric counts 0, 1, 0, 1 are counts
on F_n's own D-module. This upgrades "irreducible up to degree 40" (paper, record) and "degree 12"
(guard) to an unconditional statement. Caveat: the sub-multiset step is the standard
valuation-filtration argument, used here without a written proof in the corpus.

## 3. eq:pf on the actual function

eq:pf is not from the Mellin route; it comes from the integrand's own differential identity, and
the boundary term vanishes because Q(1)^(1/2) = 0. Re-checked here on N(D) by direct quadrature
(x = cosh v, derivatives under the integral), relative residual over the sum of the five terms,
30 digits, rho = 1/5 and 1/3, D = 0.3, 1.1, 2.7: between 0 and 4e-32 at all six points. No source
term. Irreducible by the paper's Fourier-Laplace argument, so it is N's minimal operator.
SAFE.

## 4. The triangle: is the elliptic piece really there?

The paper's operator facts are exact: L_G G = -3t^2 (checked on the function), L_G = M (theta - r)
with r = 3t^2/(4 - t^2), and M equivalent to Legendre. But "the triangle sits on the elliptic
rung" needs M to be a composition factor of G's minimal operator. Otherwise G would be
Liouvillian (algebraic and elementary only) and would sit lower still. That check was owed.

Argument. Put h = (theta - r) G. Then M h = -3t^2. G is Liouvillian iff h is, iff
O' = (theta - 2) M has a hyperexponential solution y with M y != 0 (its Liouvillian solutions form
a Galois-invariant subspace that meets ker M trivially, since M has Galois group SL2).

Exact search (`minimality_audit_ring.py`; all 21 exponent combinations with integer polynomial
degree, so no bound). Controls: the same search finds (1 - t^2/4)^(-3/2) for L_G, and nothing for
Legendre.

Result: O' **does** have one: y = t^2/(t^2 - 4), with M y = t^2. So M h = -3t^2 is solved by
h_p = 3t^2/(4 - t^2) = r, and that is (theta - r) applied to the constant -1. Exactly:
L_G(1) = Q_0(0) + t^2 Q_1(0) + t^4 Q_2(0) = 3t^2. **The source term is L_G of a constant:
L_G (G + 1) = 0.** The constant is a pure power t^0, the very kind of term the regularized Mellin
transform cannot see. So the decisive question became numerical: is h = h_p?

Numerical test on the momentum form (`minimality_audit_ring_num.py`, 15 digits; theta taken
analytically under the integral):

| t | G | h = (theta - r) G | h - h_p |
|:--|:--|:--|:--|
| 0.5 | 0.0353747404113085 | 0.0394937898604453 | -0.1605 |
| 1.3 | 0.101228313824996 | -0.134080661946601 | -2.329 |
| 2.5 | 0.161158263089377 | 1.43288715789246 | 9.766 |
| 3.5 | 0.189822786637137 | 0.925035808577962 | 5.380 |

h - h_p is O(1), not small, at all four t. It lies in ker M and is nonzero, so its minimal operator
is M, and G's minimal operator contains M. **The elliptic piece is present; the claim is SAFE.**
(15 digits is enough: the question is zero versus O(1).)

Sharpening for the record (not a defect): the "elementary" piece of G is exactly the constant -1.
G = -1 + alpha (1 - t^2/4)^(-3/2) + (image of a Legendre solution under the exact map). The
algebraic coefficient alpha is not determined here; the ratio (G + 1)/(1 - t^2/4)^(-3/2) varies
(0.940 at t = 0.5, 0.483 at t = 1.3), consistent with a nonzero elliptic part. The paper's
"(theta - 2) L_G annihilates G" is true; "L_G (G + 1) = 0" is the sharper statement, and says the
missing source in defect (a) was a constant.

## 5. The chain control: does C actually reach the K3 piece?

W = R3 R1 (theta - 2) kills C. R3 is equivalent to F_3's operator by the exact map, and that
operator is irreducible (Section 2), so R3 is irreducible. C's minimal operator contains R3 iff
v = R1 (theta - 2) C is nonzero. Evaluated with step B2's evaluator, 20 digits
(`minimality_audit_chain.py`):

| t | C | (theta - 2) C | v = R1 (theta - 2) C |
|:--|:--|:--|:--|
| 0.5 | 0.0762230487001748 | -0.0612548604169068 | 0.118274978628744 |
| 1.3 | 0.191908404097144 | -0.243802859670157 | 1.33134339062821 |
| 2.5 | 0.279759381684835 | -0.438261431342501 | 5.63722529935215 |
| 3.5 | 0.31687806634157 | -0.535138758045143 | 11.1969061569634 |

v != 0, so C reaches R3: the chain's only non-elementary piece is F_3's K3 family, as the paper
says. The Mellin route's L_C did annihilate C (residual 1e-19 at 20 digits, B2), so for the
chain defect (a) did not occur; defect (b) did (N2) and is already corrected in the paper. SAFE.

## 6. Verdict table

| Claim | Owner | Seam | Verdict | Evidence |
|:--|:--|:--|:--|:--|
| eq:pf annihilates N(D); irreducible | sec:reduction/obstruction | (P59 core) | SAFE | own-identity derivation; residual <= 4e-32 at 6 points (30 digits); FL simplicity proof |
| Modulus operator, in-module source | sec:obstruction | - | SAFE (outside scope) | classical / measured fit |
| F_n operators, n = 2..5, annihilate F_n | sec:four_density | J10 | SAFE | 3 t, 30 digits, 1e-28..1e-32; irreducible by Fuchs (new, unbounded) |
| AESZ #58, 0 quadrics, not elliptic-built | sec:four_density | J10 (F) | SAFE | as above; the degree-40 bound is no longer needed |
| n = 5 orthogonal rank five | sec:four_density | J10 (evidence) | SAFE for (a),(b) | irreducible by Fuchs; relation to n = 4 open as stated |
| Tree identity | sec:four_density | J10 (F) | SAFE (not operator-based) | delta-function argument + 3-D check |
| Three-electron K3, V3 bar = K K' | sec:three_electron | J7 (F) | SAFE (not operator-based) | geometry, Bailey identity, guarded |
| Four centres, genus 1 | sec:four_centre | J4 (F) | SAFE (not operator-based) | symbolic radicand, MD values |
| Rung rule (trees ordered; ring elliptic) | sec:f12, sec:beyond | J9 (Obs) | SAFE | rests on J10/J11, both SAFE |
| Ring: L_G G = -3t^2, factorization, elliptic rung | sec:density_ring | J11 (F) | SAFE | elliptic piece present (Section 4); source = L_G(-1) |
| Ring: no K3 piece | sec:density_ring | J11 (F) | SAFE | composition factors of orders 1, 2, 1 (Jordan-Holder) |
| Chain control: W kills C; only K3 non-elementary | sec:density_ring | J11 (F) | SAFE | v != 0 at 4 t (Section 5); R3 irreducible |

## 7. Owed, for the PI (none blocking)

1. The guards could carry the unbounded Fuchs argument for F_4 and F_5 in place of (or beside)
   the degree-12 search in `tests/test_paper59_four_density_rung.py`. Writing that guard is its
   own step, per CLAUDE.md §9.
2. Neither the ring guard nor the paper checks that G carries the elliptic piece, or that the
   chain reaches R3. Both now hold numerically here; a guard (h - h_p != 0; v != 0) is cheap.
3. Paper 59 could state L_G(G + 1) = 0 and name the elementary piece as the constant -1. The
   present wording is not wrong.
