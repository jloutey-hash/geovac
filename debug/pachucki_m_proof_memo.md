# Proof of the m != 0 r12-recursion (**) from Pachucki's w1-equation

Date 2026-11-20. Proof memo; no paper, test, CLAUDE.md, CHANGELOG or git state touched. Closes the item
"(**) itself is derived, not proved line by line" of `debug/pachucki_m_nonzero_memo.md` (Not established, l.174).
Source: K. Pachucki, Phys. Rev. A 80, 032520 (2009), arXiv:0905.1635v3. Equation numbers below are those of the
printed PDF; the earlier memo's "label 39" is PDF Eq. (42) (= Eq. (16)), "label 13" is Eq. (13), "label 61/63"
are Eqs. (64)/(66) of Sec. VII.

## Verdict: GO

Every step from Pachucki's Eq. (16) to (**) is proved below, with stated hypotheses. The one input taken without
proof is Eq. (16) itself for positive real parameters (his integration-by-parts theorem). It is not the step that
was open; it is the published starting point, and the earlier memo's transcription check (label 13) plus the exact
Laplacian controls (which do not use (16)) corroborate it.

**What is Pachucki's:** the function g (Eq. 2), its Laplace relation to f (Eq. 3), the w1-equation (16) with sigma
(12) and P (13)-(14), and the procedure "expand in w1, take coefficients, inverse Laplace transform" (Secs. IV, VII;
Eq. 64 is his exchange-case instance, obtained after dividing by t^2). **What is added here:** (i) the region and
the bounds that make g holomorphic in w1 and its Taylor coefficients the transforms of the r12-weighted masters;
(ii) the coefficient identity without division by t^2 (hence valid off the exchange locus); (iii) a boundary-free
inverse transform (integrate four times instead of differentiating), which never evaluates f_n^{(j)}(0+), so
logarithms at r -> 0 are harmless; (iv) the explicit structure of the source.

## 1. Setting and hypotheses

Exponents u2, u3, w2, w3 > 0 (u3: r1A, u2: r1B, w2: r2A, w3: r2B; t = u1 conjugate to r = R_AB; w1: r12).
Let m0 = min(u2, u3, w2, w3) and fix 0 < delta < m0/2.

Region. (**) is proved on the open set Omega+ = {u2, u3, w2, w3 > 0}; the elementary form of S_m (Sec. 5) on
Omega = Omega+ minus {u2 = u3} minus {w2 = w3}, which is exactly {q != 0}. The exchange point
u3 = w2 = a, u2 = w3 = b, 0 < a != b lies in Omega (q = -(a^2-b^2)^2, and c = 0 there).

For |w1| < delta (complex) and r >= 0 define

    f(r; w1) = r Int d^3r1/(4pi) d^3r2/(4pi) e^{-w1 r12 - u3 r1A - u2 r1B - w2 r2A - w3 r2B} / (r12 r1A r1B r2A r2B),

f_n(r) = (-1)^n d^n f/dw1^n |_{w1=0} (Pachucki Eq. 5), g(t; w1) = Int_0^inf e^{-tr} f(r; w1) dr, and
g_n = (-1)^n d^n g/dw1^n |_0, so that g = sum_n (-w1)^n g_n / n!.

## 2. Lemma 1 (bounds, smoothness, holomorphy)

**(a) Uniform linear bound.** |e^{-w1 r12}| <= e^{delta r12} <= e^{(delta/2)(r1A + r2A + r1B + r2B)} (two triangle
inequalities, averaged). Hence |f(r; w1)| <= (r/16pi^2) Int Int Phi(r1) Psi(r2) / r12 with
Phi(x) = e^{-a'|x-A| - b'|x-B|}/(|x-A||x-B|), a', b' >= m0 - delta/2 > 0, and Psi alike. By the
Hardy-Littlewood-Sobolev inequality (Lieb-Loss, Analysis, Thm 4.3; n = 3, lambda = 1, p = q = 6/5),
Int Int Phi Psi / r12 <= C ||Phi||_{6/5} ||Psi||_{6/5}. By Cauchy-Schwarz,
||Phi||_{6/5}^{6/5} <= ||h_{a'}||_2 ||h_{b'}||_2 with h_s(x) = e^{-(6/5)s|x|} |x|^{-6/5}, and
||h_s||_2^2 = 4pi Gamma(3/5) (12s/5)^{-3/5} < inf, independent of A and B. So

    |f(r; w1)| <= K r     for all r >= 0, |w1| <= delta.                                        (B1)

**(b) Prolate form; smoothness and decay.** In prolate spheroidal coordinates of each electron about A, B
(r1A = r(xi1+eta1)/2, r1B = r(xi1-eta1)/2, d^3r1/(r1A r1B) = (r/2) dxi1 deta1 dphi1, and likewise for electron 2),

    f(r; w1) = (r^2/64pi^2) Int dmu e^{-r (Q + w1 rho)} / rho,

with dmu = prod dxi deta dphi over xi >= 1, |eta| <= 1, Q = [(u3+u2)xi1 + (u3-u2)eta1 + (w2+w3)xi2 + (w2-w3)eta2]/2
and rho = r12/r, both independent of r. Bounds: Q >= m0 (xi1 + xi2), Q <= M (xi1 + xi2) (M = max exponent),
rho <= (r1A + r2A)/r <= xi1 + xi2, so Re(Q + w1 rho) >= (m0 - delta)(xi1 + xi2) >= 2(m0 - delta). Moreover
K_eps := Int dmu e^{-eps(xi1+xi2)}/rho < inf for every eps > 0: at r = 1 it is 4 times an integral of the type
bounded in (a). Differentiating under the integral sign j times in r and k times in w1 brings down
polynomials bounded by C(xi1+xi2)^{j+k}; for r >= r0 > 0 the integrand is dominated by
C (xi1+xi2)^{j+k} e^{-r0(m0-delta)(xi1+xi2)}/rho, integrable. Hence f is C^inf in r on (0, inf), holomorphic in w1 on
|w1| < delta, all mixed derivatives are obtained under the integral, and

    |d^j_r d^k_{w1} f(r; w1)| <= C_{j,k,r0} (1 + r)^{2+k} e^{-(m0 - delta) r}     (r >= r0).          (B2)

(Split e^{-r(m0-delta)(xi1+xi2)} = e^{-r(m0-delta)(xi1+xi2)/2} e^{-r(m0-delta)(xi1+xi2)/2} <= e^{-(m0-delta) r}
e^{-r0(m0-delta)(xi1+xi2)/2}.) The integrand is also entire in (u2, u3, w2, w3), with the same domination for complex
exponents whose real parts stay in a compact subset of {Re > 2 delta}; so f is jointly holomorphic in all five
exponents there (Osgood), jointly C^inf with r > 0. In particular every f_n is C^inf on (0, inf), every exponent
derivative of f_n (the "jets" d^alpha f_n) is the corresponding r1A^.. r1B^.. r2A^.. r2B^.. -weighted integral,
and r- and exponent-derivatives commute.

**(c) Coefficients.** Cauchy's estimate on |w1| = delta with (B1): |f_n(r)| <= n! K r / delta^n for all r >= 0.
So each f_n is continuous on [0, inf), O(r) at both ends of the bound, and exponentially decaying by (B2).

## 3. Lemma 2 (g, the transforms, and Eq. (16) on the disc)

(a) For real t > 0, g(t; w1) = Int_0^inf e^{-tr} f(r; w1) dr converges absolutely (B1), is holomorphic in w1 on
|w1| < delta (dominant K r e^{-tr}), and g_n(t) = Int_0^inf e^{-tr} f_n(r) dr (differentiation under the integral,
dominant from (c)). This is item (1) of the brief: **g_n = L[f_n]**, with differentiation in w1 commuting with both
integrals.

(b) For real w1 > 0, g equals Pachucki's g of Eq. (2) with u1 = t: Tonelli on the positive integrand, then
translation and rotation invariance reduce the integral over the position of B relative to A to
Int d^3R/(4pi) e^{-t|R|}/|R| (f/|R|) = Int_0^inf e^{-tr} f dr, exactly his Eq. (3) computation with w1 kept.

(c) **Input (Pachucki).** For all positive real w1, t, u2, u3, w2, w3,
sigma dg/dw1 + (1/2)(dsigma/dw1) g + P(w1, t; w2, u2; w3, u3) = 0   [Eq. (16), sigma from (12), P from (13)-(14)].
(Both sides are continuous in the parameters, so validity at generic positive values extends to all.)

(d) Extension to the disc. Fix real t > 0. sigma is a polynomial. Each Gamma(alpha1, alpha2, alpha3) in (13) is a
divided difference of the logarithm, [ln(alpha1+alpha2) - ln(alpha1+alpha3)]/(alpha2-alpha3), times
1/(alpha2+alpha3); it is holomorphic wherever Re(alpha1+alpha2), Re(alpha1+alpha3) > 0 and alpha2 + alpha3 != 0. For
the eight Gamma's of (13) with u1 = t, every alpha1+alpha2, alpha1+alpha3 and alpha2+alpha3 is a sum of at least two
of the positive numbers t, u2, u3, w2, w3, plus at most w1 (e.g. Gamma(u2+w1, u3, t+w2): u2+u3+w1, u2+w1+t+w2,
u3+t+w2). With |w1| < delta < m0 all have positive real part. So the left side of (16) is holomorphic on
|w1| < delta and vanishes on (0, delta); by the identity theorem **(16) holds on the whole disc**, for every t > 0.

## 4. Lemma 3 (the coefficient identity, no division)

By the algebra of (12) with u1 = t (checked symbolically, Sec. 9),
sigma = sigma0 + sigma2 w1^2 + t^2 w1^4, sigma0 = c - q t^2, sigma2 = t^4 - Sig t^2 + e (constants as in the
earlier memo, l.26-28). Insert g = sum (-w1)^n g_n/n! (convergent on the disc) into (16); a convergent power series
vanishes identically iff all coefficients do. Multiplying the w1^m coefficient by (-1)^{m+1} m!:

    sigma0 g_{m+1} + m^2 sigma2 g_{m-1} + m(m-1)^2(m-2) t^2 g_{m-3} + (-1)^{m+1} P_m(t) = 0,   t > 0,     (T_m)

P_m = d^m P/dw1^m |_{w1=0}, g_k := 0 for k < 0 (the factors m^2 = m(m-1)+m and m(m-1)^2(m-2) =
m(m-1)(m-2)(m-3) + 2m(m-1)(m-2) collect the sigma2 and t^2 w1^4 contributions). Unlike Pachucki's Eq. (64) route,
nothing is divided by t^2, so (T_m) holds for all exponents in Omega+, not only where c = 0.

## 5. Lemma 4 (structure of the source)

On Omega the Gamma's split termwise (all checked symbolically, Sec. 9):
* In the first two terms of (13) the t-dependent denominator cancels: t w1[(t+w2)^2 - u3^2] Gamma(u2+w1, u3, t+w2)
  = -t w1 [ln(u2+u3+w1) - ln(t+u2+w2+w1)], and likewise for the second.
* The other six have t-free denominators, (w3 - w1 - w2)(w3 + w1 + w2), (w2 - w1 - w3)(...), (u3 - u2 - w1)(...),
  (u2 - u3 - w1)(...), (w1 - w2 - w3)(...), (w1 - u2 - u3)(...), nonzero near w1 = 0 on Omega (this is where
  u2 != u3, w2 != w3 is used).

Hence near w1 = 0, P = sum_{k=1}^4 A_k(t, w1) ln(t + c_k + eps_k w1) + B(t, w1), with
(c_k, eps_k) in {(u2+w3, 0), (u3+w2, 0), (u2+w2, 1), (u3+w3, 1)}, and A_k, B polynomials in t of degree <= 2 whose
coefficients are holomorphic in w1 (B contains the t-free logarithms). Leibniz then gives

    P_m(t) = sum_k alpha_k(t) ln(t + c_k) + sum_k sum_{i=1}^m beta_{k,i}(t) / (t + c_k)^i + pi_m(t),

with alpha_k, beta_{k,i}, pi_m polynomials in t of degree <= 2 (eps_k^i (-1)^{i-1}(i-1)! absorbed into beta).

## 6. Lemma 5 (a Laplace dictionary with no boundary values)

Let I phi(r) = Int_0^r phi. For phi in L^1_loc[0, inf) with polynomial growth: L[I^j phi] = t^{-j} L[phi] (t > 0;
Fubini), and I^j phi is continuous on [0, inf) for j >= 1. If phi is in C^inf(0, inf), then
d^{j+s}/dr^{j+s} I^j phi = phi^{(s)} on (0, inf) (fundamental theorem of calculus at each r > 0, where phi is
continuous; then classical derivatives). Three pairs, all for t > 0 and c > 0:

* t^{-l} = L[r^{l-1}/(l-1)!]  (l >= 1);
* (t+c)^{-i} = L[r^{i-1} e^{-cr}/(i-1)!];
* ln(t+c)/t = L[phi_c], phi_c(r) = ln c + E_1(cr) (E_1(x) = -Ei(-x)); phi_c is log-singular at 0, hence in
  L^1_loc, bounded at infinity, smooth on (0, inf), with phi_c'(r) = -e^{-cr}/r.

**Lerch** (Widder, The Laplace Transform, Ch. II Thm 6.3): if phi is continuous on [0, inf), Int e^{-t0 r}|phi| <
inf, and L[phi](t) = 0 for all real t > t0, then phi = 0 on [0, inf).

## 7. Theorem: (**) holds pointwise for r > 0

Divide (T_m) by t^4. Term by term, every summand is now the transform of a function continuous on [0, inf) and of
polynomial growth:

* (c - q t^2) t^{-4} g_{m+1} = L[c I^4 f_{m+1} - q I^2 f_{m+1}];
* (t^4 - Sig t^2 + e) t^{-4} g_{m-1} = L[f_{m-1} - Sig I^2 f_{m-1} + e I^4 f_{m-1}];
* t^2 t^{-4} g_{m-3} = L[I^2 f_{m-3}];
* P_m t^{-4} = L[H_m] with H_m = sum_k sum_{j<=2} alpha_{k,j} I^{3-j} phi_{c_k} + (inverse of the proper rational
  part, via partial fractions in t^{-l} and (t+c_k)^{-l}: a polynomial of degree <= 3 plus polynomials times
  e^{-c_k r}) + (inverse of pi_m t^{-4}: a polynomial of degree <= 3). Here alpha_k(t) = sum_j alpha_{k,j} t^j and
  3 - j >= 1, so H_m is continuous.

(T_m)/t^4 states that the transform of

    E_m := c I^4 f_{m+1} - q I^2 f_{m+1} + m^2 [f_{m-1} - Sig I^2 f_{m-1} + e I^4 f_{m-1}]
           + m(m-1)^2(m-2) I^2 f_{m-3} + (-1)^{m+1} H_m

vanishes for every t > 0. E_m is continuous on [0, inf) with polynomial growth, so by Lerch **E_m = 0 on
[0, inf)**. On (0, inf) every term is C^inf (Lemma 1, Lemma 5); apply d^4/dr^4:

    c f_{m+1} - q f_{m+1}'' + m^2 (f_{m-1}'''' - Sig f_{m-1}'' + e f_{m-1}) + m(m-1)^2(m-2) f_{m-3}''
        + (-1)^{m+1} H_m'''' = 0,

which is (**) with **S_m := (-1)^{m+1} H_m''''**. Explicitly, by Lemma 5, on r > 0:

    S_m(r) = (-1)^{m+1} [ sum_k sum_{j<=2} alpha_{k,j} (d/dr)^j (-e^{-c_k r}/r)
                          + sum_k sum_i L^{-1}[ proper part of beta_{k,i}(t)/(t+c_k)^i ] ],

the second bracket a sum of r^l e^{-c_k r}. This is exactly the formal termwise inverse of (-1)^{m+1} P_m with
t^j ln(t+c) -> (d/dr)^j(-e^{-cr}/r) and every polynomial in t dropped (the polynomial parts of H_m have degree
<= 3 and die under d^4): it is the earlier memo's S_m, elementary, with no Ei. Item (3) and item (4) of the brief
are this paragraph and the Lerch step.

**Exponent jets.** (**) is an identity of functions jointly C^inf in (r, u2, u3, w2, w3) on (0, inf) x Omega
(Lemma 1(b) for the f's; S_m is built from finitely many elementary functions whose coefficients are rational in
the exponents with denominators from Sec. 5, nonzero on Omega). It may therefore be differentiated in the exponents
at the exchange point, r- and exponent-derivatives commuting: the Leibniz step of the earlier memo is justified, and
d^alpha f_n is the r1A/r1B/r2A/r2B-weighted master (Pachucki Eq. 6 by the same differentiation under the integral).

## 8. Behaviour at r -> 0, and where the "boundary terms" went

The proof uses only (B1) at r -> 0: f_n = O(r), continuous at 0. It never evaluates f_n^{(j)}(0+). This matters:
at the exchange point Pachucki's master is f(r) = Int_r^inf F^{(-1)}(r')/r' dr' (Eq. 62), so f' = -F^{(-1)}(r)/r,
and with Ei(-x) = gamma + ln x + O(x) in F (Eqs. 63, 66) the higher derivatives of the f_n generally carry
r^k ln r terms; f_n'' or f_n'''' may diverge at 0 (Paper 58 records the exchange integral's logarithm at R^11 ln R).
In the classical formula t^j L[phi] = L[phi^{(j)}] + sum t^{j-1-i} phi^{(i)}(0+) such divergences would make the
"boundary terms" undefined. In the integrated form they never appear: the polynomials in t of the identity (the
r = 0-supported distributions delta^{(i)}) become polynomials in r of degree <= 3 inside E_m and are removed by
d^4/dr^4 on the open half-line. The same holds for the source: ln(t+c) is not the transform of a function
(-e^{-cr}/r is not integrable at 0), but ln(t+c)/t^{1+s} is, and the formal r > 0 inverse is recovered exactly.

## 9. Corollary (decay; the earlier memo's uniqueness step)

By (B2) every f_n and every jet decays like a polynomial times e^{-(m0-delta) r}. At the exchange point c = 0, so the homogeneous
equation of (**) is q y'' = 0 with solutions A + B r; a decaying solution of (**) is therefore unique, which is the
selection the earlier memo's solver makes. (Proved here; the solver's own exactness is not re-examined.)

## 10. Checks run for this memo

`pm_proof_checks.py` (scratchpad, not kept): sigma of Eq. (12) at u1 = t minus sigma0 + sigma2 w1^2 + t^2 w1^4 = 0
(sympy, generic exponents); the w1^m coefficient of (16) for m = 0..6 equals (T_m) with symbolic g_n; in P of (13)
with each logarithm replaced by a symbol, the t-dependent log arguments are exactly the four t + c_k (+ w1), every
log coefficient has a t-free denominator and numerator of degree 2 in t, and there is no log-free remainder;
L[ln c + E_1(cr)] = ln(t+c)/t (mpmath quadrature, residuals 0 and 1e-16, the latter from float inputs).

## 11. What this does not cover

* Eq. (16) itself (Pachucki's IBP derivation, the momentum representation (7)-(8), the vanishing of surface terms).
  It is used as published. Its independent corroboration is outside the proof: the Laplacian controls (a) of the
  earlier memo (exact, built from nabla^2 r12 = 2/r12, not from (16)), and the 47-48 digit agreement with
  exchange_hp (Neumann expansion).
* The correctness of the code that implements S_m and solves (**) (`debug/pachucki_m_core.py`): the proof fixes what
  S_m is; that the code computes it is certified by the earlier memo's checks, not here.
* The elementary form of S_m on q = 0 (u2 = u3 or w2 = w3). (**) holds there too (Sec. 7 does not use Sec. 5),
  but S_m is then a limit of the termwise form; not needed for exchange at a != b.
* a = b (q = 0) and anything after (**) (span, operators, "exactly 9").

## 12. Weakest steps

1. **Lemma 4's termwise split** is read off (13) by inspection and confirmed by computer algebra on a
   transcription of (13); a mistranscribed (13) would change S_m but not the argument (any P of this Gamma shape
   gives an elementary S_m). The earlier memo's label-13 check and control (a) bear on the transcription.
2. **The input (16)**, as above: published, corroborated, not reproved.

