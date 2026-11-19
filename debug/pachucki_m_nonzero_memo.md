# Pachucki exchange bound for m != 0: r12-weighted masters without the 1/p^2

Date 2026-11-18. Derivation memo; no paper, test, record, CLAUDE.md, CHANGELOG or git state touched. Nothing
SEALED was opened. Continues `debug/pachucki_extension_memo.md` (section "m != 0 ... plan and obstacle").
Source: K. Pachucki, Phys. Rev. A 80, 032520 (2009), arXiv:0905.1635 (TeX source re-fetched; labels below are the
arXiv labels). The 2012 follow-up was not needed and not used.

## Verdicts

| case | verdict | one line |
|---|---|---|
| 2p_{+-1,A} 1s_B (2p_pi) | **GO** | controls (a), (b), (c) hold; derived operator = committed (9, 25) over Q |
| 3d_{+-2,A} 1s_B (3d_delta) | **GO** | controls (a), (b), (c) hold (f_4 rung included); derived operator = committed (9, 35) over Q |
| real 2p_x 2p_x, sigma = 2 part | **follows from 3d_delta** (same verdict) | the density is the 3d_delta density; = (3/2)(b/a)^2 x 3d_delta; hp check 48 digits at (1,2) |

The stated expectation (nine-function span; operators (9, 25) and (9, 35) equal to the committed ones) is met in
both cases. Nothing was tuned: the first complete run of each case gave the result reported.

## Route that resolved the 0/0 (and why)

Not the memo's plan (label 46 along a transverse line). The 1/p^2 of labels 45/46 is not in the r12 equation
itself; it is produced by the elimination step (label 44) that removes the r-derivatives of f_{n+1}. I kept the
r-derivatives instead. Label 39 (eq_2, first order in the r12 exponent w1) with u1 = t:

    sigma dg/dw1 + (1/2)(dsigma/dw1) g + P(w1, t; w2, u2; w3, u3) = 0,
    sigma = sigma0 + sigma2 w1^2 + t^2 w1^4,       sigma0 = c - q t^2,     sigma2 = t^4 - Sig t^2 + e,
    c = (u2^2-u3^2+w2^2-w3^2)(u2^2 w2^2-u3^2 w3^2) (= p^2 q),   q = (u2^2-u3^2)(w2^2-w3^2),
    Sig = u2^2+u3^2+w2^2+w3^2,   e = (u2^2-w3^2)(u3^2-w2^2).

With g = sum_n (-w1)^n g_n/n! (label 40), the w1^m coefficient is
sigma0 g_{m+1} + m^2 sigma2 g_{m-1} + m(m-1)^2(m-2) t^2 g_{m-3} + (-1)^{m+1} d^m P/dw1^m|_0 = 0, and the inverse
Laplace transform at r > 0 (t -> d/dr) gives, **for all exponents**,

    q f_{m+1}'' = c f_{m+1} + m^2 (d^4 - Sig d^2 + e) f_{m-1} + m (m-1)^2 (m-2) f_{m-3}'' + S_m(r)        (**)

i.e. f_{m+1}'' - p^2 f_{m+1} = (known)/q. The only division is by q, and q = -(a^2-b^2)^2 != 0 at the exchange
point. On the exchange locus (**) is label 61 differentiated twice (label 61 divides by t^2 first, which is only
possible where c = 0; that is why it "holds only on the locus"). Why this route: no 0/0 at any stage, no extra jet
order (the transverse-line plan needed f_alpha to order |alpha| + 2 + 2k, ~3000 multi-indices), and the source is
simpler than label 46's: P's t-dependent logarithms are ln(t + c_k + eps_k w1) with
c_k in {u2+w3, u3+w2, u2+w2, u3+w3} and coefficients polynomial in t (degree <= 2), so
S_m = (-1)^{m+1} L^{-1}[d^m P/dw1^m] is **elementary**: r-derivatives of e^{-c_k r}/r and r^{j} e^{-c_k r}. No Ei and
no F_i^+ enter the source; every Ei comes from the lower f's and from integrating twice.

Exponent derivatives: Leibniz on (**) at the exchange point, lowest multi-index first. Each step is u'' = known,
solved exactly in the decaying span {r^k f0, r^k e^{lam r} Ei(mu r), r^k e^{lam r} LC^j} by a triangular solver
(`MCalc.solve_D2`, modelled on the reviewed `solve_L0`; the homogeneous solutions A + B r do not decay, so the
solution is unique); the residual is re-checked as an exact identity after every solve. Two rungs are needed:
f_2 (m = 1) for 2p_pi, f_2 and f_4 (m = 3) for 3d_delta.

Angular reduction (z = (r_A^2 + R^2 - r_B^2)/(2R), X = x1 x2 + y1 y2, rho^2 = r_A^2 - z^2):
Re (x1+iy1)(x2-iy2) = X = (r1A^2 + r2A^2 - r12^2)/2 - z1 z2; Re (x1+iy1)^2 (x2-iy2)^2 = 2X^2 - rho1^2 rho2^2.
A monomial r12^{2k} r1A^n1 r1B^n2 r2A^n3 r2B^n4 / r12 (n even) maps to (4 pi)^2/R d^alpha f_{2k}, alpha = n + 1.
2p_pi: f_2 at (1,1,1,1) only, f_0 on the 144-index set of 2p_sigma. 3d_delta: f_4 at (1,1,1,1), f_2 on 144
indices, f_0 on 576. Normalisations: 2ppi = -(a^5 b^3/(2 pi^2)) int int X psi psi/r12 (the sign is the
Condon-Shortley phase of 2p_{-1}; 2ppi = -(2p_x 1s | 2p_x 1s)); 3ddel = (a^7 b^3/(12 pi^2)) int int (2X^2 - rho1^2
rho2^2) psi psi / r12. Both conventions were fixed by the definitions, then confirmed (not fitted) by exchange_hp.

## Transcription checks (before building on anything)

| check | result |
|---|---|
| label 13 verbatim minus the engine's four-logarithm decomposition is a polynomial in t (third t-derivative, generic exponents and w1) | exactly 0 (P itself: 0.089) |
| S_1 at the exchange point = -(label 63's W_1^{(-2)})'' (two separately transcribed equations) | exact, Q at (1,2), (4/5,3/2) and over Q(t) |
| f_2 at the exchange point = label 61 (n = 1), 2(u^2+w^2) f - f'' + W_1^{(-2)} over (u^2-w^2)^2 | exact, same three fields |
| independent quadrature (prolate spheroidal, phi integral closed form / Gauss, 4-D tensor Gauss, N = 40; no Pachucki equation, no Neumann series), (1,2) R = 7/5 and (4/5,3/2) R = 23/10 | f_2: 6.8e-6, 6.0e-6; d^{1111} f_2: 3.4e-6, 3.3e-6; f_4: 2.7e-8, 2.2e-8; d^{1111} f_4: 1.1e-8, 1.0e-8 (relative; limited by the r12 = 0 cone point) |

## Control (a): no r12 weight, exact

nabla_1^2 r12 = 2/r12, so int int psi psi / r12 = (1/2) int int (nabla^2 psi)(1) psi(2) r12 with
nabla^2 psi = [a^2 + b^2 - 2a/r_A - 2b/r_B + ab(r_A/r_B + r_B/r_A - R^2/(r_A r_B))] psi. In jets:

    16 f_{1111} = 8 [(a^2+b^2) f2_{1111} + 2a f2_{0111} + 2b f2_{1011} + ab (f2_{2011} + f2_{0211}) - ab R^2 f2_{0011}]

The left side is the existing 1s closed form (existing engine); the right side uses six f_2 jets from the new route,
including transverse derivatives (alpha_1 != alpha_3) and alpha_i = 0, 2. **The difference is identically zero** as
an element of the function algebra: over Q at (1,2) and (4/5,3/2), and over Q(t) (all t). One rung up
(nabla^2 r12^3 = 12 r12): the same combination of f_4 jets equals 12 f2_{1111}, identically, over Q at both pairs
(and over Q(t)). This control cannot pass with a wrong source, a wrong sign of c, q, Sig, e, or a wrong m-dependence.

## 2p_pi: results

* Over Q at (1,2), 3 s: K in the nine-function span, all nine functions present, Laurent range R^-3..R^5; rank of
  K..K^(10) is 9; **exactly one operator of order <= 9, degree <= 25, none of degree <= 24, none of order <= 8 up to
  degree 37**; coefficient degrees [23, 24, 25, ...]. It equals the committed operator over Q (scalar -1), compared
  both with `ep_heldout_2ppi_1_2_EPS.json` (`operator_Q_crt`) and with `OP_P_1_2` read from
  `tests/test_paper58_exchange_l1.py`.
* Control (b), closed form (dps 60; only quadrature = the 1-D f0 integral) against `exchange_hp` (tau_max 24, dps 50):

| (a, b), R | closed form | abs diff | digits | hp tail |
|---|---|---|---|---|
| (1,2), 1.4 | -0.013353601891216017682930492464386025458452519549795 | 4.7e-51 | 48 | 1e-65 |
| (1,2), 3.015 | -0.0012735770745720617968198177812207045324603277653535 | 3.5e-51 | 47 | 4e-52 |
| (4/5,3/2), 2.3 | -0.0072708870648988060628723129702682786806276558388417 | 1.8e-51 | 48 | 4e-63 |

  (4/5, 3/2) was not used in building or debugging the derivation; it has a != 1, so it tests the power of a.
* Control (c), Q(t) (a = 1, b = t), 374 s: every residual (160 jets: 144 f_0, 16 f_2) identically 0 in Q(t); K in
  the nine-function span and the span closed under d/dR through K^(10); rank 9 at (t,R) = (3, 7/5), hence formal
  rank 9; **denominators of every coefficient of every intermediate jet factor into t, t-1, t+1 only, and every key
  of every jet decays for all t > 1** (`pachucki_m_symcheck_2ppi.json`; the review memo's requirement 1). K's
  coefficients: max pole (t-1)^7. Specialised at t = 2 the Q(t) closed form equals the Q run exactly; at t = 15/8
  with K(R;a,b) = a K(aR;1,b/a) it equals the direct (4/5,3/2) run at R = 2.3 to 57 digits.
* Pinned entry: the coefficient of f0 in the 2ppi value is
  2 a^7 b^5 [225 - 180(a^2+b^2)R^2 - (102a^4+336a^2b^2+102b^4)R^4 - (24a^6+176a^4b^2+176a^2b^4+24b^6)R^6
  + (a^8+16a^6b^2+36a^4b^4+16a^2b^6+b^8)R^8] / (R^3 (b^2-a^2)^6).
* Bonus, a -> b: fed to the existing `pachucki_ext_limit.py`, poles (t-1)^-1..-7 cancel, no outside integral appears,
  the limit lies in the four-function span and its operator equals the committed (1,1) operator (4, 9) over Q.

## 3d_delta: results

* Over Q at (1,2), 48 s (576 f_0, 144 f_2, 16 f_4 jets, all residuals identically 0): K in the nine-function span,
  all nine present, Laurent range R^-5..R^7; rank of K..K^(10) is 9; **exactly one operator of order <= 9, degree <= 35,
  none of degree <= 34, none of order <= 8 up to degree 47**; coefficient degrees [33, 34, 35, ...]. It equals the
  committed operator over Q (scalar -1): `eq_heldout_3ddel_1_2_EQS.json` and `OP_3DDEL_1_2` in
  `tests/test_paper58_exchange_d_pp.py`.
* Control (b) against `exchange_hp` (tau_max 24, dps 50):

| (a, b), R | closed form | abs diff | digits | hp tail |
|---|---|---|---|---|
| (1,2), 1.4 | 0.0013461743233320120177136736257341354157956824513622 | 2.2e-51 | 47 | 2e-63 |
| (1,2), 3.015 | 0.00016319617726070146318353340800017014763211710773854 | 4.6e-52 | 47 | 2e-50 |
| (4/5,3/2), 2.3 | 0.00085209572345054007851997145804439878674346695593029 | 2.9e-52 | 48 | 5e-61 |

* Control (c), Q(t): run in three resumable invocations of at most 3000 s each (`pachucki_m_symres.py`; 3115 s +
  3099 s + 285 s of solving, plus a 544 s source cache), 736 jets (576 f_0, 144 f_2, 16 f_4), every residual
  identically 0 in Q(t). **Denominators of every intermediate jet factor into t, t-1, t+1 only, and every key of every
  jet decays for all t > 1** (`pachucki_m_symcheck_3ddel.json`). K is in the nine-function span over Q(t), the span is
  closed under d/dR through K^(10), and the rank at (t,R) = (3, 7/5) is 9. K's coefficients have poles up to (t-1)^9.
  At t = 2 the result equals the Q run exactly. At t = 15/8, with scaling, it equals the direct (4/5,3/2) run at
  R = 2.3 to 56 digits.
* Pinned entry: the coefficient of f0 in the 3ddel value is
  a^9 b^5 [99225 - 28350(a^2+b^2)R^2 + (15795a^4+53460a^2b^2+15795b^4)R^4 + 7920(a^6+b^6)R^6
  + 55080(a^4b^2+a^2b^4)R^6 + ...] / (12 R^5 (b^2-a^2)^8). The full polynomial is in
  `pachucki_m_case_3ddel_sym_res.json` (`nine_coords[0]`).
* a -> b limit (existing `pachucki_ext_limit.py`): the poles (t-1)^-1 to (t-1)^-9 cancel, no outside integral
  appears, the limit lies in the four-function span, and its operator equals the committed (1,1) operator (4, 10)
  over Q (`eq_heldout_3ddel_1_1_EQS.json`).

## Real 2p_x 2p_x exchange, sigma = 2 part

x^2 = rho^2 (1 + cos 2phi)/2, so (2p_x 2p_x | 2p_x 2p_x) = pppi0 + (1/2) ppdel, where pppi0 is the m = 0 part
already covered by the extension memo and ppdel = (2p+1_A 2p-1_B | 2p-1_A 2p+1_B). The ppdel density is
(x+iy)^2 e^{-a r_A - b r_B}, the 3d_delta density, so ppdel = 3 (b/a)^2 x 3ddel (the committed Derivation 3), an
identity in a, b, R. Hence the sigma = 2 part is (3/2)(b/a)^2 x 3ddel and inherits the 3d_delta verdict, and with it
the full real 2p_x 2p_x exchange integral is in the nine-function span. Check: 12 x (3ddel closed form) at (1,2),
R = 1.4 against exchange_hp for ppdel (tau 24, dps 50): 12 x 3ddel = 0.016154091879984144212564083508809624989548189416347, exchange_hp ppdel 0.0161540918799841442125640835088096249895481894163, difference 4.7e-50 (48 digits). No separate operator was derived (same function up
to a constant). Not checked: the decomposition pppi0 + ppdel/2 against a real-orbital evaluator (it is the
elementary phi average; exchange_hp takes complex orbitals).

## Closed-form structure

All nine functions enter in both cases: f0, Ei(-(a+b)R), e^{-(a+b)R}, e^{-2aR}, e^{-2bR}, e^{2aR}Ei(-2(a+b)R),
e^{2bR}Ei(-2(a+b)R), psi_a, psi_b (psi as in the extension memo, LC tied). Laurent ranges: 2p_pi R^-3..R^5,
3d_delta R^-5..R^7 at (1,2). Denominators: t, t-1, t+1 (2p_pi; in (a,b): a, b-a, b+a). The r12-weighted masters add
no function: f2_{1111} and f4_{1111} are themselves in the nine-function span.

## Surprises

1. The obstacle was an artefact of labels 45/46, not of the problem: the second-order form (**) of the r12
   recursion is regular at p = 0, for every exponent, and its source is elementary.
2. Control (a) passed on the first run, exactly; the only bug met was a sign in my own check of label 63 (I compared
   +W where the convention is S = -W).
3. sympy 1.14's `PolyElement._gcd_ZZ` is heugcd with no fallback and raised `HeuristicGCDFailed` in the Q(t) runs
   (derivative chain of K; f_4 solve). `pachucki_m_core.py` patches it with an exact python-flint gcd fallback. The
   earlier extension runs did not hit it; any future Q(t) run on larger sets can.
4. sympy `FracElement` cannot be pickled (RuntimeError in `PolyRing.__getstate__`); Q(t) state is stored as strings.

## Not established

* "Exactly 9" for all a < b: as before, needs the nine functions independent over C(R); what is shown is the bound,
  formal coordinate rank 9, and at (1,2) a unique operator at the committed degree with none of order <= 8 up to the
  tested degree.
* Operator equality is shown at (1,2) only (the only committed unequal pair) and in the a -> b limit, for both cases.
* (**) itself is derived, not proved line by line here: the termwise inverse Laplace transform at r > 0 (dropping
  the terms supported at r = 0) and analyticity of g in w1 at 0 are taken from Pachucki's own procedure. It is
  certified by the exact agreement with label 61 on the locus, the exact Laplacian controls off the locus, the
  independent quadrature, and the 47-48 digit agreement with exchange_hp.
* General l, m: only the three cases above. The route is generic (any polynomial in r12^2), cost grows with the
  f_0 set.

## What a guard must pin (for a separate guard pass)

* Control (a) as an exact identity over Q(t) (rejects any wrong source/sign in (**): e.g. S_1 -> -S_1, c -> -c,
  Sig dropped, m^2 -> m).
* f_2(exchange point) = label 61 exactly (rejects a mistranscribed label 13 term).
* 2p_pi at (1,2): unique (9, 25), none at degree 24, equal to `OP_P_1_2`; value at R = 1.4
  -0.0133536018912160176829304924643860254584525195498 (rejects sigma/pi mix-up, the sign of the r12^2 term in X, a
  shifted alpha).
* (4/5, 3/2), R = 2.3: -0.00727088706489880606287231297026827868062765583884 (rejects a wrong power of a).
* The f0 coefficient above (rejects a wrong normalisation or scaling).
* Q(t): denominators of all intermediates in {t, t-1, t+1}, all keys decaying (rejects a hidden exceptional t).
* 3d_delta at (1,2): unique (9, 35), equal to `OP_3DDEL_1_2`; value at R = 1.4
  0.00134617432333201201771367362573413541579568245136 (rejects 2X^2 - rho^2 rho^2 -> X^2, a wrong f_4 coefficient
  m(m-1)^2(m-2) = 12, a^7/12 -> a^7/18).

## Files

Scripts (`debug/`): `pachucki_m_core.py` (route (**), source S_m, solver, gcd patch), `pachucki_m_stage1.py`
(transcription checks, controls (a)), `pachucki_m_cases.py` (cases, span, rank, operator, numerics; Q or Q(t)),
`pachucki_m_symres.py` (resumable Q(t) solve, string-pickled), `pachucki_m_symres_sources.py` (source cache), `pachucki_m_symcheck.py` (intermediate denominators
and decay), `pachucki_m_hpref.py` (exchange_hp references), `pachucki_m_quadcheck.py` (independent quadrature),
`pachucki_m_verify.py` (controls (b), (c) specialisation), `pachucki_m_opcompare.py`, `pachucki_m_opcompare_tests.py`.
Data (`debug/data/`): `pachucki_m_stage1_{1_2,4-5_3-2,sym}.json`, `pachucki_m_case_{2ppi,3ddel}_{1_2,4-5_3-2}.json`,
`pachucki_m_case_{2ppi,3ddel}_sym_res.json`, `pachucki_m_K_*.pkl`, `pachucki_m_symres_*.pkl`,
`pachucki_m_symcheck_*.json`, `pachucki_m_hpref_*_T24.json`, `pachucki_m_quadcheck_*.json`,
`pachucki_m_verify_{2ppi,3ddel,ppdel}.json`, `pachucki_m_opcompare_tests.json`, `pachucki_m_sym_*_forlimit.json`,
`pachucki_ext_limit_m_{2ppi,3ddel}.json` (written by the existing limit script), logs `pachucki_m_*.log`.
