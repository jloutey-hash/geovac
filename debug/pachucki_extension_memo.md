# Pachucki exchange bound: general exponents, the a -> b limit, and l > 0

Date 2026-11-18. Derivation memo; no paper, record, test or git state touched.
Engine: `debug/pachucki_ext_core.py` (two changes from the reviewed `tests/_pachucki_exchange.py`:
(i) the sources d^alpha F, d^alpha p^2 are built as truncated Taylor jets in the four exponents, for any
downward-closed set of multi-indices, not only subsets; (ii) (r f')' = src is solved by a direct triangular
solver, with no generic rref, so the coefficient field can be Q or Q(t)). `debug/pachucki_ext_check1s.py`: identical sources and f_S to the helper at (1,2), (1,3/2). Drivers: `pachucki_ext_symbolic.py` (Part A),
`pachucki_ext_limit.py` (+ `_limit_num.py`) for a -> b, `pachucki_ext_lgt0.py` (over Q) and
`pachucki_ext_lgt0_sym.py` (over Q(t)) for Part B, `pachucki_ext_hpref.py` for exchange_hp references.
Data: `debug/data/pachucki_ext_*.json`.

## Part A: general a < b, done symbolically

**Method.** Scaling r -> r/a gives K(R; a, b) = a K(aR; 1, b/a) and f0(R; a, b) = a^-2 f0(aR; 1, b/a), so
one parameter t = b/a suffices. I ran the 15 derivative equations over the field Q(t) (sympy FracField, exact),
not at sample points: 3 s. Every residual (r f_S')' + F_S - r sum (d_U p^2) f_{S-U} is identically 0 as an
element of Q(t)[R, 1/R] x keys; no free unknowns.

**Results (all identically in t).**
* K lies in the Q(t)[R, 1/R]-span of the nine functions (the two ties, Ei((a-b)R) = -2 x Ei(2(a-b)R) and
  LC-coefficient = - Ei(2(a-b)R)-coefficient, and their b-twins, hold as identities in t).
* The nine-function span is closed under d/dR over Q(t)(R) (exact check).
* Rank of K, ..., K^(9) in those coordinates is 9 at (t, R) = (3, 7/5); since rank can only drop under
  specialisation and there are nine coordinates, the formal rank over Q(t)(R) is exactly 9.
* Denominators factor as (t-1)^i (t+1)^j only, i <= 5; in (a, b): (b-a)^5 (a+b)^5. So the statement holds for
  all real 0 < a < b (an identity of rational functions; uniqueness by decay holds pointwise for t > 1).
* Specialised to t = 2 and 3/2 the symbolic closed form equals the helper's K exactly (so it inherits the
  50-digit agreement with exchange_hp recorded there).
* Closed form in (a, b): `pachucki_ext_symbolic_1s.json`, `closed_form_coefficients_ab`; e.g. the master integral
  enters with 16 a^5 b^5 [9 + 6(a^2+b^2) R^2 + (a^4 + 4a^2 b^2 + b^4) R^4] / (R (b^2-a^2)^4).
* Extra pairs (1,3), (2,3), (1,5/2), (3/4,1), run separately over Q: nine-function span, rank 9, and a unique
  operator of order 9 and degree 16 (none of degree 15) at each pair, as at (1,2) and (1,3/2).

**The a -> b limit (`pachucki_ext_limit.py`).** With t = 1 + eps every function is expanded in eps: the
elementary and E1 pieces by Taylor's theorem; each psi bracket via Ei(y) = gamma + ln|y| + sum y^k/(k k!), where
ln eps cancels inside the bracket (psi_a -> -e^{-2R}(gamma + ln R) + O(eps)); and the master integral by
f0(R; 1+eps) = sum eps^k I_k(R), I_k = int_R^oo FM1_k(s)/s ds, FM1_k the eps^k coefficient of label 60. Each I_k
was reduced by parts (each reduction rule checked by differentiation) to the class
{R^j e^{-2R}, R^j e^{-2R}(gamma+ln R), E1(2R), R^j e^{2R}E1(4R)} plus three integrals outside it
(int e^{-2s}(gamma+ln s)/s, int e^{2s}E1(4s)/s, int E1(2s)/s). Results:
* all poles eps^-1 .. eps^-6 cancel exactly;
* no I_k contains any of the three outside integrals (so every eps-coefficient of the master integral is
  already in the four-function class);
* the eps^0 term equals Sugiura's formula exactly (symbolic difference 0), in the form
  (1/5)[-e^{-2R}(-25/8 + 23R/4 + 3R^2 + R^3/3) + (6/R)(S^2(gamma+ln R) - S'^2 E1(4R) + 2SS' E1(2R))],
  S = e^{-R}(1+R+R^2/3), S' = e^R(1-R+R^2/3); it also equals exchange_hp at R = 1.4 to 35 digits;
* its minimal operator (order 4, degree 8) equals the committed (1,1) operator `ep_heldout_1s_1_1_EPT.json`
  over Q up to a scalar.

So the nine functions degenerate, as a -> b, to the four-function span e^{-2p}, e^{-2p}(gamma + ln p), E1(2p),
e^{2p}E1(4p), and the order bound 4 at equal exponents is recovered as the limit of the derived closed form.
(A re-derivation of a classical result; its value is the continuity check at a = b.)

## Part B: l > 0

**The map (m = 0 on the axis).** Each orbital product is N_A N_B P(r_A^2, r_B^2; R) e^{-a r_A - b r_B} per
electron, with z_A = (r_A^2 + R^2 - r_B^2)/(2R), z_B = z_A - R. Since r^n e^{-ur} = (-1)^{n+1} d_u^{n+1}[e^{-ur}/r],
a monomial r1A^n1 r1B^n2 r2A^n3 r2B^n4 maps to (-1)^{sum n} d^alpha f with alpha = (n1+1, n2+1, n3+1, n4+1) in
the order (u3, u2, w2, w3), and K = N_A^2 N_B^2 (4 pi)^2 / R x sum (R-Laurent coefficient) d^alpha f. All n are
even here, so all signs are +. No z-derivative is needed: z enters only through the distance polynomial. Cases
and normalisations: 2p_sigma: P = z_A, N^2 = a^5 b^3/pi^2; 3d0: P = 3z_A^2 - r_A^2, a^7 b^3/(18 pi^2);
2p0-2p0: P = z_A z_B, a^5 b^5/pi^2; 2p+1-2p+1 (density x^2+y^2, m cancels): P = r_A^2 - z_A^2,
a^5 b^5/(4 pi^2). The equation (*) for a general multi-index carries multinomial binomials:
(r f_alpha')' = -F_alpha + r sum_{0 != U <= alpha} binom(alpha, U) (d_U p^2) f_{alpha-U}.

Results at (a, b) = (1, 2) over Q (`pachucki_ext_lgt0_<case>_1_2.json`), every residual identically 0, no
free unknowns:

| case | multi-indices | solve | span | rank | derived operator vs record |
|---|---|---|---|---|---|
| 1s (control) | 16 | 0.1 s | nine | 9 | (9,16) unique; values = record to 50 digits |
| 2p_sigma_A 1s_B | 144 | 1.9 s | nine | 9 | (9,34) unique, none at degree 33; = EP's operator_Q_crt over Q (scalar -1) |
| 2p0_A 2p0_B | 576 | 27 s | nine | 9 | (9,53) unique; = EQ's operator_Q_crt over Q |
| 2p+1_A 2p+1_B | 576 | 27 s | nine | 9 | (9,39) unique; = EQ's operator_Q_crt over Q |
| 3d0_A 1s_B | 576 | 27 s | nine | 9 | no record at (1,2); derived (9,53) unique, none at degree 52, coefficient degrees [51,52,53,...] |

Numeric check of the map (closed form, f0 by one 1-D quadrature at dps 50, against
`geovac.qfd_core.exchange_hp` at tau_max 24, dps 50, R = 1.4; `pachucki_ext_hpref_<case>_1_2_1.4_T24.json`,
tail terms <= 1e-54):

| case | closed form | exchange_hp | agreement |
|---|---|---|---|
| 2p_sigma | 0.220068065310325593968393583117085511012997133616 | identical | all 48 printed digits |
| 2p0-2p0 | 0.0552386507929406934123306194991034219264309529881 | ...529845 | 47 digits (diff 3.6e-48) |
| 2p+1-2p+1 | 0.142646548127204270123078769771284851563274015302 | identical | all 48 digits |
| 3d0 | 0.0669403480023838429687222783879717905171316950222 | ...950218 | 48 digits (diff 4e-50) |

The 2p_sigma check (which timed out before) took 1618 s. This validates the map and normalisations.

Equal exponents through the limit, Q(t) runs (`pachucki_ext_sym_<case>.json`, `pachucki_ext_limit_<case>.json`):

| case | Q(t) solve | denominators | max pole in t-1 | poles cancel | outside integrals | limit's operator vs record (1,1) |
|---|---|---|---|---|---|---|
| 1s | 3 s | t-1, t+1 | 5 | yes | none | (4,8) = EPT operator over Q |
| 2p_sigma_A 1s_B | 230 s | t-1, t, t+1 | 7 | yes | none | (4,16) = EP operator over Q |
| 2p0-2p0 | 3553 s | t-1, t, t+1 | 9 | yes | none | (4,24) = EQ operator over Q |
| 3d0_A 1s_B | 3607 s | t-1, t, t+1 | 9 | yes | none | (4,24) = EQ operator over Q |
| 2p+1-2p+1 | 3700 s | t-1, t, t+1 | 9 | yes | none | (4,16) = EQ operator over Q |

Every m = 0 case: valid identically in t (the extra denominator t = b/a is harmless), and the a -> b limit lies in
the four-function span and reproduces the committed equal-exponent operator exactly.

**m != 0 (2p+-1 1s exchange, 3d+-2): not done; plan and obstacle.** The real part of (x1+iy1)(x2-iy2) is
x1x2 + y1y2 = (r1A^2 + r2A^2 - r12^2)/2 - z1A z2A, so these integrals need r12-weighted masters (Pachucki's f_n,
n = 2 for 2p+-1, up to n = 4 for 3d+-2). The 2009 paper already has the r12 recursion (arXiv labels 45, 46; label 46
is f_2 explicitly in f, f', f'', exponentials and F_i^+). Obstacle: label 45 divides by
(u2^2-u3^2+w2^2-w3^2)(u2^2 w2^2-u3^2 w3^2) and label 46 by p^2; both vanish to second order at the exchange point
and identically along the exchange locus (u3 = w2, u2 = w3), i.e. along exactly the directions we differentiate in;
the Hessian of p^2 there has rank 2. So f_2 at the exchange point is 0/0. Plan: evaluate label 46 along a transverse
line theta0 + eps v (p^2 ~ eps^2), as a jet in eps and the four exponents, cancelling 1/p^2 against jets of f and
F^+ that the present engine already produces (cost: f_alpha to total order ~ |alpha| + 2 + 2k, ~3000 multi-indices,
est. ~1 h, not run); cross-check against an f_2 equation from the w1-equation (label 39; the 2012 paper's Eq. 13 gives
the r-ODE in the r12 exponent - not transcribed); then the same span/rank/operator pipeline. Since label 46 adds
only f, f', f'', exponentials and F^+ (same Ei arguments), the nine-function span is the expected outcome, consistent
with EP 2p_pi (9,25) and EQ 3d_delta (9,35), but not established.

## Proposed Paper 58 sentences and tiers

1. **[SYMBOLIC PROOF], all real 0 < a < b (1s).** "Differentiating Pachucki's second-order equation for the
two-centre master integral once in each exponent at the exchange point gives fifteen inhomogeneous equations; with
the exponent ratio t = b/a kept as an indeterminate, each has a solution in the span over Q(t)[R, 1/R] of eleven
decaying functions, verified as an identity of rational functions in t, and the exchange integral lies in the span
of the nine functions listed, with coefficients in Q(a, b)(R) whose denominators are powers of b - a and b + a. The
order is therefore at most 9 for every 0 < a < b." (Exact algebra over Q(t), uniqueness by decay; "exactly 9"
still needs the independence lemma or stays [MEASURED] from the guess scan, as before.)
2. **[SYMBOLIC PROOF], limit.** "As b -> a the poles in b - a cancel and the closed form tends to Sugiura's
equal-exponent formula, in the span of e^{-2p}, e^{-2p}(gamma + ln p), E1(2p), e^{2p}E1(4p); the order bound 4 at
equal exponents is recovered as this limit." (Exact series in eps; Sugiura is prior art; the limit is a consistency
check, credit stays with Sugiura/Ruedenberg.)
3. **[SYMBOLIC PROOF] at listed cases, identically in a < b; operators compared at (1,2) and (1,1).** "The same
derivation, with the l > 0 orbital products written as polynomials in r_A^2, r_B^2 and R and mapped to higher
exponent derivatives of the master integral, places (2p_sigma_A 1s_B|2p_sigma_A 1s_B), (3d_0 1s|3d_0 1s),
(2p_0 2p_0|2p_0 2p_0) and (2p_{+1} 2p_{+1}|2p_{+1} 2p_{+1}) (density x^2+y^2, m cancels) in the same nine-function span for all
0 < a < b, so their order is at most 9; at (a,b) = (1,2) the derived operators equal the guessed ones over Q
(9, 34; 9, 53; 9, 39; the 3d_0 case, unrecorded at (1,2), gives a unique (9, 53)), and their b -> a limits reproduce the guessed equal-exponent operators (4, 16; 4, 24; 4, 24; 4, 16)."
Numeric: the closed forms match exchange_hp at R = 1.4 to 47-48 digits ([MEASURED] cross-check, independent route).
4. **Not established:** 2p_{+-1} 1s exchange (EP 2p_pi) and 3d_{+-2} (EQ 3d_delta): need r12-weighted masters;
no claim beyond "plan identified". Do not extend sentence 3 to them.

## What a regression guard must pin

For a separate guard pass (each with the wrong answer it must reject):
* Over Q(t): all 15 residuals identically 0 (rejects a mistranscribed jet/source, e.g. a wrong a+- or a dropped
  binomial); K's two tie identities hold as rational-function identities (rejects ties true only at sampled
  pairs); denominators of the K coefficients factor into (t-1), (t+1) only for 1s (rejects a hidden pole
  t0 > 1, e.g. at t = 3); specialisation at t = 2, 3/2 equals tests/_pachucki_exchange.py's K (rejects drift
  between the reviewed helper and the jet engine).
* The coefficient of f0 in K equals 16 a^5 b^5 [9 + 6(a^2+b^2)R^2 + (a^4+4a^2b^2+b^4)R^4] / (R (b^2-a^2)^4)
  (a pinned closed-form entry; rejects a wrong normalisation 16 a^3 b^3, or wrong scaling a^{p+1} vs a^{p+3}).
* Limit: poles eps^-1..eps^-5 cancel, J1 = J2 = J3 coefficients 0, limit - Sugiura = 0 symbolically; and the
  limit at R = 1.4 equals 0.32329114155307318238571155893367021 (rejects a sign flip in label 60's
  e^{2bR}-e^{2aR} term, which was the one bug met here and leaves poles uncancelled).
* Part B at (1,2): unique (9, deg) operator with none at deg-1, equal over Q to the committed operator_Q_crt for
  2psig (34), ppsig (53), pppi0 (39) - and the closed-form values at R = 1.4 above to >= 45 digits (rejects a
  wrong map, e.g. odd powers / sign (-1)^{|n|}, or alpha shifted by one).
* Limits for 2psig/ppsig/3dsig/pppi0: operator equals the committed (1,1) operator (4,16)/(4,24)/(4,24)/(4,16).
