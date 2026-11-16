# Record AY: (AA|BC) with m != 0 d pairs on A (inventory row 8, the last A-side item)

*Registered 2026-11-16, before computing, under the standing autonomous-loop authorization (PI, 2026-11-14) and the
PI's direction to work down the owed items (2026-11-16: "m != 0 d functions on A"). Engine and controls were built by
a sub-agent barred from the scored quantities (P1 comparisons, beta_eff at L >= 12, constant sets beyond L <= 8); the
PM read its report, derived the mechanism, computed the amplitude toys and wrote the predictions and the scoring
driver. A fresh reviewer read the draft and the driver cold (re-derived the mechanism and both ratio signs, held;
applied: MUST-FIX the scoring inputs added to the frozen set and the commit; the P4 bands re-centred on a 9/2-admixture
estimate that reproduces AS's measured PY ratios to 1.5%; P0's zero guard; a truncation report for c; wording).*

## Object

As in records AB and AS: B at the origin, C = (0, 0, D1), A = D2 (sin g, 0, cos g), x = cos g, D1 < D2; the BC pair is
1s_B 1s_C. A-frame: n = A/|A| (z'), e' = (cos g, 0, -sin g) (x', in the plane, perpendicular to B-A), e_y (y', out of
plane). Real normalized Slater 3d pairs on A, N_d (w M w) e^{-zA s}, tr M^2 = 3/2 (N_d as AS's DA):

| integral | w M w |
|:--|:--|
| DNX | sqrt3 (n.w)(e'.w) |
| DNY | sqrt3 (n.w)(e_y.w) |
| DXY | sqrt3 (e'.w)(e_y.w) |
| DXX | (sqrt3/2)((e'.w)^2 - (e_y.w)^2) |

with AS's DA (3d-sigma along B-A) recomputed as a control. Together with DA these are the five real 3d functions on A.

**Split.** I(x) = sum_L [f0_L P_L(x) + f2_L P_L^2(x) + f4_L P_L^4(x)], P_L^m = (1 - x^2)^{m/2} d^m P_L/dx^m (no
Condon-Shortley phase). L is the BC density's harmonic degree about B; m is the azimuthal order of the squared A
density about n (A's frame). f{m}_L uses the radial kernel functions G_L^{(j)}(D2), j <= 4, of that L only. Odd m
vanish (every squared density is even in x' and in y'). DNX and DNY carry m in {0, 2}, DXY and DXX m in {0, 4}. By
construction f0(DNY) = f0(DNX), f2(DNY) = -f2(DNX), f0(DXY) = f0(DXX), f4(DXY) = -f4(DXX) (checked exactly, below), so
there are two independent new integrals. This is not AS's lab-frame split. **Total** as a plain Legendre series:
c = f0 + leg[(1 - x^2) D^2 f2] + leg[(1 - x^2)^2 D^4 f4], an exact polynomial identity at truncation. c_j uses f at L >= j
only, c_j uses the f tail above j, so c is read at L <= 221 from a run to L = 400.

Engine `debug/aay_exact.py` (imports `aaa_exact.py`, `aapc_exact.py`, `aabc_exact.py`, `na_exact.py`, unchanged). For
any traceless M with tr M^2 = 3/2, the potential of the squared d density is (d M d)^2 U4(s) + (d M^2 d) U2(s) +
tr(M^2) U0(s) with universal radial U4, U2, U0 in the engine family (U4 = (4/9) x AS's DA U4).

## Controls (measured; `debug/aay_ap0.py`, `debug/aay_controls.py`, `debug/data/aay_controls.json`, `aay_part_*.json`)

| id | what | result |
|:--|:--|:--|
| AY0 | exact channels at L = 0, 1, 2, 3, 8 (m = 0, 2, 4) vs direct quadrature with no U kernels (Gauss-law potentials, A-frame harmonic projection, BC density moments), 60 digits, both cases, the four new integrals and DA | worst relative 5.0e-50 |
| AY0 | hard-coded kernels vs a fresh sympy Fourier derivation; the kernel identity vs the Gauss-law potential for all five M, at random points | equal exactly; 0 at 40 digits |
| AY0 | jump terms at r = D2 (L <= 8); fourth and second D2-derivatives vs 5- and 3-point differences (400 digits) | zero; <= 5.5e-80 |
| AY0 | DA and PA vs `aaa_exact` key for key (L <= 12); PY, PX' totals in this split vs record AS's totals | equal; <= 8e-28 at L <= 120 |
| AY0 | total_coeffs re-expansion vs direct evaluation of the split | 2e-81 |
| SR | sum rule, exact over Q, key for key, L <= 30, both cases: f0 of DA + DNX + DNY + DXY + DXX = 5 S (S the spherically averaged pair, own kernel); f2(DNX) + f2(DNY) = 0; f4(DXY) + f4(DXX) = 0; f0(DNX) = f0(DNY); f0(DXY) = f0(DXX) | all exact |
| R2 | Gaussian fits n = 30 -> 34 at A1 x = 1/2, -1/3 and A2 x = 0, 4/5 (16 values); fit floors at n = 34 against exact one-axis cases at R = 2 (DNX = DNY vs (3d1 3d1 \| 1s 1s), DXX = DXY vs (3d2 3d2 \| 1s 1s), DA vs 3d0) | spread <= 1.2e-12; floors 1.9e-13, 2.0e-13, 2.6e-13 |
| R1 | prolate quadrature 120^2 x 144 vs R2 n = 34 | <= 3.9e-13 |
| T | toy beta_eff(220) on coefficient / rho^L, rho = 1/(x* + sqrt(x*^2 - 1)), x* = (D1^2 + D2^2)/(2 D1 D2): P_L of (x* - x)^alpha; P_L^2 of (1 - x^2)(x* - x)^alpha; P_L^4 of (1 - x^2)^2 (x* - x)^alpha | A1: P_L 5/2 -3.070, 7/2 -4.119, 9/2 -5.180; P_L^2 5/2 -5.180, 7/2 -6.255, 9/2 -7.342; P_L^4 5/2 -7.342, 7/2 -8.443, 9/2 -9.556. A2: P_L -3.055, -4.094, -5.143; P_L^2 -5.143, -6.203, -7.274; P_L^4 -7.274, -8.354, -9.446 |
| TA | amplitude toys at L = 218..222 (`python debug/aay_run.py toys`, `debug/data/aay_toys.json`): ta (P_L, 7/2), tb2 (P_L^2), tb4 (P_L^4), tc (P_L, 9/2) | ta, tb2, tc equal record AS's `amp_toys` to every printed digit; x*^2 - 1 = 0.141427 (A1), 0.237656 (A2); (x*^2 - 1)^2 = 0.020002, 0.056480 |
| P | precision ladder 1200 ... 3600 digits, L <= 400, self-differences only (no value stored or printed) | coefficients at L = 219..221 stable to 1e-10 from 1200 digits; L = 398..400 and c at 219..221 from 2000 (A1) / 2400 (A2); 3000 vs 3600 <= 1e-992 everywhere |

Timings (exact over Q, L <= 400): the kernel functions are shared by all five integrals, 2.5 h (A1) and 1.7 h (A2);
assembly 3-4 s per integral.

**Already known (disclosed):** a smoke test at L <= 8 (mod 2^61 - 1) showed every channel's constant set equal to
AB's 21. A smoke run of the scoring driver at L <= 12 (60/80 digits) checked only that it runs; its output file was
deleted unread. Nothing about the new coefficients beyond L = 12 has been computed or seen.

## Mechanism (derived by the PM; the reason DNY, DXY differ from DNX, DXX)

As in AS: the non-analytic part at x* comes from the cusp of the BC density at C against the odd powers of the A
potential. For a d pair, the degree-4 multipole of the squared density gives |a|^7 h_4(a), a = A - C, the 7/2 branch;
degree 2 reaches |a|^9 and degree 0 |a|^11 (AS's control X). On the null cone |a|^2 = 0, the degree-4 harmonic part of
(w M w)^2 equals (a M a)^2, since the harmonic projection differs from the polynomial by multiples of |a|^2. Here
a = D2 n - D1 z, so a_y = e_y . a = 0 for every x, n . a = (D2^2 - D1^2)/(2 D2) at x* (nonzero), and e' . a = D1 sin g
(nonzero at x*). Hence:

- **DNY** ~ (n.a)^2 (e_y.a)^2 = 0 and **DXY** ~ (e'.a)^2 (e_y.a)^2 = 0: **the totals lead with (x* - x)^{9/2}**.
- **DNX** ~ (n.a)^2 (e'.a)^2 and **DXX** proportional to (e'.a)^4, both nonzero (up to normalization): **the totals keep 7/2**.

**The parts.** Write y'^2 = rho'^2 (1 - cos 2 phi)/2 and x'^2 y'^2 = rho'^4 (1 - cos 4 phi)/8 about n. For DNY the m = 0
part of the squared density is (n.w)^2 (x'^2 + y'^2)/2 and the m = 2 part -(n.w)^2 (x'^2 - y'^2)/2; at the null point
(a_y' = 0) these are +(n.a)^2 a_x'^2 / 2 and -(n.a)^2 a_x'^2 / 2. So f0 and f2 each carry 7/2 and cancel in the total.
For DXY the m = 0 and m = 4 parts are +a_x'^4/8 and -a_x'^4/8: f0 and f4 each carry 7/2 and cancel. As in AS, the
cancellation should hold at every B-gradient order, since every harmonic built from the tensor e_y e_y . . . reduces on
the null cone to terms in a_y or |a|^2; this extension is argued, not derived term by term here.

**Amplitude form.** If f0 ~ s ta (leading), the m-part's function is -s (x* - x)^{7/2} near x*, i.e.
(1 - x^2)^{m/2} times -s (1 - x*^2)^{-m/2} (x* - x)^{7/2}, so f_m ~ -s (1 - x*^2)^{-m/2} tb_m. The ratio
R_m = (f0 / ta) / (f_m / tb_m) at leading order is

    DNY (m = 2):  R_2 = -(1 - x*^2)   = +(x*^2 - 1)      [AS's PY had the same value in its lab split]
    DXY (m = 4):  R_4 = -(1 - x*^2)^2 = -(x*^2 - 1)^2    (negative)

DNX and DXX have R equal to minus these, identically (their f_m are negated). At L = 220 the 9/2 terms present in one
part and absent in the other move R away from its limit. **Estimate (the cold reviewer's, checked by the PM):** f0's amplitude carries a_x'^m,
proportional to (1 - x^2)^{m/2}, which cancels out of the m-part's h but not out of f0; that leaves an f0-only 9/2 term
with relative coefficient k = m x*/(x*^2 - 1), so R/limit ~ 1 + k |tc/ta| with tc/ta = -0.00818 (A1), -0.01048 (A2) at
L = 220. For AS's PY (m = 2) this gives 1.124 / 1.098 against the measured 1.108 / 1.087. Here: DNY (m = 2) about 1.12 /
1.10; DXY (m = 4) about 1.25 / 1.20 (squaring AS's measured ratios gives 1.228 / 1.182). The bands below are set around
these.

## Predictions (scored at A1 and A2)

- **AY0** passed (above). Frozen at this commit and unchanged at scoring: `debug/aay_exact.py`, `aaa_exact.py`,
  `aapc_exact.py`, `aabc_exact.py`, `na_exact.py`, `aabc_run.py`, the driver `debug/aay_run.py`, the control drivers
  `debug/aay_ap0.py`, `aay_controls.py`, and the scoring inputs `debug/data/aay_controls.json` (referee values),
  `aay_toys.json` (P4 denominators), `aabc_p2_<p>.json` (AB's 21), `aaa_controls.json`, `aay_part_*.json`.
- **P0 (precision gate).** `python debug/aay_run.py p1p3 400 3000 3600 <case>`, one exact build per case evaluated at
  3000 and 3600 digits: every part's coefficients at L = 219..221 and 398..400, every total c at 219..221, and all
  summed values agree between the two precisions to relative 1e-10 (on the 25 and 30 printed digits), and
  `total_reexpansion_check` is below 1e-30 at all four points for the four new integrals. Failure is VOID.
- **P1 (values).** The series summed to L = 400 match R2 (n = 34) within **1e-10** at both points of both cases, for
  DNX, DNY, DXY, DXX (16 values). Probability 0.85. (DA reported as a control.)
- **P2 (constants).** Mod 2^61 - 1 and 4611686018427387847, to L = 1000: for each new integral, f0 (from L = 2), f2 or f4
  (from L = 2, resp. L = 4, its first nonzero) and their union each use exactly AB's 21 constants at the same case, the
  set fixed in L; nothing outside AB's candidate set; no Euler gamma. Probability 0.75.
- **P3 (exponents)**, beta_eff at L = 220 on coefficient / rho^L:
  - DNX: f0 in [-4.65, -3.6], f2 in [-6.8, -5.7] (7/2 in each basis); **total c in [-4.65, -3.6]** (7/2).
  - DNY: **total c in [-5.7, -4.65]** (9/2). (Its f0, f2 equal DNX's up to sign.)
  - DXX: f0 in [-4.65, -3.6], f4 in [-9.0, -7.85]; **total c in [-4.65, -3.6]** (7/2).
  - DXY: **total c in [-5.7, -4.65]** (9/2).
  Probability 0.7 (the two 9/2 totals: 0.6).
- **P4 (the cancellation, amplitude form)**, at L = 220:
  - DNY: R_2 in **[0.95, 1.30] x (x*^2 - 1)**, positive.
  - DXY: R_4 in **[-1.45, -1.00] x (x*^2 - 1)^2**, negative.
  Reported, not scored: R at 219, 221; R for DNX, DXX (minus the above by construction); c_220 / tc_220 for DNY, DXY
  and c_220 / ta_220 for DNX, DXX; c at 219..221 recomputed from f cut at L = 380 (truncation). Probability 0.65.

**Reported, not scored:** beta_eff at L = 20, ..., 380 for every part (c to L = 220); signs; timings; DA's beta and
values (control; AS's DA was scored there).

## Decision gate

- **VOID** if AY0, P0 or P1 fails, or a frozen file changed.
- **GO** if P1, P2, P3 and P4 hold.
- **STOP** if the total c of DNY or DXY has beta_eff(220) above -4.65 at either case (no cancellation: the mechanism is
  wrong), or any f0 or total c is outside [-5.7, -3.0], f2 outside [-7.9, -5.2] or f4 outside [-10.1, -7.3] (a
  leading branch other than 7/2 or 9/2), or a P4 ratio has the wrong sign or a magnitude outside [0.5, 2] times its
  limit.
- **BORDERLINE** otherwise (including P2 failing with P1, P3 and P4 holding).

## Caps

Two geometries, two primes, one orientation of the A-frame axes (x' in the plane); 3d functions only, against a 1s
pair at B and C. Constants assumed linearly independent over Qbar. The extension of the cancellation beyond the
leading B-gradient order is argued, not derived. Exchange integrals and l > 2 on A are not covered.

## Phase 2 (after this record is committed and the commit has been verified)

`python debug/aay_run.py p1p3 400 3000 3600 A1` and `... A2` (in parallel; `debug/data/aay_p1p3_<case>.json`);
`python debug/aay_run.py p2 2305843009213693951` and `... p2 4611686018427387847` (`debug/data/aay_p2_<p>.json`).
