# Pre-registered prediction: two-centre integrals over p and d orbitals stay in the 1s class (R8)

**Registered 2026-11-14** (chemistry period inventory, row 8, part 1: rows 2, 3a, 3b at l <= 2), before the panel
below is computed. **Commit status:** draft; to be committed on its own, before any computation it scores.
**Rests on:** Paper 63 sec:dictionary (Lemma orbital, Lemma integral, Proposition rung: every exponential-type orbital
integral lies in the D-module of its graph's 1s integral; for l > 0 "neither the bound nor attainment is shown");
Paper 58 sec:qfd and its class table (two-centre classes: (AA|BB) elementary; hybrid {E1, ln}, Ei for a_d > a_c,
gamma + ln R at a_d = a_c; exchange {E1, ln, gamma}); the v6.7.26 R-operator panel (`debug/two_centre_r_operators.py`,
data `debug/data/two_centre_r_operators.json`, 30 cases, 1s plus one 2p0 hybrid family).
**Guardrail check (Papers 8-9, FCI-M, Track DF):** not in scope; no molecular encoding is proposed.

## Question

The dictionary says a p or d orbital is a derivative of a 1s Yukawa in its exponent and its centre. Derivatives stay in
the D-module, but the dictionary does not bound the class on a fixed slice. Does going from 1s to p and d orbitals ever
"climb" in the two-centre layer: a new function type, a genuine (non-apparent) finite singular point in R, or a
logarithm at R = 0 in a class that has none for 1s?

## Object and instruments

Hydrogenic orbitals (centre, Z, n, l, m), rate a = Z/n, complex Y_lm, A at the origin, B at R zhat.

- **Row 2 (one-electron), new engine** `debug/row8_one_electron.py`: S = <a|b>, VA = <a|1/r_A|b>, VB = <a|1/r_B|b>,
  T = <a|-1/2 nabla^2|b>, exact in prolate spheroidal coordinates (Mulliken A_i, B_j; no quadrature). T by two routes
  (radial Laplacian; T|b> = E_b|b> + Z_b/r|b>).
- **Rows 3a, 3b (two-electron), existing engines**, unchanged: `geovac.two_center_eri.aabb_closed_form`,
  `hybrid_closed_form`.
- **Analysis, existing and unchanged:** `debug/two_centre_r_operators.py` `analyse()`: atom decomposition
  (e^{cR} x {1, ln R, E1(kR), Ei(kR)} x Q(R), constants gamma and logs of rationals split), minimal operator in R over
  Q[R] from exact data, symbolic and 60-digit verification, R = 0 indicial roots and the ln R coefficient of the local
  expansion, finite singular points with exponents, rates at R = infinity.
- **Driver (scoring):** `debug/row8_panel.py` (written with this record; its panel is the list below).

## Panel (fixed here; nothing is added after computing)

Notation `nl m`, Z on each centre in brackets. Rate relations for the hybrid: a_c = Z_A/n_c, a_d = Z_B/n_d.

**Row 2, two-centre pairs (a on A, b on B), each giving S, T, VA, VB (24 integrals):**
E1 (2p0 [3] | 1s [1]); E2 (2p1 [2] | 2p1 [2]); E3 (3d0 [3] | 2p0 [2]); E4 (3d2 [3] | 3d2 [3/2]);
E5 (3p1 [3] | 3d1 [2]); E6 (2s [2] | 3d0 [3]).
**Row 2, one-centre density on A, potential at B (VB(a, a'), 6 integrals):**
F1 (2p0, 2p0) [3]; F2 (2p1, 2p1) [2]; F3 (3d0, 3s) [3]; F4 (3d1, 2p1) [3]; F5 (3d2, 3d2) [3]; F6 (2p0, 2s) [2].

**Row 3a, (ab|cd)_A,A,B,B (8 quartets):**
C1 (3d0 3d0 | 1s 1s) [3; 1]; C2 (2p1 2p1 | 2p0 2p0) [2; 2]; C3 (3d1 3d1 | 3d2 3d2) [3; 3/2];
C4 (2p0 3d0 | 1s 2p0) [3; 2]; C5 (3d2 3d0 | 3d0 3d2) [3; 3]; C6 (2p1 2p0 | 3d0 3d1) [2; 3];
C7 (2s 3p0 | 2p0 3s) [3; 2]; C8 (3d2 3d2 | 3d2 3d2) [3; 3].

**Row 3b, (ab|cd)_A,A,A,B (10 quartets):**
H1 (3d0 3d0 | 3d0 1s) [3; 1] eq; H2 (2p1 2p1 | 2p1 2p1) [2; 2] eq; H3 (3d0 3d0 | 1s 1s) [3; 1] a_d < a_c;
H4 (2p0 2p0 | 2p0 3d0) [2; 3] eq; H5 (3d1 3d1 | 3d1 2p1) [3; 3] a_d > a_c; H6 (3d2 3d2 | 3s 2s) [3; 1] a_d < a_c;
H7 (2p0 3d0 | 2s 1s) [2; 1] eq; H8 (3p1 3d1 | 3d0 2p0) [3; 4] a_d > a_c; H9 (2s 2s | 3d0 3d0) [2; 3] a_d > a_c;
H10 (3d2 3d0 | 3d-1 3d1) [3; 3] eq.

(Every quartet satisfies the M rule m_b - m_a + m_d - m_c = 0, checked by inspection of the labels. Quartets that are nonetheless identically zero are reported as zero and excluded from scoring; the driver records them.)

## Controls measured before registering

Scratch drivers (session scratchpad; the panel driver repeats them as its first block):

| Control | what | measured |
|:--|:--|:--|
| K1 new-function refusal | `decompose` on e^{-R}(1+1/R) + Li2(e^{-R}); on e^{-R} + E1(R)^2 | both REFUSED (AssertionError) |
| K2 genuine finite singular point | `finite_singular` on the operator of e^{-R}/(R-3); of e^{-R}(R-3)^2 | R = 3 exponent -1, apparent False; exponent 2, apparent True |
| K3 ln R at R = 0, positive | 1s exchange (AB\|AB), zeta = 1 | lnR0 True; order 4; finite point apparent (exponents 0, 1, 2, 4) |
| K4 ln R at R = 0, planted | e^{-R} ln R + e^{-3R}/R | lnR0 True, order 3 |
| K5 one-electron engine vs qfd_core (1s, Z = 1 on A, 2 on B) | S, VA, VB (density on A), T, symbolic difference | 0, 0, 0, 0 |
| K6 one-electron engine vs independent quadrature (scipy dblquad, sph_harm + genlaguerre), (2p1 [3/2] \| 3d1 [5/2]) at R = 7/5 (not a panel case) | S, VA, VB | 1.0e-15, 6.4e-16, 3.5e-15 |
| K7 T routes, same pair | T_lap - T_eig, symbolic | 0 |
| K8 feasibility, two-electron engines + `analyse` at l = 2 | (3d0 3d0\|1s 1s) [3;1]; (2p1 2p-1\|2p0 2p0); hybrid (3d0 3d0\|3d0 1s) [3;1]; (2p1 2p1\|2p1 2p1) [2;2] | 4.0 s, 4.1 s, 32 s, 7.5 s; all verified symbolically |

**Seen before registering (disclosed):** K8's four cases were analysed in full. Two of them are panel cases (C1 = the
first; H2 = the fourth). Their readings: (AA|BB) order 2, rates {-2, 0}, R = 0 exponents {-5, 0}, no ln R at 0,
finite points apparent; hybrids order 4, rates {-3, -1, 1}, no ln R at 0, finite points apparent (exponents 0, 1, 2, 4).
Together with the three 2p0 hybrids of the v6.7.26 panel (orders 4, 5, 5 against 3, 4, 4 distinct rates), these seven
l > 0 readings are what P6's bound was set from. C1 and H2 are scored anyway; P6 is flagged as fitted on them.

## Predictions

- **R8-0 (instrument; run first).** K1-K4 reproduce. Every row-2 S, VA and VB of the panel agrees with the independent
  quadrature at R = 7/5 to <= 1e-11 (absolute), and T_lap = T_eig symbolically for every pair. Every two-electron
  closed form passes `analyse`'s symbolic verification. Falsifier: any failure. Probability 0.95.

- **P1 (row 2 class).** Every one-electron integral of the panel decomposes into atoms of kind "1" only (e^{cR} times a
  rational function of R times a constant algebraic normalization): no E1, Ei or ln R. Probability 0.98.
  *Disclosed:* this is close to guaranteed by the engine's construction (A_i, B_j are exponential-polynomial). Its live
  content is that the 1/r_C cancellation leaves a polynomial for every l <= 2 (no residual 1/(xi +- eta)), which the
  engine asserts.

- **P2 (row 3a class).** Every (AA|BB) quartet of the panel is elementary: atoms of kind "1" only, no gamma, no log of a
  number. Probability 0.97.

- **P3 (row 3b class).** (P3a) Every hybrid quartet's atoms lie in {1, E1, Ei, ln R}, with no square or product of
  transcendental functions (refused by `decompose`, K1) and with constants among gamma and logs of rationals.
  Probability 0.93. (P3b) ln R and gamma appear in a hybrid quartet exactly when it is an eq case (a_d = a_c), and Ei
  exactly when a_d > a_c. Probability 0.75.

- **P4 (no logarithm at R = 0 outside exchange).** For every case of rows 2, 3a and 3b the local expansion at R = 0
  carries no ln R (`analyse`'s lnR coefficient identically zero). For 1s this is the v6.7.26 reading (only exchange
  has one). Probability 0.8.

- **P5 (singular points).** For every case: R = 0 is a regular singular point (Fuchsian) or regular; every finite
  nonzero singular point is apparent (simple roots: `apparent` True; multiple roots: numeric exponents distinct
  non-negative integers to 1e-8); the rates at R = infinity are among the atom rates. Probability 0.9.

- **P6 (order bound; fitted on seven readings, see above).** order_QR <= (number of distinct atom rates) + 1 for every
  case, and for rows 2 and 3a order_QR = number of distinct atom rates. Probability 0.6.

## Decision gate

- **VOID** if R8-0 fails.
- **GO** if P2, P3a, P4 and P5 hold: no climb at l <= 2 in the two-centre layer, in function type, singular points or
  local logarithms.
- **BORDERLINE** if exactly one of P2, P3a, P4, P5 is missed on at most two cases with the mechanism identified.
- **STOP** otherwise.
- P1, P3b and P6 are scored and reported; they do not move the gate.

## Novelty

That two-centre integrals over Slater or hydrogenic orbitals of any l lie in elementary functions, E1, ln and gamma is
classical (Roothaan, Ruedenberg; Lesiuk-Moszynski 2014 for arbitrary l with exact evaluation; credited in Paper 58).
P1-P3 are therefore checks, not discoveries. New here, if they hold: the R-operator, singular-point and R = 0
statements (P4-P6) at l <= 2, the absence of a local logarithm outside exchange, and the inventory's row-2 closed forms
at l <= 2 with a guard.

## Caps

Hydrogenic orbitals only (rate Z/n with the Laguerre polynomial), n <= 3, l <= 2, the panel above. Exchange (AB|AB) at
l > 0 is not covered (no closed form at unequal exponents; record part 2). The three-centre (AA|BC) with p orbitals and
the Paper 63 fibre question for l > 0 are part 2, a separate record. The R-operator is read at fixed exponents; nothing
is claimed about the exponent dependence.

## Results (2026-11-14)

Committed as `7784d751` before the panel ran (the header's "draft" line predates the commit and was left as
registered). Driver `debug/row8_panel.py`, unchanged; data `debug/data/row8_panel.json`, log
`debug/data/row8_panel.log`. 48 cases, none identically zero.

| Prediction | verdict | reading |
|:--|:--|:--|
| R8-0 instrument | HELD | controls K1-K4 reproduce; 30 row-2 S/VA/VB vs independent quadrature at R = 7/5, max 1.2e-14; T_lap = T_eig for all 6 pairs; all 48 closed forms verified symbolically |
| P1 row 2 class | HELD | all 30 one-electron integrals: kind "1" only (constants: normalization surds only) |
| P2 row 3a class | HELD | all 8 (AA\|BB): kind "1" only, no gamma, no log of a number |
| P3a row 3b class | HELD | all 10 hybrids inside {1, E1, Ei, ln R}; constants gamma, ln 2, 3, 5, 6, 7, 11 and normalization surds |
| P3b ln/gamma at eq; Ei iff a_d > a_c | MISSED (1 case) | ln R and gamma in exactly the five eq cases (H1, H2, H4, H7, H10), held; Ei in H5 and H8, but not in H9 |
| P4 no ln R at R = 0 | HELD | 48/48 |
| P5 singular points | HELD | R = 0 Fuchsian 48/48; every finite nonzero singular point apparent 48/48; rates at infinity among the atom rates 48/48 |
| P6 order bound (fitted) | HELD | rows 2 and 3a: order = number of distinct rates, 38/38 (orders 1-3); row 3b: orders 3-5, at most rates + 1, 10/10 |

**Gate: GO** (P2, P3a, P4, P5 held). In the two-centre layer, p and d orbitals through n = 3 do not climb: no new
function type, no genuine finite singular point in R, no logarithm at R = 0 outside exchange.

**P3b's miss has a known mechanism, which the prediction failed to use.** H9 is (2s 2s | 3d0 3d0_B): the one-centre
pair (ab) is s-type, so its potential is a monopole, V_0 = Q/r - e^{-b r}(poly + c/r), and its product with c_A d_B is
a one-electron integral with a shifted rate, elementary. E1 and Ei enter the hybrid only through the r^-(L+1) tails of
the L > 0 multipoles. `geovac/qfd_core.py`'s docstring already states it ("s-type pair: {exp} only"), and the 1s
hybrids of the v6.7.26 panel are elementary for the same reason. The correct statement: E1 (and Ei when a_d > a_c, ln R
and gamma when a_d = a_c) appear iff the (ab) pair carries a multipole L >= 1. All ten panel hybrids satisfy it.

**Reported:** row 3b order = rates + 1 exactly where E1 atoms appear (H1-H8, H10); H9 (elementary) has order = rates.
Every (AA|BB) has R = 0 exponents including a negative one (the 1/R tail).

**Homing (2026-11-14, after scoring).** Engine promoted to `geovac/two_center_one_electron.py` (byte-identical
output; `debug/row8_one_electron.py` now keeps only the quadrature referee). Guard `tests/test_paper58_row8_l2.py`,
written as a separate pass, fire-tested 11/11. One thing the guard pass found that the scoring did not ask: an
operator whose integral carries E1 has exactly one logarithmic solution at R = 0 (H2 here; HYB33/HYB34 of v6.7.26
too), so "no ln R at R = 0" is a property of the integral, not of its operator. P4 was stated for the integral and
stands. Paper 58 `sec:census` carries the l <= 2 statements.
