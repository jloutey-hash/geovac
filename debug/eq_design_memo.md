# Record EQ design: two-centre exchange beyond EP (3d on A with 1s on B; 2p on both centres)

*Design stage, 2026-11-17. No scored quantity computed: at the scored cases only booleans,
|differences| against the exact reference, valuation flags and timings were recorded. The PM writes
and registers the record.*

## Object

f(R) = (ab|cd), a, c on A and b, d on B, nodeless hydrogenic orbitals (n = l + 1; rate = Z/n), Neumann
sum in prolate coordinates. Six densities were built; four are proposed for scoring:

| case | integral | sigma, H | density polynomial (with xi^2 - eta^2) | prefactor x R^9 |
|:--|:--|:--|:--|:--|
| **3dsig** | (3d0_A 1s_B \| 3d0_A 1s_B) | 0, 0 | (xi^2-eta^2)[2(1+xi eta)^2 - (xi^2-1)(1-eta^2)] | a^7 b^3 / 2304 |
| **3ddel** | (3d+2_A 1s_B \| 3d-2_A 1s_B) | 2, 2 | xi^2 - eta^2 | a^7 b^3 / 1536 |
| **ppsig** | (2p0_A 2p0_B \| 2p0_A 2p0_B) | 0, 0 | (xi^2-eta^2)(xi^2 eta^2 - 1) | a^5 b^5 / 128 |
| **pppi0** | (2p+1_A 2p+1_B \| 2p+1_A 2p+1_B) | 0, 1 | xi^2 - eta^2 | a^5 b^5 / 512 |
| ppdel | (2p+1_A 2p-1_B \| 2p-1_A 2p+1_B) | 2, 2 | xi^2 - eta^2 | a^5 b^5 / 512 |
| 3dpi (spare) | (3d+1_A 1s_B \| 3d-1_A 1s_B) | 1, 1 | (xi^2-eta^2)(1+xi eta) | a^7 b^3 / 384 |

Real (2p_x,A 2p_x,B | 2p_x,A 2p_x,B) exchange is a fixed linear combination of pppi0 and ppdel: its
density rho^2 cos^2 phi = rho^2/2 + rho^2 cos(2 phi)/2 splits into a sigma = 0 piece and a sigma = 2 piece. **ppdel is not scored: it equals 3 (b/a)^2 x 3ddel
identically** (same density polynomial, sigma, H and exponential; only the normalisations differ;
control c1i), so its operator is 3ddel's. 3ddel is the d case with sigma >= 1 and exercises the new
sigma = 2 generator path; 3dpi (sigma = 1, the EP P path with the 2psig polynomial) is built and
validated but left out to hold cost.

Scored pairs (a, b): **(1, 1)** for all four, **(1, 2)** for 3ddel, ppsig and pppi0 — **seven cases**
(a = rate of the A orbital). **3dsig at (1, 2) is dropped on cost** (estimate 10-11 h for ladder +
heldout + numeric, against about 3-4.5 h for the other unequal cases; see Cost); it can be added as an
optional eighth case at that budget, each step well inside 48 h. Reference:
`exchange_hp(*eq_gen.ref_orbitals(case, a, b), R)`.

## Instrument

`debug/eq_gen.py` imports `g_pieces_gen` and helpers from the frozen `ep_gen` and the tables and
Legendre data from the frozen `exch_gen` (neither edited nor monkeypatched); `generate` is ep_gen's,
copied with three changes. (i) Density polynomials and prefactors are **derived** from Cartesian solid
harmonics (zA = 1 + xi eta, zB = xi eta - 1, rho^2 = (xi^2-1)(1-eta^2), units R/2), with
N^2 = (2 zeta)^{2l+3}/((2l+2)! A_lm), Condon-Shortley phases and (-1)^sigma:
pref = (-1)^sigma ph N_a N_b N_c N_d 4 pi^2 / 2^rexp, rexp = 5 + 2(l_A + l_B). A different route from the
reference's `two_center_spheroidal_product`; it reproduces ep_gen's typed c, rexp and pref exactly (C0).
(ii) Any sigma: Q0^(i) = N_i/(y^2-1)^i (N_1 = -1, N_{i+1} = N_i'(y^2-1) - 2iyN_i) gives
C_j = y^j[(y^2-1)^H V^(sigma) - sum_i C(sigma,i) P^(sigma-i) N_i (y^2-1)^{H-i}], polynomial since
H = (|m_a|+|m_b|+sigma)/2 >= sigma. (iii) eta moments for any (sigma, H) by parts. Truncation: the same
proved per-tau valuation bound, here 2 tau - 4 for all four cases (226-228 tau terms at N = 450), checked
on every term.

`debug/eq_run.py` is `ep_run.py` (frozen, copied) with: generator eq_gen, outputs `eq_*`, reference
orbitals from `eq_gen.ref_orbitals`, the exch_run degree cap rebound **DMAX 60 -> 100** (d_box(r) =
min(DMAX, (E(N-HOLD-r)-MARGIN)//(r+1) - 1) still applies: at E = 3, d <= 100 for r <= 10, 94 at r = 11,
80 at r = 13, 74 at r = 14), and a `gate` reproducing three frozen operators (1s EXS, 2ppi (1,1) EPS,
2psig (1,2) EPS) on both scored primes. Ladder, k0 rule, CAP 96, heldout and numeric are EP's. Why
DMAX 100: at DMAX 60 a degree above 60 is reached only through the gcrd of a higher-order hit (C5d),
which runs out near r = 14; the margin is validated by C4/G0 (known operators unchanged) and C5c/C5d
(a planted order-9, degree-70 operator is the first hit and equals the DMAX-60 gcrd). Margin, not necessity. Prime budget: that plant first
reconstructs at 12 primes (EP's real cases: k0 4-8); CAP 96 leaves wide room.

## Derivation (before guessing)

1. **Constant class.** Unchanged from row 3c and EP: sigma = 2 adds Q0' and Q0'' terms, rational in y
   with poles (y^2-1)^{-i}, i <= 2, cleared by (y^2-1)^H, H = 2; so every tau-term is again a
   Laurent polynomial times {1, gamma + ln p, ln 2, e^{2p}E1(2p), e^{4p}E1(4p)} x e^{-2p} x beta(q)^2.
   Prior art (Lesiuk-Moszynski); not scored.
2. **Equal exponents: order <= 4 for all four cases, theorem-level.** At q = 0, beta_{tau,k}(0) is a
   moment of eta^k (1-eta^2)^{H-sigma} x (1-eta^2)^sigma D^sigma P_tau and vanishes for
   tau > k_max + 2H - sigma = 4 in every case (3dsig, ppsig: k_max 4, H 0; 3ddel: 2 + 4 - 2;
   pppi0: 2 + 2). Measured tau counts 5, 3, 5, 5 (c3). So f lies in the rank-4 C(R)-module spanned by
   e^{-2p}, e^{-2p}(gamma + ln p), E1(2p), e^{2p}E1(4p), closed under d/dR: order <= 4 (exactly 4 absent
   an accidental relation; EP: 4). Rates in {0, 2a, -2a}.
3. **The delta identity.** ppdel = 3 (b/a)^2 3ddel as functions of R (proved from the densities;
   c1i confirms it on the series). So any operator statement for 3ddel holds for the sigma = 2 part of
   real p_x p_x exchange.
4. **A <-> B symmetry of the PP family.** Reflection z -> R - z swaps the centres and maps
   2p0_A(a) -> -2p0_B(a), 2p+-1 to themselves: f_PP(a, b) = f_PP(b, a) for ppsig, pppi0, ppdel
   (C7; D is not, negative control). So PP at (2,1) is PP at (1,2), and its rates are symmetric in a, b.
5. **Unequal exponents: no proof.** The tau-sum is infinite; the exponential skeleton is the 1s one,
   so the heuristic is EP's: order 9, the 1s rate multiset, higher degree. EP held this exactly.
6. **R = 0.** Onset heuristic (EX4b, held in all four EP cases): Sum(l_i + 1) - 1 + 6 =
   **R^13 ln R** for all four cases (3+1+3+1 = 2+2+2+2 = 8). Scaling f_{la,lb}(R) = l f_{a,b}(lR)
   fixes the coefficient's total degree at 14; the EP split (normalisation x a^2 b^2) gives a^9 b^5
   for D (ratio (1,2)/(1,1) = 2^5 = 32) and a^7 b^7 for PP (ratio 2^7 = 128; symmetric, as item 4
   requires). Weaker than in EP: ppsig's (xi^2 eta^2 - 1) or pppi0's H = 1 could cancel the leading term.
7. **Degree.** No bound. EP's degree rose with density size (1s 8/16, 2ppi 9/25, 2psig 16/34). Each new
   apparent point adds about the order (EP C5b +18 for +-i; (1+R^2)^3 also +18, C5c'; five real points
   +45, C5c), so degrees above 60 at (1,2) are not excluded: hence the DMAX 100 margin.

## Controls (`debug/data/eq_controls.json`, `eq_controls_c*.json`, `eq_gate_EQG0.json`, `eq_*_EQT.json`)

| id | what | result |
|:--|:--|:--|
| C0 | eq_gen vs frozen exch_gen (1s) and ep_gen (2psig, 2ppi), mod 67108859 to R^450 at (1,1), (1,2); Q to R^60 at (1,2) for 1s, 2ppi; builder vs ep_gen's typed c, rexp, (sigma, H), pref | identical, bit for bit, all ten; builder equal on every field |
| C1 | (1,1): Q series N = 90 vs `exchange_hp` (dps 50), R = 1/2, 6/5, all six densities | every difference 2e-65 to 1.7e-62 (3ddel, ppdel: <= 5e-64); sign, prefactor and the sigma = 2 path confirmed |
| C1f | (1,2), (2,1): Q series N = 60 at R = 1/2 vs reference at tau 16 (tau 12 vs 16 differ by <= 1.5e-34, 3dsig (2,1)) | 3dsig 1.2e-51, 1.8e-50; 3ddel 2.0e-52, 3.1e-51; ppsig 1.0e-50 (both); pppi0 1.9e-51 (both); all pass 1e-35 |
| C1i | ppdel = 3(b/a)^2 3ddel, series mod p, N = 200, at (1,1), (1,2), (2,1) | true at all three |
| C7 | PP symmetry f(1,2) = f(2,1), mod p, N = 200; 3dsig as negative control | true for ppsig, pppi0, ppdel; false for 3dsig |
| C2 | scored cases: Q series reduced mod 33554393 = mod-p series (N = 60); two mod-p runs identical; valuation | all true, eight cases (3dsig (1,2) included) |
| C3 | N = 450 mod 67108859: valuation flag, tau count, wall time (series discarded unopened) | all pass; tau 5, 3, 5, 5 at (1,1) (3dsig, 3ddel, ppsig, pppi0); 228, 226, 228, 228 at (1,2); times in Cost |
| C4 | positive control: frozen `exch_run.readings` on eq_gen 2ppi (1,1) and 2psig (1,2), DMAX 60 and 100 | order 4 degree 9 and order 9 degree 34, first hits (4,9), (9,34) under both boxes; coefficients equal the frozen EPS operators mod p |
| C5 | normalised heights of the certified EP operators | 2psig (1,1) 29/10 bits, 3 primes; (1,2) 92/60, 8; 2ppi (1,1) 16/6, 2; (1,2) 66/34, 6 (EP used 8-12 accepted primes) |
| C5c | widened box: planted (1+R)(2+R)(3+R)(4+R)(5+R) f_2ppi(1,2) (known order 9, five new apparent points), DMAX 60 vs 100, then a CRT ladder under DMAX 100 | DMAX 100: first hit (9, 70), order 9 degree 70 on every prime; reconstruction fails at k = 2..11 primes, first succeeds at k0 = 12 (151/35 bits), same operator at 13, 14 |
| C5d | same planted series, one prime: minimal operator under DMAX 60 vs 100 | DMAX 60 first hit (10, 38), nullity 2, its gcrd = order 9 degree 70; DMAX 100 first hit (9, 70); coefficients equal |
| G0 | `eq_run gate` (suffix EQG0): 1s (1,2) vs EXS, 2ppi (1,1) and 2psig (1,2) vs EPS, both scored primes, DMAX 100 | every field equal, all six runs; valuation passes |
| C6 | `eq_run` ladder, heldout, numeric end to end on 2ppi (1,1), suffix EQT | 8 primes accepted, k0 = 4, op(4) = op(8); rates {0, 2, -2 x2} on both primes; order 4 degree 9; exact Q residual zero; third prime nullity 1 and matches; numeric_pass (4.5e-93 at R = 1.4, 2.3e-90 at 3.015) |

**Box reading (C5c, C5d).** The frozen box already reaches the planted degree-70 operator through the gcrd of an
r = 10 hit; DMAX 100 only makes it the first hit. The PM may keep DMAX 60 at no loss on the controls.
(A first plant, (1+R^2)^3 f_2ppi(1,2), reached only degree 43 and was replaced: C5c'.)

## Cost (measured under 15-16 concurrent jobs on 16 cores; serial runs are faster)

Generator per prime at N = 450, (1,2): 3dsig 1013 s, 3ddel 581 s, ppsig 478 s, pppi0 298 s (EP: S 401-488 s,
P 175 s). Equal pairs: 5 s. Exact-Q generator at (1,2), N = 20/30/40: 3dsig 18/44/116 s (about N^2.7),
3ddel 5/12/24, ppsig 9/22/47, pppi0 5/12/24. 3dsig (1,2) re-timed at 708 s with 7 jobs. Reference `exchange_hp` at (1,2), R = 7/5, tau 16,
dps 80: 2psig 651 s, 3dsig 794, 3ddel 189, ppsig 273, pppi0 175 s. EP's actual steps for S (1,2):
ladder 1.4 h (12 primes), heldout 0.5 h, numeric 4.6 h (reference tau 8 -> 16 -> 32 at two R).

Estimates (12-16 primes; the planted degree-70 control needed 12, EP's (1,2) ladders accepted 10-12, so budget up to ~20; heldout N_Q up to ~180; numeric scaled from EP S by the reference ratio):

| case (1,2) | ladder | heldout | numeric | total |
|:--|:--|:--|:--|:--|
| pppi0 | 1-1.3 h | 0.3 h | 1.3 h | ~3 h |
| 3ddel | 1.5-2 h | 0.4 h | 1.4 h | ~3.5 h |
| ppsig | 1.5-1.8 h | 0.5 h | 2.0 h | ~4.3 h |
| 3dsig (dropped) | 2.5-3.2 h | 1-1.5 h | 5.8 h | ~10-11 h |

Equal pairs: minutes each. Seven cases: about 11 h CPU, about 4.5 h wall with the three unequal cases in
parallel; every step far inside the 48 h clock. 3dsig (1,2) breaks the 6 h rule on its numeric step alone;
add it (about 11 h) only if the sigma = 0 d integral at unequal exponents is wanted.

## Proposed clauses (PM to set; probabilities mine)

- **EQ0 (gate).** `eq_run gate 1s --suffix EQS` reproduces, under DMAX 100 and on both scored primes, the
  frozen 1s (1,2) EXS operator, the 2ppi (1,1) EPS operator and the 2psig (1,2) EPS operator in order,
  degree, coefficient degrees and coefficients mod p; valuation passes in every output; frozen files
  hash-match. Failure VOID. (sigma = 2 is covered by C1/C1f and EQ1's numeric bar, not by the gate.)
- **EQ1 (existence, all seven cases).** As EP1: ladder accepts on both scored primes with equal shape,
  k0 <= 92 and op(k0) = op(k0 + 4); exact Q residual zero; third prime nullity 1 and matches;
  `numeric_pass`. Instrument parts failing after an operator is found: VOID. Probability 0.75.
- **EQ2a (equal exponents).** Order exactly 4 for all four cases at (1,1) (<= 4 derived; above 4 is VOID).
  0.75.
- **EQ2b (unequal).** Order in [9, 13] at (1,2) for 3ddel, ppsig, pppi0 (0.8); exactly 9 reported (0.5).
- **EQ3 (no climb).** Every finite nonzero point `apparent` or `exponents_only: apparent-type`, both primes,
  all seven. 0.7.
- **EQ4 (R = 0).** Fuchsian, a log solution, f has a nonzero ln R component, all seven. 0.8.
- **EQ5 (rates).** Rank <= 1; every rate rational, equal on both primes, in {+-(a+b), +-2a, +-2b, 0, +-(a-b)};
  equal-exponent set {0, 2, -2}. 0.6.

**Reported, not gating (probabilities in advance):** exactly 9 at (1,2) in all three (0.5); the 1s rate
multiset {-3, 4, 2, 0^2, -2^2, -4^2} at (1,2) (0.45); no a - b rate (0.55); onset R^13 ln R in all seven
(0.35; per family: D 0.5, PP 0.4); coefficient ratio (1,2)/(1,1) = 32 for 3ddel (0.4) and 128
for ppsig and pppi0 (0.3); degree above the EP degree of the same sigma at the same pair (3dsig > 16 at (1,1),
3ddel and pppi0 > 9 / 25, ppsig > 16 / 34; 0.65); degree above 60 in at least one (1,2) case (0.35).
Joint GO probability, gating clauses treated as independent: about 0.15-0.25.

## Gate

- **VOID** (overrides): EQ0 fails; scored primes disagree in shape; an EQ1 instrument part fails after an
  operator is accepted; reconstruction never stabilises within CAP 96; order above 4 at (1,1) (contradicts
  the proved bound); a step unfinished 48 h after launch.
- **STOP**: for some case no operator in the box (r <= 14, d <= d_box(r), DMAX 100) on both scored primes
  (reading: "not in the box"), or a finite nonzero point `genuine` on both scored primes (a climb).
- **GO**: EQ1-EQ5 hold in all seven cases. **BORDERLINE** otherwise (including order below
  4 at (1,1)).

Scoring order (cheapest first, so an instrument failure shows early): gate; pppi0 (1,1), 3ddel (1,1),
ppsig (1,1), 3dsig (1,1); pppi0 (1,2), 3ddel (1,2), ppsig (1,2); each `ladder`, `heldout`,
`numeric`, suffix EQS.

## Caps

Nodeless orbitals only (n = l + 1); l <= 2 on A with 1s on B, l = 1 on both; sigma <= 2; two pairs,
a = 1 at both (so the power of a in the onset coefficient is untested, as in EP; a (2,1) pair would test it
for the D family at the cost of four more unequal cases); PP's (2,1) equals (1,2) by symmetry; ppdel and
3dpi not scored (identity; cost); 3dsig only at equal exponents (cost); mod-p singular readings; gcrd not certified minimal; constants formal.

Files: `debug/eq_gen.py`, `debug/eq_controls.py`, `debug/eq_run.py`, `debug/eq_c5d.py`; outputs
`debug/data/eq_*`.
