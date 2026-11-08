# P3e-2W to P3e-5: results (scored against docs/predictions/2026-10-19_p3e2_widened.md, commit bb3ffe11)

Run 2026-10-19/20. Nothing tracked was edited; nothing committed. Scripts `debug/p3e2w_*.py`,
data `debug/data/p3e2w_*.json|log`. Predictions, spans, triples, precisions, the counting rule and
the controls were used as written; every deviation in how they were carried out is listed in Sec. 6.

## 1. Verdicts

| ID | Verdict | Class (README) | One line |
|:--|:--|:--|:--|
| P3e-3 | **HELD** | CHECK | G(t;1,1) is annihilated by W = R3 R1 (theta+2); R3 is exactly equivalent to F_3's operator; R1 and theta+2 are elementary; no split (planted split found); singular at t = 1. |
| P3e-2W (a) | **HELD where testable** (19 of 22 cells; SE void at T2, T3, T4) | CHECK | No relation for G meets the counting rule at T1-T4 in any span whose controls pass. |
| P3e-2W (b) | **HELD** | CHECK | Swap at {1,3,7}: no relation among G(1;3,7), G(3;1,7), the reduced K3 set and T3. |
| P3e-4 | **FAILED** at S1 and at S3 (scored separately) | MISS, MISS | No relation meets the counting rule in any threshold span with passing controls (SB at S1 void). |
| P3e-5 | **HELD** | CHECK | Every Gctrl reduction, at all seven triples, uses degree-3 elements only; no G relation exists to test. Fallback rerun not triggered. |

If the PM prefers to score a prediction with untestable cells as PARTLY, P3e-2W(a) is the only
candidate: the falsifier did not fire anywhere, but SE could not be tested at three of the four
generic triples (Sec. 4).

## 2. P3e-3: the slice G(t;1,1), Slater weights (run first)

Object: G(t;1,1) = int_0^inf x K1(tx) J(x)^2 dx, J = 3L - 3xK0 - x^2 K1, t on the central density.

- **Operator (Mellin route, density_loop_B_ode code imported unmodified).** theta J = x^3 K2 gives
  theta^3 J = (z-5) theta J + 6 theta^2 J; symmetric square; Phi(s) = A(2-s) E(s), E = M[K1],
  E(s+2) = (s^2-1)E(s). Self-check reproduces 4 theta^3 - t^2(theta+1)^3. Result L_G, order 7,
  t-degree 6, leading coefficient -(t-1)(t+1)(t-2)(t+2) (`p3e2w_slice_A.json`).
  L_G G at t = 0.5, 0.9, 1.3, 2.5, 3.5: <= 1.8e-18 at 20 digits, <= 8.6e-31 at 32 digits
  (no source term; unlike the triangle, nothing is invisible to the Mellin transform here).
- **Exact factorization.** L_G = N2 R3 R1 (theta+2), every division with zero remainder:
  - theta+2: common right factor, solution t^-2.
  - R1 = (1-t^2) theta - (1+7t^2), solution t/(1-t^2)^4: the unique algebraic solution of
    L_G/(theta+2), found by an exponent-bounded ansatz (control: the chain's R1 solution
    t^3/(1-t^2)^2 recovered by the same code).
  - R3 (order 3, monic): [60(4t^4+22t^2+7), 4(46t^4+232t^2+91), 12(3t^4+11t^2-7)]/((t^2-4)(2t^2+7)),
    apparent singularity 2t^2+7. The rational map
    A = a0 + a1 theta + a2 theta^2,
    a0 = (4t^12-82t^10+745t^8-5924t^6-3328t^4+4416t^2-4096)/(t^3 (t^2-4)^5),
    a1 = (8t^12-168t^10+1621t^8-14416t^6-18384t^4+16064t^2-10240)/(t^3 (t^2-4)^5),
    a2 = 2(2t^14-47t^12+517t^10-5012t^8-16400t^6+16160t^4-21504t^2+16384)/(t^5 (t^2-4)^5)
    sends every solution of L3 = 4theta^3 - t^2(theta+1)^3 to a solution of M5 = L_G/(R1 (theta+2));
    M5 o A = 0 modulo L3 checked exactly (no series). R3 = annihilator of A(Sol L3), so R3 is
    equivalent to F_3's operator.
  - N2 (order 2, irregular at infinity): [(24t^8+12t^6-398t^4-1211t^2-980)/(2t^2+7)^2,
    6(3t^2+28)/(2t^2+7), (t^2-4)].
  - mod-p screens before the exact step: f3 found as a sub of M5 only once the box reached
    (t^2-4)^5, numerator degree 24 (none at the D5 default box); f3 never a quotient; the elliptic
    families 2F1(1/2,1/2;1;t^2) and 2F1(1/2,1/2;1;t^2/4) neither sub nor quotient of M6 or M5.
- **N2 is not needed.** W = R3 R1 (theta+2) applied to G (under the integral):

  | t | G | (den W) G, 20 dps | (den W) G, 32 dps | control (one R3 coeff x 199/198) | R1(theta+2) G |
  |:--|:--|:--|:--|:--|:--|
  | 0.5 | 63.71736087 | -1.96e-19 | -1.16e-31 | -1839 | 191.19 |
  | 0.9 | 11.71310974 | 7.5e-20 | 3.8e-32 | 1100 | 88.59 |
  | 1.3 | 3.646420117 | 1.8e-19 | 7.9e-32 | 4232 | 51.07 |
  | 2.5 | 0.3771185738 | -2.0e-19 | -1.1e-31 | 18151 | 16.70 |
  | 3.5 | 0.1084941668 | 2.7e-19 | 1.6e-31 | 35498 | 8.898 |

  R1(theta+2)G is nonzero, so G does carry the K3 piece. The composition factors of the annihilator
  W are {theta+2, R1, R3}: two elementary, one equivalent to F_3's (irreducible) operator. No
  elliptic factor and nothing else.
- **No split.** A split (K3 also a sub of W) needs B = b0 + b1 theta + b2 theta^2 with W(Bf) = 0 on
  Sol(L3). Searched once, mod p, at the maximal box t^-12 (t^2-1)^-6 (t^2-4)^-6 (2t^2+7)^-3 and
  numerators of degree <= 40 (159 unknowns, 644 rows). Exponent bounds from leading terms:
  pole <= 5 at t = 0, +-1, +-2, <= 2 at the apparent points, rational degree <= -1 at infinity, all
  inside the box. W: nullity 0. Planted split E2 o R3 (E2 = R1 (theta+2)): nullity 1, found by the
  same call. Chain W (density loop D12): nullity 0 (consistent with the record).
  Caveat: the exponent bound is the leading-term argument, not a written proof; the search is mod p.
- **t = 1.** W's leading coefficient carries (1 - t^2); W's exponents at t = 1 are {-3,0,1,2,3}.
  The -3 comes from the R1 piece, a genuine pole, so the singularity is not apparent.

Result: G(t;1,1) = c t^-2 + (an elementary solution of R1(theta+2)) + (a K3 function, image of
Sol(L3) under A, integrated twice against elementary kernels). Same pattern as Paper 59 sec:density_ring's Delta^-3 chain.

## 3. Grids and C0

- Grid: x = exp(t - e^{-t}), coarse h = 1/128, fine h/2, one pass, working precision 520,
  x_max = ((520+10) ln 10 + 50)/min mass. 3652 nodes (3372 at T4); 8-22 s per triple.
  max |h - h/2| over all 18 integrands: T1 6.4e-519, T2 6.9e-519, T3 2.6e-520, T4 4.5e-520,
  S1 5.0e-520, S3 2.3e-520, all < 1e-480 (the D = 460 criterion, which implies 340 and 400).
  T1 at 560 vs 520: agree to ~2e-519. **T1 vs debug/data/three_electron_p3e_values_320.json:**
  every key agrees to <= 6e-320 absolute (the cache's own precision), so the 150-digit control
  passes. Bailey check |V - S_KK| <= 1.5e-520 at the generic triples (7e-520 at S1).
- **C0 control (1;1,1):** PSLQ of {G, c30, c32, pi^3, pi^2, pi, 1} gives (216, -80, 6600, -81, 0, 0, 0)
  at 100 digits (cap 1e11) and at 150 digits (cap 1e18). |h - h/2| = 1.7e-139 and 4.5e-189.
  **PASS.**

## 4. Controls by triple and span (pre-step reductions first, without G)

Pre-step removals (counting rule applied to basis-internal relations):
- SA: S_EK and S_EE at every generic triple. Example at T1: 3V - 19D1 + 26D2 + 16S_KE - 16S_EK = 0.
  The Galois-adapted Bailey products are partly derivatives of V.
- SD: Uc at every triple, by Euler homogeneity; e.g. at T1 2U - Ua - 2Ub - 4Uc = 0. At T1, Ub also,
  via 2U - 15Ua + 15Ub = 0 (the geometric progression).
- S1 (b = c): duplicate R, RR and LR entries and kappa EK' (= KE' at k = k').
- Swap: W, D3, S23 and S13, leaving {V, D1, D2, S12} + T3 (rank 4).

| Triple | SA | SB | SC | SD | SE |
|:--|:--|:--|:--|:--|:--|
| T1 | pass | pass | pass | pass | pass |
| T2 | pass | pass | pass | pass | **VOID** |
| T3 | pass | pass | pass | pass | **VOID** |
| T4 | pass | pass | pass | pass | **VOID** |
| S1 | pass | **VOID** | pass | pass | pass |
| S3 | pass | pass | pass | pass | pass |
| Swap | (single span) pass | | | | |

"pass" = structural Gctrl relation identical at 340 and 400, supported on the K3 block (DegB at
S1, S3), and the planted relation recovered exactly at 340 and 400. Gctrl reductions (stable,
also at 460):
- T1: 363825 Gctrl = 11050792 V - 27272112 D1 - 52735824 D2 + 87003024 S12. Denominator 363825 reproduced.
- T2: 8085825 Gctrl = 3324335 V + 37576217 D1 - 32345691 D2 - 156240189 S12.
- T3: 946041525 Gctrl = -1397087125 V + 11915647527 D1 + 3198557411 D2 - 17252396217 S12.
- T4: 29604575 Gctrl = -14049804456 V + 90827118000 D1 + 118528848288 D2 - 423024584544 S12.
- S1: 16 Gctrl = -29 kKK' + 174 kKE' - 174 kEE'.
- S3: 11340 Gctrl = -8615 kKK' + 37330 kKE' + 21961 kEK' - 59291 kEE'.
- Swap: both Gctrl's reduce on {V, D1, D2, S12} (T3's with D1 and D2 swapped, as it must be).

Why the voids:
- **SE at T2-T4.** n = 30, so the cap at 340 digits is 10^(11-3) = 1e8. That is below the planted
  height (1e9) and below the Gctrl heights (1.6e8 to 4e11). Both controls fail at 340 and pass at
  400 (the planted one). This was foreseeable from the file before any G was computed. It was not
  repaired.
- **SB at S1.** The reduced span has 7 elements, so an 8-element planted relation cannot be built.

**Negative control** (G(1;2,4) in T2's SA and SB): no relation meets the counting rule. The vectors
returned at the cap differ between 340 and 400: in SA, none at 340 and a dense 4e17 vector at 400;
in SB, dense vectors of 2e23 and 3e27. **PASS**, and the stability rule is calibrated.

## 5. G outcomes (counting rule)

In every cell, valid or void, PSLQ either returned nothing or returned dense vectors at the cap
that differ between 340 and 400. No relation was identical at 340 and 400, so no run at 460 was
triggered and nothing needed the adaptive 200-digit residual check. Largest coefficient of the
returned vectors (340 | 400):

| | SA | SB | SC | SD | SE |
|:--|:--|:--|:--|:--|:--|
| T1 | none / 4e17 | 1e26 / 6e29 | 8e12 / 4e15 | 1e16 / 3e19 | none / 9e11 |
| T2 | none / 9e16 | 2e23 / 2e27 | none / 3e14 | none / 9e16 | (void) none / none |
| T3 | none / 1e17 | 6e23 / 2e27 | none / 6e14 | none / 6e16 | (void) none / none |
| T4 | none / 1e17 | 2e23 / 1e27 | 9e11 / 6e14 | none / 4e16 | (void) none / none |
| S1 | 2e28 / 2e33 | (void) 2e36 / 2e42 | 8e23 / 7e27 | 4e23 / 1e27 | 2e26 / 1e30 |
| S3 | 4e19 / 1e24 | 2e28 / 1e33 | 8e14 / 5e18 | 2e16 / 8e18 | none / 2e14 |
| Swap | 4e22 / 3e25 (no relation) | | | | |

(Full vectors: `debug/data/p3e2w_pslq_<ID>.json`; every dense hit sits at its cap.)

So:
- P3e-2W(a) holds in all 19 valid generic cells, and (b) holds.
- P3e-4 fails at S1 (SA, SC, SD, SE valid) and at S3 (all five valid). The file's caveat applies:
  a null at S3 does not exclude a collapse into L-values of non-CM forms; at S1 (CM) the
  Chowla-Selberg argument puts them inside the degenerate block, so the S1 null is the more
  informative of the two.
- P3e-5 holds on the controls; there is no G relation for it to apply to.

## 6. Deviations from the pre-registration (how, not what)

1. **Function evaluation.** K0, K1 and L were computed by a fixed-point power series (gmpy2
   integers, guard bits for the e^{2X} cancellation; L by termwise integration of K0's series) and
   by the Hankel asymptotic series at large argument (L = pi/2 - Ki, with Ki's own asymptotic).
   They were not computed by mpmath besselk/struvel. Reason: mpmath's integer-order besselk costs
   1.5-3.5 s per call at 430 digits for 20 < X < 600, which would put a triple at about 25 h, not
   the memo's 1.5 h. Validated against mpmath at 16 arguments from 1e-30 to 900, at 120 and 430
   digits (all digits agree). The T1 sums agree with the existing 320-digit cache to its precision.
   The per-node precision is lowered by the decay of K(ax), which every integrand carries; this is
   an absolute-error budget of 10^-(D+10).
2. **One grid per triple at working precision 520, truncated to 340/400/460,** instead of a grid
   per D. The criterion |h - h/2| < 10^-(D+20) was checked at D = 460, which implies the other two.
   The tail cut used D = 520, which is more conservative. C0 was computed at working precisions
   140 and 190 for the 100- and 150-digit controls.
3. **Order.** The grids and the PSLQ stage ran while the later P3e-3 steps (exact gauge, N2 test,
   split search) were still running. P3e-3 does not use them. No G value was inspected before the
   controls, apart from the pre-registered T1 grid comparison.
4. **Interpretations, fixed before any G PSLQ:**
   - Lambda's residue moduli are k_ab and k_ac, the R-block curves, as in the triple table for T2
     and T4. The T1 Pr and Lambda reproduce the file exactly.
   - Bailey moduli: c = max mass and p <= q the other two, k+- = [sqrt((c+p)^2-q^2) +- sqrt((c-p)^2-q^2)]/(2c),
     kappa = pi/(2c). V = S_KK was verified numerically at all six triples, and the moduli match the
     file's stated values.
   - sqrt(d) is the squarefree part of (c-p)^2 - q^2 (5 at T1, 3 at T2 and T3); at T4, 2Re and 2Im.
   - The swap's "first and second derivatives" are {D1, D2, D3, S12, S13, S23, W}, with
     W = int x^2 K0K0K0 carrying the pure second derivatives.
5. **Rules the file left implicit:**
   - Internal relations are removed only if they meet the full counting rule; the latest nonzero
     element in span order is removed.
   - The planted RNG is random.Random(20261019), seeded fresh for each span: 8 indices sampled, each
     k_i = +-randint(1, 1e9).
   - The structural control passes on a stable 340/400 relation supported on K3 (or DegB); the 460
     run is recorded too.
6. **P3e-3 tooling.** The density_loop_D* modules were imported unmodified, with new drivers.
   - The no-split search used one maximal-box mod-p call per operator instead of D12's grid of
     exact calls. A maximal box contains every smaller one.
   - Step H's first back end (sympy factor on huge rational constants) stalled and was replaced by
     H2: same gauge vector, cancel-only arithmetic.

## 7. Reading

- The Slater-weight slice behaves like the Delta^-3 chain: its only non-elementary piece is the K3
  family of F_3, it does not split off, and the extra order-2 piece of the Mellin operator is unused.
  The transfer the scoping memo flagged as assumed (risk 6) is shown on this slice, which supports
  reading the generic null as "a K3 period integrated once".
- The widened bases catch every control and no G at four generic triples and both thresholds; the
  threshold collapse (the one possible finding) did not appear.
- Process note: with n = 30 the cap 10^(floor(340/n)-3) = 1e8 is below a 1e9 planted height, so any
  span over ~26 entries voids itself at 340 digits.
