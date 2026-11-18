# Pre-registered: the contact amplitude of (AA|BC) at unequal arms, and a Yukawa line that leaves it unchanged (record L52)

**Status:** registered 2026-11-17, before any scored computation (PI: do the search recommendations; the 5/2 test dissolved in design, so the PM registered this amplitude record in its place). Cold review: REGISTRABLE after fixes; all must-fix and should-fix items applied. Design: `debug/l52_design_memo.md`. Driver: `debug/l52_run.py` (gate fire test 23/23).
**Scope of a GO, stated first.** Function level, not operator level: the clauses are about the amplitudes of
f = (AA|BC) and of f_kappa - f at the contact x*, never about which exponents an annihilating operator has (AB
P3a's question, open). P1 tests the **transfer** of an already twice-confirmed mechanism (the Landau-pinch segment
average: (XY|XZ) water records, matched to 3.2e-5; AE2b's equal-arm amplitudes, 0.04%) to (AA|BC) at unequal arms; it
is not a new mechanism. The kappa-independence of v3 behind P2 follows from (-lap + kappa^2) V = 4 pi rho_A alone
(one line, step 1), so a P2a GO tests the Yukawa assembly and the segment formula at order 9/2, not a deep claim.
**Owner:** Paper 59 `sec:scope` ((AA|BC) paragraph; the Landau-pinch average paragraph for (XY|XZ)); inventory row 6.
**Rests on:** record AB (P3b HELD: no 3/2 or 5/2 branch, leading 7/2), AE2b (equal-arm amplitudes), record
`2026-10-28_water_landau_gaussian_fresh.md` (segment average for (XY|XZ)).
**Guardrail check (Papers 8-9, FCI-M, Track DF):** not in scope.

## Object

(AA|BC) = int rho_A(r1) phi_B(r2) phi_C(r2) / r12 (normalized 1s), AB's frame: B origin, C = (0,0,D1),
A = D2 (sin g, 0, cos g), x = cos g, f(x) = sum_L F_L P_L(x). Contact x* = (D1^2+D2^2)/(2 D1 D2), where
|A - C|^2 = 2 D1 D2 (x* - x) vanishes. Second object: the same integral with the Yukawa line e^{-kappa r12}/r12,
f_kappa, and the difference D = f_kappa - f.

Reduction (engine `debug/aabc_exact.py`, unchanged): f = NA - Y_mu - zA W_mu (AB), and
f_kappa = mu^4 a (Y_kappa - Y_mu) - mu^3 / (2 (mu^2 - kappa^2)) W_mu, a = (mu^2 - kappa^2)^-2, mu = 2 zA
(Fourier: V_kappa^ = 4 pi mu^4 / ((k^2 + kappa^2)(k^2 + mu^2)^2), partial fractions). V_kappa checked symbolically:
it solves (-lap + kappa^2) V = 4 pi rho_A (residual 0) and tends to AB's V_A at kappa -> 0 (control K1).

## Derivation (before computing)

1. **Singular terms.** Near A, V(s) = E(s^2) + sum_{p odd} v_p s^p. From (-lap + kappa^2) V = 4 pi rho_A with rho_A
   bounded: v_{-1} = 0 (no point charge), v_1 = 0 (would need rho ~ 1/s), v_3 = -4 pi rho_1 / 12 = 2 zA^4/3 for every
   kappa; v_5 = zA^4 (kappa^2 + 8 zA^2)/45. Near C, rho_BC = smooth + sum_{q odd} N (-zC)^q / q! e^{-zB r_B} |r - C|^q.
2. **Pinch (Feynman parameters, 3-D).** int |s|^p |s - a|^q g(A + s) d^3s has non-analytic part
   C_pq |a|^{p+q+3} [ int_0^1 u^{(p+1)/2} (1-u)^{(q+1)/2} g(A + u a) du + O(|a|^2) ],
   C_pq = pi^{3/2} Gamma(-(p+q+3)/2) / (Gamma(-p/2) Gamma(-q/2)). At x*, |A + u(C-A) - B|^2 = (1-u) D2^2 + u D1^2.
3. **f.** Lowest pair (3,1): f = analytic + Phi7 (x* - x)^{7/2} + O((x* - x)^{9/2}),
   Phi7 = v3 N (-zC) C_31 (2 D1 D2)^{7/2} int_0^1 u^2 (1-u) e^{-zB sqrt((1-u) D2^2 + u D1^2)} du,  C_31 = -2 pi/35.
   (The 3/2 and 5/2 coefficients carry v_{-1} or v_1: zero identically.)
4. **D.** v3 is kappa-free, so every term built on v3 (the (3,1) pair to all orders, the (3,3) pair) cancels in D;
   the lowest surviving pair is (5,1) with Delta v5 = kappa^2 zA^4/45:
   D = analytic + DPhi9 (x* - x)^{9/2} + O((x* - x)^{11/2}), DPhi9 = Delta v5 N (-zC) C_51 (2 D1 D2)^{9/2}
   int_0^1 u^3 (1-u) e^{-zB R_B(u)} du, C_51 = -2 pi/63.
5. **Legendre side.** T^gamma_L = Legendre coefficients of (x* - x)^gamma (Rodrigues, positive integrand;
   `l52_local.T`, checked against direct quadrature and AE2b's pred_F to 1e-15). Prediction:
   R_L = F_L / (Phi7 T^{7/2}_L) -> 1 and R^D_L = D_L / (DPhi9 T^{9/2}_L) -> 1, read by a four-term fit
   r0 + c1/L + c2/L^2 + c3/L^3 on L = 120, 140, ..., 240.

## Cases (A1, A2: AB's geometries, amplitude never computed; A3 fresh, no F_L computed at registration)

| Case | D1 | D2 | zB | zC | zA | kappa | x* | Phi7 (predicted) | DPhi9 (predicted) |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| A1 | 9/5 | 13/5 | 13/10 | 1 | 9/10 | 7/10 | 125/117 | 0.48406711402945841531 | 0.02636189906870907033 |
| A2 | 3/2 | 12/5 | 1 | 6/5 | 4/5 | -- | 89/80 | 0.31156574196371525736 | -- |
| A3 | 2 | 3 | 11/10 | 9/10 | 6/5 | 1 | 13/12 | 2.4426119167955476068 | 0.34932701023832986156 |

(`debug/l52_predict.py` -> `debug/data/l52_predictions.json`.) Toy decay values beta_eff(220): 7/2 -4.119 (A1),
-4.094 (A2), -4.108 (A3); 9/2 -5.180 (A1), -5.164 (A3).

## Controls (measured; `debug/data/l52_controls_*.json`, `debug/data/l52_ctrl_*.json`)

| Control | what | measured |
|:--|:--|:--|
| K1 sym | odd coefficients of V (Coulomb, Yukawa); V_kappa against its Helmholtz equation | v_-1 = v_1 = 0 identically; v3 = 2 zeta^4/3 both; Yukawa(kappa->0) - Coulomb = 0; Helmholtz residual 0 |
| K2 const | C_pq * Beta((p+3)/2, (q+3)/2) vs the corpus constant K(q,p) of `aaeq2_local.py` | equal to 20 digits for (-1,1), (3,1), (5,1), (3,3), (-1,3) |
| K3 equal arms (AE2b's closed form, E3) | segment formula at D1 = D2, 1s 7/2 amplitude | 1.3356449576164494005e-6 vs AE2b 1.3356449576164494005e-6 (rel 4e-20); amp * T_L vs AE2b pred_F at L = 220, 400, 580: rel <= 1e-15 |
| K4 NA (record NA's 3/2), exact F_L from `na_exact`, 2000 digits, L <= 240 | R_L = F_L / (Phi3 T^{3/2}_L), Phi3 = segment formula p = -1, weight (1-u) | N1: R_240 = 0.99616, extrapolated 0.99999987; N2: 0.99508, 0.99999991 |
| K5 instrument, 7/2 toy (x*-x)^{7/2} + 3(..)^{9/2} - 5(..)^{11/2}, A1's x* (`l52_controls.py toy`) | extrapolated ratio | 1.000012 (R_240 = 0.977) |
| K6 instrument, 9/2 toys with nuisances (5, -10) and (-8, 20), A1 and A3 x* (`l52_controls.py toy92`) | extrapolated ratio | 1.000053, 0.999895, 1.000056, 0.999890 (R_240 0.948..1.084); a three-term fit gave up to 1.7e-3 and was replaced before any scored data existed |
| K7 Y at mass kappa, per L (`l52_run.py ctrl_perL`; Y_kappa alone, D not formed), A1 | engine (1200 digits) vs angular + Bessel quadrature (50 digits) | L = 0, 1, 4, 12: 1.3e-51, 2.7e-51, 1.3e-51, 4.7e-45 |
| K9 assembled f_kappa per L (`l52_run.py ctrl_fk`), A1, 30-digit quadrature of closed-form V_kappa projected on P_L about B, against rho_L by quadrature | engine (Y_kappa, Y_mu, dY_mu combined as in Phase 2) vs quadrature | L = 0: 0.03815115818873339393958243, rel 2.0e-31; L = 1: 0.02019853900114325788620231, 2.0e-31; L = 4: 4.705346535503054890364366e-4, 3.9e-31; L = 12: 1.888270464650091236791601e-7, 1.4e-30. 12 to 53 s per L |
| K8 gate fire test | `l52_run.py fire` | 23/23 (first version 11/12: a missing case read BORDERLINE, fixed; the P0 extensions add 10 plants, all VOID as intended) |

Disclosed contact with scored objects: K9 forms f_kappa at A1 for L <= 12 only (scored L are 120..240); K7 forms
Y_kappa alone. No amplitude of f or D has been evaluated; nothing at A3 beyond the predicted numbers. AB's
beta_f(220) at A1, A2 (-4.106, -4.083, measured) sits near the 7/2 toys; the exponent is not scored again.

## Predictions

- **P0 (instrument; VOID if any part fails).**
  (i) AB's stored values from L <= 240: A1 at x = 1/2, 0.21286449876124374093484119906, to 1e-25; A2 at x = -1/3,
  0.223136849487580864384315474643, to 1e-20.
  (ii) F_L and D_L at 2000 and 2600 digits agree to 1e-12 relative at L = 120..240.
  (iii) Y at mass kappa per L, L = 0, 1, 4, 12, against the independent quadrature, to 1e-25, at A1 and A3.
  (iv) Assembled f_kappa per L, L = 0, 1, 4, 12, against the direct quadrature of closed-form V_kappa (K9's route,
  30 digits), to 1e-20, at A1 and A3. A P2a STOP therefore cannot be an assembly or engine error.
  (v) Coulomb beta_eff(220) within 0.1 of the 7/2 toy at every case, A3's (-4.108) required: a singularity nearer
  than x* at A3 reads VOID, not STOP.
  Probability 0.92 (A3's (iii)-(v) not yet measured; A1's measured).
- **P1 (7/2 amplitude, unequal arms).** |r0 - 1| <= 1e-3 at A1, A2 and A3. K4 puts the one-electron formula at 1e-7
  and K5 the instrument at 1e-5; the residual risk is two-electron subleading terms too large for a four-term fit (K6
  bounds that at 1e-4 for nuisances up to 20). Probability 0.88.
- **P2a (Yukawa line: same 7/2 amplitude; difference leads with 9/2).** |r0^D - 1| <= 3e-3 at A1 (kappa = 7/10)
  and A3 (kappa = 1). With assembly errors moved to VOID by P0(iv), the residual risk is a 9/2 piece missed in step 4.
  Probability 0.80 unconditionally; 0.85 given P1 (same segment formula and instrument).
- **P2b (reported only).** beta_eff(220) of D within 0.1 of the 9/2 toy. (The f half was dropped: AB already
  measured it.) Reported at scoring with R_240 and c1 for f and for D.

## Decision gate

VOID if P0 fails. Otherwise GO if P1 and P2a hold at every case; BORDERLINE if all five ratios are present and each
is within 5x its tolerance; STOP otherwise (a missing case is STOP). P2b never moves the verdict.
**P(GO) = 0.92 x 0.88 x 0.85 = 0.69,** conditioning P2a on P1 because both rest on the segment formula and the same
fit; treated as independent the product would be 0.92 x 0.88 x 0.80 = 0.65.

## Caps

Function level only, not operator level (MO's operator exponents are not addressed). 1s functions; three geometries,
two Yukawa masses; leading amplitude only; the segment formula's O(|a|^2) remainder is not bounded; the
Darboux/Lehmann reading of the coefficients is assumed, as in every prior record; P1 is a transfer test of a
confirmed mechanism; P2's v3 kappa-independence is a one-line consequence of the Helmholtz equation.

## Frozen (at registration)

sha256 of the working-tree files at registration (git may store text files with LF line endings, so `git show`
copies can hash differently with identical content).

| File | sha256 |
|:--|:--|
| `debug/l52_local.py` | b29f5a0ffbb454487caaf6fd403c3f0411375ec6ef33978b2972666bff9461d4 |
| `debug/l52_predict.py` | 6c06db7f0face33d24ea415109847c996d4e44bf0deff794dee00e8e0d10124c |
| `debug/l52_run.py` | 6df234cf754a3492619789add76bcec0736748406b8c3d57d9fafa86ca756be4 |
| `debug/l52_controls.py` | 572e6a78660d1a89322672c3578f7bc21bcf3e28c8ca6fa9e4aa4ae92a7dd844 |
| `debug/aabc_exact.py` | 218b7283f525d25ed8a87df61f0749ca6f576b20abcc4bde1272de32436b2d30 |
| `debug/na_exact.py` | 80acd793e8537b0e8a731ce7d0c47df2ed6ac5ab6564d4f7bbc14c0598be1b61 |
| `debug/aaeq2_local.py` | c6dd1c1775293e4e10be15092734891b3ee00011c4130791f8f7cd79ee6fed15 |
| `debug/data/l52_predictions.json` | 2a62d39003cf4c5d84dee5630e074c51eb71124468ba846c5a698ad3679b00e0 |
| `debug/data/l52_controls_sym.json` | 96fc88b7764ef97da40a39c49a79e2cf86418fc2a9f43a9fa8af345ecca7aa51 |
| `debug/data/l52_controls_const.json` | bf566aa75d965d1c2d8f2e3e408c27a6bd3376363267db96f6de456a2c29ece6 |
| `debug/data/l52_controls_equal.json` | eeb2289a91dbb7c51d01bbbbd4df69e30d7ac101ec2beec3dec54f5627b8026c |
| `debug/data/l52_controls_na.json` | ce7ed9c976dcd0a8e2bb71dc0fd7a769fa41704900ac0f8b0b86f3b777837b6f |
| `debug/data/l52_controls_toy.json` | 3c65cd5684859815fa4801753dcb6123181a47cb0c4027cdfac64882a653762b |
| `debug/data/l52_controls_toy92.json` | c7c097357d2e6f924fb79a8203e3d26ff59be3debe0e5e57f1a2fb660f9b0b34 |
| `debug/data/l52_ctrl_perL_A1.json` | 7bffa20fe594a4b903bb57d83125a9640cb4cae3ca6dd3df2638ab826882bdf0 |
| `debug/data/l52_ctrl_fk_A1_0_1_30.json` | 9d69a98565c762a0355833beeeefd4268fc165c322141a4e7c1cd7da1d509072 |
| `debug/data/l52_ctrl_fk_A1_4_12_30.json` | 8a19110708c0aeffecbbe54a5ddf37068b03f3315cd26a38ae2c620bf458ef19 |

## Phase 2 (after registration; at most 2 heavy processes)

```
python debug/l52_run.py fire
python debug/l52_run.py perL A1 ; python debug/l52_run.py perL A3
python debug/l52_run.py fk A1 ;   python debug/l52_run.py fk A3
python debug/l52_run.py run A1 2000 ; python debug/l52_run.py run A1 2600
python debug/l52_run.py run A2 2000 ; python debug/l52_run.py run A2 2600
python debug/l52_run.py run A3 2000 ; python debug/l52_run.py run A3 2600
python debug/l52_run.py score            # -> debug/data/l52_score.json (verdict, clauses, R_240, c1, beta per case)
```
Estimated 35-55 min total (AB: about 4.5 min per case at L = 240, 2000 digits; the Yukawa piece adds about a third;
perL and fk about 2 min per case).
