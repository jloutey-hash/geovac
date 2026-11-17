# Record EP design: two-centre exchange with an l > 0 orbital (row 3c, "exchange beyond 1s")

*Design stage, 2026-11-16/17. No scored quantity computed: at the scored cases only booleans,
|differences| against the exact reference, valuation flags and timings were recorded. The PM writes
and registers the record.*

## Object

f(R) = (ab|ab), a on A, b = 1s on B, Neumann sum in prolate coordinates, hydrogenic orbitals
(`geovac.qfd_core` convention, rate = Z/n). Two families:

- **S** = (2pσ_A 1s_B | 2pσ_A 1s_B): σ = 0; density polynomial (ξ²−η²)(1+ξη), so kernels X_{j1j2},
  j ≤ 3, and η-moments k ≤ 3.
- **P** = (2p₊₁,A 1s_B | 2p₋₁,A 1s_B) (the radial function of real 2p_x exchange; the e^{±i(φ1−φ2)}
  part of cos φ1 cos φ2): σ = 1, H = 1, kernels j ≤ 2. Uses P_τ¹, Q_τ¹; the (ξ²−1) of the density
  cancels the pole of Q_τ'.

Scored pairs (a = A-orbital rate, b = 1s rate): **(1, 1)** and **(1, 2)** for each family
(reference call `exchange_hp(2a, orb_A, (1,0,0), b, orb_A', (1,0,0), R)`).

## Instrument

`debug/ep_gen.py` generalises the frozen `exch_gen.py` by import (Tables, Legendre data, the
bivariate evaluators reused unchanged): g(u) = ∫_{u/2}^{u−1}[F_{j1}(u−y)G_{j2}(y) + F_{j2}(u−y)G_{j1}(y)] dy
with G = B·Q0 − C; for σ = 1, C = y^j[(y²−1)^H V_τ' + (y²−1)^{H−1}P_τ]. The η-moments for σ = H = 1
are by parts from the closed Legendre moments. The extension cost was small (one function body
copied with general F, B, C; the frozen assertion `top <= 2τ+6` is the only reason it could not be
called as is). Prefactors derived from the normalised orbitals: S a⁵b³/32 R⁷, P a⁵b³/64 R⁷ (sign
confirmed by control C1). Truncation: a proved valuation bound per τ-term,
v(τ) = min[rexp + v_β(k1) + v_β(k2) − (j1+j2+4H−2σ+1)], v_β(k) = max(0, τ+σ−k−2H) (orthogonality of
D^σP_τ under (1−η²)^σ); v(τ) = 2τ − 2 for both families (2τ for 1s), so τ ≤ (N+2)/2. Checked on
every term built (pass/fail).

`debug/ep_run.py` (written, run only on the 1s control): rebinds the frozen `exch_run` module's
`generate`, `outname`, `_run_jsons` (file not edited) and so uses `readings`, `crt_operator`,
`mode_heldout`, `frobenius` unchanged (r ≤ 14, d ≤ 60, HOLD 40, MARGIN 50, N 450; heldout exact Q
residual to max(120, r+d+80), P3 = 16777213 with 60 held-out orders). New: `ladder` mode — the two
scored primes 67108859, 33554393 first, then descending primes < 2²⁶; acceptance = same (order,
degree, per-coefficient degrees) as on 67108859, all rows annihilated, valuation passed; after each
acceptance from k = 4, CRT + reconstruction; stop when k0 and k0 + 4 both reconstruct the same
operator; cap 96 candidates. `numeric` mode: Frobenius vs `exchange_hp` on the case's orbitals.

## Derivation (before guessing)

1. **Constant class.** By construction, identical to row 3c: every τ-term is a Laurent polynomial
   times {1, γ + ln p, ln 2, e^{2p}E1(2p), e^{4p}E1(4p)} times e^{−2p} times β(q)², and l > 0 only
   changes the polynomials (S) or replaces P_τQ_τ by P_τ¹Q_τ¹, whose (ξ²−1)-cleared form is again
   polynomial × {1, Q0} (P). So components ⊂ {1, γ, ln 2, ln(primes of κ)} × (ln R)^{≤1}, weight one.
   This is the prior-art class (Lesiuk–Moszynski 2014 II; the corpus guard
   `tests/test_paper58_exchange_fast.py` covers σ = 1 values); not a scored claim.
2. **Equal exponents: order ≤ 4, theorem-level.** At q = 0 the τ-sum terminates (τ ≤ 3 for S and P,
   control C3), so f = r1 e^{−2p} + r2 e^{−2p}(γ + ln p) + r3 E1(2p) + r4 e^{2p}E1(4p) (+ constants
   absorbed) with r_i Laurent polynomials in R. This C(R)-module is closed under d/dR and has rank 4,
   so the minimal operator has order ≤ 4 for any polynomial weights, l included. It is exactly 4 iff
   no accidental relation among the r_i; 1s gave 4. Rates ⊂ {0, 2ζ, −2ζ}.
3. **Unequal exponents: no proof.** The τ-sum is infinite, the coefficient functions are not
   rational, and EX's order 9 is empirical. Heuristic: the exponential skeleton is the same as for
   1s, so the same nine formal solutions {−(a+b), 2a, 2b, 0², −2a², −2b²} are expected; l > 0 should
   raise the degree, not the order. The rate set {s(a+b)+t(a−b)} is forced for any l.
4. **R = 0.** The ln R part enters only through γ + ln(κR), so ln R and γ components coincide up to
   sign. Onset heuristic (EX4b's, the "coalescence of the four leading non-smooth parts"): an orbital
   r^l Y e^{−ar} has leading non-smooth part −a r·r^l Y, degree l + 1; onset = Σ(l_i + 1) − 1 + 6 =
   **R^11 ln R** for S and P (9 for 1s, which held). Scaling f_{λa,λb}(R) = λ f_{a,b}(λR) fixes the
   coefficient's total degree at 12; the heuristic split is normalisation a⁵b³ × cusp factors a²b²,
   i.e. c·a⁷b⁵ with c rational (ratio (1,2)/(1,1) = 32). Weak: the l > 0 pieces z·1 of the 2p
   orbital are smooth but not constant, and nothing rules out a lower onset from mixed terms.
5. **Does l > 0 raise the order?** Equal exponents: no (≤ 4, proved above). Unequal: not by the
   skeleton count; a rise would mean the infinite τ-sum activates a new companion (e.g. ln at
   −(a+b), or the a − b atoms that 1s left empty). Degree: the planted control C5b shows a degree
   rise alone (16 → 34) costs almost no height.

## Controls (`debug/data/ep_controls.json`; `ep_ladder/heldout/numeric_1s_1_1_EPT.json`)

| id | what | result |
|:--|:--|:--|
| C0 | ep_gen "1s" vs frozen exch_gen, mod 67108859 to R^450 at (1,1), (1,2), (1,3/2); Q to R^60 at (1,3/2) | identical, bit for bit, all four |
| C1 | ep_gen S, P series (Q, N = 90) vs `exchange_hp` (dps 50, finite τ) at (1,1), R = 1/2, 6/5 | S 8.1e-65, 8.2e-63; P 3.7e-63, 3.2e-63 (sign and prefactor confirmed) |
| C1f | same at (1,2) and (2,1): Q series N = 60, R = 1/2, reference τ_max 16 (τ 12 vs 16 differ by ≤ 1.3e-37) | S 2.0e-54, 7.9e-54; P 7.3e-55, 2.9e-54; all pass 1e-35 |
| C2 | Q series reduced mod 33554393 = mod-p series (N = 60); two mod-p runs identical; valuation, all four scored cases | all true |
| C3 | scored cases, N = 450 mod 67108859: valuation flag, timing (series discarded unopened) | all pass; τ count 4, 227 (S); 3, 226 (P); 3.9 s, 488 s (S); 4.2 s, 175 s (P) (5 jobs sharing the CPU) |
| C4 | positive control: frozen `exch_run.readings` on ep_gen "1s" | (1,1) order 4, degree 8, first hit (4,8); (1,2) order 9, degree 16, first hit (9,16); both equal the frozen EXS operators mod p coefficient for coefficient |
| C5 | normalised heights of the certified 1s operators (num/den bits; primes needed) | (1,1) 12/6, 2; (1,2) 44/22, 4; (1,3/2) 52/26, 5 (EX's 4-prime failure reproduced) |
| C5b | planted (1+R²)·f_1s(1,2): known order 9 | order 9, degree 34 on every prime; reconstruction fails at k ≤ 4, k0 = 5, the same operator at every k = 5..16; heights 55/23 bits |
| C6 | end-to-end `ep_run` on 1s (1,1), suffix EPT: ladder, heldout, numeric | 8 primes accepted, k0 = 4, k0+4 agrees; order 4 degree 8; exact Q residual zero on -4..116; P3 nullity 1, matches; Frobenius vs exchange_hp 3.7e-92 (R = 1.4), 1.2e-88 (R = 3.015) |

**Prime budget.** Heights depend on the exponent arithmetic far more than on degree (C5b: degree
16 → 34 moves the numerator 44 → 55 bits; (1,3/2) needed more than (1,2)). S and P at (1,2) carry
a⁵b³ and higher-degree polynomials; expect 5–10 primes. The adaptive ladder makes the count an
output: accept from k = 4, require k0 and k0 + 4 to agree; cap 96. Use only pairs with small
denominators ((1,2)); do not use (1,3/2) here. EX's VOID (4 fixed primes) cannot recur unless k0 > 92.

**Cost.** Generator per prime at N = 450 for the unequal pair: S about 5–8 min, P about 3 min (C3,
under load); scan comparable to 1s (C4: 169 s at (1,2)). With 12 primes per unequal case, 1 to 2 h
each. The exact-Q generator is the slow step: S at (1,2) took 5.6, 20, 45 s to N = 20, 30, 40 and
~13 min to N = 60 (with the reference), roughly N^3.5, so heldout's N_Q = max(120, r+d+80) costs
about 2–4 h for S (1,2) and numeric's N_Q = 120 the same again (P about a third). Budget the 48 h VOID
clock per case, or run cases in parallel. (An N = 90 C1 at the unequal pairs was stopped after ~3 h;
C1f replaces it.)

## Proposed clauses (PM to set; probabilities mine)

- **EP0 (gate).** `ep_run run 1s 1 2` on both scored primes reproduces the frozen EXS (1,2) operator
  mod p; C0 and C4 rerun identical; valuation passes in every output; frozen files hash-match.
  Failure VOID.
- **EP1 (existence, all four cases).** Ladder finds an accepted operator on both scored primes with
  equal shape, k0 ≤ 92 and k0 + 4 agree; heldout exact residual zero, P3 nullity 1 and matches;
  numeric |Frobenius − exchange_hp(dps 80)| < 1e-30 at R = 1.4, 3.015, dps 50 vs 80 < 1e-35.
  Instrument parts failing after an operator is found: VOID. Probability 0.8.
- **EP2a (equal exponents, order).** Order exactly 4 for S and P at (1,1) (≤ 4 is derived; = 4 is
  the prediction). 0.8.
- **EP2b (unequal, order).** Order 9 for S and P at (1,2), same as 1s. 0.55. Sub-clause **EP2b'**:
  order in [9, 14]. 0.85.
- **EP2c (degree rises).** Degree above the 1s degree at the same pair (> 8 at (1,1), > 16 at (1,2))
  in every case. 0.75.
- **EP3 (no climb at finite R).** Every finite nonzero point apparent or exponents_only, both primes,
  all cases. 0.75.
- **EP4 (R = 0).** Fuchsian, log solution present, f has a nonzero ln R component, all cases. 0.8.
  **EP4b**: onset R^11 ln R in all four cases. 0.45. **EP4c**: coefficient ratio (1,2)/(1,1) = 32 in
  each family (a⁷b⁵ at a = 1). 0.4.
- **EP5 (rates).** Rank ≤ 1 at infinity; every rate in {s(a+b)+t(a−b)}; equal: set {0, 2, −2};
  unequal: the 1s multiset {−3, 4, 2, 0², −2², −4²} at (1,2). 0.6 (set) / 0.5 (multiset).
  **EP5b**: no a − b rate (±1 at (1,2)). 0.6.

## Gate

- **VOID**: EP0 fails; scored primes disagree in shape; any EP1 instrument part fails after an
  operator is found; reconstruction never stabilises within the cap; a step unfinished 48 h after launch.
- **STOP**: no operator accepted in the box on both scored primes for some case (box, not class, is
  the reading); or a finite nonzero point genuine on both primes (a climb).
- **GO**: EP1–EP5 hold at both families (EP4b/c, EP5b, the EP5 multiset and EP2b exactness may miss,
  reported). **BORDERLINE** otherwise.

Scoring order: S (1,1), P (1,1), S (1,2), P (1,2), each `ladder`, `heldout`, `numeric`, suffix EPS.

## Caps

l = 1 on A only, 1s on B; σ ≤ 1 in the generator (d orbitals need D^σ for σ = 2, a small extension);
two pairs; mod-p singular readings; gcrd not certified minimal; constants formal.

Files: `debug/ep_gen.py`, `debug/ep_controls.py`, `debug/ep_run.py`; outputs `debug/data/ep_*`.
