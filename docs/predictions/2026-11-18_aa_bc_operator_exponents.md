# Pre-registered: the contact exponents of the bond-angle operator move with the A density's odd terms (record OP52)

**Status:** registered 2026-11-18, before any scored computation (PM-directed follow-up to record L52 and the operator-level question AB P3a left open). Cold review: REGISTRABLE after fixes (2 must-fix, 5 should-fix, 4 nits), all applied. Nothing scored; NA - Y has never been formed at A1 or A4. Driver `debug/op52_run.py` (gate fire test 35/35); design `debug/op_exponents_52_memo.md`.

## Object

f'(x) = int V'(|r - A|) rho_BC(r) d^3r, AB's frame (B origin, C = (0, 0, D1), A at distance D2, x = cos of the
angle at B), with V'(s) = (1 - e^{-mu s})/s, the potential of the normalized density rho_A = mu^2 e^{-mu s}/(4 pi s).
In AB's engine (`debug/aabc_exact.py`, unchanged) this is **f' = NA - Y_mu** (kernels NA and Y only; no dY), mu = 2 zA.
Legendre coefficients L <= 1000 modulo primes 1073741441 and 1073741381 (MO's), and over Q to L = 240 for the
decay clause.

| case | D1 | D2 | zB | zC | zA (mu = 2 zA) | x* | D1 endpoint constants | D2 endpoint constants | status |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| A1 | 9/5 | 13/5 | 13/10 | 1 | 9/10 | 125/117 | e^{-117/50}, e^{-189/50}, e^{-513/50} | e^{-209/50}, e^{-389/50} | AB's geometry; f' never formed |
| A4 | 8/5 | 5/2 | 13/10 | 4/5 | 7/10 | 881/800 | e^{-52/25}, e^{-167/50}, e^{-391/50} | e^{-397/100}, e^{-653/100} | fresh: nothing ever computed |

Endpoint constants: D1 = e^{-zB D1} {1, e^{-mu(D2-D1)}, e^{-mu(D2+D1)}}, D2 = e^{+-zC D1 - (zB+zC) D2}. A4 was
chosen with no rate degeneracy: the inner rates -zB +- zC +- mu are 0.9, -1.9, -0.7, -3.5. None vanishes, none has
|b| = 1 (a log would vanish), none equals zB, zC or mu, and zB + zC - mu = 7/10 > 0.

## Derivation (before computing)

1. **Counting (derived).** Near x* every solution is u + (x* - x)^{1/2} v with u and v analytic. The half-integer
   exponents are 1/2 + ord(v) over a valuation-adapted basis of the space H of v's. The minimal operator of the
   generic sum has H = the sum of the pieces' H's.
2. **Endpoint bookkeeping (derived from the engine).** The radial integral about B splits at D1 (the sphere through
   C) and D2 (the sphere through A). D1-endpoint pieces carry e^{-zB D1} times the kernel's value at D1: NA gives
   no exponential, and Y gives the e^{+-mu D1} e^{-mu D2} halves of i_L(mu D1) k_L(mu D2). D2-endpoint pieces carry
   the density at D2, e^{-(zB+zC)D2} e^{+-zC D1}, times the kernel's jump there.
3. **mu-free jump (derived).** The jump of the Yukawa kernel at D2 is a Wronskian of i_L and k_L. Its exponentials
   cancel, so the D2 pieces see the whole kernel NA - Y, i.e. V'. (Measured consistent at A1 and A2 for f: the D2
   constants are fed by NA, Y and dY together; `op52_subsets_*` kernel split.)
4. **Odd terms (derived; sympy).** The odd powers make a branch (p meets C's cusp q = 1: (x* - x)^{(p+4)/2}):

   | potential | s^-1 | s^1 | s^3 |
   |:--|:--|:--|:--|
   | 1/s (NA) | 1 | 0 | 0 |
   | e^{-mu s}/s (Y alone) | 1 | mu^2/2 | mu^4/24 |
   | V' = 1/s - e^{-mu s}/s | 0 | **-mu^2/2** | -mu^4/24 |

   Contrast f (MO): the screened potential there was S = (1/s + zA) e^{-mu s}, whose s^1 term vanishes (v_1 = 0).
   Here it does not.
5. **Matched assumption M1 (under test): one branch dimension per (endpoint, potential).** In f this held for the
   Y/dY halves at D1 (one dimension, not two) and for the e^{+-zC D1} halves at D2. Here Y enters without dY, which
   is a new combination.
6. **Matched assumption M2 (under test): a collision of equal valuations resolves at the next order (+1).** In f:
   D1's V_A branch and D2's V_A branch (both valuation 3) gave 9/2. In NA: two point-charge branches (valuation 1)
   gave 5/2.

7. **The collision step, read two ways (both give 7/2 here).** In f the D1 collision of 1/s and S (both valuation 1)
   resolved at valuation 3, i.e. +2. That is not a failure of "+1": the difference potential V_A has no s^1 term
   (v_1(S) = 0). Reading (a): a collision at one endpoint resolves at the valuation of the difference potential.
   Reading (b), "+1": a collision between different endpoints (D1's V_A branch against D2's in f, giving 9/2). In f',
   (a) gives D1's 5/2 (V' has an s^1 term) and (b) gives the D1-D2 collision at valuation 2 resolving at 3 (7/2).
   The two readings do not differ on this object; OP2 tests (a), OP1 tests (b).

**Prediction from 1-7.** D1: span{1/s, e^{-mu s}/s} gives valuation 1 (3/2), and their difference V' gives
valuation (1+3)/2 = 2 (5/2). So the D1 block is {3/2, 5/2}. D2: V', valuation 2, so the D2 block is {5/2}. Union:
valuations {1, 2, 2}, and the collision at 2 resolves at 3. **Full operator: {3/2, 5/2, 7/2}; 9/2 absent.**
Function level (derived, steps 2-4 with the L52 segment mechanism): f' = analytic + Phi5 (x* - x)^{5/2} + ...,
with Phi5 proportional to v_1 = -mu^2/2, which is not zero. So f' leads with 5/2 and the Legendre coefficients decay
as rho^L L^{-3}.

## Controls (measured; files under `debug/data/`)

| id | what | measured |
|:--|:--|:--|
| K1 (AA\|BC) subsets, A1 | generic sums of MO's cached f pieces: D1 trio / D2 pair / all five | p 1073741441: {3/2, 7/2} / {7/2} / {3/2, 7/2, 9/2} (`op52_subsets_A1_1073741441.json`); p 1073741381: identical (`op52_subsets_A1_1073741381.json`) |
| K2 (AA\|BC) subsets, A2 | same | p 1073741381: {3/2, 7/2} / {7/2} / {3/2, 7/2, 9/2} (`op52_subsets_A2_1073741381.json`); p 1073741441: identical (`op52_subsets_A2_1073741441.json`) |
| K3 kernel split | which kernel feeds each endpoint constant (engine, L <= 12) | A1 and A2: D1 NA piece NA only; D1 mu pieces Y and dY; D2 pair NA, Y and dY (same files) |
| K4 NA at N3 (5/2 present) | one-electron, D1 piece / D2 pair / D1 + pair / all six | p 1073741441: {3/2} / {3/2} / {3/2, 5/2} / {3/2, 5/2} (`op52_na_1073741441.json`; prediction written before running, held); p 1073741381: identical, orders 5/6/7/8 (`op52_na_1073741381.json`) |
| K5 per-piece exponents, A1 f | the five x* carriers alone (orders 5, 6, 6, 10, 10) | 3/2, 3/2, 3/2, 7/2, 7/2; the 16 other pieces have no x* point (`op52_pieces_A1_1073741441.json`) |
| K6 absent exponents, known answers (toys, both primes) | Legendre coefficients of (x_t - x)^{k-1/2} (rational recurrence from sum t^L P_L, t = 1/3) | k = 3: {5/2}; k = 4: {7/2} (3/2 and 5/2 absent, read absent); {k = 2, 4, 5} generic: {3/2} (rank one: rational multiples of one function, lowest exponent only) (`op52_ctrl_toy_<p>.json`, pass on both primes) |
| K6b rank-one integrals | int_0^t c(s)(1 - 2xs + s^2)^{-1/2} ds with c of orders 1, 3, 4 at s = t | {3/2} only, both primes. **Disclosed:** first run with the expectation {3/2, 7/2, 9/2}, which was wrong: integration by parts reduces every such integral's x_t branch to a rational multiple of sqrt(1 - 2xt + t^2), so the module has rank one. The expectation was corrected, not the instrument |
| K6c rank two, independent of every engine (both primes) | order-2 Fuchsian operator with solutions (x_t - x)^{3/2} and (x_t - x)^k (x_u - x)^{1/2}, x_t = 5/3, x_u = -3 (sympy, verified to annihilate both); translated exactly into a Legendre-coefficient recurrence ((1-x^2) P_m' and x P_m three-term rules, (1-x^2)^2 L applied to P_m); one solution from random initial values, L < 500, through `generic_operator` / `xs_exps` | k = 5/2: recovered order 2, degree 4, x_t exponents {3/2, 5/2}; k = 7/2: order 2, degree 4, {3/2, 7/2}, **5/2 read absent**; x_u = -3 read class M {0, 1/2}. Identical on both primes (`op52_ctrl_toy_<p>.json`, sets `rank2_52`, `rank2_72`) |
| K7 decay toys | beta_eff(220) of (x* - x)^g Legendre coefficients (`op52_ctrl_betatoys.json`) | A1: g = 3/2 -2.0337, 5/2 **-3.0697**, 7/2 -4.1186; A4: -2.0277, **-3.0577**, -4.0986 |
| K8 decay instrument on f (record AB, L52) | beta_eff(220) of f vs the 7/2 toy | A1 -4.106 vs -4.119; A2 -4.083 vs -4.094 (AB P3b, 2000 and 2600 digits identical) |
| K9 gate fire test | `python debug/op52_run.py fire` | 35/35: 15 VOID, 11 STOP, 7 BORDERLINE, 2 GO (one a reported-only OP4/OP5 miss that must read GO) (`op52_ctrl_fire.json`) |
| K10 smoke | f' assembly code path at A2 (unscored), L <= 8, one prime | 21 constants; the five endpoint keys present and distinct. No operator and no exponent read |

**Disclosed contact with the object.** NA - Y has never been formed at A1 or A4. At A1 the per-kernel pieces have
been computed: AB's P3b reports beta_eff of NA alone (-2.029) and Y alone (-2.024); `debug/op52_groups.py` held all
three kernels' exact pieces at L <= 200 in memory (only f's constant groups were saved); `debug/op52_subsets.py` did
a kernel split at L <= 12. Nothing has been computed at A4. K10 formed f' at A2 (not a scored geometry) to L = 8.
**A1 is the geometry the rule was matched on** (MO's f and its subsets), so A4 is the clean test; both gate.

## Predictions (every probability is joint over both cases and both primes; a prime disagreement is VOID)

- **OP0 (instrument; VOID on failure).** K1, K2, K4, K6 and K6c pass on both primes; K9 passes; the frozen hashes
  match; no exception. At each case on each prime: a minimal operator is found (R <= 24); the five endpoint keys are
  present and distinct; the two weight vectors give the same exponent set at x* for the full operator and for both
  blocks; and the two primes agree on the full operator's (order, degree), its x* half-integer set and both blocks.
  The Q-side F_L at 2000 and 2600 digits agree to 1e-12 relative, L = 120..240. Probability 0.93.
- **OP1 (gates): the full operator's half-integer exponents at x* are exactly {3/2, 5/2, 7/2}**, 9/2 absent, at A1
  and A4. One clause; partial outcomes map to BORDERLINE or STOP through the gate. Probability 0.65.
- **OP2 (gates): blocks.** The generic sum of the three D1 pieces has {3/2, 5/2}; that of the two D2 pieces has
  {5/2}; both cases. Probability 0.75; 0.85 given OP1.
- **OP3 (gates): f' leads with 5/2.** |beta_eff(220) - toy_{5/2}| <= 0.10 at both cases (toy -3.070 at A1, -3.058 at
  A4; the 3/2 and 7/2 toys lie about 1.0 away on either side). **Implication:** f' is in the span of the pieces, so
  its branch has valuation 2 and 5/2 must be an operator exponent. OP3 therefore implies "5/2 present", and
  P(5/2 present at both cases) = 0.92 >= P(OP3) = 0.90. Given OP1 and OP2, OP3 is nearly implied (the remaining risks
  are a nearer singularity of f' and slow subleading terms): P(OP3 | OP1, OP2) = 0.97. Subleading size, from the L52
  segment formula: Phi7/Phi5 = 0.60 (A1) and 0.31 (A4). The reviewer's combined-toy computation shifts beta_eff(220)
  by about 0.004, so the 0.10 tolerance is safe. K8 puts f's distance from its own toy at 0.012 and 0.011.
- **OP4 (reported).** Which pieces carry x* (predicted: exactly the five endpoint pieces), orders, the full
  (order, degree), integer exponents at x*, other singular points, and whether beta_eff falls monotonically toward -3
  over L = 40..220.
- **OP5 (reported, not gating; added after review, the gate is not widened).** The 5/2 amplitude:
  Phi5 = v_1 N (-zC) C_11 (2 D1 D2)^{5/2} int_0^1 u (1-u) e^{-zB R_B(u)} du, with v_1 = -mu^2/2
  (`l52_local.amplitude(1, 1, ...)`). **Phi5 = -0.81014798223614941827 (A1), -0.2287660499486816037 (A4).**
  R_L = F_L / (Phi5 T^{5/2}_L) on L = 120, 140, ..., 240, using L52's normalization of F_L and its four-term fit
  r0 + c1/L + c2/L^2 + c3/L^3. Tolerance |r0 - 1| <= 1e-3. Probability 0.80.

**What each other outcome would mean.**
- {3/2, 5/2, 9/2}: the D1-D2 collision resolves two orders late (reading (b) fails). BORDERLINE.
- {3/2, 5/2, 7/2, 9/2} or more: an extra branch, e.g. Y's e^{+-mu D1} halves separate when dY is absent (M1 fails;
  OP2's D1 block would show three). BORDERLINE.
- {3/2, 5/2} only: D2's branch lies in D1's span to all orders. BORDERLINE.
- 5/2 absent at either case (including 5/2 present only at A4): the exponent set is not tied to the potential's odd
  terms, so the rule is wrong at its core. STOP.
- A D2 block other than {5/2}: the D2 jump is not V' alone (step 3 fails). BORDERLINE.
- beta near -2 or -4: the function-level derivation fails (a surviving point charge, or a cancelled s^1 term). STOP.

**Paper consequence per outcome (the edit is the PI's call).**
- GO: Paper 59's (AA|BC) paragraph gains one clause at [ARGUMENT + MEASURED]: the counting rule, with both objects
  (f: {3/2, 7/2, 9/2}; f': {3/2, 5/2, 7/2}) and NA's {3/2, 5/2} as instances.
- BORDERLINE: the clause is limited to what held, with the missed part named.
- STOP: no rule in the paper; the (AA|BC) per-piece table may be stated as [MEASURED] only.

## Decision gate (`op52_run.gate`; VOID > STOP > BORDERLINE > GO)

- **VOID:** any OP0 item fails (control, hash, exception, operator not found, keys, weight disagreement for the full
  operator or either block, prime disagreement, precision).
- **STOP:** a case or prime missing; or 3/2 or 5/2 absent from the full operator at either case; or OP3 misses.
- **GO:** OP1, OP2 and OP3 hold at both cases.
- **BORDERLINE:** otherwise, i.e. 3/2 and 5/2 present and OP3 holds, but the third exponent or a block differs.
OP4 and OP5 never move the verdict.

**P(GO) = 0.93 x 0.65 x 0.85 x 0.97 = 0.50.** P(VOID) about 0.07.

## Caps

1s functions only (plus the 1/s-weighted density on A); two geometries, two primes. Exponents are read mod p as
candidate indicial roots of a minimal operator found by guessing (gcrd of annihilators), not certified minimal and not
certified over Q. One fresh object. The constants are assumed linearly independent over Q-bar (as AB). M1 and M2 are
tested, not proved. The gating decay clause reads beta at L = 220 only; the amplitude (OP5) is reported.

## Frozen (working-tree sha256 at registration, from `python debug/op52_run.py freeze` -> `debug/data/op52_frozen.json`)

| File | sha256 |
|:--|:--|
| `debug/aabc_exact.py` | 218b7283f525d25ed8a87df61f0749ca6f576b20abcc4bde1272de32436b2d30 |
| `debug/na_exact.py` | 80acd793e8537b0e8a731ce7d0c47df2ed6ac5ab6564d4f7bbc14c0598be1b61 |
| `debug/intclass_lib.py` | b1655d6450fec0c1ed5f51d29e89e4d0cdb222da7eadf3d570286e1bff7467a6 |
| `debug/s1_controls.py` | 1c50016620c5f6291ac735028d9c87a658d8a9c226f28f08a1c667eba9436687 |
| `debug/na_modguess.py` | 84691c87a4a280b6823438af574c0fff6af9bf3f7ba7b76d93f2bbe866847d7d |
| `debug/na_modanalyze.py` | de92dc3acdcc7082001e32b50ef44db8c3bd75d990f615603b5f97f6ea178005 |
| `debug/l52_local.py` | b29f5a0ffbb454487caaf6fd403c3f0411375ec6ef33978b2972666bff9461d4 |
| `debug/l52_predict.py` | 6c06db7f0face33d24ea415109847c996d4e44bf0deff794dee00e8e0d10124c |
| `debug/op52_subsets.py` | ea760b51004502f8acea87dc6f2c1af8fb8da58a942b1014b5ebbd0612420a50 |
| `debug/op52_na.py` | 56579da7cc78773430d4cbe745d024f27d108f6425a3766198c505829be503ab |
| `debug/op52_pieces.py` | 064ec9dc7b512f3fb41b00709dbfbb44d580a07f2db3bcbb0d030d533d4517fc |
| `debug/data/op52_subsets_A1_1073741441.json` | 7c091560389bc20ee362836867b95dcb895b8407e3f257eccd06e16cc5305948 |
| `debug/data/op52_subsets_A1_1073741381.json` | f9e7be48b7baa536d366b6981b8264992a10aadf7b4037e73baacc0b6a0b1cb8 |
| `debug/data/op52_subsets_A2_1073741441.json` | c57369eb67759eb0d659c6c2f4867996fb19365107f25347960362fd50186091 |
| `debug/data/op52_subsets_A2_1073741381.json` | 2f2a60082ad2d6b4999d4a60387f05ede8a659bd7d04f341ae9f8baded390d32 |
| `debug/data/op52_na_1073741441.json` | 5bd15f46076681ed0e52d1eb509d09bd0cdf0696a39c0c271b4a786023f4ab50 |
| `debug/data/op52_na_1073741381.json` | 874e1233ca0d512d4520217763a46b45cc18eafd2486e2510766f7556a24a1bc |
| `debug/data/op52_ctrl_toy_1073741441.json` | a8932041cce44781dcd02bef46b5b5f5fafa37eb9a77fedd9da3bcfe022885ff |
| `debug/data/op52_ctrl_toy_1073741381.json` | a28bd819816bfa2f932250a126fb2914e95a32c4738b8c94977884c17d17d94e |
| `debug/data/op52_ctrl_betatoys.json` | dad4dffca871fe379f81f4e128cb7ced0f88639fdb6df17ad7d841633cbd9351 |
| `debug/data/op52_ctrl_fire.json` | e07141541564fc2c206b7fc78f915d9172e12f634224ac54e3c0d751506bce41 |
| `debug/data/s1_coeffs_A1_1073741441.json` | 8a7176193bf8324a94c212581b147efbb6950ed26df10c9b8edc12a988a8a288 |
| `debug/data/s1_coeffs_A1_1073741381.json` | f7faacfba11b473ce3e73eea25313c89ed13b31ecc6afa301072910cfce48c8b |
| `debug/data/s1_coeffs_A2_1073741441.json` | 8fc1a69159a6bdbb0c1e2b5fc5ed05828d5cd6b9b1511fe58014fbd16e350036 |
| `debug/data/s1_coeffs_A2_1073741381.json` | a73f47c7c30a5df9b2e3b8f952d7f47950412888a4e1e5dac517d51ebf5aaf19 |
| `debug/op52_run.py` | 0cd002500e564415074bd97f702a6876467721c5c8a9159cc4ceaafcfefd59b1 |

The four `s1_coeffs_*` files are MO's cached f sequences behind K1 and K2 (provenance).

Environment: Python 3.14.0, python-flint 0.9.0, mpmath 1.3.0 (as MO). The same hashes are in `debug/data/op52_frozen.json`, which `score` checks.

## Phase 2 (after registration and commit; one heavy process at a time, or two in parallel)

```
python debug/op52_run.py fire
python debug/op52_run.py run A1 1073741441 --registered ; python debug/op52_run.py run A1 1073741381 --registered
python debug/op52_run.py run A4 1073741441 --registered ; python debug/op52_run.py run A4 1073741381 --registered
python debug/op52_run.py beta A1 2000 --registered ; python debug/op52_run.py beta A1 2600 --registered
python debug/op52_run.py beta A4 2000 --registered ; python debug/op52_run.py beta A4 2600 --registered
python debug/op52_run.py score --registered       # -> debug/data/op52_score.json
```

**Cost.** Each run is about 7-9 min of coefficients (two of AB's three kernels; MO's three took 10.5-12.5 min),
plus about 8 min of operators (full operator and both blocks under 2 weights, 21 per-piece operators for OP4):
about 17 min each, 68 min for four. Each beta run is about 3-4 min of exact split plus evaluation: about 15 min for four. **Total about
85 min serial, about 45 min with two processes.**

## Results (scored 2026-11-18; `python debug/op52_run.py score --registered` -> `debug/data/op52_score.json`)

**Verdict: GO.** OP0 holds; OP1, OP2 and OP3 hold at A1 and at the fresh A4, on both primes and both weight
vectors. No VOID flag. The reported clauses OP4 and OP5 hold as well.

**OP0 (instrument).** All controls pass (`controls_pass`), including the independent rank-2 control K6c; the 26
frozen files hash-match (`hash_ok`); no run carries an error; keys present and distinct in all four runs; the two
primes and the two weight vectors agree on every reading below; 2000 against 2600 digits give identical beta values.
Gate fire test 35/35.

| reading | predicted | A1, p1 | A1, p2 | A4, p1 | A4, p2 |
|:--|:--|:--|:--|:--|:--|
| full operator (order, degree) | (reported) | (14, 109) | (14, 109) | (14, 109) | (14, 109) |
| OP1: half-integer exponents at x*, full operator | {3/2, 5/2, 7/2}, 9/2 absent | {3/2, 5/2, 7/2} | same | same | same |
| OP2: D1 block | {3/2, 5/2} | {3/2, 5/2} | same | same | same |
| OP2: D2 block | {5/2} | {5/2} | same | same | same |

p1 = 1073741441, p2 = 1073741381; every cell is the same under both weight vectors.

| | A1 | A4 |
|:--|:--|:--|
| OP3: beta_eff(220) of f' against the 5/2 toy (tolerance 0.10) | -3.0607 against -3.0697 | -3.0504 against -3.0577 |
| OP4 (reported): beta falls monotonically toward -3 | yes | yes |
| OP5 (reported): 5/2 amplitude, predicted Phi5 | -0.810147982236 | -0.228766049949 |
| OP5: measured / predicted, extrapolated r0 (tolerance 1e-3) | 0.99999926 | 0.99999969 |
| OP5: raw ratio at L = 240, fitted c1 | 0.99213, -1.82 | 0.99354, -1.50 |

**Reading.**
- The counting rule read from the (AA|BC) per-piece operators predicted the contact exponents of a different
  object before it was computed, and every gating clause held at both geometries: with a density on A whose
  potential has an s^1 term, the minimal bond-angle operator shows 3/2, 5/2 and 7/2 at x* and no 9/2; the near-end
  block shows 3/2 and 5/2, the far-end block 5/2 alone; and f' itself leads with 5/2. The rule now accounts for
  three objects: the two it was matched on ((AA|BC), {3/2, 7/2, 9/2}, no 5/2; one-electron NA, {3/2, 5/2}) and one
  predicted (f', fresh only at A4; A1 is the geometry the rule was matched on, and reads the same).
- This is consistent with the reading that 5/2 is absent from the (AA|BC) operator because the potential of a
  bounded density has no s^1 term (derived in L52); f' supplies an s^1 term and 5/2 appears. f' differs from f in
  more than that term (Y enters without dY, and the far-end jump is V' rather than S), so this is a reading, not
  an isolated cause.
- OP5: the segment-average amplitude formula of record L52 holds at p = 1, between L52's gating p = 3 and its
  p = -1 control: |r0 - 1| is 7.4e-7 and 3.1e-7 with the registered four-term fit, and about 1e-9 with five terms.
- Caps as registered: 1s functions; two geometries; exponents are candidate indicial roots read mod p on two
  primes, not certified over Q; the operator is not certified minimal; the 21 constants are assumed linearly
  independent; the gating decay clause reads beta at L = 220 only; the rule's two matched steps (one branch
  dimension per endpoint and potential; the collision step) are supported by this prediction, not derived, and
  the two readings of the collision step coincide on this object, so the data do not tell them apart; one fresh
  object.
- Paper consequence as registered for GO: Paper 59 gains one clause at [ARGUMENT + MEASURED] with the rule and the
  three objects.

**OP4 (reported), in full.** x* is a singular point of exactly the five endpoint pieces: the three D1 pieces have
orders 5, 6, 6, each with {3/2}; the two D2 pieces have order 10, each with {5/2}. The integer exponents of the full
operator at x* are 0..10 (all 14 indicial roots are accounted for: 0..10, 3/2, 5/2, 7/2). The D1 block has order 11
and the D2 block order 10. The full operator has the same order as f's (14) with three half-integer branches in
both; Y without dY did not add a branch dimension. Its degree is 109 against f's 111; that difference is not
interpreted. The other singular points of the operator were not recorded.

Cold review of the scoring: 0 MUST-FIX, 5 SHOULD-FIX (the Reading stated a cause the data do not isolate; "three
objects" now names the two matched ones; four registered caps restored; OP4 reported in full; the OP5 figure
relabelled), 4 NIT; applied. The reviewer recomputed one case with its own weights and its own indicial reader (all
14 roots, 9/2 not among them), re-read every gating clause from the raw runs, rechecked the 26 hashes and the file
times against the registration commit, refitted OP5 with 3 to 6 terms, and confirmed the GO. Its caveat: the
coefficients and the operator guesser are the same frozen engine; only the weights and the indicial reader were
independent.

**Added by the guard pass (2026-11-18), not in the scoring above.** The full operator found for the random
combination annihilates 19 of the 21 pieces individually. The two it does not are non-endpoint pure exponentials,
regular at x*, the same two at both cases and both primes. (A generic combination of pieces that differ by a
rational factor has a lower-order annihilator that depends on the weights.) All five endpoint pieces are annihilated,
so the x* exponents do not depend on the weights; the guard pins this.
