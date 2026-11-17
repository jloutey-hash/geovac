# Pre-registered: on the matched slice the four-centre 1s integral's L-recurrence is found within a fixed budget, needs depth-2 exponential integrals, and keeps 7/2 at both contacts (record FX)

**Registered:** 2026-11-17, before any scored computation. This is the owed item "row 7 exact recurrence"
(inventory row 7). It was blocked at a constant-class step until the depth-2 reducer was built.
**Owner:** Paper 63 `sec:four_centre` (record FC GO: numeric decay, no operator); inventory row 7, operator column
MISSING.
**Design:** `debug/fx_design_memo.md` §1–5 (object, the matched condition, the representation), §7 (reducer
built, FX0) and §8 (canonical basis, driver, clauses).
Cold review before commit: 2 MUST-FIX, 8 SHOULD-FIX, 7 NIT, all applied. The MUST-FIXes: (1) a driver bug left FX1's L-independence unable to fail, because a late-appearing key was counted present everywhere; fixed and fire-tested. (2) FX0 named a rerun no command performed and had no pass criterion; it is now scored mechanically from named fields. The reviewer checked the canonical reduction (disjoint identity sets, an order-independent normal form), the matched geometry, the x* values and the control numbers, and found no leak.

## Object

Take the four-centre 1s integral (XY|ZW) on FC's slice, at the fingerprint exponents zX = 13/10, zY = 1:
- X sits at the origin and Y at D1 z;
- Z = a n and W = b n, with x = n·z.

Write F_L for its Legendre coefficient in x. The **matched slice** is the condition zW/zZ = √(a/b), with a/b a
ratio of squares. Within FC's slice it is the only exponent ratio for which ρ2_L is elementary: e-terms plus E1-terms
whose coefficient is R_L ∝ (b − a)², checked to 1e-39. Off it, ρ2_L is a Bessel-K form, closed but not elementary. There F_L reduces exactly to a finite set of constants with
rational coefficients:

G_L = F_L/π² = Σ_c q_{L,c} c.

The constant set includes depth-2 iterated exponential integrals, which is new against the AB/NA class. Two
cases are scored:

| case | (a, b, D1) | (zZ, zW) | contacts x* = (D1² + c²)/(2 D1 c) |
|:--|:--|:--|:--|
| M1 | (1, 4, 6/5) | (6/5, 3/5) | x_Z* = 61/60 (dominant), x_W* = 109/60 |
| M2r | (4/9, 1, 6/5) | (3/2, 1) | x_Z* = 829/540, x_W* = 61/60 (dominant) |

On M2r, zW/zZ = 2/3 = √(4/9), so it is on the matched slice. The design's original M2 had coinciding contacts and
was replaced. The design memo fixed M2r's geometry but not its exponents. The PM confirmed (3/2, 1), M2's own pair,
before registration; it is the matched pair at a = 4/9, b = 1. The cold reviewer checked this independently.

**Existence is not the question.** F(x) is an integral of holonomic functions, so it is holonomic, and a recurrence
exists by closure. FX3 tests that the recurrence is **found within the budget**: R ≤ 32 on L ≤ 1300 (degree ≤ 34 at
R = 32). A STOP on FX3 means the operator is beyond the budget, not that none exists.

**Terms.**
- A **depth-2 key** is a key containing a named base B with one or two E1 factors. As iterated integrals, the E1·E1
  bases are depth 3.
- FX1's **alphabet** is mechanised from the same engine structure (§3's rule), so its containment clause mostly
  checks the engine against its own rule. It does not constrain the e^{e0} prefactors, and it always admits γ.

## Instrument

- **`debug/fx_reducer.py`.** An integration-by-parts reducer for e^{λr} X(r) E(r), where E is a product of zero
  to two factors E1(cr + d). It works panel by panel, with partial fractions over the poles and a Hadamard finite
  part at r = 0. It runs over fmpq and nmod.
- **`debug/fx_exact.py`.** Produces the exact G_L.
- **`debug/fx_indep.py`.** Reduces to a **canonical basis** modulo the product-rule identity
  Σ_f m_f e^{−d_f} B(−c_f, F − f, q_f) = −[∏E1]_A^B. The proof is in the module docstring. Without this step, a
  nonzero coefficient on a named depth-2 base would be necessary but not sufficient for "depth-2 content
  survives", because the bases depend on the order of integration and are linked by this identity.
- **`debug/fx_run.py`.** The scoring driver: the mod-p generic sequence, NA's bisection on order, MO's gcrd for
  the minimal operator, and AB's `local()` exponents. It refuses to score without `--registered`.

Settings:
- LMAX = 1500, with a guess window of L ≤ 1300 and 200 rows held out;
- RMAX = 32, seeds 2026 and 7;
- primes 1073741441 and 1073741381 (AK's).

## Controls (measured; `debug/data/fx_gate_controls.json`, `fx_indep_controls.json`, `fx_runctl_*.json`)

| id | result |
|:--|:--|
| V1 | Exact G_L vs numeric `fc_engine`, L ≤ 30: M1 8.7e-29 / 7.2e-32, M2r 3.0e-38 / 3.4e-42 (bar 1e-25). Evaluation at 320 vs 400 digits agrees to 3e-180 (M1), 8e-174 (M2r) |
| V2 | a = b with zZ ≠ zW (K4A1 = AB's A1) reproduces `aabc_exact` exactly, L ≤ 40 |
| V3 | mod p on two primes (2305843009213693951, 4611686018427387847; not the scoring primes, which are covered by A4 and the positive controls) equals fmpq reduced mod p, L ≤ 40, for M1, M2r, K4A1 |
| V4 | NA's (R, Q0) = (8, 31) reproduced on V3's two primes |
| A1 | Product-rule identity: 80 sets, all panel types, to 3.6e-60 |
| A2 | Planted reducible combination: 26 → 0 depth-2 terms. Non-reducible combination keeps one. mod p equal |
| A3 | Raw vs canonical at 350 digits: value preserved on C1/C2 (≤ 3e-320) and on M1/M2r (booleans only) |
| A4, A5 | Canonical mod p = canonical fmpq mod p; a = b gives no named base |
| B1–B4 | LLL (300–550 digits) recovers the identity and a planted relation. No relation among the 77 canonical constants of C1 at height ≤ 1e4. Adding back one eliminated base is detected |
| alphabet | FX1's candidate set, calibrated on the unscored matched geometries C1 and C2 (L ≤ 60). One rule bug was found and fixed (an E1 at the r = 0 endpoint). After the fix: 0 keys outside |
| posctrl NA | search shape (8, 31), minimal (8, 25); singular points identical to the frozen S1 record; 200 held-out rows verified |
| posctrl AB | search (17, 48), minimal (14, 111); points identical to MO's W0; 7/2 at x* |
| posctrl K4A1 | through the FX engine: (17, 48), (14, 111) = AB's; 200 rows verified |
| e2e K4A1 | end-to-end score + summary after the FX1 fix (`fx_runctl_e2e_k4a1.json`): FX1 holds (19 keys, no interior zeros); FX2 false at all 1199 L (correct for the a = b class); FX3 (14, 111); FX4 holds; FX5 None; gate BORDERLINE (correct for this negative case), `as_expected` |
| FX1 fire test | planted keys through the real accumulator and classifier (`fx_runctl_firetest_fx1.json`). Caught as misses: a late key (initial run 2..9) and a vanishing key (terminal run 90..100). Classified correctly and not counted as misses: an isolated zero (interior), an L = 0, 1-only key (out of scope), a zero on one prime only (not absent). An isolated L without depth 2 fails FX2. The pre-fix rule misses the late key |

**Disclosure.** On the unscored matched geometries C1 and C2, for L ≤ 60, two things were seen during calibration:
the canonical key set did not depend on L, and depth-2 keys were nonzero at every L. This bears directly on FX1 and
FX2, whose probabilities below were raised from the design's 0.7 and 0.85 after seeing it.

## Predictions (each clause on both M1 and M2r)

- **FX0 (gate; VOID on failure).** `fx_run.summary` scores this clause mechanically, as `fx0_check` (34 checks).
  The control steps of Phase 2 are rerun with the frozen files, and FX0 then requires their pass fields:
  - `fx_gate_controls.json`:
    - V1 `pass_1e_25` (M1, M2r; levels A and B);
    - V2 `all_equal_L_le_40`;
    - V3 equal on each case and both primes;
    - V4 `reproduced` on both primes.
  - `fx_indep_controls.json`:
    - A1 `pass_1e_50`;
    - A2 `pos_reduced`, `neg_survives` and `modp_equals_fmpq_mod_p`;
    - A3 `value_preserved` and `rewrite_exercised`;
    - A4 all true; A5 both true;
    - B1 `equals_identity`; B2 `involves_planted`;
    - B3 and B4neg `relation_found` false at both precisions;
    - B4pos `involves_eliminated`.
  - `fx_runctl_ctrl_alphabet.json`: `n_outside` = 0 and `L_independent`, for C1 and C2.
  - The positive-control files on both primes:
    - NA: reproduces frozen, 3/2 at x*, verified;
    - AB: reproduces frozen, 7/2 at x*;
    - K4A1: equals AB A1, verified, 7/2 at x*.
  - `fx_runctl_e2e_k4a1.json`: `as_expected`.

  In addition, the natural factor sets must be L-independent for L = 2..6, and every annihilator must vanish on its
  held-out rows. The reruns overwrite the committed control JSONs; the committed copies (git) are the reference.
  The checker was fire-tested by flipping one field in each of five control files: each flip was caught
  (`fx_runctl_firetest_fx0.json`).
- **FX1 (constant set; P = 0.8).** Holds when both of the following hold:
  - In the canonical basis, no key of the key set (the keys nonzero mod either prime at some 2 ≤ L ≤ 1500) is
    absent on an initial or terminal range of L.
  - Every atom lies in `alphabet(case)`.

  Sporadic isolated zeros of a coefficient are reported, not counted as a miss.
- **FX2 (class; P = 0.9).** For every 2 ≤ L ≤ 1500, some canonical depth-2 key is nonzero on either prime.
  - This holds modulo the product-rule family. B4 found no further relation at height ≤ 1e4 on C1.
  - `fx_run.json` reports the uncapped count of L without a depth-2 key, and whether that is every L or isolated L.
  - **What a miss means depends on where it happens.**
    - No depth-2 key at every L: F_L lies in the depth-0 span modulo the family, and the matched slice stays in the
      AB/NA class. That would be a finding.
    - Misses only at isolated L: sporadic zeros, from which no class reading follows.
  - **Cap.** B4's LLL basis left out π², Li₂ special values and products of logs. A half-line base with a Li₂-type
    value would count as depth-2 by name. This is unlikely, since no panel is [0, ∞), but it is untested.
- **FX3 (operator; P = 0.5).** For seeds 2026 and 7 on both primes, all of the following hold:
  - NA's bisection finds an operator at order R0 ≤ 32 on the window L ≤ 1300;
  - every annihilator at orders R0..R0+2 vanishes on all rows to L = 1500;
  - the minimal operator (the gcrd) has the same (order, degree) in all four runs.

  The search shape and the minimal shape are reported, not scored. A plausible range was 12–30. An annihilator that
  fails its held-out rows is VOID (FX0), not an FX3 miss.
  - Reported, not scored: the full set of half-integer indicial roots at each contact. At AB's x* the candidates
    included 3/2, 7/2 and 9/2, so "7/2 is a root" alone discriminates weakly.
- **FX4 (contacts; P = 0.8 given FX3).** In all four runs, x_Z* and x_W* are factors of the minimal operator's
  leading coefficient, and 7/2 is an indicial root at both (`local()`). Whether 7/2 is the leading branch cannot be
  decided mod p, because a root's branch can cancel (AB's 3/2). That question stays with FC's numerical GO.
- **FX5 (FC's assumption; reported, outside the gate).** This asks for genuine, weight-fixed singular points of the
  minimal operator other than {±1, x_Z*, x_W*}, located by CRT over the two primes at height ≤ 1e6. Three outcomes:

  | outcome | meaning | probability |
  |:--|:--|:--|
  | False | some such point lies in (1, x_dom*) | 0.2 |
  | None | some genuine factor is non-linear, or not reconstructable | 0.4 |
  | True | otherwise | 0.4 |

  Strict containment in {±1, x_Z*, x_W*} is reported only, and is expected to fail: NA and AB both carry extra
  genuine log points. In both cases the dominant contact is 61/60, so the interval (1, 61/60) is very narrow, which is
  why P(False) is low.

## Decision gate (per case; overall precedence VOID > STOP > BORDERLINE > GO: the record is VOID if either case is VOID, else STOP if either is STOP, else GO only if both are GO)

- **VOID:** any of the following:
  - FX0 fails;
  - any instrument exception during scoring;
  - a step has not finished 48 h after the launch timestamp recorded in `fx_run_<case>_<p>.json`.

  Primes and seeds are not substituted after launch. The design memo's "6 h cap" is a cost budget, not a gate.
- **GO:** FX1, FX2, FX3 and FX4 hold.
- **STOP:** either of
  - no operator up to R = 32 on either prime (that is, on neither prime);
  - a contact present without 7/2 among its exponents.
- **BORDERLINE:** otherwise. This includes an FX2 miss.

## Caps

- Two geometries on one slice (X–Z–W collinear, matched exponents). Read in x = cos(bond angle), 1s only.
- Operators and singular points are read mod p. The gcrd is not certified minimal, and the constants are formal.
- FX2's "depth-2 survives" holds modulo the product-rule family. Independence beyond it is tested only by LLL at
  bounded height on an unscored geometry.
- FC's leading-branch claim is not rescored here.

## Frozen (this commit; unchanged at scoring)

| file | sha256 |
|:--|:--|
| `debug/fx_run.py` | `63ba74df1653c7824f644aca20e8ebf1ff03c225c9f1e5506c283b5ed06fc5c6` |
| `debug/fx_indep.py` | `aa178b8a6ac1ec37aa34db07c6ffed9d982faf9f2eccaa2f666895e3cc472110` |
| `debug/fx_reducer.py` | `39afb54f67858e859b8db9fe274d9022ed5dfed187655ba49c0da9e8e070d246` |
| `debug/fx_exact.py` | `c585a9262a72b3d572697d27245ce54a3463d33ff38af8dec7fe58ebc7ea3134` |
| `debug/fx_gate_controls.py` | `ff72b0d8e2d18696e4a660b82990e65605d92bc172b6199afbdbe3a530678727` |
| `debug/fx_design_memo.md` | `e9b033818ff68ab5ca0f7645f0c7b655568b509dd31d12a304bedfd061ab5ab6` |
| `debug/fx_engine.py` | `b6131b84a8ed51036eacc14d3f6531aad23dd9fc12baee44fcc61b8877c616bb` |
| `debug/na_exact.py` | `80acd793e8537b0e8a731ce7d0c47df2ed6ac5ab6564d4f7bbc14c0598be1b61` |
| `debug/na_modguess.py` | `84691c87a4a280b6823438af574c0fff6af9bf3f7ba7b76d93f2bbe866847d7d` |
| `debug/na_modanalyze.py` | `de92dc3acdcc7082001e32b50ef44db8c3bd75d990f615603b5f97f6ea178005` |
| `debug/intclass_lib.py` | `b1655d6450fec0c1ed5f51d29e89e4d0cdb222da7eadf3d570286e1bff7467a6` |
| `debug/s1_controls.py` | `1c50016620c5f6291ac735028d9c87a658d8a9c226f28f08a1c667eba9436687` |
| `debug/aabc_exact.py` | `218b7283f525d25ed8a87df61f0749ca6f576b20abcc4bde1272de32436b2d30` |
| `debug/fc_engine.py` | `9fef5cf60eeec833283f2211050af56e38d163b71b4fcf6a6bbe1744f56707aa` |
| `debug/data/s1_controls.json` | `adf21b1ee456af41cfd3df5871a9a8e92d314b2365d0e0290a3473d23cbbfd8a` |
| `debug/data/s1_coeffs_A1_1073741441.json` | `8a7176193bf8324a94c212581b147efbb6950ed26df10c9b8edc12a988a8a288` |
| `debug/data/s1_coeffs_A1_1073741381.json` | `f7faacfba11b473ce3e73eea25313c89ed13b31ecc6afa301072910cfce48c8b` |
| `debug/data/fx_gate_controls.json` | `903b291cd22936cdabb63f73d4e7205e62fa83db397d76141e02a81493241c3e` |
| `debug/data/fx_indep_controls.json` | `0b8ded4c46d1502be0f9b5a766bcf9cdc643984ce5dcfb05d56621196f8eee1b` |
| `debug/data/fx_runctl_ctrl_alphabet.json` | `a73b21f18b3f6e52fa929a7e76d0285d256fb3bd3b2dcf979d00601ff1d3ecc0` |
| `debug/data/fx_runctl_posctrl_na_1073741441.json` | `4670de7b692b3f18ebcf9ec2ff6cfb0b766b39ec84afe07d321ecacb6635a008` |
| `debug/data/fx_runctl_posctrl_na_1073741381.json` | `93f76371db884139905da4dbb9b9fa064b493910fc8fc396af4445ead806e5b4` |
| `debug/data/fx_runctl_posctrl_ab_1073741441.json` | `388b14d58fe68d5bfe71d77acf43e82fe6aa522d4a44ab653049553b4a8da50c` |
| `debug/data/fx_runctl_posctrl_ab_1073741381.json` | `8fc82e2ec7a7ef0feacddce68133235b6586722fbbb5e47e64119615df81427a` |
| `debug/data/fx_runctl_posctrl_k4a1_1073741441.json` | `78c7e919cf763247c4db91d928bd5b07e4d5fcd0b26908bf1f08b742c03273e6` |
| `debug/data/fx_runctl_posctrl_k4a1_1073741381.json` | `0b06993cb57b6ccb3056ac010fadfbca2d27724a6bdda56aa741baedcc83c0cb` |
| `debug/data/fx_runctl_e2e_k4a1.json` | `83b21b4a8bbb8c5cfc8a50b4ae8707a656e9d6a6931e44f795b42c086f61604c` |
| `debug/data/fx_runctl_e2e_K4A1_1073741441.json` | `a1a7923256a78d8a20511f808b8c1706839b386d5974668d05a8f0232a6a259b` |
| `debug/data/fx_runctl_e2e_K4A1_1073741381.json` | `f929ef872f50214b26525af5165a0e152a80cdd68af037888453374e05d948e8` |
| `debug/data/fx_runctl_firetest_fx1.json` | `4dd0057aa7ceec8c345467944fc9b58a8673a01486f122da80a8c4851ff11275` |
| `debug/data/fx_runctl_firetest_fx0.json` | `c21c2acc4d8ddca9c6db544b4eed805ff225e31098b885233b01a082a6fc8bf6` |

Settings: LMAX 1500, HOLD 200, RMAX 32, MARGIN 30, EXTRA 2; seeds 2026, 7; primes 1073741441, 1073741381.
Environment: Python 3.14.0, python-flint 0.9.0, mpmath 1.3.0, numpy 2.4.6.

## Phase 2 (after this record is committed and the commit verified)

```
python debug/fx_gate_controls.py V1 V2 V3 V4
python debug/fx_indep.py controls
python debug/fx_run.py ctrl_alphabet
python debug/fx_run.py posctrl {na,ab,k4a1} {1073741441,1073741381}
python debug/fx_run.py e2e_run 1073741441 ; python debug/fx_run.py e2e_run 1073741381 ; python debug/fx_run.py e2e_summary
python debug/fx_run.py score {M1,M2r} {1073741441,1073741381} --registered     (4 jobs in parallel)
python debug/fx_run.py summary                                                  (-> debug/data/fx_run.json)
```
