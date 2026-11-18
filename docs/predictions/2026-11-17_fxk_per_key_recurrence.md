# On the matched slice, the bond-angle operators of the four-centre 1s integral's canonical keys reach both Landau contacts within FX's budget, with 7/2 an indicial root at each (record FXK)

**Status:** registered 2026-11-17, before scoring. PI-approved as a follow-up to record FX. Nothing scored; the only contact with the scored geometries is the L <= 300 pilot on M1 at one prime, disclosed under Controls. Gate revised after the C1 calibration (unscored geometry, controls measured first) and after a cold review (REGISTRABLE after fixes; all applied, gate fire test 17/17).

**Scope of a GO (read this before the clauses).** A GO shows that each contact is a 7/2-type local exponent of
some key's D-module. It does **not** show that f_k or F carries a nonzero (1-x)^{7/2} branch on the physical sheet,
nor that this branch leads. The only Paper 63 sentence a GO licenses is: **"a key operator has 7/2 among its
exponents at each contact."** The sentence "the contacts keep the 7/2 branch", said of F, is **not** licensed by FXK.
**Owner:** Paper 63 `sec:four_centre` (record FC GO numerically; FX STOP on FX3); inventory row 7, operator column.
**Design:** `debug/fxk_design_memo.md`. **Driver:** `debug/fxk_run.py`.

## Object (reused from FX, unchanged)

FX's matched-slice integral, cases M1 and M2r, and its canonical decomposition G_L = F_L/pi^2 = sum_k c_k(L) key_k
(67 keys at M1, 76 at M2r; FX1 and FX2 held: key set L-independent, inside the alphabet, depth-2 at every L).
For each key k, f_k(x) = sum_L c_k(L) P_L(x). F(x) = pi^2 sum_k key_k f_k(x), with the keys treated as formal
constants (FX's cap).

**FX facts reused, not re-scored:**
- the coefficient sequences themselves: FX's saved `debug/data/fx_seq_<case>_<p>.json`, L <= 1500, both primes
  (sha256 in the Frozen table; written by FX's scored run, under FX0);
- FX1 and FX2 (key set, alphabet, depth 2);
- FX3's negative: the generic combination has no operator at R <= 32 on L <= 1300.

**The question FX left.** The generic combination's minimal operator has at most the order of the LCLM of the key
operators, and at least the largest key order. FXK asks how large the key operators are, and reads the contacts
off the ones that carry them.

## Instrument

`debug/fxk_run.py` imports FX's frozen `fx_run` (minimal_windowed, analyse, generic, landau, coeffs) and changes
nothing in it.

- **Per key:** `fx_run.minimal_windowed` on c_k: NA's bisection on the window L <= 1300, R <= 32, margin 30; every
  annihilator at R0..R0+2 verified on all rows to L = 1500 (200 held out); minimal operator = their gcrd (MO's
  construction). This is exactly FX3's search applied to one key instead of the generic combination.
- **Structure (FXK4, reported):** distinct operators (content-free and monic); the F_p-rank of the key x L matrix
  and the proportional classes; pairwise gcrd orders (if <= 12 distinct operators); the LCLM order, read as the
  rank of the stacked truncated power-series solution bases at two random regular points.
- **Contacts:** for every distinct operator, at x_Z*, x_W*, +1 and -1: present as a factor of the leading
  coefficient, its class (`s1_controls.points`), and the candidate indicial roots (`na_modanalyze.local`), exactly as
  FX4 does. An operator is **present** at a contact if the contact is a factor of its leading coefficient, in any
  class (M, L, A, ...). That is FX4's condition. The gate reads the half-integer candidate roots, never the class. In the control rows below, "carrier" means an operator present with class M; in AB, K4A1 and C1, every present operator was class M.
- `score` refuses to run without `--registered`. It raises (giving VOID) if a sequence file's sha256 differs
  from `SEQ_SHA`, or if `fx_run.py`, `na_modanalyze.py` or `s1_controls.py` differs from the sha256 pinned in
  `DEP_SHA`. Those values are FX's registered values.

Settings: LMAX 1500, HOLD 200, RMAX 32, MARGIN 30, EXTRA 2; primes 1073741441, 1073741381; series points from seeds
11 and 12. The weights' seeds play no role per key.

## Controls (measured; `debug/data/fxk_ctrl_*.json`)

| id | result |
|:--|:--|
| series | LCLM-order instrument on five toy operator sets with known answers (distinct exponentials 2; a repeated operator 1; a contained solution space 2; Legendre with an exponential 3; Legendre alone 2): all correct, both primes, both series points agree |
| pos AB | AB A1 frozen pieces (21, L <= 1000, window = all rows, as FX's AB control), both primes: all 21 found and verified; key shapes equal MO's per-piece shapes exactly; 7 distinct operators (4,13) (4,26) (5,13) (5,16) (5,28) (6,36) (10,76); F_p-rank 13; LCLM order 15 on both series points; generic minimal (14, 111) = AB's, 14 <= 15, and adding the generic operator does not raise the LCLM order; x* carried by three operators, half-integer roots 3/2 (two carriers) and 7/2 (the (10,76) carrier only) |
| pos K4A1 | K4A1 (a = b) through FX's engine and canonical basis, 19 keys, L <= 1200, window 1000 + 200 held out, both primes: all found and verified; same 7 distinct operators as AB (identical fingerprints), same carriers and roots; generic (14, 111), LCLM 15 |
| fire: contact | the contact read at x* + 1/7 instead of x*: no carrier, no 7/2 (caught, all four runs) |
| fire: LCLM | the largest key operator, (10, 76), dropped: LCLM of the rest 13; adding the generic operator raises it to 15 (caught, all four runs). An operator set missing part of F's solution space is detected |
| fire: gate | 17 synthetic run pairs through the final gate (`fxk_ctrl_firetest_gate.json`, `all_pass` true). They give GO for: a single carrier; an incomplete key set; two carriers; present class L with 7/2 (a consistency case only, since `points()` classes a half-integer root as M); a non-carrier key whose shape differs between primes (reported, not gated). They give STOP for: a complete set without 7/2; one contact without 7/2; no contact present; present class L without 7/2; present class A without 7/2. They give BORDERLINE for one contact only. They give VOID for: xZ present on prime 1 only and xW on prime 2 only; a carrier's shape differing between primes; a carrier's roots differing between primes; a carrier found on one prime only; an unverified annihilator; an exception |
| pilot M1 | **unscored, one prime (1073741441), L <= 300, window 250 + 50 held out.** Recomputed sequences equal the first 301 terms of FX's saved file. F_p-rank 31 of 67 keys (31 proportional classes). 2 of 67 keys have an operator in the pilot window: both (6, 21), one operator, verified on 50 rows. For the other 65 there is no operator with R <= 4 and degree <= 35, R <= 6 and degree <= 26, R <= 10 and degree <= 16 (the window's ceiling). No singular point or indicial root was read on M1 |
| **calib C1** | **FX's unscored matched geometry C1**: (a, b, D1) = (1, 9/4, 7/5), zZ, zW = 3/2, 1 (a < D1 < b, as M1); contacts x_Z* = 37/35 (dominant), x_W* = 2809/2520. Sequences through FX's engine to L = 1500 on both primes (about 5700 s each; saved as `fxk_seq_calib_C1_<p>.json`), then the full scored pipeline (window 1300 + 200 held out, R <= 32; search about 2230 s). Identical on both primes: 77 keys, F_p-rank 32 (32 proportional classes). **14 keys found (10 of 32 classes), 63 not**, including all 16 depth-2 keys. Found shapes: (6,21) x 3 keys, (10,37), (10,44) x 6, (11,44) x 2, (15,77) x 2; 8 distinct operators. Pairwise gcrd orders are 6 or 10: the operators share right factors. LCLM of the found operators 22 (both series points); with 63 keys unfound this bounds nothing about F. **Each contact is carried by exactly one operator, a (15,77) from a pure-exponential key (prefactors e^{-307/100}, e^{-739/200}; search shape (16,52)). Its only half-integer candidate root there is 7/2, at both contacts, on both primes.** No other operator is present at either contact, in any class. +-1: class L in all 8 operators. Through the final gate: GO, with no VOID flag (`fxk_ctrl_calib_C1_summary.json`) |

**Readings from the controls, which set the probabilities.**
1. In AB's class the per-key route recovers everything FX4 needs, and 7/2 at x* sits only in the **largest** key
   operator, (10, 76). The small carriers have 3/2 only.
2. **What a root does not show.** A candidate root 7/2 of the indicial polynomial says that the key's D-module
   has a local solution of type (x - x*)^{7/2} there. It does not say that f_k itself, on the physical sheet, has
   a nonzero component along that solution. Nor does it say the component leads: AB's minimal operator also has
   3/2 and 9/2 at x*, and 9/2 is plausibly an integer shift of 7/2. That gap separates FXK3 from FC's numerical
   leading-branch GO, and FXK does not close it.
3. The generic minimal order (14) is **below** the LCLM order (15). So the LCLM order is an upper bound for F's
   operator order, not its value. A plausible mechanism: two order-5 operators each add a one-dimensional
   solution, and these merge in a generic combination. That mechanism is not measured.
4. The largest AB key needs about (11)(78) = 858 ansatz columns. AB's generic search shape (17, 48) needs 900.
   The hardest key was nearly as hard as the generic combination. FX's generic at M1 was not found up to about
   1250 columns at any R <= 32. So it is plausible that M1's largest keys, which are the likely 7/2 carriers, are
   also beyond the window.
5. **C1 (reading 4 tested).** On an M1-like geometry the key set is far from complete within the window: 10 of 32
   classes, and no depth-2 key. But the contacts are reached anyway, each by a (15,77) pure-exponential key with
   7/2 as its only half-integer root. Here the hard keys are not the 7/2 carriers, unlike AB. Neither contact had a
   3/2-only carrier.
6. M1 has 31 independent sequences, C1 32, AB 13. The pilot found (6,21) at M1, a shape C1 also has. The pilot's ceiling (about 220 columns) is far below the
   scored window's (about 1250), so it cannot price the large keys. It shows only that the smallest M1 operator
   found so far, (6, 21), is larger than AB's smallest, (4, 13).

## Predictions (each clause on both M1 and M2r; revised after the C1 calibration)

**Change from the first draft (accepted by the cold review).** The first draft put FXK1 (a complete key set)
inside GO. C1 found 14 of 77 keys, which would put P(GO) near 0.03. FXK1 is now reported, and GO rests on FXK2 and
FXK3. STOP(b) follows FX's own condition: the contact is present, in any class, in some found operator, and no
operator present there has 7/2. That is FX4's own STOP condition ("a contact present without 7/2").

**The change rested on the unscored geometry C1, not on the M1 pilot.** The pilot read no contact and no exponent,
and its window could not reach operators of the carriers' size.

**Which clause carries the information.** FXK2 is the informative clause: does the per-key route reach the
contacts within FX's budget at M1 and M2r? FXK3 is a low-information confirmation, given FC's numerical 7/2, AB's
7/2 carrier and C1's clean 7/2.

- **FXK0 (instrument; VOID on failure).** The following must hold:
  - `fxk_ctrl_series_<p>`, `fxk_ctrl_pos_ab_<p>` and `fxk_ctrl_pos_k4a1_<p>` have `pass` true on both primes,
    including both fire tests;
  - `fxk_ctrl_firetest_gate` has `all_pass`;
  - the four sequence files match the frozen sha256;
  - every annihilator vanishes on all held-out rows, with at least 200 of them;
  - no instrument exception occurs.

  `fxk_run.summary` reads these mechanically.
- **FXK1 (coverage; reported, outside the gate).** Reports the keys and proportional classes that have a minimal
  operator at R <= 32 on the window, on both primes, with the same (order, degree) on both. Expected from C1: a
  minority of the classes (C1: 10 of 32) and no depth-2 key, so the set is not complete (P(complete) about 0.05).
- **FXK2 (contacts reached; P = 0.75 at M1, 0.6 at M2r).** x_Z* and x_W* are each present in some found key
  operator (any class), consistently on both primes.
  - M1's figure comes from C1, where both contacts were reached by (15,77) operators found with room to spare
    (search shape (16,52): 918 of about 1250 columns). C1 is **n = 1**. M1's dominant contact, 61/60, sits closer to
    +1 than C1's 37/35.
  - M2r's geometry (D1 > b) has no measured analogue, since C2 was not run; hence the lower figure.
  - FX4's premise, that the contacts are singular for F, is not re-scored. A contact present in F's operator must
    appear in some key operator (MO3's carriage theorem), but that key operator need not be found.
- **FXK3 (7/2; P = 0.9 given FXK2, each case).** At each contact, some operator present there has 7/2 among its
  half-integer candidate roots, with the same root set on both primes. In C1, 7/2 was the only half-integer root
  of the sole present operator at both contacts. In AB, 7/2 was in one of three present operators.
- **FXK4 (structure; reported, outside the gate).** Reported: per-key (order, degree); the distinct operators; the
  pairwise gcrd orders; the LCLM order, which bounds F's operator order from above, with the largest key order
  below it (if FXK1 holds; otherwise only the LCLM of the found operators, which bounds nothing); classes at +1
  and -1; the half-integer roots and the lowest one per carrier. Expected from C1: an LCLM of the found operators
  of about 20-25, and a largest found order of about 15.

## Decision gate (per case; overall VOID > STOP > BORDERLINE > GO; GO only if both cases are GO)

- **VOID:** any of the following.
  - FXK0 fails.
  - An instrument exception, including a sha256 mismatch of a sequence file or a pinned dependency.
  - **A cross-prime inconsistency at a contact.** For each contact, every key whose operator is present there on
    either prime must be found on both primes, with the same (order, degree), present on both, and with the same
    half-integer root set. Any mismatch is VOID: an instrument problem, not physics. This covers 7/2 found on one
    prime only, and a contact present on one prime only. Shape differences of keys present at no contact are
    reported under FXK1, not gated.
  - A step not finished 48 h after launch. The PM applies this by hand; the driver records no launch time.

  Disagreement of the LCLM's two series points (seeds 11, 12) is reported under FXK4, not gated.
- **STOP:** either of
  - (a) neither contact is present in any found key operator on either prime, so the per-key route reaches none
    of the contact structure within FX's budget;
  - (b) some contact is present, in any class, in a found operator, and no operator present there has 7/2 among
    its half-integer candidate roots. This is FX4's own STOP condition, read through the found keys. While keys
    are unfound it prices rather than refutes, because the 7/2 could sit in an unfound key (AB's reading 1). The
    Reading must say so.
- **GO:** FXK2 and FXK3 hold, and there is no VOID. Scope: as stated in the Status block.
- **BORDERLINE:** otherwise; in practice, one contact reached with 7/2 and the other not reached.

**Probabilities, with the arithmetic.**

P(VOID) per case is 0.03. The cross-prime and instrument path never fired in the 6 control and calibration runs.

| | M1 | M2r |
|:--|:--|:--|
| P(FXK2): both contacts reached | 0.75 | 0.60 |
| P(no contact reached) | 0.10 | 0.20 |
| P(one contact only) | 0.15 | 0.20 |
| P(FXK3 given FXK2) | 0.9 | 0.9 |
| GO = 0.97 x P(FXK2) x 0.9 | 0.65 | 0.52 |
| STOP = 0.97 x [P(none) + P(FXK2) x 0.1] | 0.17 | 0.25 |
| BORDERLINE = 0.97 x P(one only) | 0.15 | 0.19 |
| VOID | 0.03 | 0.03 |

The record needs both cases to GO. Independent cases would give 0.65 x 0.52 = 0.34; perfect correlation would give
0.52. The cases share an engine, a key structure and a budget, but differ in geometry (D1 inside versus outside
[a, b]). I take partial correlation, so the record-level GO is about 0.4.

The other record-level outcomes:
- STOP: about 0.33. That lies between the larger single case (0.25) and independence (1 - 0.83 x 0.75 = 0.38).
- VOID: about 0.06.
- BORDERLINE: about 0.21, the rest.

## Caps

- The caps of FX are inherited: one slice, two geometries, 1s only, mod-p readings, formal constants, and a gcrd
  not certified minimal.
- "Found" means found in the window. An unfound key is beyond the budget, not without an operator (holonomic closure).
- The LCLM order is an upper bound for F's minimal order only when every key is found, and AB shows that the two
  can differ. At the coverage C1 leads us to expect, FXK4's LCLM bounds nothing about F.
- A GO shows that each contact is a 7/2-type local exponent of some key's D-module. It does not show that f_k or F
  carries a nonzero (1-x)^{7/2} branch on the physical sheet, nor that the branch leads (reading 2). Whether 7/2
  leads stays with FC's numerical GO. Paper 63 may say only "a key operator has 7/2 among its exponents at each
  contact".
- The sequences are FX's scored sequences. FXK adds no new geometry, and its contact reading is fresh only in that
  FX never reached an operator.

## Frozen (at registration)

| file | sha256 |
|:--|:--|
| `debug/fxk_run.py` | `e89086ad6b3d5577a4acc27de119a53fc6f5f0a465e9d9ce66623e1080b5e47d` |
| `debug/fx_run.py`, `debug/na_modanalyze.py`, `debug/s1_controls.py` | pinned in `fxk_run.DEP_SHA`; FX's values (`63ba74df…5c6`, `de92dc3a…005`, `1c500166…7a6`) |
| `debug/data/fxk_ctrl_series_1073741441.json` | `c02f5eaa46bedc2c64d618ff9f1db5dd8b2a493a4ca17668ae8f68549dc59ba9` |
| `debug/data/fxk_ctrl_series_1073741381.json` | `c02f5eaa46bedc2c64d618ff9f1db5dd8b2a493a4ca17668ae8f68549dc59ba9` |
| `debug/data/fxk_ctrl_pos_ab_1073741441.json` | `ae585751b932be3c2eb18ae81639257922ec57cc18d1850db4958b2cac6f7531` |
| `debug/data/fxk_ctrl_pos_ab_1073741381.json` | `4e3abe5f106027aa9b6b226c89d73f5f62e2a5f87b698a4294e960e93895356d` |
| `debug/data/fxk_ctrl_pos_k4a1_1073741441.json` | `8a962e9b0811d10d63905339eb3bfca86447209b4af39c37c70175149f84920c` |
| `debug/data/fxk_ctrl_pos_k4a1_1073741381.json` | `065ff8f38e49c9fb98d2e904eae7e51c84d53c8e9a04a9733acf77b6794831ff` |
| `debug/data/fxk_ctrl_firetest_gate.json` | `40e27f6b38817c18fe0c3f4bdb0394858ea378e8b49b4873edbb00604922e37d` |
| `debug/data/fxk_ctrl_pilot_M1_1073741441_L300.json` | `05bf732f77e2da1f3c35932268d921472659bf0326659c8fe4835ab85f3cf445` |
| `debug/data/fxk_seq_pilot_M1_1073741441_L300.json` | `57502ae3ebb018e4851838d50bdd6c5e23155ebe3efd95d33f8ef20ced39fe9c` |
| `debug/data/fxk_seq_calib_C1_1073741441.json` | `73415996067f37cd4ded6aaa3f3727d147dc482c5ff21e2f4987979c2d698b9f` |
| `debug/data/fxk_seq_calib_C1_1073741381.json` | `cc86bc736dda36352a93d20029fc7aba6c1a3fa31794e04b1824c091c4ac98c0` |
| `debug/data/fxk_ctrl_calib_C1_1073741441.json` | `4a78fed0c0d2eac35ee0f13f2ff27d66107adba5acb99ad00db502706d189eb5` |
| `debug/data/fxk_ctrl_calib_C1_1073741381.json` | `973cdd4934b4b586b17dd182cb18ecc40eb53f1a1f2e56d49ad6132f0b1cf466` |
| `debug/data/fxk_ctrl_calib_C1_summary.json` | `b68d9f55fcd7c10634d23804d0cf14cd313530ee98ccff48f1310007e63e739e` |
| `debug/data/fx_seq_M1_1073741441.json` | `edcf845b3186660750e59e595724967ddb97e6903906ecbedf8d9d2ec4cbc3e1` |
| `debug/data/fx_seq_M1_1073741381.json` | `faf2f135acc4f6b7166af44a24734360d36f39dd0448cb0fa58743b8acc0526e` |
| `debug/data/fx_seq_M2r_1073741441.json` | `55b8dfdf882c2a06ff765eb395d95244491f5129063172fa56fe31caa00e8576` |
| `debug/data/fx_seq_M2r_1073741381.json` | `e6374838388e38404cdcf9af7a90a19237a4d50de39d418b06e7f557239d489c` |
| FX's frozen files | unchanged; sha256 as in FX's record |

## Phase 2 (after registration)

```
python debug/fxk_run.py ctrl_series {1073741441,1073741381}
python debug/fxk_run.py ctrl_pos {ab,k4a1} {1073741441,1073741381}
python debug/fxk_run.py firetest_gate
python debug/fxk_run.py score {M1,M2r} {1073741441,1073741381} --registered     (4 jobs in parallel)
python debug/fxk_run.py summary                                                  (-> debug/data/fxk_run.json)
```

Budget: unchanged from FX (LMAX 1500, window 1300, R <= 32). C1 shows the contact carriers inside it with room to
spare. Since the sequences are reused, a longer L means a recompute at about 4x the cost, and C1 does not justify
that.

Cost: the controls take about 4 min each. The C1 calibration is a frozen control and is not rerun in Phase 2. A scored job needs no coefficient recompute (it reads the saved
sequences). Measured on C1, the per-key search over 77 keys took about 2230 s per prime (about 20-40 s per key). That is about
35-40 min per scored job, and under 1 h wall with 4 parallel jobs. If the sequence files are lost, the fallback recompute is about 7000 s per job.
