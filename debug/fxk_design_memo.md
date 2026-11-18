# Record FXK design memo (per-key operators after FX's STOP), 2026-11-17

**Status:** instrument built, controls measured, draft record written (`docs/predictions/2026-11-17_fxk_per_key_recurrence.md`). Not registered;
nothing scored. Frozen FX files untouched. Git not touched.

## 1. What FX left, and the per-key question

FX scored STOP on FX3 (`docs/predictions/2026-11-17_four_centre_recurrence.md`). The generic combination
sum_k w_k c_k(L) of the canonical coefficient sequences has no operator at R <= 32 on L <= 1300. FXK applies FX's
own search (`fx_run.minimal_windowed`, unchanged) to each key's sequence c_k separately. It then reads two things:
- F's structure through the key operators: distinct operators, gcrd table, and the LCLM order;
- FX4's question (7/2 at both Landau contacts) from the indicial roots of the key operators that carry each contact.

## 2. A correction to the brief: the sequences were saved

The brief said the coefficient sequences were not saved. They were (FX's `score()` writes
`data/fx_seq_<case>_<p>.json`), and they are now committed (e9604c49), with hashes matching `SEQ_SHA`.
- The scored FXK run needs no 7000 s recompute. It reads these files, and `score` checks their sha256 against values
  hard-coded in `fxk_run.SEQ_SHA`.
- The pilot's recompute of M1 to L = 300 matches the saved prefix.

## 3. Instrument (`debug/fxk_run.py`)

- **Per key.** `fx_run.minimal_windowed` on each in-scope key. The settings are FX's: window L <= 1300, R <= 32,
  margin 30, annihilators at R0..R0+2, every one verified on all rows to L = 1500, minimal operator = gcrd.
  `content_reduce` then gives a fingerprint per operator, so keys sharing an operator are grouped.
- **Sequence span.** The F_p-rank of the key x L matrix, and the classes of keys with proportional sequences.
- **LCLM order** (new code, mod p). This is the dimension of the sum of the key operators' solution spaces. Each
  distinct operator's power-series solution basis is computed at a random regular point (by the recurrence of
  sum c_i(x) d^i, to sum-of-orders + 40 terms), the bases are stacked, and the rank is taken. Two random points must
  agree. The LCLM annihilates F(x) = sum_k key_k f_k(x) and the generic combination, so it is an upper bound for
  their minimal order.
- **Pairwise gcrd orders** use `intclass_lib.gcrd` (skipped above 12 distinct operators).
- **Contacts.** `fx_run.analyse` (FX4's routine) runs on every distinct operator at x_Z*, x_W*, +1 and -1. A
  carrier of a contact is an operator in which that point is present with class M. The output records the union of
  the carriers' half-integer candidate roots, whether some or every carrier has 7/2, and each carrier's lowest
  half-integer root.
- **Gate** (`case_summary`; final after C1 and the cold review; the draft record's gate section owns the wording):
  - VOID: an exception (including a sha256 mismatch), a held-out failure, or a cross-prime inconsistency of any key
    present at a contact (shape, presence or half-integer roots);
  - STOP: no contact present in any found operator, or a contact present (any class) with no 7/2 among the roots of
    the operators present there;
  - GO: FXK2 and FXK3;
  - BORDERLINE: otherwise.
- **Refusals.** `score` refuses without `--registered`. `pilot` refuses L > 300 on a scored case.

## 4. Controls, measured (`debug/data/fxk_ctrl_*.json`)

**Series instrument (`fxk_ctrl_series_<p>`).** Five toy operator sets with known LCLM orders (1, 2, 2, 2, 3) all come
out correct, on both primes and at both series points.

**Positive controls AB and K4A1 (`fxk_ctrl_pos_{ab,k4a1}_<p>`), both primes.**
- AB: A1's 21 frozen pieces, L <= 1000, the whole range as the window.
- K4A1: FX's engine and canonical basis, 19 keys, window 1000 with 200 rows held out.

Results:
- Every key is found and verified.
- The AB shapes equal MO's per-piece shapes exactly.
- Both give the same 7 distinct operators, with identical fingerprints across the two engines: (4,13), (4,26) x 10,
  (5,13), (5,16), (5,28), (6,36) x 2 and (10,76) x 2.
- F_p-rank 13.
- LCLM order 15.
- Generic minimal (14, 111), which is AB's.

**Contacts in the controls.** x* is carried by three operators:
- (5,16), half-integer root 3/2;
- (6,36), 3/2;
- (10,76), 7/2.

So the per-key route does recover 7/2 at x*, but only from the largest key operator. The generic minimal operator
also has 9/2 at x* (MO); with formal keys the 7/2 coefficients cannot cancel, and 9/2 is plausibly an integer shift
of 7/2.

**LCLM versus generic.** The generic order (14) is below the LCLM order (15). Two key operators, (5,28) and (5,13),
do not right-divide the generic minimal operator. Each shares an order-4 gcrd with the (4,26) family. A plausible
reading is that each adds a one-dimensional solution, and the two merge in a generic sum; this is not measured.
Either way, the LCLM order is an upper bound only, and the draft says so.

**Fire tests** (all caught, all four runs):
- **Wrong contact.** Reading at x* + 1/7 gives no carrier and no 7/2.
- **Missing operator.** Drop (10,76); the LCLM of the rest is 13. Adding the generic operator raises it to 15. An
  operator set that misses part of F's solution space is detected by `lclm_closed_under_generic`. The unperturbed
  set passes it (15 -> 15).
- **Gate** (`fxk_ctrl_firetest_gate.json`). Seventeen synthetic run pairs, including the cross-prime VOID cases, each
  reach the intended branch of the final gate.

`fxk_run.py summary` reads FXK0 as holding: 0 failed checks.

## 5. Pilot at M1 (`fxk_ctrl_pilot_M1_1073741441_L300.json`; unscored)

Setup: one prime (1073741441), L <= 300 recomputed (212 s, with the EQ jobs sharing the CPU), window 250 with 50 held
out, R <= 32. The pilot read no singular point and no indicial root, so FXK2 and FXK3 stay unseen.

- **Sequence span.** F_p-rank 31 of 67 keys. Every dependency is a proportionality, giving 31 classes; AB has 13.
- **Found.** 2 of 67 keys: E1(9/5) and E1(18/5) under the prefactor e^{-6/5}. Both are (6, 21) with one shared
  operator, verified on 50 rows.
- **Not found.** For the other 65 keys there is no operator within the window's ceiling. The ceiling is degree 35 at
  R = 4, 26 at R = 6, 16 at R = 10 and 4 at R = 32: about 220 ansatz columns, against about 1250 for the scored
  window. So the pilot cannot price the large keys, and it does not show that they exceed the budget.

**AB's worry** was that its 7/2 carrier was its hardest key. Section 6 tests it on C1.

## 6. C1 calibration (PM-directed; `fxk_ctrl_calib_C1_<p>.json`, both primes)

C1 is FX's unscored matched geometry: (a, b, D1) = (1, 9/4, 7/5), zZ, zW = 3/2, 1. Like M1 it has a < D1 < b. Its
contacts are 37/35 (dominant) and 2809/2520. The run computed sequences through FX's engine to L = 1500, saved as
`fxk_seq_calib_C1_<p>.json` (no FX file is touched), then ran the full scored pipeline. Cost: about 5700 s for the
coefficients and 2230 s for the search, per prime. The results are identical on both primes:

- 77 keys, F_p-rank 32, in 32 proportional classes (M1: 31).
- **14 keys found, covering 10 of 32 classes; 63 not found.** The unfound include all 16 depth-2 keys.
- Shapes: (6,21) x 3 (also the pilot's M1 shape), (10,37), (10,44) x 6, (11,44) x 2, (15,77) x 2. There are 8
  distinct operators. Every pairwise gcrd has order 6 or 10, so the operators share right factors.
- The LCLM of the found operators is 22. With most keys unfound it bounds nothing about F.
- **Each contact is carried by exactly one operator.** Both carriers are (15,77) operators from pure-exponential
  keys (search shape (16,52), 918 of about 1250 columns). At its contact, each carrier's only half-integer
  candidate root is **7/2**, on both primes. There is no 3/2-only carrier.
- +-1 are class L in every operator.

**Consequences for the draft:**
1. A complete key set is out of reach: P(FXK1) about 0.05. Under the first draft's gate (GO needs FXK1), P(GO) would
   be about 0.03, which makes a poor record.
2. The contacts are within reach, and 7/2 reads cleanly. I therefore revised the gate:
   - FXK1 is now reported, not scored.
   - GO = FXK2 (both contacts reached) and FXK3 (7/2 at each). STOP and VOID are as in the gate summary in
     section 3. The cold review accepted this change.
   - Run through the final gate, C1 reads GO (`fxk_ctrl_calib_C1_summary.json`).
3. Probabilities (the arithmetic is in the draft). C1 is n = 1.
   - FXK2: 0.75 at M1, 0.6 at M2r (no measured analogue for M2r).
   - FXK3 given FXK2: 0.9. This clause is low-information; FXK2 carries the information.
   - VOID: 0.03 per case.
   - Record level, with the cases partially correlated: GO about 0.4, STOP about 0.33, BORDERLINE about 0.21, VOID
     about 0.06.
4. Budget: unchanged (LMAX 1500, window 1300, R <= 32). C1's carriers fit with room to spare. The search cost
   measured on C1 is about 37 min per job.

**The cap to keep in view.** A GO shows that each contact is a 7/2-type local exponent of some key's D-module. It
does not show that f_k or F carries a nonzero (1-x)^{7/2} branch on the physical sheet, nor that the branch leads.
The gap between a candidate root and the physical sheet is FC's, not FXK's. Paper 63 may say only "a key operator
has 7/2 among its exponents at each contact".

## 7. The peeling lever (not built)

**What it is.** C1's found operators share right factors: every pair has a gcrd of order 6 or 10. The idea is to take
a found operator L_s, write it in the (x, T) basis the guesser uses, and apply it to a hard key's sequence. That
leaves the sequence of L_s f_k. If L_s right-divides f_k's operator, this sequence satisfies the left quotient, whose
order is R - ord L_s. Its operator may then fit the window. The full operator is recovered by composition, M o L_s,
followed by `analyse`. The exponent reading must be taken on the composed operator: on the quotient, the exponents
are shifted.

**Estimate.**
- Peeling an order-6 factor from a (15,77)-like operator leaves about (9, 60): 600 columns, against 1264.
- It could therefore reach keys up to roughly twice the present column budget, **if** they share a found right
  factor.
- It does not change the revised gate's likely outcome. On C1 the contact carriers were already found without
  peeling.
- What it would change is coverage (FXK1, reported). I estimate it would add a few classes, not complete the set.
  The 16 depth-2 keys are the hardest, and nothing shows that they share the exponential families' right factors;
  that is not measured.
- I do not recommend building it before FXK is scored. It is a separate instrument with its own controls: an AB
  control (peel (4,26) from (10,76) and recover it) and a fire test (a non-dividing factor must fail).

## 8. Recommendation (revised after C1)

**Register with the revised gate and FX's budget unchanged.** The reasons:
- C1, an unscored M1-like geometry, shows the per-key route reaching both contacts within the budget, with 7/2
  cleanly.
- The first draft's GO would have hinged on a complete key set, which C1 shows to be out of reach.
- The cost is under an hour of wall time, with no recompute, provided the PM commits FX's four `fx_seq` files.

If the reviewer rejects the gate change, the honest alternative is to register the first draft with P(GO) of about
0.03 and say so. I recommend against that: the record would be a near-certain BORDERLINE that tests little. M2r is the
weaker case; running C2 (D1 > b) as a second calibration would firm up its 0.7, at about 2.5 h of CPU.

## 9. Files

- `debug/fxk_run.py`: driver. Its sha256 at writing is
  filled at commit by the PM (cold-review revision: pinned dependencies, cross-prime VOID, final gate).
- `docs/predictions/2026-11-17_fxk_per_key_recurrence.md`: draft record.
- `debug/data/fxk_ctrl_series_{1073741441,1073741381}.json`
- `debug/data/fxk_ctrl_pos_{ab,k4a1}_{1073741441,1073741381}.json`, with `.log` files
- `debug/data/fxk_ctrl_firetest_gate.json`
- `debug/data/fxk_ctrl_pilot_M1_1073741441_L300.json`, `fxk_ctrl_pilot_M1.log`
- `debug/data/fxk_seq_pilot_M1_1073741441_L300.json`: the pilot's sequences.
- `debug/data/fxk_seq_calib_C1_{1073741441,1073741381}.json`, `fxk_ctrl_calib_C1_{1073741441,1073741381}.json` (+ `.log`),
  `fxk_ctrl_calib_C1_summary.json`: the C1 calibration.
