# Record AE2b: AE2 re-run at higher precision (precision-only amendment)

*Registered 2026-11-15, before computing. PI-approved after record AE2 came back VOID as registered
(`docs/predictions/2026-11-15_aa_bc_equal_arms_2.md`, Results): its precision gate P0 failed because converting the
exact rational coefficients to floating point loses digits as L grows, and 2000 digits were not enough at
L = 579..581 (2600 held there except E3's PC and PD from L = 581). AE2's VOID stands; this record does not rescore it. A fresh reviewer read the draft cold before commit: no
MUST-FIX; its findings were applied (all clauses at one precision; P0 also certifies the full sums; the precision
estimate reworded; the engine hashes recorded; the disclosure completed).*

## What changes, and what does not

**Changed (precision only):**
- Working precision 3600 and 4400 digits, replacing 2000 and 2600. The exact coefficients are computed once per
  integral and converted at both precisions in the same process (`debug/aaeq2_run.py p1b <case> <name> 600 3600 4400`);
  each precision's output is the same as AE2's `p1` at that precision.
- **P0 widened:** F_579, F_580, F_581 **and F_598, F_599, F_600** of every integral and case agree between 3600 and
  4400 digits to relative 1e-10, **and the P1 sums (to L = 600, both points of both cases, 30 digits printed) agree
  between 3600 and 4400 digits to relative 1e-10** (this certifies every L <= 600 at once, not six sampled L).
  Failure is VOID. (AE2 checked only 579..581, while P1 sums to L = 600.)
- **All clauses (P1, P3, P4) are scored at 3600 digits.**
- Engine files at registration (unchanged since AE's scoring commit 95131691; `git diff` empty), git blob hashes:
  `debug/aapc_exact.py` 4f219a4a, `debug/aabc_exact.py` 40d93108, `debug/na_exact.py` 4841b06f. Re-checked at scoring;
  a change is VOID.

**Unchanged:** cases E3 and E4, the three integrals, the engine (unchanged since AE's scoring commit), AE2-0
(already passed, 24/24; not re-run), the referees and tolerances, every P1/P3/P4 band and probability, the predicted
amplitudes in `debug/data/aaeq2_controls.json`, the "reported, not scored" list, the decision gate, and the caps.
All are as written in the AE2 record, which is the registration this record re-runs.

## Disclosure: this re-run is not blind

AE2's runs produced valid coefficients at L = 220 and 400 (both precisions agree there), and at L = 580 at 2600
digits for S at both cases, PC and PD at E4, and (580 only) PC and PD at E3. The PM has seen them: every decay
exponent within 0.02 of its toy value, PC at E3 at -11.97 / -11.98, and every amplitude ratio between 0.996 and
0.9998 (table in the AE2 record), and F(PC)/F(PD) at E3 falling 3.30-fold from L = 220 to 400 against
(400/220)^2 = 3.31. beta(220) = -11.969 for PC at E3 already makes the drift clause very likely to pass. The predictions were fixed before those were computed and are not changed here.
What this record adds that AE2 did not deliver is a valid P0 gate, the P1 value check against the referees (no AE2
sum was valid), and the L = 580 clauses with a precision certificate.

## Precision estimate (measured from AE2)

The total loss is roughly proportional to L, about 4.5 L digits: 2000 digits are valid at L = 400 and fail by
L = 579; 2600 digits fail from L = 581 at E3 (PC, PD) and near L = 585 elsewhere (2600 / 580 = 4.5). L = 600 needs
roughly 2700 digits even if the loss grows like L log L; 3600 leaves about 800 digits of margin and 4400 is the
check. The gate's sum-agreement clause covers a non-monotone loss at intermediate L.

## Results (scored 2026-11-15)

**Verdict: GO.** Engine files unchanged at scoring (`git diff 95131691` empty). All clauses at 3600 digits.

- **P0 held.** F_579..581 and F_598..600 of all six integrals, and all twelve P1 sums, agree between 3600 and 4400
  digits in every printed digit (25 significant for F, 30 for the sums).
- **P1 held (12/12).** Differences from R2 (n = 34): E3 4.2e-11 to 1.8e-10 (tolerance 1e-9); E4 6.7e-13 to 1.6e-11
  (tolerance 2e-10). R1 (reported): E3 within 4.8e-12, E4 within 7.2e-9.
- **P3 held (6/6).**

  | | beta(220) | beta(400) | beta(580) | band |
  |:--|:--|:--|:--|:--|
  | E3 S | -7.980 | -7.989 | -7.993 | [-8.25, -7.75] |
  | E3 PC | -11.969 | -11.984 | -11.989 | [-12.3, -11.7] at all three; drift 0.020 < 0.05 |
  | E3 PD | -9.974 | -9.987 | -9.991 | [-10.25, -9.75] |
  | E4 S | -7.979 | -7.989 | -7.993 | [-8.25, -7.75] |
  | E4 PC | -9.974 | -9.987 | -9.991 | [-10.25, -9.75] |
  | E4 PD | -9.977 | -9.987 | -9.991 | [-10.25, -9.75] |

  Toy at 580: -7.993, -9.992, -11.990.
- **P4 held (6/6).** F_580 / pred: E3 S 0.99975, PD 0.99960; E4 S 0.99974, PC 0.99961, PD 0.99982 (band [0.98, 1.02]);
  E3 PC 0.99943 against the (1 - x)^{11/2} amplitude from the local expansion (band [0.9, 1.1]). Ratios rise toward 1
  with L (L = 220: 0.996 to 0.999).

**Reported, not scored.** F(PC)/F(PD) at E3: -2.147e-4 (L = 220), -6.511e-5 (400), -3.099e-5 (580); the 220 -> 580
drop is 6.93 against (580/220)^2 = 6.95, the predicted L^-2.

**Reading.** At equal arms the decay exponents, the closed-form amplitudes with the corrected sign rule, and the
predicted cancellation of PC's (1 - x)^{9/2} term at zB D = 9 all hold at two geometries AE never computed. Caveat
as registered: the re-run was not blind at L = 220 / 400 (disclosure above).
