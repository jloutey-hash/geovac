# Pre-registered prediction: a deterministic tie-break for greedy QWC grouping (QT)

**Registered 2026-11-13**, before any quantity below is computed. PI-directed ("let's do your recommendations",
after v6.7.23 flagged that Paper 20's balanced QWC counts move by 3-12 groups under ~1e-10 coefficient changes).
**Commit status:** committed on its own, before any computation it scores.
**Rests on:** `geovac/measurement_grouping.qwc_groups` (greedy sorted insertion, Verteletskyi et al. 2020).

## Why

`qwc_groups` sorts terms by descending |coefficient| and inserts greedily. Coefficients that are equal by symmetry
differ only by float round-off, so their visiting order -- and the greedy group count -- is set by round-off and by
the order in which the builder happened to accumulate terms. Published counts are reproducible from today's code,
but any change to accumulation order (the v6.7.23 dict route, a scaled transform) can move them.

## The change (to be scored, not assumed)

Sort key: |coefficient| rounded to 10 significant digits (descending), then the Pauli string (ascending). Values
equal to 10 digits are treated as ties and broken by the string, so the order no longer depends on round-off.

## Predictions (scored by `debug/qwc_tiebreak_audit.py`)

- **QT1 (invariance).** Under the new key, every audited system's QWC count is identical for the shipped operator
  and for the same operator with each coefficient multiplied by (1 + 1e-13 u), u uniform in [-1, 1] (seeded).
  Under the old key, at least one of the five Paper 20 balanced systems HCl, H2S, HBr, H2Se, AsH3 changes.
- **QT2 (size of the move).** No published QWC count changes by more than 3% between the old and the new key.
- **QT3 (composed rows).** Every composed library row that publishes 64 groups (LiH, BeH2, H2O, HF, NH3, CH4,
  LiF, CO, N2, F2, NaCl) still gives 64.
- **QT4 (He exponent).** The three-point He QWC exponent (n_max = 2, 3, 4) moves by less than 0.02 from 4.049.

## Result (scored 2026-11-13, after the record was committed)

Driver `debug/qwc_tiebreak_audit.py`, output `debug/data/qwc_tiebreak_audit.json` (43 systems: the 11 composed
library rows, 12 + 3 balanced, coupled LiH, 11 relativistic, He n_max = 2, 3, three STO-3G baselines); He
n_max = 4 by `debug/library_he_fix.py 4 library --qwc`.

- **QT1: MISSED for the registered key; HELD for its replacement.** The 10-significant-digit key failed on
  composed LiH (64 shipped, 69 perturbed): its coefficients are exact binary fractions such as 0.12158203125,
  which sit exactly halfway at the 10th digit, so the rounding split them by round-off. The order shipped instead
  (descending |coefficient|, every run within 1e-9 relative of its largest treated as a tie, ties by Pauli
  string) gives the same count shipped and perturbed on all 43 systems. The old order moves on 38 of them,
  including all five named (HCl 1173 -> 1218, H2S 1264 -> 1342, HBr 1176 -> 1219, H2Se 1264 -> 1342,
  AsH3 1468 -> 1487), and on every composed row (64 -> 76).
- **QT2: MISSED.** Three counts move by more than 3%: PH3 1427 -> 1486 (+4.1%), HCl 1173 -> 1212 (+3.3%),
  MgH2 903 -> 875 (-3.1%). Largest elsewhere: CaH2 -2.8%, HBr +2.6%, BeH2 balanced +2.5%.
- **QT3: HELD.** All 11 composed rows stay at 64.
- **QT4: HELD.** He n_max = 4: 86,251 -> 86,222 groups (-0.03%). The three-point exponent over
  64 / 5,558 / 86,222 is 4.039 (R^2 0.9974, max |log residual| 0.214), from 4.049: a move of 0.009. The
  n_max = 4 invariance run (another ~1 h) was not done; n_max = 2, 3 are invariant.

**Not re-measured.** LiH STO-3G (179, Paper 20 / synthesis / Paper 19): its build script is not on file. The
OpenFermion-cached LiH STO-3G at 1.45 A gives 631 terms but 175 groups under the old order (178 new), so it is
not the source. Paper 19's prolate Route C (2,810) has no script either.

**Verdict:** the deterministic order is in (QT1 held for the shipped order, QT3 and QT4 held); the registered
key was wrong (QT1 as registered) and the moves were larger than predicted (QT2). Every published count above
is updated to the new order in Papers 14, 19, 20 and the group-4 synthesis (CHANGELOG v6.7.24).
