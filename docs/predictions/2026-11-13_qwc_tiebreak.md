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
