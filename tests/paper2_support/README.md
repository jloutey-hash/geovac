# Paper 2 support artifacts (durability migration, 2026-07-04)

`phase2_combinatorial.py` — the Phase-2 combinatorial numerology sweep behind
Paper 2's p-value headline (p = 5.2e-9 over 1.92e9 candidate formulas).
Migrated verbatim from `debug/alpha_audit/` during the group5 certifying
`/qa` run so the cited artifact lives in the permanent record (the same
pattern as `tests/wilson_rule_b_support/`).

Status: **archived-measured, driver preserved**. The full 1.92e9-formula
sweep is far beyond any regression-test budget, so the p-value number is
NOT regression-gated; this script is the reproducibility artifact. A
reduced deterministic-slice pin test is a named follow-up
(docs/qa/group5.done.md).

Not collected by pytest (no `test_` prefix).
