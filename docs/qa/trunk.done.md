# Trunk — `/qa` profile

> **Inherits the shared criteria in [`docs/qa/criteria.md`](criteria.md).** This
> file supplies only trunk-specific scope + deltas. C1–C13, the verdict rule,
> the review-dimensions map, and the hard rules live in `criteria.md`.

> **STATUS: FROZEN — certified PASS** (`/qa trunk` run #4, 8/8 seeds, v4.16.2).

**Scope:** Papers **0, 1, 7** (group3 foundations) + **32, 38** (group1
operator-algebras) + the **group3 foundations synthesis**. The trunk is the
foundation every branch depends on; it is QA'd first so a finding at a root
re-prices everything above it.

**Deterministic `--gate`:** `trunk` (the default scope for the check scripts).

## Branch deltas (the only non-inherited content)

- **C7 (trunk-dependent status).** The trunk *is* the foundation, so C7 reduces
  to WH1 self-consistency: Paper 38 / WH1 is PROVEN **scoped to the van
  Suijlekom state-space GH distance** (translation-seminorm metrization), with
  no residual "Latrémolière propinquity" overclaim where the proved object is
  the state-space GH distance.
- **C8 (headline honesty), per-paper.** κ = −1/16 is an **Observation**
  (coincidence, no bridge), not a derivation; **4/π** is **derived
  (numerics-pinned)**, not asserted as full symbolic proof; the **Forced-Count**
  moduli chain is stated at its **full-axiom** count, not the matter-sector
  subcount.
- **No branch-specific C14+.**

## Change log
- 2026-06-14 — created (co-authored PM + PI) as the first pre-registered `/qa`
  target.
- 2026-06-14 — added the all-dimensions-mandatory rule after run #3 found
  runs #1–2 had exercised only the claims + citation dimensions; the code
  (C1–C2) and synthesis (C9) dimensions surfaced 2 real synthesis defects
  (κ "derivable", K "conjectural"). Run #4 = **PASS** (8/8 seeds).
- 2026-06-16 — **slimmed to a profile**; C1–C13 + rules moved verbatim to
  `docs/qa/criteria.md` (no criterion changed). PASS status unaffected.
