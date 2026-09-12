"""Finish the paper_60 DoD delegation honestly: fix the retired value, delegate
what has an owner, and DECLARE the rest as debt rather than pretending it is done.

A first pass delegated W1, W6, the branch criterion and C8.3/C8.4.  A wider sweep
of the criteria sections found more literals.  They split two ways:

  * W3 still asserts the RETIRED K^0.84.  That is the live defect and must go.
  * ~8 others (Q^3.33/Q^1.19, N^1.97, N^1.85, the tab:resource row, the water
    numbers, H2+/H2 toy figures) are currently CORRECT but unregistered -- no
    C21 key and no C17 family owns them.  Only two C17 families exist for this
    paper: paper60-atomic-sublinear-exponent and paper60-molecular-lambda-exponent.

Delegating an unregistered number to a key that does not exist would leave the
DoD pointing at nothing -- strictly worse than the literal, and exactly the
dangling-pointer failure Sec.9 warns about for debug/ citations.  So: fix the
retired one, delegate the one with an owner, and DECLARE the remainder as dated
debt so the exposure is visible rather than silent.  Registering them is real
work under Sec.15 rule 3 (never register a value you have not measured or cited)
and is recorded as owed, not smuggled in here.
"""
import io

DOD = "docs/qa/paper_60.done.md"
d = io.open(DOD, encoding="utf-8").read()

# --- the live defect: W3 still asserts the retired exponent
OLD_W3 = """- **W3 — sublinear-axis conflation [framing-zombie].** The atomic K^0.84 (config-count
  LCU-λ) and the plane-wave "sublinear in basis size N" (Babbush 2019, Toffoli-in-N) are"""
NEW_W3 = """- **W3 — sublinear-axis conflation [framing-zombie].** The atomic exponent of
  `eq:sublinear` (config-count LCU-λ, C21 `p60_onenorm_exponent`) and the plane-wave
  "sublinear in basis size N" (Babbush 2019, Toffoli-in-N) are"""
assert OLD_W3 in d, "W3 locus not found"
d = d.replace(OLD_W3, NEW_W3, 1)

# --- delegate the one that has an owner
OLD_C8_10 = """    at R_eq) has standard block-encoding λ∼n_orb^2.2 — polynomial, no sublinear behaviour;"""
NEW_C8_10 = """    at R_eq) has a standard POLYNOMIAL block-encoding λ — exponent owned by C17 family
    `paper60-molecular-lambda-exponent` — with no sublinear behaviour;"""
assert OLD_C8_10 in d, "C8.10 locus not found"
d = d.replace(OLD_C8_10, NEW_C8_10, 1)

# --- declare the remainder as dated debt
ANCHOR = "## Paper-60-specific watch-notes (the risk surface — ranked)"
DEBT = """## Un-delegated literals — DECLARED DEBT (2026-09-11)

The criteria sections below still write these numbers as literals, because **no
C21 registry key and no C17 family owns them yet** (this paper has exactly two
C17 families: `paper60-atomic-sublinear-exponent`, `paper60-molecular-lambda-exponent`).
Delegating to a key that does not exist would leave the DoD pointing at nothing —
worse than the literal. They are recorded here so the exposure is **visible and
dated** rather than silent:

| locus | literal | status |
|---|---|---|
| branch criterion 4 | Gaussian ratio-dependence (`N^6` / `N^1.4`) | correct; unregistered |
| W5 | water A1 `cond∼N^1.97` | correct; unregistered |
| W8 | H2+ / H2 toy-validation figures | correct; unregistered |
| C8.1 | naive L2 inflation `Q^3.33` vs `Q^1.19` | correct; unregistered |
| C8.6 | SW `cond(S)∼N^1.85` (vs L2 `N^1.70`) | correct; unregistered |
| C8.7 | the `tab:resource` row | correct; unregistered |
| C8.9 | water gerade-lever-fails numbers | correct; unregistered |

**Every one is currently CORRECT** — none is a retired value. The risk is
structural, not present: each is a literal that will rot the next time its
measurement moves, exactly as `K^0.84` did in W1. **Owed before the next FULL
certifying run:** register each under C21 or a C17 family (Sec.15 rule 3 applies —
measure or cite, never guess), then delegate. Until then a reviewer treats them
as literals and checks them against the body.

"""
assert ANCHOR in d and "DECLARED DEBT" not in d
d = d.replace(ANCHOR, DEBT + ANCHOR, 1)

io.open(DOD, "w", encoding="utf-8").write(d)
print("paper_60.done.md: W3 retired value removed; C8.10 delegated; 7 un-delegated literals declared as dated debt")
