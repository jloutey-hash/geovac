"""Record the owed-items pass (v5.11.6). Idempotent."""
from __future__ import annotations

import sys

CL = "CHANGELOG.md"
CL_ANCHOR = "## [v5.11.5] - 2026-09-12\n"
CL_MARKER = "## [v5.11.6]"

ENTRY = """## [v5.11.6] - 2026-09-12

**The owed items from the `/qa paper_60` FULL run, cleared -- and the recheck found one more wrong number.** PI-directed.

### The recheck found a second stale rung, and a mixed-basis chain

The FULL run confirmed only the K=164 endpoint of the He convergence chain was stale. The other two rungs were flagged as candidates. Recomputed on a converged grid (box 500, 40 000 points), all four at a common `n_max = 10`:

| rung | l_max | K | E (Ha) | paper printed |
|:--|--:|--:|--:|:--|
| 1s^2 | - | 1 | -2.847651 | -2.847 OK |
| s | 0 | 55 | **-2.874468** | **-2.873 WRONG** |
| +p | 1 | 100 | -2.894672 | -2.894 OK |
| spdf | 3 | 164 | -2.896432 | -2.896 OK (fixed earlier) |

The printed `-2.873` is the **n_max = 4 (K=10)** value, so **the chain silently mixed basis sizes** -- against this paper's own branch criterion that every quantity name its basis-growth family. Corrected, all rungs now carry their `K`, the common `n_max` is stated, the truncation convention is declared, and the s value is **registered** (`p60_he_chain_s_k55`). *Grid note worth keeping:* an apparent box drift (-2.8744 -> -2.8739 over boxes 300..1200) is a GRID artifact -- at 40 000 points the value is stable to 2e-5 across the same boxes. Measuring box sensitivity at fixed, insufficient resolution would have produced a third wrong "correction".

### The full-shell family: measured, corrected, and now tested

The abstract's "substantive finding" had **no test, no registry key and no C17 family**. Re-measured (n = 3..10): `||T'||_1 ~ K^1.0745` confirmed superlinear (local slopes 1.057-1.076), total exponent rising 0.849 -> 0.911 and never reaching 1. **The window label was wrong**: `0.867` is the K = 35 -> 56 rung and lies OUTSIDE the stated K=56--220 window, whose first interior slope is `0.879` (global fit `0.893`). Prose corrected; new `tests/test_paper60_full_shell_family.py` (3 tests, slow) asserts the superlinear block, the rising-but-bounded total, and the endpoint windows. Both new guards fire-tested.

### Attributions closed -- and the "missing" reference existed

Three verified bibitems added (`wulfman_takahata1967`, `red_weatherford2004`, `goscinski2002`), each checked at source. Two notes:

- The citation reviewer had reported Wulfman & Takahata as **unlocatable after three searches**, recommending the paper drop that prior-art count and re-price its novelty concession from three to two. **The work exists** (*J. Chem. Phys.* **47**(2), 488-498 (1967)) and its abstract names the Lie algebras of E4, R5 and O(4,1) exactly as attributed. The concession stands at three. The paper also had **Red and Weatherford's author order reversed**.
- The **"Bernstein floor" reframes rather than resolves**: its published statement was already in this paper's bibliography as `gslw2019` Theorem 73 (their Corollary 67 covers `x^{-c}` for every `c > 0`, so `c = 1/2` is included, and states the `delta` dependence is optimal *by* Theorem 73). So the fix was a pinpoint cite, not a new reference. The dangling back-reference was separate and real: "quoted above" had no antecedent, the only other occurrence being below it.

### Instrument fix: the staleness banner measured one file of two

`check_cert_staleness.py` globbed only the paper for single-paper targets, while the trunk and group branches both add their synthesis. C9 is a GATING dimension, so **every single-paper certification on record under-reported drift**. Fixed; `paper_60` now reports 2 changed files where it reported 1.

### The durable half: a reading rule, not another pattern

New **Summary-Surface Reading Rule** in CLAUDE.md Sec. 9. When a claim changes, reread the abstract, conclusion, Scope paragraphs and Acknowledgments in the same edit; read the whole paper once per session that touches claims; a paper's synthesis moves with the paper.

*Why a reading rule.* The corpus has three mechanisms for stale claims -- the C16 phrase registry, `cited_by` dependents, and the Sec. 13.8 `rests on:` edges -- and **all three are document-granular**. Every defect the FULL run found was **locus-granular, inside one file**, so `cited_by` correctly reported no dependents because there were none; and no phrase registry can catch a paraphrase, which is what a summary is by construction.

*Why it is affordable, measured:* abstract + conclusion ~6.2k tokens, whole paper ~29.5k, against ~1.5M for the review pass that found these. **Reading the paper costs about 2% of reviewing it.** Honest limit recorded in the rule itself: roughly half that run's findings were in summary surfaces; the rest were body-text self-contradictions, which are a different failure caught by the reviewers' internal-consistency mandate. The phrase registries stay as a backstop.

### A self-inflicted defect, caught in regression

The new full-shell test called `SV.set_grid` with **no restoring fixture**, so it leaked a module global and broke `test_paper60_split_is_box_sensitive_and_ordering_is_not` -- a test that deliberately pins a box artifact and therefore reads the global it is handed. Passed alone, failed in suite. Fixed with the same fixture the resource-ladder file uses. Worth recording because it is the same class this whole arc is about: a change whose effect on its neighbours was not checked.

### Gates

Deterministic layer 14/14 in scope `paper_60`; group2 compiles; 141 tests pass with `--slow`.

"""


def main() -> int:
    with open(CL, encoding="utf-8") as fh:
        t = fh.read()
    if CL_MARKER in t:
        print("ALREADY APPLIED")
        return 1
    if t.count(CL_ANCHOR) != 1:
        print(f"anchor count={t.count(CL_ANCHOR)}; ABORT")
        return 2
    with open(CL, "w", encoding="utf-8") as fh:
        fh.write(t.replace(CL_ANCHOR, ENTRY + CL_ANCHOR))
    print("applied: v5.11.6 CHANGELOG entry")
    return 0


if __name__ == "__main__":
    sys.exit(main())
