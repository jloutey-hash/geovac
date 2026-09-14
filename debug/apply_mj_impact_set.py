"""Record the Monkhorst-Jeziorski impact set in the sprint memo.

Enumerated from the ARGUMENT before the verdict is known, so the sweep is not
reconstructed from memory afterwards -- which is the failure mode CLAUDE.md
Sec. 9's retraction->dependents rule exists to stop.

Idempotent.
"""
from __future__ import annotations

import sys

M = "debug/sprint_contraction_seam_memo.md"
MARKER = "## 7e. Monkhorst-Jeziorski impact set"

NEW = """
## 7e. Monkhorst-Jeziorski impact set, enumerated BEFORE the verdict

Built while the read was running, so that whichever way it goes the sweep is
already scoped. Four loci carry the exposure, and only two are genuinely at
risk.

| locus | what it says | exposure |
|:--|:--|:--|
| `memory/avery_method_and_prior_art_gaps.md` items 2-3 | "no prior art for secular-matrix norm growth"; "conditioning is a blind spot in the whole Avery canon" | **HIGH** — the direct target; already SUSPENDED pending the read; auto-loads every session |
| Paper 60 `sec:obstruction` (~L186) | shared-scale Coulomb Sturmians "grow linearly dependent" at common `k` | **MEDIUM** — a measured `L^2` fact; needs a scoping clause only if M-J's claim is about a DIFFERENT inner product |
| Paper 60 `sec:quantum` (~L673) | "no existing quantum algorithm combining (i) Sturmian basis, (ii) quantum eigenvalue routine, (iii) isoenergetic inversion" | **LOW** — scoped to *quantum algorithms*; a 1979 classical paper cannot reach it |
| group2 synthesis L737-740 | "measures rather than assumes its cost"; "the sharpest ... basis-independent conditioning" | **LOW** — measurements and a comparative, not novelty claims |

**The resolution the corpus should test first**, because it is the one its own
recent work predicts: GeoVac's linear-dependence statement is about the **L^2**
overlap, while the momentum-space method's natural inner product is the
**V_0-weighted** one — and v5.10.13/v5.10.15 established that the SW matrix IS
`V_0`, with the L^2 metric cancelling identically in the metric-free posing. If
M-J's "no linear dependence" is a statement in the momentum-space metric, both
claims can be true at once and the corpus already owns the reason.

*But that resolution is incomplete as it stands*, and the gap is the
interesting part: GeoVac ALSO measures the SW metric itself — i.e. `V_0` — to be
ill-conditioned at two centres (`cond ~ n^2`). So if M-J are making a
`V_0`-metric claim, either they are single-centre, or their construction avoids
the `p = 0` degeneracy in a way this corpus has not identified. **That second
branch would be a lever rather than a correction**, which is why the read is
worth its cost regardless of the attribution verdict.
"""

ANCHOR = "\n---\n\n## 8. Follow-on items (2026-09-12, PI-directed)\n"


def main() -> int:
    with open(M, encoding="utf-8") as fh:
        t = fh.read()
    if MARKER in t:
        print("ALREADY APPLIED")
        return 1
    if t.count(ANCHOR) != 1:
        print(f"anchor count={t.count(ANCHOR)}")
        return 2
    with open(M, "w", encoding="utf-8") as fh:
        fh.write(t.replace(ANCHOR, "\n" + MARKER + NEW + ANCHOR))
    print("applied: M-J impact set recorded")
    return 0


if __name__ == "__main__":
    sys.exit(main())
