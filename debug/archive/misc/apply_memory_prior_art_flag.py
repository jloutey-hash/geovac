"""Flag the Monkhorst-Jeziorski threat in the auto-loaded Avery memory file.

Items 2 and 3 of that file ("there is NO prior art for secular-matrix norm
growth"; "conditioning is a blind spot in the whole Avery canon") are
load-bearing for Paper 60's novelty framing, are loaded into every session, and
have survived three literature scans.  A 1979 JCP paper titled "No linear
dependence or many-center integral problems in momentum space quantum
chemistry" bears on both by its title alone, and nobody has opened it.

It is flagged rather than acted on -- the honest state is "unread", and writing
either conclusion would be guessing.

Idempotent.
"""
from __future__ import annotations

import os
import sys

M = os.path.join(os.environ.get("USERPROFILE", os.path.expanduser("~")),
                 ".claude", "projects",
                 "C--Users-jlout-Desktop-Project-Geometric", "memory",
                 "avery_method_and_prior_art_gaps.md")

MARKER = "Monkhorst"

NEW = """

**8. UNREAD SOURCE THAT BEARS DIRECTLY ON POINTS 2 AND 3 (flagged 2026-09-12).**
**Monkhorst & Jeziorski, "No linear dependence or many-center integral problems
in momentum space quantum chemistry", J. Chem. Phys. 71, 5268 (1979).** Not
opened. Surfaced as by-catch by the primaries scan
(`debug/lit_scan/primaries_sw1965_toeplitz_pencil_memo.md`). A 1979 title of
that exact shape speaks to *linear dependence* in *momentum-space* quantum
chemistry, which is points 2 and 3 above almost word for word.

**Do not repeat "there is no prior art for the conditioning analysis" or
"conditioning is a blind spot in the whole Avery canon" until it has been
read.** Those claims are load-bearing for Paper 60's novelty framing, they are
loaded into every session from this file, and they have now survived three
scans that never looked at this paper. Points 2 and 3 are NOT retracted --
nothing has been read that contradicts them -- they are suspended pending a
read. Related: [[polyatomic_state_of_play]], [[branch_qa_sweep_phase]].

Also owed from the same scan, lower stakes: the Serra-Capizzano paper this
corpus once declined to cite as "LAA 270 (1998)" is the **wrong reference** --
that is Tyrtyshnikov-Zamarashkin on multilevel Szego distribution. The intended
ones are LAA **267** (1997) 139-161 and LAA **282** (1998) 161-183.
"""


def main() -> int:
    if not os.path.exists(M):
        print(f"memory file not found: {M}")
        return 2
    with open(M, encoding="utf-8") as fh:
        t = fh.read()
    if MARKER in t:
        print("ALREADY APPLIED")
        return 1
    with open(M, "a", encoding="utf-8") as fh:
        fh.write(NEW)
    print("appended: point 8 (Monkhorst-Jeziorski flag)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
