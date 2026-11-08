"""OWED -- correct the full-shell window mislabel.

/qa paper_60 FULL 2026-09-12, code-B M2.  Re-measured by the PM before editing
(box max(80, 5n^2), 12000 pts, n = 3..10, l_max = n-1):

    K      |M|_1     local slope
     35              0.8492
     56    68.218    0.8670   <- this rung is K=35 -> 56
     84    97.407    0.8785   <- the first slope wholly inside K=56..220
    ...
    220   231.614    0.9111

So 0.867 is the K = 35 -> 56 rung and lies OUTSIDE the stated K=56--220 window;
the first slope inside that window is 0.879, and a global fit over it is 0.893.
The companion exponent is confirmed: full-shell ||T'||_1 fits K^1.0745 over
K >= 56 (local slopes 1.074--1.076), so the superlinear-block claim stands.

The fix widens the stated window to K=35--220 so both printed values keep their
meaning, and states the in-window figures alongside.  This is the paper's own
branch-defining criterion applied to itself: every quantity names its
evaluation domain.

Idempotent.
"""
from __future__ import annotations

import sys

P = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
MARKER = "across the rungs from"

OLD = ("\\emph{total} exponent rises but stays below $1$ ($0.867\\to0.911$ over\n"
       "$K=56$--$220$).")
NEW = ("\\emph{total} exponent rises but stays below $1$ (local slopes\n"
       "$0.867\\to0.911$ across the rungs from $K=35$ to $K=220$;\\ the first\n"
       "slope lying wholly inside $K=56$--$220$ is $0.879$, and a global fit\n"
       "over that window is $0.893$).")


def main() -> int:
    with open(P, encoding="utf-8") as fh:
        t = fh.read()
    if MARKER in t:
        print("ALREADY APPLIED")
        return 1
    if t.count(OLD) != 1:
        print(f"anchor count={t.count(OLD)}; ABORT")
        return 2
    with open(P, "w", encoding="utf-8") as fh:
        fh.write(t.replace(OLD, NEW))
    print("applied: full-shell window corrected to K=35--220, "
          "in-window figures stated")
    return 0


if __name__ == "__main__":
    sys.exit(main())
