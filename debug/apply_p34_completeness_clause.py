"""Reconcile candidate (e) with the Sprint-3 "structurally complete" verdict.

That verdict is scoped -- it reads "structurally complete: the master Mellin
engine ... is fully accounted for through M1+M2+M3 at the existing §III
entries".  Candidate (e)'s content is M2 throughout, so what is missing is a
SLOT, not a mechanism, and the two statements do not conflict.  Saying so in
the entry stops a future claims-reviewer from reading a contradiction into it
(and stops a future author from "fixing" the completeness sentence).

Idempotent.
"""
from __future__ import annotations

import sys

P = "papers/group6_precision_observations/paper_34_projection_taxonomy.tex"
MARKER = "a slot rather than a mechanism"

ANCHOR = """Paper~60 \\S~molecular;\\ CHANGELOG
v5.11.4.  \\emph{Caution if promoted:}"""

NEW = """Paper~60 \\S~molecular;\\ CHANGELOG
v5.11.4.  This does not disturb the Sprint~3 completeness verdict above,
which is scoped to the master Mellin engine's \\emph{accounting}:\\ the
content here is M2 on both sides, so what candidate~(e) would add is a
slot rather than a mechanism.  \\emph{Caution if promoted:}"""


def main() -> int:
    with open(P, encoding="utf-8") as fh:
        t = fh.read()
    if MARKER in t:
        print("ALREADY APPLIED")
        return 1
    if t.count(ANCHOR) != 1:
        print(f"anchor count={t.count(ANCHOR)}; aborting")
        return 2
    t = t.replace(ANCHOR, NEW)
    with open(P, "w", encoding="utf-8") as fh:
        fh.write(t)
    print("applied: completeness-scope clause added to candidate (e)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
