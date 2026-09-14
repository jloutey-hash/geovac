"""Repair the sentence split by apply_p60_mcentre_scope.py.

The anchor ended mid-sentence ("... basis size.  For"), so the inserted
paragraph stranded a dangling "For" and left "water's $A_1$ block" starting a
paragraph mid-sentence.  Move the "For" to the head of the resumed sentence.
"""
from __future__ import annotations

import sys

P = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"

HEAD_BAD = "or with basis size.  For\n"
HEAD_OK = "or with basis size.\n"
TAIL_BAD = "is not claimed.\n\nwater's $A_1$ block"
TAIL_OK = "is not claimed.\n\nFor water's $A_1$ block"


def main() -> int:
    with open(P, encoding="utf-8") as fh:
        t = fh.read()
    if TAIL_OK in t and HEAD_BAD not in t:
        print("ALREADY REPAIRED")
        return 1
    for label, s in (("head", HEAD_BAD), ("tail", TAIL_BAD)):
        if t.count(s) != 1:
            print(f"{label} anchor count={t.count(s)}; aborting")
            return 2
    t = t.replace(HEAD_BAD, HEAD_OK).replace(TAIL_BAD, TAIL_OK)
    with open(P, "w", encoding="utf-8") as fh:
        fh.write(t)
    print("repaired: dangling 'For' moved to the resumed sentence")
    return 0


if __name__ == "__main__":
    sys.exit(main())
