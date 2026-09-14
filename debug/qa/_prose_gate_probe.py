"""Prove check_prose_continuity actually examines the corpus.

Zero advisory findings across three scopes is the same shape as the gates this
corpus has repeatedly caught examining nothing (C19 scoped away, C17 zero
families, C21 zero annotations).  Two checks:

  1. Coverage: how many paragraphs does the analyser actually judge per file?
  2. Discrimination: re-plant the REAL 2026-09-12 defect into a copy of
     Paper 60 and confirm the gate fires on it.
"""
from __future__ import annotations

import os
import shutil
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_prose_continuity import analyse, paragraphs  # noqa: E402

P60 = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"

REPAIRED = ("or with basis size.\n\n\\textbf{[MEASURED]} \\emph{The null space "
            "is geometry-independent;")
SPLICED = ("or with basis size.  For\n\n\\textbf{[MEASURED]} \\emph{The null "
           "space is geometry-independent;")
TAIL_OK = "is not claimed.\n\nFor water's $A_1$ block"
TAIL_BAD = "is not claimed.\n\nwater's $A_1$ block"


def main() -> int:
    with open(P60, encoding="utf-8") as fh:
        text = fh.read()

    n_par = sum(1 for _ in paragraphs(text))
    print(f"1. COVERAGE: {n_par} paragraph blocks parsed from Paper 60")
    if n_par < 50:
        print("   FAIL -- parser is not seeing the document")
        return 1

    u0, l0 = analyse(P60)
    print(f"   current state: {len(u0)} unterminated, {len(l0)} lowercase-open")

    with tempfile.TemporaryDirectory() as d:
        planted = os.path.join(d, "planted.tex")
        bad = text.replace(REPAIRED, SPLICED).replace(TAIL_OK, TAIL_BAD)
        if bad == text:
            print("2. FAIL -- could not re-plant the defect (anchors moved)")
            return 1
        with open(planted, "w", encoding="utf-8") as fh:
            fh.write(bad)
        u2, l2 = analyse(planted)
        print(f"2. DISCRIMINATION: with the real defect re-planted -> "
              f"{len(u2)} unterminated, {len(l2)} lowercase-open")
        for ln, txt in u2:
            print(f"     UNTERMINATED L{ln}: ...{txt!r}")
        for ln, txt in l2:
            print(f"     LOWERCASE-OPEN L{ln}: {txt!r}")
        if not (u2 or l2):
            print("   FAIL -- the gate does NOT catch the defect it was "
                  "written for")
            return 1
    print("\nPROBE: PASS -- parses the document and fires on the real defect")
    return 0


if __name__ == "__main__":
    sys.exit(main())
