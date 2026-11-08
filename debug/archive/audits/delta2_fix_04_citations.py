"""DELTA #2 -- the citation dimension (5 SMALL, 0 LARGE).

The two NOMINATED fixes were both verified correct at the primary source and
need no change:
  * the GSLW recast -- Theorem 73 IS "Lower bound for eigenvalue transformation"
    and bounds T = applications of U; GSLW themselves pair it with Corollary 67
    in their own proof paragraph, so the chain is the source's, not ours;
  * the Bernstein removal is complete (0 occurrences remain).
And `wulfman_takahata1967` is real -- JCP 47(2), 488-498 (1967) matches
character for character, so the FULL run's drop recommendation was wrong and
overturning it was right.  Recorded because a reviewer finding was overturned
by the PM and has now been independently confirmed a second time.

WHAT IS FIXED HERE

F2 -- the Bernstein removal left the inverse-closedness mechanism asserted with
      NO cite at its locus, and `grochenig_leinert2006` -- which is exactly that
      mechanism, "Symmetry and inverse-closedness of matrix algebras and
      functional calculus for infinite matrices" -- sits in the bibliography
      never cited (verified: 1 occurrence in the whole file, the bibitem).
      This is the precise shape the delta was told to hunt: a mechanism still
      asserted, its attribution removed with the wrong name.
F3 -- `lowdin1950` carries NO title, the only such bibitem in the file; and
      `rokob2008`'s venue is given only as a Festschrift, which a reader cannot
      resolve.  Both corrected from the verified record.
F4 -- the KMS surrender sentence is flat ("it is the KMS asymptotic") while the
      paper establishes 55 lines later that our symbol does NOT meet the
      smoothness hypothesis under which the constant is proved.  Fixed with a
      NAVIGATIONAL pointer only.  Deliberately NOT with the reviewer's suggested
      "what the measurement confirms is its conclusion outside its proven
      hypotheses" -- that phrasing moves credit toward this paper, and the
      reviewer said so itself.  A surrender is not the place to take some back.
F5 -- the `gslw2019` bibitem leads with the 12-page STOC proceedings while
      Theorem 73 and Corollary 67 exist only in the arXiv version.  Both
      resolve (the arXiv ID is in the same bibitem), so this is navigational.

Idempotent.
"""
from __future__ import annotations

import sys

P = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"

EDITS = [
    (P, "F2-grochenig-cite", "inverse-closedness~\\cite{grochenig_leinert2006}",
     "operator is boundedly invertible};\\ inverse-closedness then makes the",
     "operator is boundedly invertible};\\ inverse-closedness~\\cite{grochenig_leinert2006}\nthen makes the"),

    (P, "F3-lowdin-title", "On the non-orthogonality problem",
     "P.-O.~L\\\"owdin, \\textit{J.\\ Chem.\\ Phys.}\\ \\textbf{18}, 365 (1950).",
     "P.-O.~L\\\"owdin, ``On the non-orthogonality problem connected with the use of\natomic wave functions in the theory of molecules and crystals,''\n"
     "\\textit{J.\\ Chem.\\ Phys.}\\ \\textbf{18}, 365 (1950)."),

    (P, "F3-rokob-venue", "Collect.\\ Czech.\\ Chem.\\ Commun.",
     "properties of L\\\"owdin's orthogonalization schemes,'' in\n\\textit{Zahradn\\'ik Festschrift} (2008).",
     "properties of L\\\"owdin's orthogonalization schemes,''\n"
     "\\textit{Collect.\\ Czech.\\ Chem.\\ Commun.}\\ \\textbf{73}, 937 (2008)\n"
     "(Zahradn\\'ik Festschrift issue)."),

    (P, "F5-gslw-numbering", "theorem numbering follows the arXiv version",
     "arithmetics,'' in \\textit{Proc.\\ 51st ACM STOC} (2019), p.~193; arXiv:1806.01838.",
     "arithmetics,'' in \\textit{Proc.\\ 51st ACM STOC} (2019), p.~193; arXiv:1806.01838.\n"
     "(Theorem and corollary numbering here follows the arXiv version;\\ the STOC\n"
     "proceedings text is the 12-page abridgement.)"),

    (P, "F4-kms-scope-pointer", "whose smoothness hypothesis our symbol does not meet",
     "\\textbf{[PRIOR ART]} Eq.~\\eqref{eq:sigma_law} is not new \\emph{as an\nasymptotic law}, and we claim only the identification.",
     "\\textbf{[PRIOR ART]} Eq.~\\eqref{eq:sigma_law} is not new \\emph{as an\nasymptotic law}, and we claim only the identification (whose smoothness\nhypothesis our symbol does not meet --- see below;\\ that is a caveat on the\nproof we can cite, not a claim on our side of the line)."),
]


def main() -> int:
    loaded: dict[str, str] = {}
    applied, skipped, missed = [], [], []
    for path, name, marker, old, new in EDITS:
        if path not in loaded:
            with open(path, encoding="utf-8") as fh:
                loaded[path] = fh.read()
        t = loaded[path]
        if marker in t:
            skipped.append(name)
            continue
        if t.count(old) != 1:
            missed.append((name, t.count(old)))
            continue
        loaded[path] = t.replace(old, new)
        applied.append(name)
    for path, t in loaded.items():
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(t)
    for n in applied:
        print(f"  ok    {n}")
    for n in skipped:
        print(f"  skip  {n} (already applied)")
    for n, c in missed:
        print(f"  MISS  {n}: anchor count={c}")
    print(f"applied {len(applied)}, skipped {len(skipped)}, MISSED {len(missed)}")
    return 3 if missed else 0


if __name__ == "__main__":
    sys.exit(main())
