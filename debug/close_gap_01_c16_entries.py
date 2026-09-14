"""KNOWN-GAP 1/5 -- the C16 entries owed for the withdrawn readings.

Hard rule (qa.md): "whenever a /qa run retires/withdraws a claim, ADD its phrase
to the REGISTRY ... AND declare the claim's cited_by dependents."  The
2026-09-12 FULL run retired six readings across headlines 15-27 and the
completeness-critic found that NO C16 entry fires on any of them.  My own DoD
extension flagged two as owed; it is six.

Step 1 of 2.  This script adds the standardized withdrawal marker to the ONE
locus where the corrected text quotes the retired wording verbatim -- the
"What this does not license" paragraph, which restates the removability
corollary in order to deny it.  Without the marker, any honest pattern for that
entry fires on the correction, which is the failure mode the registry notes
call out: "a guard that fires on the correction is worse than no guard -- it
makes the right answer unwritable."

Verified before writing: every OTHER retired phrase is already absent from the
paper (0 occurrences of "price of one-cent", "price of complete", "growth law
is derived", "Proposition D", "multiplication by the translation"), so their
entries can use exempt_if_nearby = (?!) and never exempt anything.

Idempotent.
"""
from __future__ import annotations

import sys

P = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
MARKER = "retracted 2026-09-12: p60-removability-corollary"

OLD = """\\emph{What this does not license.}  An earlier version of this paragraph read
the provenance as \\emph{predicting} removability --- the first price being ``a
property of the matrix, which a preconditioner reaches'' and the second ``a
property of the symbol, which no congruence can touch''.  Both halves are
wrong."""

NEW = """\\emph{What this does not license.}  An earlier version of this paragraph read
the provenance as \\emph{predicting} removability --- the first price being ``a
property of the matrix, which a preconditioner reaches'' and the second ``a
property of the symbol, which no congruence can touch''
[retracted 2026-09-12: p60-removability-corollary].  Both halves are
wrong."""


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
    print("applied: withdrawal marker added at the denial locus")
    return 0


if __name__ == "__main__":
    sys.exit(main())
