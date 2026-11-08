"""KNOWN-GAP 1/5, step 3 -- markers at the three loci the new entries found.

The gate-first order worked exactly as qa.md says it should: register the
class, run it, let it enumerate the loci.  It named three, and ALL THREE are
denials -- text that quotes the retired wording precisely in order to withdraw
it.  None is a re-assertion.

The correct fix is the standardized marker at each locus, NOT a weaker pattern.
Weakening the pattern so it misses a denial would also make it miss the
re-assertion, and the registry's own note records why that is the worse error:
"a guard that fires on the correction is worse than no guard -- it makes the
right answer unwritable, and the cheapest escape is to reword the correction
until the gate goes quiet."

Also flips `p60-frames-completeness` from exempt_if_nearby = (?!) to its marker,
because legitimate denials of that reading demonstrably exist (two of the three
loci).  `p60-removability-corollary` already carried its marker.

Idempotent.
"""
from __future__ import annotations

import sys

EDITS = [
    ("docs/claim_test_matrix.md",
     "Withdraws the 2026-09-11 frames reading** (completeness of the one-centre "
     "set forces lam_min -> 0)",
     "Withdraws the 2026-09-11 frames reading** (completeness of the one-centre "
     "set forces lam_min -> 0) [retracted 2026-09-12: p60-frames-completeness]"),

    ("docs/claim_test_matrix.md",
     'WITHDRAWN** ("a truncation-side price is matrix-level and reachable; a '
     'continuum-side price is symbol-level and is not")',
     'WITHDRAWN** [retracted 2026-09-12: p60-removability-corollary] ("a '
     'truncation-side price is matrix-level and reachable; a continuum-side '
     'price is symbol-level and is not")'),

    ("tests/test_paper60_one_direction.py",
     "completeness of the ONE-CENTRE set forces lam_min -> 0 (if g lies in the closed",
     "completeness of the ONE-CENTRE set forces lam_min -> 0 (if g lies in the closed\n"
     "[retracted 2026-09-12: p60-frames-completeness]"),
]

REG = "debug/qa/check_retracted_terms.py"
REG_OLD = '''                   r"|completeness of the one-cent(?:er|re) set[^.\\n]{0,40}forces",
        "exempt_if_nearby": r"(?!)",'''
REG_NEW = '''                   r"|completeness of the one-cent(?:er|re) set[^.\\n]{0,40}forces",
        "exempt_if_nearby": r"\\[retracted \\d{4}-\\d\\d-\\d\\d: p60-frames-completeness\\]",'''


def main() -> int:
    applied = 0
    for path, old, new in EDITS:
        with open(path, encoding="utf-8") as fh:
            t = fh.read()
        if new[-60:] in t:
            print(f"  skip {path} (marker present)")
            continue
        if t.count(old) != 1:
            print(f"  MISS {path}: count={t.count(old)}")
            return 2
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(t.replace(old, new))
        applied += 1
        print(f"  ok   {path}")

    with open(REG, encoding="utf-8") as fh:
        r = fh.read()
    if "p60-frames-completeness\\]" in r:
        print("  skip registry exemption (already flipped)")
    elif r.count(REG_OLD) == 1:
        with open(REG, "w", encoding="utf-8") as fh:
            fh.write(r.replace(REG_OLD, REG_NEW))
        applied += 1
        print("  ok   registry: frames entry now exempts on its marker")
    else:
        print(f"  MISS registry exemption: count={r.count(REG_OLD)}")
        return 3
    print(f"applied {applied}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
