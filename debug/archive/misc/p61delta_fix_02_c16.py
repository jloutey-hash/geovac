"""Paper 61 DELTA -- add the C16 entry that makes the K(1/2)=lemniscate
conflation gateable. Separate pass from the fix (guard-writing rule).

Root cause the claim-impact reviewer identified: retired claim #6 (K(1/2) named
the lemniscate constant) was a paper-prose-only fix and got NO
check_retracted_terms.py entry, so it was never swept corpus-wide and survived
in a tracked driver through two prior passes. This closes that.

DISCRIMINATION (proven by the fire-test that follows this applier):
  FIRES on the conflation:  "K(1/2), lemniscate constant" (K(1/2)=varpi labelled
      as the lemniscate constant, with no sqrt2 distinction).
  SILENT on correct usage:  the paper's "the classical lemniscate constant is
      sqrt2 varpi", the drivers' "lemniscate constant = Gamma(1/4)^2/(2 sqrt(2
      pi))", the done.md record "...(they differ by exactly sqrt2...)", and the
      lit-memo "no lemniscate constant" -- all carry sqrt2 / 'differ' / lack the
      adjacent K(1/2) label.

Idempotent.
"""
from __future__ import annotations

import sys

REG = "debug/qa/check_retracted_terms.py"
ANCHOR = '        "id": "p61-every-cm-fibre-universal",'

ENTRY = '''        "id": "p61-k12-is-lemniscate",
        "scope": "paper_61 group3 group2",
        "severity": "fail",
        "retired": "2026-09-07 (/qa paper_61, claims dimension); corpus-swept "
                   "2026-09-13 (DELTA). K(1/2) = varpi = Gamma(1/4)^2/(4 sqrt pi) "
                   "was labelled 'the lemniscate constant'. FALSE: the classical "
                   "lemniscate constant is sqrt(2)*varpi = Gamma(1/4)^2/(2 sqrt(2 "
                   "pi)), a DIFFERENT number (off by exactly sqrt 2). The SYMBOL "
                   "K(1/2) is the corpus convention; the NAME was wrong. Retired "
                   "claim #6 got no C16 entry (paper-prose-only fix), so it "
                   "survived in the tracked driver debug/routeC_pslq_fit.py:6 "
                   "through two passes -- the locus-by-locus-remediation class. "
                   "This entry makes the conflation gateable.",
        "pattern": r"K\\(\\s*1\\s*/?\\s*2\\s*\\)[,;]?\\s*(?:the\\s+)?lemniscate\\s+const",
        "exempt_if_nearby": r"differ|√2|\\\\sqrt|sqrt\\s*\\(?\\s*2|NOT\\s+the\\s+lemniscate|corpus\\s+convention",
        "cited_by": {
            "papers/group3_foundations/paper_61_bessel_moment_periods.tex":
                "reviewed 2026-09-13 -- owner; L113-115 states varpi = K(1/2) is "
                "the corpus symbol and the classical lemniscate constant is sqrt2 "
                "varpi (carries the sqrt2 distinction)",
            "debug/routeC_pslq_fit.py":
                "reviewed 2026-09-13 -- the survivor; comment corrected to '= "
                "K(1/2); this is NOT the lemniscate constant, which is "
                "sqrt(2)*varpi'",
        },
        "files": [
            "papers/group3_foundations/paper_61_bessel_moment_periods.tex",
            "papers/group2_quantum_chemistry/paper_59_elliptic_bessel_moment.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
            "docs/qa/paper_61.done.md",
            "debug/routeC_pslq_fit.py",
            "debug/routeC_cosmic_galois_rung3c.py",
        ],
    },
    {
'''


def main() -> int:
    with open(REG, encoding="utf-8") as fh:
        t = fh.read()
    if '"p61-k12-is-lemniscate"' in t:
        print("already applied")
        return 0
    if t.count(ANCHOR) != 1:
        print(f"  MISS anchor count={t.count(ANCHOR)}")
        return 3
    # insert a new dict right before the every-cm-fibre entry's opening.
    # the anchor is the "id" line; the dict opens with "{\n" on the line above.
    # Replace the "{\n        \"id\": every-cm" boundary: prepend our entry dict.
    needle = '    {\n' + ANCHOR
    if t.count(needle) != 1:
        print(f"  MISS needle count={t.count(needle)}")
        return 3
    replacement = '    {\n' + ENTRY + '        "id": "p61-every-cm-fibre-universal",'
    with open(REG, "w", encoding="utf-8") as fh:
        fh.write(t.replace(needle, replacement, 1))
    print("  ok    p61-k12-is-lemniscate entry inserted")
    return 0


if __name__ == "__main__":
    sys.exit(main())
