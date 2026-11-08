"""KNOWN-GAP 1/5, step 5 -- fix the p60-prop-d-as-new pattern.

The discrimination proof caught it: the pattern was written `Proposition~?D`,
which requires LaTeX's non-breaking tilde and therefore misses the plain-text
form "Proposition D" -- the most likely way the label would come back, and the
form a docstring or a CHANGELOG entry would use.

Fires 1/2 before, and would have shipped as a guarded class.  This is the second
time today a check written in the same breath as its subject encoded the
author's assumption instead of testing it.

Idempotent.
"""
from __future__ import annotations

import sys

R = "debug/qa/check_retracted_terms.py"
OLD = '        "pattern": r"Proposition~?D\\b"\n'
NEW = '        "pattern": r"Proposition[~\\s]?D\\b"\n'


def main() -> int:
    with open(R, encoding="utf-8") as fh:
        t = fh.read()
    if NEW in t:
        print("ALREADY APPLIED")
        return 1
    if t.count(OLD) != 1:
        print(f"anchor count={t.count(OLD)}; ABORT")
        return 2
    with open(R, "w", encoding="utf-8") as fh:
        fh.write(t.replace(OLD, NEW))
    print("applied: Proposition D pattern now matches the plain-text form")
    return 0


if __name__ == "__main__":
    sys.exit(main())
