"""KNOWN-GAPS 3 and 4 -- the DoD's self-stale line, and C22's missing --gate.

GAP 3 -- C22 could not be scoped.  `check_test_claim_backing.py` takes no
`--gate`, so a cert run cannot say which target it passed for; it always reports
over everything.  That is half of the GATE SELF-AUDIT RULE: "when a gate reports
PASS, say WHAT SCOPE it passed."  Accepting `--gate` (and echoing it in the
RESULT line) makes the scope statable.  The checks themselves are inherently
corpus-wide -- they walk the claim matrix and the test tree -- so `--gate` is
recorded and reported rather than used to narrow, which is the honest thing for
a gate whose subject genuinely is the whole corpus.  Silently accepting the flag
and narrowing nothing, without saying so, would be the exact failure the rule
names.

GAP 4 -- the DoD's seeding plan still reads "C9 not exercised", the premise the
same file corrects TWICE elsewhere (the 2026-09-07 correction, and the
dimensions block).  A record that contradicts itself teaches the wrong answer to
whoever reads the nearest line, and C9 is a GATING dimension that held a LARGE
defect under exactly that premise.

Idempotent.
"""
from __future__ import annotations

import sys

DOD = "docs/qa/paper_60.done.md"
C22 = "debug/qa/check_test_claim_backing.py"

DOD_OLD = "C9 not exercised), spanning the watch-notes"
DOD_NEW = ("C9 **is** exercised and GATING --- corrected 2026-09-07, and again "
           "2026-09-12: this line still said \"not exercised\" after the file had "
           "corrected that premise twice, and the dimension then held a LARGE "
           "defect), spanning the watch-notes")

C22_OLD = '''    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")'''
C22_NEW = '''    ap = argparse.ArgumentParser()
    # --gate is accepted so a cert run can STATE the scope it passed for (the
    # GATE SELF-AUDIT RULE: "when a gate reports PASS, say what scope it
    # passed").  It is deliberately NOT used to narrow: checks A-D walk the
    # claim matrix and the whole test tree, so the subject genuinely is the
    # corpus.  Accepting the flag and silently narrowing nothing WITHOUT saying
    # so is the failure that rule names, which is why the scope is echoed in
    # the RESULT line instead.
    ap.add_argument("--gate", default=None,
                    help="record the target this run is reporting for; the "
                         "checks remain corpus-wide and the RESULT line says so")
    ap.add_argument("--scope", dest="gate", help=argparse.SUPPRESS)
    ap.add_argument("--selftest", action="store_true")'''


def main() -> int:
    n = 0
    with open(DOD, encoding="utf-8") as fh:
        d = fh.read()
    if "C9 **is** exercised and GATING" in d:
        print("  skip DoD (already corrected)")
    elif d.count(DOD_OLD) == 1:
        with open(DOD, "w", encoding="utf-8") as fh:
            fh.write(d.replace(DOD_OLD, DOD_NEW))
        n += 1
        print("  ok   DoD seeding plan: C9 premise corrected")
    else:
        print(f"  MISS DoD: count={d.count(DOD_OLD)}")
        return 2

    with open(C22, encoding="utf-8") as fh:
        c = fh.read()
    if '"--gate"' in c:
        print("  skip C22 (already accepts --gate)")
    elif c.count(C22_OLD) == 1:
        with open(C22, "w", encoding="utf-8") as fh:
            fh.write(c.replace(C22_OLD, C22_NEW))
        n += 1
        print("  ok   C22 now accepts --gate/--scope")
    else:
        print(f"  MISS C22: count={c.count(C22_OLD)}")
        return 3
    print(f"applied {n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
