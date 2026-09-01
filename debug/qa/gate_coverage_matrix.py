#!/usr/bin/env python
"""Gate x target coverage matrix -- proof that every gate examines something.

WHY
---
The 2026-08-31 trunk run found C19 reporting PASS having examined ZERO trunk
papers, and noted that seven of the ten gates print no file count at all, so
their coverage could not be confirmed from output for ANY target.  A gate that
scopes its verdict away is indistinguishable, from the outside, from a working
one -- so "the deterministic layer is green" was an unverified claim across the
whole sweep, not just for C19.

This script closes that: it runs every gate against every pre-registered
target, parses the coverage figure out of the gate's own output, and FAILS if
any cell reports zero.  It is the standing answer to "did the gate actually
look at anything?", and it should be run before the deterministic layer of any
/qa run is recorded as green.

C10 (check_compiles) is EXCLUDED by default: it invokes pdflatex, which
regenerates the tracked PDFs and dirties the working tree.  Run it separately
with --with-compiles when that is acceptable.

Usage:
    python debug/qa/gate_coverage_matrix.py
    python debug/qa/gate_coverage_matrix.py --with-compiles
"""
from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
import qa_scopes  # noqa: E402

TARGETS = list(qa_scopes.SCOPES)

# gate id -> (script, criterion, regex capturing the coverage count)
GATES = [
    ("C5",  "check_k_label.py",              r"all (\d+) paper\(s\)"),
    ("C11", "check_internal_titles.py",      r"(\d+) paper\(s\) in scope"),
    ("C13", "check_paper_test_refs.py",      r"(\d+) paper\(s\) in scope"),
    ("C14", "check_file_refs.py",            r"(\d+) paper\(s\) in scope"),
    ("C15", "check_inline_arxiv.py",         r"(\d+) paper\(s\) in scope"),
    ("C16", "check_retracted_terms.py",      r"(\d+)/\d+ entries"),
    ("C17", "check_headline_numbers.py",     r"(\d+)/\d+ families"),
    ("C18", "check_duration_language.py",    r"(\d+) paper\(s\) in scope"),
    ("C19", "check_latex_escapes.py",        r"in (\d+) paper\(s\)"),
    ("C20", "check_inline_attributions.py",  r"over (\d+) paper\(s\)"),
    ("C21", "check_numeric_consistency.py",  None),   # count via qa_scopes
]
COMPILES = ("C10", "check_compiles.py", r"(\d+) paper\(s\) in scope")

# Cells that are legitimately empty: the gate is correctly scoped, and the
# registry simply holds no entry for that target.  These are NOT bugs -- but
# they are NOT green either.  "Absence is not compliance" (the /qa
# completeness-critic rule): a criterion the target does not exercise is an
# UNMEASURED criterion, not a passed one, and must be reported that way on the
# cert scorecard rather than counted as a clean gate.
#
# Add a cell here only after confirming the emptiness is a property of the
# registry, not of the scope resolution.
ACKNOWLEDGED_EMPTY = {
    ("paper_59", "C16"): "no retracted claim has a locus in Paper 59",
    ("paper_60", "C16"): "no retracted claim has a locus in Paper 60",
}


def run(script: str, target: str):
    proc = subprocess.run(
        [sys.executable, os.path.join(HERE, script), "--gate", target],
        cwd=ROOT, capture_output=True, text=True, encoding="utf-8",
        errors="replace")
    return proc.returncode, (proc.stdout or "") + (proc.stderr or "")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--with-compiles", action="store_true",
                    help="also run C10 (regenerates tracked PDFs)")
    args = ap.parse_args()

    gates = list(GATES)
    if args.with_compiles:
        gates.append(COMPILES)

    width = max(len(t) for t in TARGETS) + 1
    header = "target".ljust(width) + "".join(g[0].rjust(7) for g in gates)
    print(header)
    print("-" * len(header))

    zero_cells, error_cells, unexercised = [], [], []
    for target in TARGETS:
        row = target.ljust(width)
        for gid, script, pat in gates:
            rc, out = run(script, target)
            if pat is None:
                n = len(qa_scopes.resolve(target)[0])
            else:
                m = re.search(pat, out)
                n = int(m.group(1)) if m else None
            if n is None:
                row += "  ERR!".rjust(7)
                error_cells.append((target, gid))
            elif n == 0 and (target, gid) in ACKNOWLEDGED_EMPTY:
                row += "   -".rjust(7)
                unexercised.append((target, gid))
            elif n == 0:
                row += f" {n}!".rjust(7)
                zero_cells.append((target, gid))
            else:
                row += str(n).rjust(7)
        print(row)

    print()
    if unexercised:
        print(f"UNEXERCISED ({len(unexercised)}) -- correctly scoped, but the "
              f"registry holds nothing for this target.  Report these as "
              f"UNMEASURED on the cert scorecard; do NOT record them as green:")
        for t, g in unexercised:
            print(f"   {g} on '{t}'  ({ACKNOWLEDGED_EMPTY[(t, g)]})")
        print()
    if error_cells:
        print(f"UNPARSEABLE COVERAGE ({len(error_cells)}) -- the gate prints no "
              f"figure this script can read; its coverage is unverifiable:")
        for t, g in error_cells:
            print(f"   {g} on '{t}'")
    if zero_cells:
        print(f"ZERO-COVERAGE CELLS ({len(zero_cells)}) -- gate examined nothing "
              f"while still returning a verdict:")
        for t, g in zero_cells:
            print(f"   {g} on '{t}'")

    if error_cells or zero_cells:
        print("\nRESULT: FAIL (a gate's verdict covers nothing, or cannot be "
              "shown to cover anything)")
        return 1
    print(f"RESULT: PASS (every one of {len(gates)} gate(s) examines a "
          f"non-empty, named scope on all {len(TARGETS)} target(s))")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
