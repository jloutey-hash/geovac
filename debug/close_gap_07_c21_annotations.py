"""KNOWN-GAP 2/5 -- annotate the three registry keys cited from no .tex.

The completeness-critic found that `p60_weighted_collapse_control`,
`p60_window_richardson_pi2` and `p60_window_rms_richardson_pi` -- all registered
on 2026-09-12 -- appear in NO `.tex` at all.  A registered key nothing cites is
the mirror of an unregistered literal: C21 recomputes and stores it, and then
has nothing in the paper to check it against, so the paper's copy can drift
freely while the gate reports PASS.

Wrapping a literal in \\gvq is a no-op by construction (the macro renders the
literal only), so this changes no rendered text -- verified by the C10 rebuild
and by C21's own annotation check.

NOT attempted here, and recorded as standing debt rather than silently dropped:
the critic also counted ~141 unannotated decimal literals in the molecular half
beyond these three.  Annotating those is a corpus-wide maintenance path, not
part of this gap; the three registered-but-uncited keys are the actual defect,
because they are the ones where the registry believes it has coverage.

Idempotent.
"""
from __future__ import annotations

import sys

P = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
MARKER = "gvq{p60_window_richardson_pi2}"

EDITS = [
    ("(Richardson $9.86949$ against $\\pi^2=9.86960$)",
     "(Richardson $\\gvq{p60_window_richardson_pi2}{9.86949}$ against "
     "$\\pi^2=9.86960$)"),
    ("(Richardson $3.14158$)",
     "(Richardson $\\gvq{p60_window_rms_richardson_pi}{3.14158}$)"),
    ("constant to $0.828$ while leaving the exponent at $-1.97$",
     "constant to $\\gvq{p60_weighted_collapse_control}{0.828}$ while leaving "
     "the exponent at $-1.97$"),
]


def main() -> int:
    with open(P, encoding="utf-8") as fh:
        t = fh.read()
    if MARKER in t:
        print("ALREADY APPLIED")
        return 1
    for i, (old, _) in enumerate(EDITS, 1):
        if t.count(old) != 1:
            print(f"  edit {i} anchor count={t.count(old)}; ABORT")
            return 2
    for old, new in EDITS:
        t = t.replace(old, new)
    with open(P, "w", encoding="utf-8") as fh:
        fh.write(t)
    print("applied: 3 registered keys now annotated at their loci")
    return 0


if __name__ == "__main__":
    sys.exit(main())
