"""group2 Batch-4 remediation: COMPLETE the FCI-atoms energy re-measurement that
Batch 3 (v5.11.16) left locus-incomplete. Batch 3 fixed the abstract, Table I,
and the convergence-detail table, but the correction did not propagate to the
paper's conclusion, its He monotone sequence, a Li-range sentence, the
hybrid-vs-exact-h1 comparison, and the graph-native comparison paragraph. Caught
by the Batch-4 C9-synthesis and completeness-critic reviewers.

All values re-measured this session (paper's own conventions; He hybrid via
slater_full+hybrid-h1 direct CI at the 5995-determinant convention; Li/Be
exact-h1; He exact-h1 n4=1.99%):
  He hybrid: n2 0.56->0.50, n3 0.45->0.37, n4 0.38->0.29, n5 0.35->0.26
  Li exact-h1: n2 5.04->5.00 ... n5 1.07->1.03 (n4 1.10->1.06)
  Be exact-h1 n4: 0.90->0.71
  He exact-h1 n4: 2.08->1.99

The L976 graph-native paragraph needs a REFRAME, not a number swap: with the grid
(hybrid) He now at 0.26% (Table I) and graph-native at 0.25%, they AGREE to
0.01pp -- the old "0.35% grid vs 0.25% analytical, analytical more accurate"
gap was the wrong-sign-q grid bug, not a real accuracy difference. Post-fix the
two routes confirm each other.

Write-first; LaTeX raw strings; idempotent (marker = new value / new phrase).
"""
from __future__ import annotations
import sys

P = "papers/group2_quantum_chemistry/paper_fci_atoms.tex"

EDITS = [
    # He hybrid monotone convergence sequence
    ("he-seq", r"0.50\% ($n_{\max}=2$) $\to$ 0.37\%",
     r"""0.56\% ($n_{\max}=2$) $\to$ 0.45\% ($n_{\max}=3$) $\to$ 0.38\%
($n_{\max}=4$) $\to$ 0.35\% ($n_{\max}=5$).""",
     r"""0.50\% ($n_{\max}=2$) $\to$ 0.37\% ($n_{\max}=3$) $\to$ 0.29\%
($n_{\max}=4$) $\to$ 0.26\% ($n_{\max}=5$)."""),

    # Li exact-h1 convergence range
    ("li-range", r"converging monotonically from 5.00\% to 1.03\%.",
     r"converging monotonically from 5.04\% to 1.10\%.",
     r"converging monotonically from 5.00\% to 1.03\%."),

    # He hybrid vs exact-h1 at n_max=4 (both re-measured; make it a clean n4 comparison)
    ("he-vs-exact", r"lower errors (0.29\% versus 1.99\% at $n_{\max}=4$ with exact $h_1$)",
     r"lower errors (0.35\% versus 2.08\% at $n_{\max}=4$ with exact $h_1$)",
     r"lower errors (0.29\% versus 1.99\% at $n_{\max}=4$ with exact $h_1$)"),

    # Remaining-basis-error paragraph
    ("rem-error", r"The helium error of 0.26\% at $n_{\max}=5$, lithium error of",
     r"""The helium error of 0.35\% at $n_{\max}=5$, lithium error of
1.07\% at $n_{\max}=5$, and beryllium error of 0.90\% at""",
     r"""The helium error of 0.26\% at $n_{\max}=5$, lithium error of
1.03\% at $n_{\max}=5$, and beryllium error of 0.71\% at"""),

    # Conclusion
    ("conclusion", r"0.26\% accuracy for helium, 1.03\% for lithium, and 0.71\% for",
     r"0.35\% accuracy for helium, 1.07\% for lithium, and 0.90\% for",
     r"0.26\% accuracy for helium, 1.03\% for lithium, and 0.71\% for"),

    # graph-native comparison paragraph -- REFRAME (post-fix the two routes agree)
    ("graph-native-reframe", r"essentially matching the\ngrid-based hybrid-$h_1$ value",
     r"""graph-native achieves 0.25\% error vs.\ the 0.35\% reported in
the original paper, confirming that the analytical integrals are more
accurate than the grid-based numerics (the error decreases to
0.23\% at $n_{\max} = 6$).  At $n_{\max} = 7$: 0.19\%
with 1,218 configurations.""",
     r"""graph-native achieves 0.25\% error, essentially matching the
grid-based hybrid-$h_1$ value (0.26\%, Table~\ref{tab:comparison}) to
within 0.01 percentage points; the analytical integrals confirm the
grid-based numerics rather than improving on them---the earlier
apparent gap was a grid-integration artifact, since corrected by the
exact-rule ERI evaluator.  The graph-native error decreases to
0.23\% at $n_{\max} = 6$, and at $n_{\max} = 7$ reaches 0.19\%
with 1,218 configurations."""),
]


def main() -> int:
    with open(P, encoding="utf-8") as fh:
        t = fh.read()
    applied, skipped, missed = [], [], []
    for name, marker, old, new in EDITS:
        if marker in t:
            skipped.append(name); continue
        if t.count(old) != 1:
            missed.append((name, t.count(old))); continue
        t = t.replace(old, new); applied.append(name)
    with open(P, "w", encoding="utf-8") as fh:
        fh.write(t)
    for n in applied: print(f"  ok    {n}")
    for n in skipped: print(f"  skip  {n} (already applied)")
    for n, c in missed: print(f"  MISS  {n}: anchor count={c}")
    print(f"applied {len(applied)}, skipped {len(skipped)}, MISSED {len(missed)}")
    return 3 if missed else 0


if __name__ == "__main__":
    sys.exit(main())
