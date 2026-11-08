"""group2 Batch-4: the Li convergence sequence + range in the figure discussion
(L383-388) -- a parallel to the He sequence, missed by the reviewers' line lists,
caught by the systematic re-read. Li exact-h1 measured: n2 5.00, n3 1.12, n4 1.06.
Range at n>=3 (shown points n3,n4): 1.06--1.12%. Still above the 0.61% HF limit,
so the paragraph's conclusion is unchanged.
Write-first; LaTeX raw strings; idempotent."""
from __future__ import annotations
import sys

P = "papers/group2_quantum_chemistry/paper_fci_atoms.tex"

EDITS = [
    ("li-seq", r"5.00\% ($n_{\max}=2$) $\to$ 1.12\%",
     r"""5.04\% ($n_{\max}=2$) $\to$ 1.15\% ($n_{\max}=3$) $\to$ 1.10\%
($n_{\max}=4$).""",
     r"""5.00\% ($n_{\max}=2$) $\to$ 1.12\% ($n_{\max}=3$) $\to$ 1.06\%
($n_{\max}=4$)."""),

    ("li-range-fig", r"at $n_{\max} \geq 3$ (1.06--1.12\%)",
     r"at $n_{\max} \geq 3$ (1.10--1.15\%)",
     r"at $n_{\max} \geq 3$ (1.06--1.12\%)"),
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
