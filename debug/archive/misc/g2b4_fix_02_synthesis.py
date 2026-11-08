"""group2 Batch-4 remediation: the C9 synthesis carried the retired pre-ERI-fix
FCI-atoms accuracies (Batch-3 zombie in the synthesis). Update to the re-measured
values (He 0.35->0.26, Li 1.07->1.03, Be 0.90->0.71). The 0.19%@n7 graph-native,
the determinant count, and everything else on the line are current and untouched.
Write-first; LaTeX raw strings; idempotent."""
from __future__ import annotations
import sys

P = "papers/synthesis/group2_quantum_chemistry_synthesis.tex"

EDITS = [
    ("fciatoms", r"$0.26\%$ at $n_{\max}=5$ (grid-based)",
     r"""$0.35\%$ at $n_{\max}=5$ (grid-based) and $0.19\%$ at $n_{\max}=7$ with
exact rational Slater integrals; lithium $1.07\%$ (at $n_{\max}=5$)
and beryllium $0.90\%$ (at $n_{\max}=4$, with $487{,}635$""",
     r"""$0.26\%$ at $n_{\max}=5$ (grid-based) and $0.19\%$ at $n_{\max}=7$ with
exact rational Slater integrals; lithium $1.03\%$ (at $n_{\max}=5$)
and beryllium $0.71\%$ (at $n_{\max}=4$, with $487{,}635$"""),
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
