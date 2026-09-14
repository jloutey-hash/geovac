"""group2 Batch-3 remediation: Paper 19 R_eq-drift consistency (reviewer M2).

Two inconsistencies about the same balanced-LiH R_eq drift:
  (number) abstract L79 + Step-4-adjacent L793 say the n2->n3 drift is +0.057
    bohr; the authoritative analytical value is +0.053 (3.227->3.280, stated
    explicitly at L2019: "the drift from n2 to 3 is +0.053 bohr").
  (framing) Step-4 L852 says "+0.053 bohr per n_max step" -- a CONSTANT-rate
    framing -- while the abstract says "decelerating, not a constant rate."
Reconcile all three: n2->3 = +0.053, n3->4 = +0.023, decelerating.

Write-first; LaTeX raw strings; idempotent.
"""
from __future__ import annotations
import sys

P = "papers/group2_quantum_chemistry/paper_19_coupled_composition.tex"

EDITS = [
    ("abs-drift", r"$+0.053$ then $+0.023$~bohr at the two available steps",
     r"$+0.057$ then $+0.023$~bohr at the two available steps",
     r"$+0.053$ then $+0.023$~bohr at the two available steps"),

    ("l793-drift", r"$+0.053$~bohr from $n_{\max}=2$ to 3 and $+0.023$~bohr",
     r"$+0.057$~bohr from $n_{\max}=2$ to 3 and $+0.023$~bohr",
     r"$+0.053$~bohr from $n_{\max}=2$ to 3 and $+0.023$~bohr"),

    ("step4-drift", r"of $+0.053$ then $+0.023$~bohr per step (decelerating, not a constant rate; $3\times$ smaller than PK's",
     r"of $+0.053$~bohr per $n_{\max}$ step ($3\times$ smaller than PK's",
     r"of $+0.053$ then $+0.023$~bohr per step (decelerating, not a constant rate; $3\times$ smaller than PK's"),
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
