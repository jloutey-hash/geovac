"""FULL-run remediation, batch 1: the three content/citation NITs.

F1 (claims PASS, SMALL) -- the abstract tags the third-lever sentence [MEASURED],
   but its closing clause "a direct block-encoding of the preconditioned metric
   then takes the penalty from n^3 to n" is a [RESOURCE MODEL] result in the body
   (the circuit is CITED, not compiled; C8#19). Split the tag: the cond-bounding
   is measured, the n^3->n encoding is the resource model.

C9-NIT (synthesis PASS, NIT) -- "a metric-free standard eigenproblem whose
   eigenvalues are the energies". Strictly the eigenvalues are the scaling
   parameters p_kappa = sqrt(-2E); the energy follows by one algebraic map.
   Tighten so the one-step map is not collapsed.

C4-NIT (citations PASS, NIT) -- the monkhorst_jeziorski1979 bibitem title reads
   "No linear dependence OR MANY-center integral problems"; the published title
   (verified verbatim at the author's own publication list, item #37) is
   "No Linear Dependence AND MULTI-Center Integral Problems in Momentum Space
   Quantum Chemistry". DOI resolves regardless, so this is a NIT, but the exact
   title is now known. Correct it.

Write-first; LaTeX via raw strings (no heredoc). Idempotent.
"""
from __future__ import annotations

import sys

P = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
S = "papers/synthesis/group2_quantum_chemistry_synthesis.tex"

F1_OLD = r"""unrotated frame the same band is $106\times$ \emph{worse} than no treatment at all.  A direct block-encoding of the preconditioned metric then takes the
penalty from $n^3$ to $n$."""
F1_NEW = r"""unrotated frame the same band is $106\times$ \emph{worse} than no treatment at all.  \textbf{[RESOURCE MODEL]} A direct block-encoding of the preconditioned metric then takes the
penalty from $n^3$ to $n$ (the circulant-embedded circuit is analysed, not compiled)."""

S_OLD = r"""reformulation removes the metric \emph{for atoms} entirely---a metric-free
standard eigenproblem whose eigenvalues are the
energies~\cite{loutey_paper60}."""
S_NEW = r"""reformulation removes the metric \emph{for atoms} entirely---a metric-free
standard eigenproblem whose eigenvalues are the scaling parameters
$p_\kappa=\sqrt{-2E}$, from which the energies follow
directly~\cite{loutey_paper60}."""

MJ_OLD = r"""H.~J.~Monkhorst and B.~Jeziorski, ``No linear dependence or many-center integral
problems in momentum space quantum chemistry,''"""
MJ_NEW = r"""H.~J.~Monkhorst and B.~Jeziorski, ``No linear dependence and multi-center integral
problems in momentum space quantum chemistry,''"""

EDITS = [
    (P, "F1-resource-tag", r"circulant-embedded circuit is analysed", F1_OLD, F1_NEW),
    (S, "C9-eigenvalue-nit", r"scaling parameters", S_OLD, S_NEW),
    (P, "C4-mj-title", r"and multi-center integral", MJ_OLD, MJ_NEW),
]


def main() -> int:
    loaded: dict[str, str] = {}
    applied, skipped, missed = [], [], []
    for path, name, marker, old, new in EDITS:
        if path not in loaded:
            with open(path, encoding="utf-8") as fh:
                loaded[path] = fh.read()
        t = loaded[path]
        if marker in t:
            skipped.append(name)
            continue
        if t.count(old) != 1:
            missed.append((name, t.count(old)))
            continue
        loaded[path] = t.replace(old, new)
        applied.append(name)
    for path, t in loaded.items():
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(t)
    for n in applied:
        print(f"  ok    {n}")
    for n in skipped:
        print(f"  skip  {n} (already applied)")
    for n, c in missed:
        print(f"  MISS  {n}: anchor count={c}")
    print(f"applied {len(applied)}, skipped {len(skipped)}, MISSED {len(missed)}")
    return 3 if missed else 0


if __name__ == "__main__":
    sys.exit(main())
