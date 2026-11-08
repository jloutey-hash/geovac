"""Paper 61 DELTA remediation batch 1: the driver comment + the Chowla-Selberg
bibitem/cites. (The C16 entry is a separate pass -- guard-writing rule.)

F1 (claim-impact NIT) -- `debug/routeC_pslq_fit.py:6` labels K(1/2)=varpi as
   "the lemniscate constant". The value/variable are correct; only the label is
   the retired conflation (the lemniscate constant is sqrt(2)*varpi, a different
   number). The paper fixed this in its own prose but the fix was never
   gate-backed, so it survived here. Corrected.

Citations MATERIAL -- "Chowla-Selberg" is attributed by name at three
   load-bearing loci (abstract [MEASURED], L117, L182 [MEASURED]) with NO
   bibitem and NO \\cite anywhere. The attribution is CORRECT; add the reference
   (A. Selberg and S. Chowla, J. Reine Angew. Math. 227, 86 (1967) -- author
   order and details verified at the de Gruyter/EUDML listing) and cite it at
   the two body loci, which makes the abstract attribution resolvable.

Write-first; LaTeX via raw strings. Idempotent.
"""
from __future__ import annotations

import sys

P = "papers/group3_foundations/paper_61_bessel_moment_periods.tex"
DRV = "debug/routeC_pslq_fit.py"

BIBITEM = (r"""\bibitem{chowla_selberg1967}
A.~Selberg and S.~Chowla, ``On Epstein's zeta-function,''
\textit{J.\ Reine Angew.\ Math.}\ \textbf{227}, 86 (1967).

\bibitem{beilinson_levin1994}""")

EDITS = [
    (DRV, "driver-comment", "this is NOT the lemniscate constant",
     "varpi = Gamma(1/4)^2/(4 sqrt(pi))  [K(1/2), lemniscate constant]",
     "varpi = Gamma(1/4)^2/(4 sqrt(pi))  [= K(1/2); this is NOT the lemniscate "
     "constant, which is sqrt(2)*varpi -- a DIFFERENT number]"),

    (P, "cite-bibitem", r"\bibitem{chowla_selberg1967}",
     r"\bibitem{beilinson_levin1994}",
     BIBITEM),

    (P, "cite-L117", r"Chowla--Selberg~\cite{chowla_selberg1967} periods attached",
     r"are Chowla--Selberg periods attached to the",
     r"are Chowla--Selberg~\cite{chowla_selberg1967} periods attached to the"),

    (P, "cite-L182", r"by Chowla--Selberg~\cite{chowla_selberg1967} \textbf{[MEASURED",
     r"$\Gamma$-values by Chowla--Selberg \textbf{[MEASURED, $\sim\!10^{-51}$;",
     r"$\Gamma$-values by Chowla--Selberg~\cite{chowla_selberg1967} \textbf{[MEASURED, $\sim\!10^{-51}$;"),
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
