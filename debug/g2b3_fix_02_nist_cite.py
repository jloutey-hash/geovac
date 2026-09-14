"""group2 Batch-3 remediation: the NIST_ASD misattribution in FCI-atoms.

The EXACT non-relativistic total energies E_He=-2.9037, E_Li=-7.4781,
E_Be=-14.6674 Ha are cited to \\cite{NIST_ASD} -- the NIST Atomic Spectra
Database, a spectroscopic (ionization-energy) compilation that does NOT
tabulate total non-relativistic electronic energies.  The VALUES are correct;
only the source label is wrong.  Correct primary sources (verified 2026-09-13):
  * He -2.903724: C. L. Pekeris, "Ground state of two-electron atoms," Phys.
    Rev. 112, 1649 (1958). [same key Paper 13 already uses: Pekeris1958]
  * Li -7.4781 and Be -14.6674: S. J. Chakravorty, S. R. Gwaltney, E. R.
    Davidson, F. A. Parpia, C. Froese Fischer, "Ground-state correlation
    energies for atomic ions with 3 to 18 electrons," Phys. Rev. A 47, 3649
    (1993) -- covers 3-18 electrons, so Li (3e) and Be (4e) both.
    [https://journals.aps.org/pra/abstract/10.1103/PhysRevA.47.3649]

NIST_ASD is cited ONLY at these two loci, so after repointing it is an orphan;
its bibitem is replaced by the two correct ones.

Write-first; LaTeX raw strings; idempotent.
"""
from __future__ import annotations
import sys

P = "papers/group2_quantum_chemistry/paper_fci_atoms.tex"

EDITS = [
    # He + Li: split the single misattributed cite
    ("he-li-cite", r"\cite{Pekeris1958}",
     r"""$E_{\text{He}} = -2.9037\,\text{Ha}$
and $E_{\text{Li}} = -7.4781\,\text{Ha}$~\cite{NIST_ASD}.""",
     r"""$E_{\text{He}} = -2.9037\,\text{Ha}$~\cite{Pekeris1958}
and $E_{\text{Li}} = -7.4781\,\text{Ha}$~\cite{ChakravortyDavidson1993}."""),

    # Be
    ("be-cite", r"E_{\text{Be}} = -14.6674\,\text{Ha}$~\cite{ChakravortyDavidson1993}",
     r"E_{\text{Be}} = -14.6674\,\text{Ha}$~\cite{NIST_ASD}.}",
     r"E_{\text{Be}} = -14.6674\,\text{Ha}$~\cite{ChakravortyDavidson1993}.}"),

    # bibitem: NIST_ASD -> Pekeris1958 + ChakravortyDavidson1993
    ("bibitem-swap", r"\bibitem{ChakravortyDavidson1993}",
     r"""\bibitem{NIST_ASD}
A.~Kramida, Yu.~Ralchenko, J.~Reader, and NIST ASD Team,
\textit{NIST Atomic Spectra Database} (ver.~5.10),
\url{https://physics.nist.gov/asd} (2023).""",
     r"""\bibitem{Pekeris1958}
C.~L. Pekeris,
``Ground state of two-electron atoms,''
Phys. Rev. \textbf{112}, 1649--1658 (1958).

\bibitem{ChakravortyDavidson1993}
S.~J. Chakravorty, S.~R. Gwaltney, E.~R. Davidson, F.~A. Parpia,
and C.~Froese Fischer,
``Ground-state correlation energies for atomic ions with 3 to 18
electrons,''
Phys. Rev. A \textbf{47}, 3649--3670 (1993)."""),
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
