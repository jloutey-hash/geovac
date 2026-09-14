"""Capture the Monkhorst-Jeziorski resolution in Paper 60.

The read (debug/lit_scan/monkhorst_jeziorski_1979_memo.md) answers the question
the corpus could not answer itself: how can a 1979 paper assert that
overcompleteness instabilities "do not appear at all" in momentum-space quantum
chemistry while this paper measures cond(S) ~ n^2 for the same metric?

Answer: same pencil, different extraction.  They never invert S.  The corpus
inverts because a quantum block-encoding wants a standard Hermitian
eigenproblem.

Verified by the PM before writing: Crossref record for DOI 10.1063/1.438337 --
title, journal, vol 71, issue 12, pp. 5268-5269, 1979, and the abstract
verbatim, which names "instabilities due to overcompleteness of basis sets" and
"diagonalizations of simple overlap matrices" for a "many-center one electron
system".  NOT verified: the two-page body, which is closed with no repository
copy anywhere; the mechanism paragraph is reconstructed from the lineage and
says so.

Novosadov is deliberately NOT cited here -- the attribution is second-hand and
unverified, so it lives in the memo as an owed thread rather than as an inline
attribution with no bibitem (a C20 defect).

Idempotent.
"""
from __future__ import annotations

import sys

P = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
MARKER = "out of the denominator"

ANCHOR = """\\section{The isoenergetic reformulation is metric-free for atoms}
\\label{sec:iso}
"""

NEW = r"""\textbf{[PRIOR ART]} The classical momentum-space literature met this and
stepped around it, and \emph{how} it did so prices exactly what the quantum
encoding costs.  Monkhorst and Jeziorski~\cite{monkhorst_jeziorski1979} titled a
1979 note ``No linear dependence or many-center integral problems in momentum
space quantum chemistry'', and its abstract states that for a many-center
\emph{one-electron} system the eigenvalues follow from ``diagonalizations of
simple overlap matrices'' and that ``the problems of many-center integrals and
instabilities due to overcompleteness of basis sets do not appear at all''.
Taken at face value that contradicts Eq.~\eqref{eq:sigma_law}.  It does not,
and the reconciliation is worth stating because it locates our own cost.

The degeneracy is not absent from their construction:\ it is the same pencil,
the same matrix of translation phases, and $\sigma_{\max}\to1$ at $p=0$ is a
property of the symbol, indifferent to who is looking.  What differs is the
extraction.  Letting the overlap enter \emph{only multiplicatively} --- solving
a determinantal condition in the scale as a scalar root search, or
diagonalising the overlap directly --- inverts nothing, so a near-null
direction produces a harmless spurious branch instead of amplified error;\ and
the basis is exactly orthonormal in the metric actually used, the statement this
paper measures independently as the intra-center block being the identity, so
the $L^2$ Gram matrix that carries the classical linear-dependence problem never
enters at all.  \textbf{They do not remove the degeneracy;\ they keep it out of
the denominator.}  We invert because a block-encoding wants a standard Hermitian
eigenproblem, and that requirement --- not the basis, and not the metric --- is
what converts a benign spectral feature into a resource multiplier.

The trade runs both ways and we do not claim the better side of it:\ a
determinantal re-posing removes the conditioning multiplier but reinstates the
outer nonlinear search over the scale that Eq.~\eqref{eq:secular} exists to
eliminate, and which of the two is cheaper at scale is an open resource question
this paper does not answer.  \emph{Provenance, because it bounds the paragraph:}\
the bibliographic record and the abstract are verified at source;\ the two-page
body is closed with no repository copy and is unread, so the mechanism above is
reconstructed from the surrounding lineage rather than read.

"""

BIB_ANCHOR = "\\bibitem{shibuya1965}\n"
BIB_NEW = """\\bibitem{monkhorst_jeziorski1979}
H.~J.~Monkhorst and B.~Jeziorski, ``No linear dependence or many-center integral
problems in momentum space quantum chemistry,'' \\textit{J.\\ Chem.\\ Phys.}\\
\\textbf{71}(12), 5268--5269 (1979);\\ doi:10.1063/1.438337.

"""


def main() -> int:
    with open(P, encoding="utf-8") as fh:
        t = fh.read()
    if MARKER in t:
        print("ALREADY APPLIED")
        return 1
    if t.count(ANCHOR) != 1:
        print(f"body anchor count={t.count(ANCHOR)}")
        return 2
    if t.count(BIB_ANCHOR) != 1:
        print(f"bib anchor count={t.count(BIB_ANCHOR)}")
        return 3
    t = t.replace(ANCHOR, NEW + ANCHOR)
    t = t.replace(BIB_ANCHOR, BIB_NEW + BIB_ANCHOR)
    with open(P, "w", encoding="utf-8") as fh:
        fh.write(t)
    print("applied: M-J resolution paragraph + verified bibitem")
    return 0


if __name__ == "__main__":
    sys.exit(main())
