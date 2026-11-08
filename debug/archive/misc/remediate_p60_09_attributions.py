"""OWED ITEM 2 -- the four named-but-uncited attributions.

/qa paper_60 FULL 2026-09-12, citation dimension M1-M4.  All four are works the
paper names in prose with no resolvable bibliography entry, which is invisible
to the automated attribution gate because the prose names authors without an
adjacent year token.

  1. Wulfman & Takahata -- VERIFIED BY THE PM (Crossref 10.1063/1.1711921):
     "Noninvariance Groups in Molecular Quantum Mechanics. I", J. Chem. Phys.
     47(2), 488-498 (1967).  Abstract names "the Lie algebras of E4, R5, and
     O4,1, all noninvariance groups of quantal electrostatics" -- exactly the
     content attributed.  NOTE: the citation reviewer reported this work as
     UNLOCATABLE after three searches and recommended dropping the count and
     re-pricing the paper's novelty concession from three to two.  The work
     exists; the concession stands at three.
  2. Red & Weatherford -- VERIFIED BY THE PM (Crossref 10.1002/qua.20122):
     "Derivation of a general formula for the Shibuya-Wulfman matrix," IJQC
     100, 208-213 (2004).  The paper had the AUTHOR ORDER REVERSED.
  3. Goscinski -- VERIFIED (Crossref 10.1016/S0065-3276(02)41046-5): Adv.
     Quantum Chem. 41, 51-56 (2002); the 1968 Uppsala report No. 217 is
     confirmed from published secondary literature, not read directly.
  4. The "Bernstein Theta(kappa) floor" -- the scan REFRAMED this one, and the
     reframing is the useful part: its published statement is ALREADY in this
     paper's bibliography, as gslw2019's Theorem 73 (their Corollary 67 covers
     x^{-c} for any c > 0, c = 1/2 included, and states the delta and c
     dependence is optimal BY Theorem 73; delta = 1/kappa gives Theta(kappa)).
     So the fix is a pinpoint cite, not a new reference.  The dangling
     back-reference is real and separate: "quoted above" at the first locus has
     NO antecedent -- the only other occurrence is BELOW it.

Idempotent.
"""
from __future__ import annotations

import sys

P = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
MARKER = "wulfman_takahata1967"

EDITS = [
 ("prose-attributions",
  """structure in substance;\\ Wulfman and Takahata gave the explicit
continuous-group formulation two years later;\\ and Weatherford and Red titled
papers on representing that operator in a Coulomb--Sturmian basis.""",
  """structure in substance;\\ Wulfman and Takahata gave the explicit
continuous-group formulation two years later, in terms of the Lie algebras of
$E_4$, $R_5$ and $O(4,1)$~\\cite{wulfman_takahata1967};\\ and Red and
Weatherford derived the general formula for that matrix in a Coulomb--Sturmian
basis~\\cite{red_weatherford2004}."""),

 ("goscinski-cite",
  "Avery~\\cite{avery1989,avery2006,averyphd,averymsc}, built on Goscinski's",
  "Avery~\\cite{avery1989,avery2006,averyphd,averymsc}, built on Goscinski's~\\cite{goscinski2002}"),

 ("bernstein-floor",
  """This escapes the Bernstein $\\Theta(\\kappa)$ floor quoted
above rather than contradicting it:\\ that floor constrains polynomial
approximation of $x^{-1/2}$ on $[\\kappa^{-1},1]$, and preconditioning changes""",
  """This escapes the $\\Theta(\\kappa)$ degree floor for QSVT matrix inversion
---\\ optimal in $\\delta=\\kappa^{-1}$ by Theorem~73
of~\\cite{gslw2019}, whose Corollary~67 covers $x^{-c}$ for every $c>0$ and so
includes $c=\\tfrac12$ ---\\ rather than contradicting it:\\ that floor
constrains polynomial approximation of $x^{-1/2}$ on $[\\kappa^{-1},1]$ under
the QSVT parity and boundedness constraints, and preconditioning changes"""),
]

BIB_ANCHOR = "\\bibitem{shibuya1965}\n"
BIB_NEW = """\\bibitem{wulfman_takahata1967}
C.~E.~Wulfman and Y.~Takahata, ``Noninvariance groups in molecular quantum
mechanics.\\ I,'' \\textit{J.\\ Chem.\\ Phys.}\\ \\textbf{47}(2), 488--498 (1967);\\
doi:10.1063/1.1711921.

\\bibitem{red_weatherford2004}
E.~Red and C.~A.~Weatherford, ``Derivation of a general formula for the
Shibuya--Wulfman matrix,'' \\textit{Int.\\ J.\\ Quantum Chem.}\\ \\textbf{100},
208--213 (2004);\\ doi:10.1002/qua.20122.

\\bibitem{goscinski2002}
O.~Goscinski, ``Conjugate eigenvalue problems and generalized Sturmians,''
\\textit{Adv.\\ Quantum Chem.}\\ \\textbf{41}, 51--56 (2002);\\
doi:10.1016/S0065-3276(02)41046-5.  A presentation of the previously
unpublished \\textit{Conjugate Eigenvalue Problems and the Theory of Upper and
Lower Bounds}, Preliminary Research Report No.~217, Quantum Chemistry Group,
Uppsala University (1968).

"""


def main() -> int:
    with open(P, encoding="utf-8") as fh:
        t = fh.read()
    if MARKER in t:
        print("ALREADY APPLIED")
        return 1
    bad = 0
    for name, old, _ in EDITS:
        if t.count(old) != 1:
            print(f"  {name}: anchor count={t.count(old)}")
            bad += 1
    if t.count(BIB_ANCHOR) != 1:
        print(f"  bib anchor count={t.count(BIB_ANCHOR)}")
        bad += 1
    if bad:
        print("ABORT")
        return 2
    for name, old, new in EDITS:
        t = t.replace(old, new)
        print(f"  ok  {name}")
    t = t.replace(BIB_ANCHOR, BIB_NEW + BIB_ANCHOR)
    with open(P, "w", encoding="utf-8") as fh:
        fh.write(t)
    print("applied: 3 prose fixes + 3 verified bibitems")
    return 0


if __name__ == "__main__":
    sys.exit(main())
