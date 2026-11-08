"""Scope the M-centre geometry-independence claim in Paper 60.

C23 run #2 (debug/lit_scan/contraction_seam_e3_memo.md) flagged the sentence
"spanned by a fixed vector set that does not move with geometry or with basis
size" plus "a fixed rank-(M-1) rotation removes it" as an over-claim: the null
SPACE is geometry-independent (it is always 1-perp), but the ORDERS at which
its directions open are not.  Re-derived and re-measured locally before editing
(C23 hard rule).

Mechanism, analytic:  j0(p d) = 1 - (p d)^2/6 + O(p^4), so the order-p^2
behaviour on 1-perp is carried by P D2 P with (D2)_ij = d_ij^2.  For COLLINEAR
centres d_ij^2 = (i-j)^2 h^2 = h^2 (i^2 1^T + 1 (j^2)^T - 2 x x^T), and the two
rank-one outer terms are annihilated by P on both sides, leaving
P D2 P = -2 h^2 P x x^T P -- RANK ONE.  So exactly one direction opens at order
2 and the rest at 4, 6, ..., 2(M-1).  Non-collinear geometries generically give
full rank M-1 and every direction opens at order 2.

Measured (mpmath dps=60, two independent p-ratios): collinear M=3 -> (2,4);
collinear M=4 -> (2,4,6); equilateral M=3 -> (2,2); tetrahedral M=4 -> (2,2,2);
bent water-like M=3 -> (2,2).  rank(P D2 P): 1 / 1 / 2 / 3 / 2 respectively.

Adds the verified Batenkov-Demanet-Goldman-Yomdin bibitem (arXiv abstract read
at source 2026-09-12: title, authors and the cluster-size-controls-the-exponent
content all confirmed).

Idempotent.
"""
from __future__ import annotations

import sys

PAPER = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
MARKER = "the rates at which its directions open are not"

ANCHOR = (
    "fixed vector set that does not move with geometry or with basis size.  For\n"
)

INSERT = r"""
\textbf{[MEASURED]} \emph{The null space is geometry-independent;\ the rates at
which its directions open are not, and that scopes the lever.}  Expanding
$j_0(pd)=1-(pd)^2/6+O(p^4)$, the order-$p^2$ behaviour on $\mathbf{1}^\perp$ is
carried by $PD_2P$ with $(D_2)_{ij}=d_{ij}^2$ and $P$ the projector off the
constants.  For \emph{collinear} centers
$d_{ij}^2=(i-j)^2h^2=h^2\bigl(i^2\mathbf{1}^{\!\top}+\mathbf{1}(j^2)^{\!\top}-2xx^{\!\top}\bigr)$,
and $P$ annihilates the two outer terms from both sides, leaving
$PD_2P=-2h^2Pxx^{\!\top}P$ --- \emph{rank one}.  Exactly one of the $M-1$
directions then opens at order $2$ and the rest at $4,6,\dots,2(M-1)$;\ a
non-collinear geometry generically gives $PD_2P$ full rank $M-1$ and every
direction opens at order $2$.  Measured at $\mathrm{dps}=60$ over two
independent $p$-ratios:\ collinear $M=3$ gives orders $(2,4)$ and $M=4$ gives
$(2,4,6)$, against $(2,2)$ for the equilateral triangle, $(2,2,2)$ for the
tetrahedron and $(2,2)$ for a bent water-like triangle, with
$\mathrm{rank}(PD_2P)=1,1,2,3,2$ respectively.  This is the cluster-size
phenomenon of Batenkov, Demanet, Goldman and
Yomdin~\cite{batenkov_demanet_goldman_yomdin2019}, whose exponent is controlled
by the largest cluster and whose $\ell=2$ case is the two-center problem here.
Two consequences, and the second is a limit on what is claimed below.  Water's
$A_1$ block is \emph{bent}, hence in the full-rank class --- which is why one
$\mathrm{tri}(1,2,1)$ on the rotated component suffices and the table below
works.  A \emph{linear} polyatomic is not:\ matching a zero of order $2k$
requires a polynomial with a zero of the same order, so the single
$\mathrm{tri}(1,2,1)$ cures only the one order-$2$ direction there.  The lever
is therefore established for $M=2$ and for non-collinear $M=3$;\ the collinear
case is open and is not claimed.

"""

BIB_ANCHOR = "\\bibitem{barthelme_usevich2021}\n"
BIB = r"""\bibitem{batenkov_demanet_goldman_yomdin2019}
D.~Batenkov, L.~Demanet, G.~Goldman and Y.~Yomdin, ``Conditioning of partial
nonuniform Fourier matrices with clustered nodes,'' arXiv:1809.00658 (2018;
rev.\ 2019).

"""


def main() -> int:
    with open(PAPER, encoding="utf-8") as fh:
        text = fh.read()
    if MARKER in text:
        print("ALREADY APPLIED -- marker present; nothing done.")
        return 1
    for name, anc in (("body", ANCHOR), ("bib", BIB_ANCHOR)):
        if text.count(anc) != 1:
            print(f"{name} ANCHOR not unique (count={text.count(anc)}); aborting.")
            return 2
    text = text.replace(ANCHOR, ANCHOR.rstrip("\n") + "\n" + INSERT)
    text = text.replace(BIB_ANCHOR, BIB + BIB_ANCHOR)
    with open(PAPER, "w", encoding="utf-8") as fh:
        fh.write(text)
    print("applied: M-centre scope clause + 1 verified bibitem")
    return 0


if __name__ == "__main__":
    sys.exit(main())
