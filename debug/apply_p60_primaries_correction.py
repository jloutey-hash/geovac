"""Apply the three corrections from the primaries scan to Paper 60.

All three are conservative (they give claims away or sharpen a caveat), so the
worst case of acting on a relayed read is under-claiming, never over-claiming.

  (1) THE TRANSLATION IDENTIFICATION IS NOT OURS.  Prior art on three counts:
      Shibuya-Wulfman's OWN 1965 abstract builds the molecular p0 operator from
      "a sum of unitary transformations, one for each nucleus in the molecule";
      Wulfman & Takahata 1967 gave the explicit continuous-group formulation;
      Weatherford & Red titled 2002-2004 papers on representing that
      translation operator in a Coulomb-Sturmian basis.  What survives as ours
      is the SYMBOL.  Provenance note: the RS page is 403 from here as it was
      for two prior scans, so the abstract quotation is relayed from the scan's
      direct read of the publisher abstract, recorded in
      debug/lit_scan/primaries_sw1965_toeplitz_pencil_memo.md.  The BODY
      remains unread.

  (2) THE WEIGHT-INDEPENDENCE RESULT IS PRIOR ART, and in the stronger form.
      Ahmad, Al-Aidarous, Alrehaili, Ekstroem, Furci & Serra-Capizzano,
      Numer. Algorithms 78(3) 867-893 (2018) -- Crossref-verified here: title,
      six authors, volume, issue, pages, year -- give it EXACTLY for the tau
      (DST-I) algebra, which is precisely the Toeplitz-minus-Hankel structure
      this paper works in.  The C23 run #2 verdict of ABSENT was reached
      without this source and is superseded.

  (3) THE ~1% RESIDUE IS MOSTLY THE n -> n+1 GRID CONVENTION, not purely an
      O(1/n) asymptotic tail.  Measured here: at n=160, kR=2 the residue moves
      from -0.99% to +0.26% under (n+1)^2, a ~4x reduction with a sign flip.
      An O(1/n) term persists in both, so the sentence says "dominated by",
      not "is".

Adds two bibitems, shibuya1965 (text taken verbatim from Paper 19's existing
entry) and ahmad2018 (Crossref-verified).  Idempotent.
"""
from __future__ import annotations

import sys

P = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
MARKER = "one unitary per center"

OLD1 = """with symbol $j_0(kR\\cot(\\chi/2))$;\\ equivalently that the Shibuya--Wulfman
operator is multiplication by the translation phase
$e^{i\\mathbf{p}\\cdot\\mathbf{R}}$ on the Fock sphere, whose angular average
$j_0(pR)$ becomes trivial at $p=0$."""

NEW1 = """with symbol $j_0(kR\\cot(\\chi/2))$.  It is \\emph{not} the reading of the
Shibuya--Wulfman operator as a translation:\\ that is prior art on three
counts, and earlier versions of this paper claimed it.  Shibuya and Wulfman's
own abstract already builds the molecular $p_0$ operator from the united-atom
one ``by a sum of unitary transformations, one for each nucleus in the
molecule''~\\cite{shibuya1965} --- one unitary per center, which is the
structure in substance;\\ Wulfman and Takahata gave the explicit
continuous-group formulation two years later;\\ and Weatherford and Red titled
papers on representing that operator in a Coulomb--Sturmian basis.  That the
angular average $j_0(pR)$ becomes trivial at $p=0$ is then a property of the
symbol, which is the part we are entitled to."""

OLD2 = """Relatedly, the
$\\sim\\!1\\%$ figure quoted above is the asymptotic's own $O(1/n)$ term rather
than scatter --- the relative residue halves under each doubling of $n$,
falling $0.134\\to0.010$ across $n=10\\to160$ at $kR=2$."""

NEW2 = """Relatedly, the
$\\sim\\!1\\%$ figure quoted above is structure rather than scatter --- the
relative residue halves under each doubling of $n$, falling
$0.134\\to0.010$ across $n=10\\to160$ at $kR=2$ --- and it is \\emph{dominated by
the grid convention}:\\ the $\\tau$-algebra statement places the eigenvalues at
$j\\pi/(n+1)$, and collapsing with $(n+1)^2$ in place of $n^2$ moves the residue
at $n=160$, $kR=2$ from $-0.99\\%$ to $+0.26\\%$, a fourfold reduction with a
change of sign.  A genuine $O(1/n)$ term survives in both conventions."""

OLD3 = """\\textbf{[MEASURED]} Reading the symbol as a quotient yields a prediction, and
it survives its own control.  Let the metric be deformed by any smooth positive"""

NEW3 = """\\textbf{[MEASURED + PRIOR ART]} Reading the symbol as a quotient yields a
prediction, it survives its own control, and it is \\emph{known} --- in a
stronger form than the one measured here.  Ahmad, Al-Aidarous, Alrehaili,
Ekstr\\"om, Furci and Serra-Capizzano~\\cite{ahmad2018} give the eigenvalues of
exactly this preconditioned pencil, and for matrices in the $\\tau$ (DST-I)
algebra --- which is to say Toeplitz minus Hankel, the structure of this
section --- their result is an identity rather than an asymptotic:\\ the
eigenvalues are the ratio symbol sampled on the grid, with the weight
cancelling identically.  So the measurement below confirms a theorem instead of
establishing one, and the paper claims it as such.  Let the metric be deformed
by any smooth positive"""

BIB_ANCHOR = "\\bibitem{batenkov_demanet_goldman_yomdin2019}\n"
BIB_NEW = """\\bibitem{shibuya1965}
T.~I.~Shibuya and C.~E.~Wulfman, ``Molecular orbitals in momentum space,''
\\textit{Proc.\\ Roy.\\ Soc.\\ A} \\textbf{286}, 376 (1965).

\\bibitem{ahmad2018}
F.~Ahmad, E.~S.~Al-Aidarous, D.~A.~Alrehaili, S.-E.~Ekstr\\"om, I.~Furci and
S.~Serra-Capizzano, ``Are the eigenvalues of preconditioned banded symmetric
Toeplitz matrices known in almost closed form?,'' \\textit{Numer.\\ Algorithms}
\\textbf{78}(3), 867--893 (2018);\\ doi:10.1007/s11075-017-0404-z.

"""


def main() -> int:
    with open(P, encoding="utf-8") as fh:
        t = fh.read()
    if MARKER in t:
        print("ALREADY APPLIED")
        return 1
    for name, old in (("translation", OLD1), ("residue", OLD2), ("weights", OLD3)):
        if t.count(old) != 1:
            print(f"  {name} anchor count={t.count(old)}; ABORT")
            return 2
    if t.count(BIB_ANCHOR) != 1:
        print("  bib anchor missing; ABORT")
        return 3
    t = t.replace(OLD1, NEW1).replace(OLD2, NEW2).replace(OLD3, NEW3)
    t = t.replace(BIB_ANCHOR, BIB_NEW + BIB_ANCHOR)
    with open(P, "w", encoding="utf-8") as fh:
        fh.write(t)
    print("applied: 3 corrections + 2 bibitems")
    return 0


if __name__ == "__main__":
    sys.exit(main())
