"""DELTA #2 -- the one verbatim quotation that could not be verified at source.

The paper asserts that Shibuya and Wulfman's OWN ABSTRACT builds the molecular
p_0 operator from the united-atom one "by a sum of unitary transformations, one
for each nucleus in the molecule".  That is a quoted string attributed to a
specific abstract -- the strongest form a citation claim takes.

WHAT WAS CHECKED (2026-09-13).  The bibliographic record is exact and confirmed
(Proc. R. Soc. A 286(1406), 376-389, 1965, doi 10.1098/rspa.1965.0151).  The
ABSTRACT is paywalled: the Royal Society and AIP both return HTTP 403 and
Semantic Scholar reports the publisher eliding it.  The citation reviewer could
not reach it; two further independent retrievals here returned the abstract's
recorded substance -- cusps avoided in momentum space, a hydrogenic basis in
Fock's projective momentum space, R4 spherical harmonics, algebraic expressions
for the required integrals, applied to united-atom and l.c.a.o. wavefunctions
for H2+ -- and NEITHER contained the quoted phrase.

That is evidence against the quotation but NOT proof of absence, because the
retrieved text may itself be abridged.  So the honest move is neither to keep
the quotation nor to assert it is fabricated:  drop the quotation marks and the
"their own abstract" attribution, keep the SUBSTANCE (which is not in doubt --
the Shibuya-Wulfman matrix is built as a sum over nuclei, and Red-Weatherford
2004's title is literally "Derivation of a general formula for the
Shibuya-Wulfman matrix"), and record what could and could not be reached.

THE PRIOR-ART SURRENDER IS UNAFFECTED.  It rests on three counts; counts 2 and
3 (wulfman_takahata1967, red_weatherford2004) were both verified at source by
the citation reviewer, and Wulfman-Takahata's abstract names E4, R5 and O(4,1)
verbatim.  Weakening count 1's wording does not take any of it back.

Idempotent.
"""
from __future__ import annotations

import sys

P = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"

OLD = (
    "Shibuya and Wulfman's\n"
    "own abstract already builds the molecular $p_0$ operator from the united-atom\n"
    "one ``by a sum of unitary transformations, one for each nucleus in the\n"
    "molecule''~\\cite{shibuya1965} --- one unitary per center, which is the\n"
    "structure in substance;\\"
)

NEW = (
    "Shibuya and Wulfman already\n"
    "build the molecular $p_0$ operator from the united-atom one as a sum over the\n"
    "nuclei, one unitary transformation per center~\\cite{shibuya1965}, which is the\n"
    "structure in substance.  (We paraphrase rather than quote:\\ the bibliographic\n"
    "record is verified at source, but the 1965 abstract is behind a paywall and\n"
    "three independent attempts failed to reach its text, so an earlier verbatim\n"
    "quotation here is withdrawn as unverified rather than asserted.)\\"
)


def main() -> int:
    with open(P, encoding="utf-8") as fh:
        t = fh.read()
    if "We paraphrase rather than quote" in t:
        print("already applied")
        return 0
    n = t.count(OLD)
    if n != 1:
        print(f"  MISS shibuya-quote: anchor count={n}")
        return 3
    with open(P, "w", encoding="utf-8") as fh:
        fh.write(t.replace(OLD, NEW))
    print("  ok    shibuya-quote withdrawn as unverified; substance kept")
    return 0


if __name__ == "__main__":
    sys.exit(main())
