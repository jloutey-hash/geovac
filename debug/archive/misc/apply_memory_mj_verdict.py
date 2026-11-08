"""Resolve the suspended Avery memory claims against the Monkhorst-Jeziorski read.

Point 2 SURVIVES narrowed.  Point 3 DOES NOT SURVIVE and is retracted as
phrased.  Replaces the 2026-09-12 suspension block rather than appending to it
(CLAUDE.md Sec. 13.11 rule 9: status updates replace, never append).

Idempotent.
"""
from __future__ import annotations

import os
import sys

M = os.path.join(os.environ.get("USERPROFILE", os.path.expanduser("~")),
                 ".claude", "projects",
                 "C--Users-jlout-Desktop-Project-Geometric", "memory",
                 "avery_method_and_prior_art_gaps.md")

MARKER = "RESOLVED 2026-09-12"
OLD_HEAD = "**8. UNREAD SOURCE THAT BEARS DIRECTLY ON POINTS 2 AND 3 (flagged 2026-09-12).**"

NEW = """**8. THE MONKHORST-JEZIORSKI READ -- RESOLVED 2026-09-12. Point 2 survives
narrowed; POINT 3 DOES NOT SURVIVE.**
**Monkhorst & Jeziorski, "No linear dependence or many-center integral problems
in momentum space quantum chemistry", J. Chem. Phys. 71(12), 5268-5269 (1979),
doi:10.1063/1.438337.** Abstract and bibliographic record verified at source;
the two-page body is closed with no repository copy anywhere and is UNREAD.
Full memo: `debug/lit_scan/monkhorst_jeziorski_1979_memo.md`.

- **Point 2 ("no prior art for secular-matrix norm growth") SURVIVES, narrowed.**
  No condition number, spectrum, smallest eigenvalue or growth law is in
  evidence; the reference list carries no numerical-linear-algebra source, and
  no citing paper cites it for a conditioning result. *Honest cap:* this rests
  on the abstract, the reference list and the page count, not on a read of the
  body.
- **Point 3 ("conditioning is a blind spot in the whole Avery canon") IS
  RETRACTED as phrased.** A 1979 *JCP* paper whose abstract names
  "instabilities due to overcompleteness of basis sets" is explicit published
  engagement with the topic. Do not write "blind spot", "nobody studied", or
  "no prior art for the conditioning question" again. The defensible remainder
  is narrower and is point 2: nobody *priced* it -- no growth law, no condition
  number as a function of basis size.

**The resolution, which is worth more than the attribution.** Their claim and
ours are about different objects and both are true. It is the SAME pencil --
same translation-phase matrix, and `sigma_max -> 1` at `p = 0` is a property of
the symbol. The difference is extraction: they let the overlap enter only
multiplicatively (determinantal root search in the scale, or direct
diagonalization), so **nothing is ever inverted** and a near-null direction
gives a harmless spurious branch instead of amplified error; and the basis is
exactly orthonormal in the metric actually used (our own "the intra-center block
is exactly the identity"), so the L2 Gram never enters. **They do not remove the
degeneracy; they keep it out of the denominator.** GeoVac inverts because a
quantum block-encoding wants a standard Hermitian eigenproblem -- so the
conditioning exposure is created by the ENCODING REQUIREMENT, not by the basis
or the metric. Captured in Paper 60 `sec:obstruction`.

**The lever, and its price (OPEN).** Determinantal re-posing would remove the
conditioning multiplier but reinstates the outer nonlinear energy search that
`eq:secular` exists to eliminate. Which is cheaper at scale is an unasked,
answerable resource question.

**Two owed threads from the same read.** (i) The method is credited in the later
literature to **Novosadov**, not to Monkhorst-Jeziorski (via Duchon,
Dumont-Lepage & Gazeau, JCP 76, 445 (1982)) -- unverified, and Novosadov was
still publishing on overcompleteness in 1986, an untouched Russian-language
thread. (ii) The lineage obtains three-centre **one-electron** integrals with no
three-centre evaluation, via a completeness sum over Fock-sphere states. That is
a diagnostic probe candidate and **NOT** a breach of the three-centre wall --
this corpus's wall is the two-electron ERI, and the same source concedes the
advertisement weakens at three or more centres. Related:
[[polyatomic_state_of_play]], [[composition_wall_non_commuting_projections]].

Also owed, lower stakes: the Serra-Capizzano paper this corpus once declined to
cite as "LAA 270 (1998)" is the **wrong reference** -- that is
Tyrtyshnikov-Zamarashkin. The intended ones are LAA **267** (1997) 139-161 and
LAA **282** (1998) 161-183.
"""


def main() -> int:
    if not os.path.exists(M):
        print(f"not found: {M}")
        return 2
    with open(M, encoding="utf-8") as fh:
        t = fh.read()
    if MARKER in t:
        print("ALREADY APPLIED")
        return 1
    i = t.find(OLD_HEAD)
    if i < 0:
        print("suspension block not found")
        return 3
    with open(M, "w", encoding="utf-8") as fh:
        fh.write(t[:i] + NEW)
    print("applied: point 3 retracted, point 2 narrowed, resolution recorded")
    return 0


if __name__ == "__main__":
    sys.exit(main())
