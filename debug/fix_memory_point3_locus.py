"""Mark point 3 AT ITS OWN LOCUS, not only in the resolution block below it.

The retraction was written as item 8 at the foot of the file while point 3
stayed live at line 34 -- so anything auto-loading this file reads the retracted
claim first and the retraction 170 lines later, if at all.  That is precisely
the corrected-in-the-owner-left-standing-in-the-citer class CLAUDE.md Sec. 9
exists to stop, committed inside the very edit that was fixing it.

Idempotent.
"""
from __future__ import annotations

import os
import sys

M = os.path.join(os.environ.get("USERPROFILE", os.path.expanduser("~")),
                 ".claude", "projects",
                 "C--Users-jlout-Desktop-Project-Geometric", "memory",
                 "avery_method_and_prior_art_gaps.md")

MARKER = "RETRACTED 2026-09-12"

OLD3 = """**3. Conditioning is a blind spot in the whole Avery canon.** The isoenergetic
framing is championed because fixing the energy gives every basis function the
correct asymptotic decay at the classical turning points — solving the continuum
problem without a continuum. The *linear-algebraic* consequences (condition
number, norm bounds as the basis expands) are simply not studied."""

NEW3 = """**3. ~~Conditioning is a blind spot in the whole Avery canon.~~
RETRACTED 2026-09-12 — see item 8.** Monkhorst & Jeziorski (JCP 71, 5268
(1979)) name "instabilities due to overcompleteness of basis sets" in their
abstract, which is explicit published engagement. **Do not write "blind spot",
"nobody studied", or "no prior art for the conditioning question".** What
survives is the narrower point 2: nobody *priced* it — no growth law, no
condition number as a function of basis size. The rest of the original
observation stands and is worth keeping: the isoenergetic framing is championed
because fixing the energy gives every basis function the correct asymptotic
decay at the classical turning points, and the *linear-algebraic* consequences
are not what that literature is about."""

OLD2_TAIL = """So
GeoVac's diagonal exponent is a **novel measurement, not citable prior art**, and
deriving it in closed form would be a genuinely new result."""

NEW2_TAIL = """So
GeoVac's diagonal exponent is a **novel measurement, not citable prior art**, and
deriving it in closed form would be a genuinely new result. *(Narrowed
2026-09-12, item 8: this SURVIVES the Monkhorst–Jeziorski read, but the
survival rests on that paper's abstract, reference list and page count — the
two-page body is unread.)*"""


def main() -> int:
    if not os.path.exists(M):
        print(f"not found: {M}")
        return 2
    with open(M, encoding="utf-8") as fh:
        t = fh.read()
    if MARKER in t:
        print("ALREADY APPLIED")
        return 1
    if t.count(OLD3) != 1 or t.count(OLD2_TAIL) != 1:
        print(f"anchors: p3={t.count(OLD3)} p2={t.count(OLD2_TAIL)}")
        return 3
    t = t.replace(OLD3, NEW3).replace(OLD2_TAIL, NEW2_TAIL)
    with open(M, "w", encoding="utf-8") as fh:
        fh.write(t)
    print("applied: points 2 and 3 marked at their own loci")
    return 0


if __name__ == "__main__":
    sys.exit(main())
