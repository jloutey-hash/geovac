"""DELTA remediation -- the two remaining memory-file loci.

Two failures in delta_fix_03, both mine and both avoidable:

  * F2a was SKIPPED by a false-positive idempotency check.  I compared the
    first 55 characters of the NEW text, which are identical to the OLD text's
    opening -- the third time today I have made exactly this mistake.  Fixed
    here by keying on a marker that exists ONLY in the new text.
  * F2c missed because the anchor spanned a line wrap I had guessed rather
    than read.

Idempotent, keyed on distinctive markers.
"""
from __future__ import annotations

import os
import sys

M = os.path.join(os.environ.get("USERPROFILE", os.path.expanduser("~")),
                 ".claude", "projects",
                 "C--Users-jlout-Desktop-Project-Geometric", "memory",
                 "avery_method_and_prior_art_gaps.md")

EDITS = [
 ("F2a-translation-instruction",
  "do NOT claim the translation reading",
  "(Toeplitz minus Hankel), with symbol `j_0(kR cot(chi/2))`; equivalently the SW\n"
  "operator is multiplication by `e^{ip.R}` on the Fock sphere, going trivial at\n"
  "`p = 0`. Do not re-claim the asymptotic; do claim the identification.",
  "(Toeplitz minus Hankel), with symbol `j_0(kR cot(chi/2))`. Do not re-claim the\n"
  "asymptotic; do claim the SYMBOL.\n\n"
  "**CORRECTED 2026-09-12 — do NOT claim the translation reading.** This entry\n"
  "used to fold \"equivalently the SW operator is multiplication by `e^{ip.R}` on\n"
  "the Fock sphere\" into the claimable identification, and then instructed the\n"
  "reader to claim it. That reading is **prior art on three counts**:\n"
  "Shibuya-Wulfman's own 1965 abstract builds the molecular p0 operator from \"a\n"
  "sum of unitary transformations, one for each nucleus in the molecule\";\n"
  "Wulfman & Takahata gave the explicit continuous-group formulation in 1967\n"
  "(JCP 47, 488 -- Lie algebras of E4, R5, O(4,1)); Red & Weatherford derived the\n"
  "general formula for that matrix in a Coulomb-Sturmian basis in 2004 (IJQC 100,\n"
  "208). **What survives as ours is the SYMBOL alone.**"),

 ("F2c-west-ruedenberg",
  "cannot support the",
  "West-Ruedenberg 2013 is still\nunread (HTTP 403) and is the one named source that uses an SVD/principal-angle\nconstruction.",
  "West-Ruedenberg 2013 was **DROPPED**\nfrom Paper 60 on 2026-09-12: its abstract was reached and describes localizing\n"
  "orbital transformations with **no principal angles, no SVD and no corresponding\n"
  "orbitals**, so it cannot support the attribution it carried. Amos-Hall (1961)\n"
  "and King (1967) are the verified lineage. It is STILL cited in Paper 58's twin\n"
  "sentence (L860) -- owed there, not a Paper 60 defect."),
]


def main() -> int:
    if not os.path.exists(M):
        print(f"not found: {M}")
        return 2
    with open(M, encoding="utf-8") as fh:
        t = fh.read()
    applied = 0
    for name, marker, old, new in EDITS:
        if marker in t:
            print(f"  skip {name} (already applied)")
            continue
        if t.count(old) != 1:
            print(f"  MISS {name}: count={t.count(old)}")
            return 3
        t = t.replace(old, new)
        applied += 1
        print(f"  ok   {name}")
    with open(M, "w", encoding="utf-8") as fh:
        fh.write(t)
    print(f"applied {applied}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
