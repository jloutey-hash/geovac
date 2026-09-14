"""DELTA remediation -- the non-paper loci.

F2 (LARGE) -- the RECALL layer. `memory/avery_method_and_prior_art_gaps.md`
   auto-loads into every session and carried four pre-2026-09-12 readings, one
   of them as an OPERATIVE INSTRUCTION ("do claim the identification") folding
   in the translation reading that is prior art on three counts.  Also asserts
   the WITHDRAWN frames reading as live paper content, uses the demoted
   "Proposition D" label, and calls west_ruedenberg2013 "the one named source
   that uses an SVD/principal-angle construction" -- contradicted by the
   2026-09-12 abstract read ("no principal angles, no SVD").
   This layer sits outside the repo, outside C16/C21/C22, and the same file was
   self-caught for this exact class one day earlier.

F3 (SMALL) -- `docs/walls/register.md` states the breach scope as "M = 2, 3",
   which covers the collinear case the paper explicitly declines.  The
   register's stated purpose is DISPATCH, and BeH2 and CO2 are collinear and in
   this corpus's own library.

F4 (SMALL) -- the same register still calls the block-diagonal congruence
   result "Proposition D".

F6 (SMALL) -- `docs/code_architecture.md` describes the sigma-law as "derived",
   the word the owning module's docstring explicitly forbids.

F7 (SMALL) -- the synthesis and `docs/claims_register.md` quote the accuracy
   floor as a single value (6.4 / 6.44 mHa) where the paper says "it is the
   BRACKET we quote", [6.47, 6.62].  6.4 sits BELOW the bracket's own lower
   endpoint.  An honest-scope denial reversed downstream.

M7 (SMALL) -- `docs/qa/paper_60.done.md` C17 note freezes two RETIRED values as
   goalposts (K^0.84, retired 2026-09-07; He -2.897, retired by this run) on a
   premise the same file contradicts twice.  This is the precise class that
   stopped the 2026-09-11 run at protocol step 1.

Idempotent.
"""
from __future__ import annotations

import os
import sys

MEM = os.path.join(os.environ.get("USERPROFILE", os.path.expanduser("~")),
                   ".claude", "projects",
                   "C--Users-jlout-Desktop-Project-Geometric", "memory",
                   "avery_method_and_prior_art_gaps.md")

EDITS = [
 (MEM, "F2a-translation-instruction",
  "(Toeplitz minus Hankel), with symbol `j_0(kR cot(chi/2))`; equivalently the SW\n"
  "operator is multiplication by `e^{ip.R}` on the Fock sphere, going trivial at\n"
  "`p = 0`. Do not re-claim the asymptotic; do claim the identification.",
  "(Toeplitz minus Hankel), with symbol `j_0(kR cot(chi/2))`. Do not re-claim the\n"
  "asymptotic; do claim the SYMBOL.\n\n"
  "**CORRECTED 2026-09-12 — do NOT claim the translation reading.** This entry\n"
  "used to add \"equivalently the SW operator is multiplication by `e^{ip.R}` on\n"
  "the Fock sphere\" to the claimable identification. That reading is **prior art\n"
  "on three counts**: Shibuya-Wulfman's own 1965 abstract builds the molecular\n"
  "p0 operator from \"a sum of unitary transformations, one for each nucleus in\n"
  "the molecule\"; Wulfman & Takahata gave the explicit continuous-group\n"
  "formulation in 1967 (JCP 47, 488, Lie algebras of E4/R5/O(4,1)); Red &\n"
  "Weatherford derived the general formula for the matrix in a Coulomb-Sturmian\n"
  "basis in 2004 (IJQC 100, 208). What survives as ours is the SYMBOL alone."),

 (MEM, "F2b-frames-and-propd",
  "*Two further results from the same round, both in Paper 60 now:* the\n"
  "overcompleteness wall is an elementary theorem (completeness of the ONE-centre\n"
  "set forces `lam_min -> 0`; \"translate\" is incidental), and the `l`-selection\n"
  "loss is **independent of conditioning** (Proposition D) -- so it does not relax\n"
  "as `cond(S) -> 1+`, and the wall is stronger than the paper had stated.",
  "*One further result from the same round:* the `l`-selection loss is\n"
  "**independent of conditioning** -- it does not relax as `cond(S) -> 1+`, and\n"
  "the wall is stronger than the paper had stated. **Re-attributed 2026-09-12:**\n"
  "this is Loewdin symmetry preservation specialised to the `l` grading, known\n"
  "since Slater-Koster (1954); the label \"Proposition D\" is retired and only the\n"
  "`l`-vs-`m` application is claimed.\n\n"
  "**The overcompleteness-as-price-of-completeness reading is WITHDRAWN**\n"
  "[retracted 2026-09-12: p60-frames-completeness] and never reached the paper.\n"
  "It held that completeness of the one-centre set forces `lam_min -> 0`; measured,\n"
  "the Bessel deficit PLATEAUS at 0.380/0.696/0.907 for kR = 1/2/4, flat over\n"
  "N = 16..256, so the one-centre set is far from complete IN THE MOLECULAR\n"
  "METRIC and the bound holds only vacuously. Ron-Shen is the surviving\n"
  "mechanism: the near-dependence is ONE DIRECTION."),

 (MEM, "F2c-west-ruedenberg",
  "West-Ruedenberg 2013 is still unread (HTTP 403) and is the one named source\nthat uses an SVD/principal-angle construction.",
  "West-Ruedenberg 2013 was DROPPED from Paper 60 on 2026-09-12: its abstract was\n"
  "reached and describes localizing orbital transformations with **no principal\n"
  "angles, no SVD and no corresponding orbitals**, so it cannot support the\n"
  "attribution it carried. Amos-Hall (1961) and King (1967) are the verified\n"
  "lineage. (It is still cited in Paper 58's twin sentence -- owed there.)"),

 ("docs/walls/register.md", "F3-collinear-scope",
  "Remaining scope: `s`-sector shared-scale bases at `M = 2, 3`.",
  "Remaining scope: `s`-sector shared-scale bases at `M = 2` and **non-collinear** "
  "`M = 3`; the **collinear case is open and is NOT claimed** (corrected "
  "2026-09-12 -- the unqualified \"M = 2, 3\" authorised dispatch into exactly the "
  "regime Paper 60 declines, and BeH2 and CO2 are collinear and in this corpus's "
  "own library). For collinear centres `P D2 P` is rank ONE, so only one of the "
  "`M-1` directions opens at order 2 and the rest at 4, 6, ..., 2(M-1)."),

 ("docs/code_architecture.md", "F6-derived-word",
  "(derived N^2 law; collapse pi^2/24)",
  "(N^2 law -- **Kac-Murdock-Szego PRIOR ART, not derived here**; ours is the "
  "identification of the SW metric as such a finite section; collapse pi^2/24)"),
]


def main() -> int:
    applied = 0
    loaded: dict[str, str] = {}
    for path, name, old, new in EDITS:
        if not os.path.exists(path):
            print(f"  MISS {name}: {path} not found")
            return 2
        if path not in loaded:
            with open(path, encoding="utf-8") as fh:
                loaded[path] = fh.read()
        t = loaded[path]
        if new[:55] in t:
            print(f"  skip {name} (already applied)")
            continue
        if t.count(old) != 1:
            print(f"  MISS {name}: count={t.count(old)}")
            return 3
        loaded[path] = t.replace(old, new)
        applied += 1
        print(f"  ok   {name}")
    for path, t in loaded.items():
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(t)
    print(f"applied {applied}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
