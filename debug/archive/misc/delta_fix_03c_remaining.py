"""DELTA remediation -- the three remaining non-paper loci.

delta_fix_03 aborted before writing (its F2c anchor was stale once
delta_fix_03b had applied it), so F2b, F3 and F6 never landed.  This file
carries them with markers unique to the NEW text.

  F2b -- memory: the WITHDRAWN frames reading asserted as live paper content,
         and the demoted "Proposition D" label.
  F3  -- walls register: breach scope stated as "M = 2, 3", covering the
         collinear case the paper declines.  The register's purpose is
         DISPATCH, and BeH2/CO2 are collinear and in this corpus's library.
  F6  -- code architecture: the sigma-law called "derived", the word the
         owning module's docstring explicitly forbids.

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
 (MEM, "F2b-frames-and-propd", "p60-frames-completeness",
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
  "**The overcompleteness-as-the-price-of-completeness reading is WITHDRAWN**\n"
  "[retracted 2026-09-12: p60-frames-completeness] and never reached the paper.\n"
  "It held that completeness of the one-centre set forces `lam_min -> 0`;\n"
  "measured, the Bessel deficit PLATEAUS at 0.380/0.696/0.907 for kR = 1/2/4,\n"
  "flat over N = 16..256, so the one-centre set is far from complete IN THE\n"
  "MOLECULAR METRIC and the bound holds only vacuously. Ron-Shen is the surviving\n"
  "mechanism: the near-dependence is ONE DIRECTION."),

 ("docs/walls/register.md", "F3-collinear-scope", "collinear case is open and is NOT claimed",
  "Remaining scope: `s`-sector shared-scale bases at `M = 2, 3`.",
  "Remaining scope: `s`-sector shared-scale bases at `M = 2` and **non-collinear** "
  "`M = 3`; the **collinear case is open and is NOT claimed** (corrected "
  "2026-09-12: the unqualified \"M = 2, 3\" authorised dispatch into exactly the "
  "regime Paper 60 declines, and BeH2 and CO2 are collinear and in this corpus's "
  "own library). Mechanism: for collinear centres `P D2 P` is rank ONE, so only "
  "one of the `M-1` null directions opens at order `p^2` and the rest at "
  "4, 6, ..., 2(M-1)."),

 ("docs/code_architecture.md", "F6-derived-word", "PRIOR ART, not derived here",
  "(derived N^2 law; collapse pi^2/24)",
  "(N^2 law -- **Kac-Murdock-Szego PRIOR ART, not derived here**; ours is the "
  "identification of the SW metric as such a finite section; collapse pi^2/24)"),
]


def main() -> int:
    applied = 0
    loaded: dict[str, str] = {}
    for path, name, marker, old, new in EDITS:
        if not os.path.exists(path):
            print(f"  MISS {name}: {path} not found")
            return 2
        if path not in loaded:
            with open(path, encoding="utf-8") as fh:
                loaded[path] = fh.read()
        t = loaded[path]
        if marker in t:
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
