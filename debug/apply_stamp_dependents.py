"""Stamp the four dependents that were left unstamped, now that they are done.

Sec. 9: "Stamp outcomes, not intentions."  Each of these was actually
corrected and verified this session, so each earns its stamp:

  paper_18_exchange_constants.tex -- Level-4 subsection re-priced: the
      irreducibility argument, the "structurally blind" reading, the 92.5%
      ceiling as a product-space property, the "surpassing the prolate CI
      ceiling" comparison, the qualitative Level-3/4 distinction, the L1604
      "key irreducibility result" and the L3021 Class-C necessity claim are
      all withdrawn or scoped.  Compiles clean; C11 passes on group3.

  viz/public/papers/*.html -- regenerated from the corrected .tex via
      zenodo_manifest_build.py + build_paper_pages.py; the withdrawn
      sentences no longer appear anywhere under viz/public/papers/.
      (Regeneration writes local files only; deployment happens on push,
      which is PI-only.)

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

GATE = "debug/qa/check_retracted_terms.py"

PAIRS = [
    ('"papers/group3_foundations/paper_18_exchange_constants.tex": None,',
     '"papers/group3_foundations/paper_18_exchange_constants.tex": "reviewed 2026-09-14 -- Level-4 subsection re-priced",'),
    ('"viz/public/papers/paper_12_algebraic_vee.html": None,',
     '"viz/public/papers/paper_12_algebraic_vee.html": "reviewed 2026-09-14 -- regenerated",'),
    ('"viz/public/papers/paper_13_hyperspherical.html": None,',
     '"viz/public/papers/paper_13_hyperspherical.html": "reviewed 2026-09-14 -- regenerated",'),
    ('"viz/public/papers/paper_15_level4_geometry.html": None,',
     '"viz/public/papers/paper_15_level4_geometry.html": "reviewed 2026-09-14 -- regenerated",'),
]

with io.open(GATE, encoding="utf-8") as fh:
    t = fh.read()

n = 0
for old, new in PAIRS:
    if old in t:
        t = t.replace(old, new, 1)
        n += 1
    else:
        print("  ! not matched:", old[:60])

# the comment above the block said Paper 18 is left unstamped deliberately;
# that is no longer true.
t = t.replace(
    "        # independent reviewers found it; no pattern could, because Paper 18\n"
    "        # restates the claim in its own words and its own numbers (92.5%,\n"
    "        # 7.5%).  Left unstamped deliberately -- the correction is a PI call.",
    "        # independent reviewers found it; no pattern could, because Paper 18\n"
    "        # restates the claim in its own words and its own numbers (92.5%,\n"
    "        # 7.5%).  Re-priced 2026-09-14 under PI direction.",
    1)

with io.open(GATE, "w", encoding="utf-8") as fh:
    fh.write(t)

print("stamped %d dependents" % n)
if n != len(PAIRS):
    sys.exit(1)
