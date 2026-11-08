"""B2 and B3: two doc rows the round-1 fix missed because it was locus-by-locus.

B2 -- docs/topic_to_paper_lookup.md has the cusp-diagnosis row TWICE; line 21
      was corrected, line 72 was not.
B3 -- docs/validation_benchmarks.md line 14 names the solver class ("2D
      solver"), which positively implies a variational bound, for a figure
      that is cusp-corrected.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

EDITS = [
    ("docs/topic_to_paper_lookup.md",
     "| Cusp diagnosis (7.6% gap) | 12 | Sec VII | Core |",
     "| Azimuthal-channel diagnosis of the 7.6% gap (the cusp reading is withdrawn) | 12 | Sec VII | Core |",
     "B2: duplicate cusp-diagnosis index row corrected"),
    ("docs/validation_benchmarks.md",
     "| H2 Level 4 (2D solver) | 96.0% D_e | Molecule-frame hyperspherical |",
     "| H2 Level 4 (2D solver + Schwartz cusp correction) | 96.0% D_e | Molecule-frame hyperspherical; ~95% pure-variational |",
     "B3: benchmark row states the cusp correction"),
]

applied, failed = [], []
for path, old, new, label in EDITS:
    with io.open(path, encoding="utf-8") as fh:
        t = fh.read()
    if old in t:
        with io.open(path, "w", encoding="utf-8") as fh:
            fh.write(t.replace(old, new, 1))
        applied.append(label)
    else:
        failed.append(label)

print("applied %d edits" % len(applied))
for a in applied:
    print("  +", a)
if failed:
    print("UNMATCHED:")
    for f in failed:
        print("  -", f)
    sys.exit(1)
