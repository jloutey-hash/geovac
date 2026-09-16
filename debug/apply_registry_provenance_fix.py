r"""Two provenance defects inside one registry entry (round-3 code review NIT).

1. "NON-VARIATIONAL values at some alpha" understates what was measured: six of
   nine grid points over alpha in [0.90, 1.30] fail at the declared threshold,
   including alpha = 1.00.  The paper says "a majority"; the registry said
   "some".
2. The provenance string still carried 99.15 in the threshold ladder while this
   entry's OWN alias had already been corrected to 99.14 (measured 99.1446) --
   the two disagreed inside one entry.
3. The grid floor is added, because the paper now states it and a reader of the
   registry should not have to reconstruct it.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

REG = "debug/qa/numeric_registry.py"

OLD = ('                   "NOT a stable fourth digit: at this basis cond(S)=2e16, the "\n'
       '                   "solver returns NON-VARIATIONAL values at some alpha "\n'
       '                   "(1.15 -> -4.1 Ha, 1.20 -> -8.1 Ha) and the value moves "\n'
       '                   "99.15/99.09/98.99/98.41 across thresholds "\n'
       '                   "1e-12/1e-11/1e-10/1e-8. Quote 99.1% at summary surfaces; "\n'
       '                   "the precise value only where alpha and threshold are "\n'
       '                   "stated. Envelope: 99.0-99.1%.",')

NEW = ('                   "NOT a stable fourth digit: at this basis cond(S)=2e16, the "\n'
       '                   "solver returns NON-VARIATIONAL values at a MAJORITY of "\n'
       '                   "alpha -- six of nine grid points over [0.90, 1.30] at the "\n'
       '                   "declared threshold, including alpha=1.00, the natural "\n'
       '                   "default (0.95 -> -277 Ha, 1.15 -> -4.1, 1.20 -> -8.1) -- "\n'
       '                   "and the value moves 99.14/99.09/98.99/98.41 across "\n'
       '                   "thresholds 1e-12/1e-11/1e-10/1e-8. Over the full "\n'
       '                   "(alpha, threshold) grid the weakest VARIATIONAL value is "\n'
       '                   "76.11% at (0.90, 1e-8), where the sigma-only solve gives "\n'
       '                   "91.06% -- i.e. the channels lose ground there; 95.5% is "\n'
       '                   "the floor only for alpha >= 1.10 and must not be quoted "\n'
       '                   "as a grid floor (C16: p12-grid-floor-955). Quote 99.1% at "\n'
       '                   "summary surfaces; the precise value only where alpha and "\n'
       '                   "threshold are stated. Envelope: 99.0-99.1%.",')

with io.open(REG, encoding="utf-8") as fh:
    t = fh.read()

if OLD not in t:
    print("FAILED to match")
    sys.exit(1)

with io.open(REG, "w", encoding="utf-8") as fh:
    fh.write(t.replace(OLD, NEW, 1))

print("  + p12_azimuthal_de_pct provenance corrected")
