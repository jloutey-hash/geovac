r"""Round-3 code-review LARGE-2: the stability-envelope paragraph, corrected
for the second time, plus its C16 entry.

PM VERIFICATION (I re-ran this myself rather than accepting the finding; the
scan is scratchpad/verify_large2.py, using the production module):

  (3,3), |m|<=1, N = 144        (3,3), sigma only, N = 72
  alpha  1e-11  1e-10  1e-9  1e-8      alpha  1e-8
  0.90   NONVAR 96.85  80.15 76.11     0.90   91.06
  1.05   NONVAR NONVAR NONVAR 92.09    1.05   92.10
  1.10   99.03  98.92  98.24 95.50     1.10   92.22

So 95.5% is the cell at (alpha = 1.10, 1e-8).  It is the grid minimum ONLY if
alpha is restricted to >= 1.10 -- which excludes exactly the low-alpha half
that the same paragraph has just called pathological.  Over the declared range
alpha in [0.90, 1.30] the weakest VARIATIONAL value is 76.11%, and at that cell
the sigma-only solve gives 91.06%:  opening the azimuthal channels there makes
the answer WORSE by fifteen points.  "Robust across the grid" is false at a
grid point.

This is the FOURTH consecutive round in which the largest defect was in the
previous round's remediation -- and this time the previous round was mine, from
today:  the sentence being corrected here is the one I wrote to fix round 2's
self-refuting "insensitive to the choice".  Recorded plainly in the note.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

P12 = "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex"
CHK = "debug/qa/check_retracted_terms.py"

OLD = r"""statement is therefore $99.0$--$99.1\%$, with $99.09\%$ the best
variational value at the stated parameters;\ the conclusion of this
section---that the azimuthal channels, and not growth along the
$\sigma$ axis, are what close the gap---is robust across the grid.  The
\emph{size} of the closure is not.

\textbf{[SCOPE 2026-09-14]}
[retracted 2026-09-14: p12-envelope-insensitive]  An earlier version of this
sentence called the conclusion ``insensitive to the choice'' and then quoted
the figure
that refutes it.  Every variational point on the two one-dimensional
slices above exceeds $98.4\%$, but the weakest variational value
anywhere on the full $(\alpha, \mathrm{threshold})$ grid is $95.5\%$,
which closes about $41\%$ of the $7.6$-point gap rather than
``essentially all'' of it.  The defensible statement is that the
qualitative effect is robust while the quantitative value ranges over
$95.5$--$99.1\%$ across the grid, the independent and well-conditioned
Gaussian route giving $99.10\%$."""

NEW = r"""statement is therefore $99.0$--$99.1\%$, with $99.09\%$ the best
variational value at the stated parameters.

\textbf{[SCOPE 2026-09-14, corrected twice]}
[retracted 2026-09-14: p12-envelope-insensitive]
[retracted 2026-09-14: p12-grid-floor-955]
Two earlier versions of this passage overstated how stable that number is.
The first called the conclusion ``insensitive to the choice'' and then
quoted the figure that refutes it.  The second replaced it with a range,
``$95.5$--$99.1\%$ across the grid'', and called the qualitative conclusion
``robust across the grid''.  Both are wrong, and the second is wrong in the
direction that flatters the result.

$95.5\%$ is the cell at $\alpha = 1.10$, threshold $10^{-8}$.  It is the
weakest variational value only if $\alpha$ is restricted to $\ge 1.10$---which
excludes precisely the low-$\alpha$ half this paragraph has just identified as
pathological.  Over the full declared range $\alpha \in [0.90, 1.30]$ the
weakest variational value is $76.11\%$, at $\alpha = 0.90$ and threshold
$10^{-8}$;\ and at that same cell the $\sigma$-only solve returns $91.06\%$,
so opening the azimuthal channels there makes the answer \emph{worse} by
fifteen points.  The qualitative conclusion is therefore \emph{not} robust
pointwise across the grid, and this paper does not claim that it is.

What the result rests on instead, stated exactly:\ at the tight discard
thresholds, where the generalized eigenproblem is actually being solved rather
than regularized, every surviving variational point lies between $99.03$ and
$99.14\%$;\ and the independent Gaussian-basis route, which is well
conditioned and needs no discard threshold at all, gives $99.10\%$.  The claim
is the agreement of those two.  It is not a property of the whole
$(\alpha, \mathrm{threshold})$ grid, most of which is conditioning noise."""

NEW_ENTRY = r'''    {
        "id": "p12-grid-floor-955",
        "note": "2026-09-14 (round-3 DELTA, code review).  Paper 12 said the "
                "weakest variational value 'anywhere on the full (alpha, "
                "threshold) grid is 95.5%' and that the conclusion is 'robust "
                "across the grid'.  REFUTED by re-measurement: 95.5% is the "
                "cell at alpha = 1.10, threshold 1e-8, and is the minimum only "
                "if alpha is restricted to >= 1.10 -- excluding the low-alpha "
                "half the same paragraph calls pathological.  Over the declared "
                "range alpha in [0.90, 1.30] the weakest variational value is "
                "76.11% (alpha = 0.90, 1e-8), where the sigma-only solve gives "
                "91.06% -- so the azimuthal channels LOSE fifteen points at "
                "that cell.  A restricted-evaluation artifact: a clean floor "
                "produced by deleting the part of the object that breaks it.  "
                "Recorded also as process: this was the fourth consecutive "
                "round whose largest defect sat in the previous round's "
                "remediation, and the sentence was written the same day by the "
                "PM to fix the round-2 defect next to it.",
        "pattern": r"weakest variational value\s+anywhere on the full"
                   r"|ranges over\s+\$?95\.5\$?[^0-9]{0,12}99\.1"
                   r"|robust across the grid",
        "exempt_if_nearby": r"(?!)",
        "severity": "fail",
        "scope": "group2 synthesis",
        "cited_by": {
            "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex": "reviewed 2026-09-14 -- owner, corrected in place",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex": "reviewed 2026-09-14 -- quotes the 99.0-99.1 envelope, never the grid floor",
            "docs/claim_test_matrix.md": "reviewed 2026-09-14 -- no row asserts a grid floor",
        },
        "files": [
            "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
            "docs/claim_test_matrix.md",
        ],
    },
'''

applied, failed = [], []

with io.open(P12, encoding="utf-8") as fh:
    t = fh.read()
if OLD in t:
    with io.open(P12, "w", encoding="utf-8") as fh:
        fh.write(t.replace(OLD, NEW, 1))
    applied.append("LARGE-2: envelope paragraph corrected to the measured grid")
else:
    failed.append("LARGE-2 paper block")

with io.open(CHK, encoding="utf-8") as fh:
    src = fh.read()
ANCHOR = '    {\n        "id": "p15-sigma-pi-decoupled",'
if "p12-grid-floor-955" in src:
    applied.append("registry entry already present")
elif ANCHOR in src:
    with io.open(CHK, "w", encoding="utf-8") as fh:
        fh.write(src.replace(ANCHOR, NEW_ENTRY + ANCHOR, 1))
    applied.append("C16 entry p12-grid-floor-955 added")
else:
    failed.append("registry anchor")

print("applied %d" % len(applied))
for a in applied:
    print("  +", a)
if failed:
    print("UNMATCHED:")
    for f in failed:
        print("  -", f)
    sys.exit(1)
