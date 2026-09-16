"""Paper 15 F7: the delta-channel sentence is internally garbled, and the two
sigma+pi figures (93.6 and 94.1) look contradictory but are not.

Table tab:extended_convergence at l_max = 4 reads sigma-only 87.0, sigma+pi
93.6, and the delta row 87.6 with footnote "gain relative to sigma-only at the
same l_max".  So the prose "from 93.0% to 93.6% for sigma-only; from 87.0% to
87.6% without pi" is wrong twice: 93.6 is the sigma+pi column, not a delta
result, and 93.0 appears nowhere.  Only the 87.0 -> 87.6 pair is a delta gain.

The 93.6 vs 94.1 difference is NOT an inconsistency: that table's own caption
says its sigma+pi column freezes the pi channels at their l_max = 2 values,
while tab:comparison's 94.1 is the full sigma+pi solve.  The paper never says
so where a reader meets both numbers.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

P15 = "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex"
EDITS = []


def edit(old, new, label):
    EDITS.append((old, new, label))


edit(
    r"""yields a modest $+0.65$~percentage-point gain (from 93.0\% to
93.6\% for $\sigma$-only; from 87.0\% to 87.6\% without $\pi$).""",
    r"""yields a modest $+0.65$~percentage-point gain, from $87.0\%$ to
$87.6\%$ relative to $\sigma$-only at the same $l_{\max}$ (the figure
quoted in Table~\ref{tab:extended_convergence}'s footnote).""",
    "F7a: delta-gain sentence corrected to match the table")

edit(
    r"""  $\sigma{+}\pi$ uses $m_{\max} = 1$ with $\pi$ channels frozen at
  their $l_{\max} = 2$ values.""",
    r"""  $\sigma{+}\pi$ uses $m_{\max} = 1$ with $\pi$ channels frozen at
  their $l_{\max} = 2$ values --- which is why this column's
  $l_{\max} = 4$ entry ($93.6\%$) sits below the $94.1\%$ of
  Table~\ref{tab:comparison}, where the $\pi$ channels are solved
  at full $l_{\max}$;\ the two are different calculations, not a
  discrepancy.""",
    "F7b: the 93.6-vs-94.1 difference explained where both are met")

with io.open(P15, encoding="utf-8") as fh:
    t = fh.read()

applied, failed = [], []
for old, new, label in EDITS:
    if old in t:
        t = t.replace(old, new, 1)
        applied.append(label)
    else:
        failed.append(label)

with io.open(P15, "w", encoding="utf-8") as fh:
    fh.write(t)

print("applied %d edits" % len(applied))
for a in applied:
    print("  +", a)
if failed:
    print("UNMATCHED:")
    for f in failed:
        print("  -", f)
    sys.exit(1)
