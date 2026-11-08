r"""Bring fig:convergence's caption in line with the regenerated figure.

The old caption described a single dashed reference line and a bar series whose
numbers came from a superseded solver vintage (four of five appeared nowhere in
the paper).  The figure is now drawn from tab:extended_convergence directly, with
both azimuthal sectors plotted and both Paper 12 reference values labelled by
sector.

Also adds the \ref that was missing:\ fig:convergence carried a \label and was
never cited in the text, so nothing pointed a reader at it.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

P15 = "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex"
EDITS = []


def edit(old, new, label):
    EDITS.append((old, new, label))


edit(
    r"""\caption{Percentage of exact dissociation energy recovered as a function
  of $l_{\max}$.  The dashed line marks Paper~12's $\sigma$-only 92.4\% (its
  $|m|\le1$ value is 99.1\%, not a like-for-like reference for a
  $\sigma{+}\pi$ curve).  The dramatic
  jump at $l_{\max}=2$ reflects the onset of quadrupolar nuclear
  coupling (Sec.~\ref{sec:multichannel}).}""",
    r"""\caption{Percentage of exact dissociation energy recovered as a function
  of $l_{\max}$, for both azimuthal sectors.  Bars are the $\sigma$ and
  $\sigma{+}\pi$ columns of Table~\ref{tab:extended_convergence}
  (2D solver with Schwartz cusp correction);\ the figure is drawn from
  that table and adds no data of its own.  Both Paper~12 reference
  values are shown and each is labelled with the sector it belongs to:\
  its $\sigma$-only $92.4\%$ and its $|m|\le1$ $99.1\%$.  Neither is a
  like-for-like reference for these curves---the truncation and the
  solver class differ---and no ordering between the two coordinate
  systems is asserted here or anywhere in this paper.  The jump at
  $l_{\max}=2$ ($+42.3$ points, the largest step in the table) reflects
  the onset of quadrupolar nuclear coupling
  (Sec.~\ref{sec:multichannel}).}""",
    "caption rewritten for the two-sector figure, with the no-ordering scope")

# The figure had a \label and no \ref anywhere in the document.
edit(
    r"""\paragraph{Extended $l_{\max}$ convergence study.}
Table~\ref{tab:extended_convergence} extends the convergence analysis""",
    r"""\paragraph{Extended $l_{\max}$ convergence study.}
Table~\ref{tab:extended_convergence} (plotted in
Fig.~\ref{fig:convergence}) extends the convergence analysis""",
    "fig:convergence is now actually referenced from the text")

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

print("applied %d of %d" % (len(applied), len(EDITS)))
for a in applied:
    print("  +", a)
if failed:
    print("UNMATCHED:")
    for f in failed:
        print("  -", f)
    sys.exit(1)
