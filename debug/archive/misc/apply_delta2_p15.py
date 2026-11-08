"""The three Paper 15 edits whose anchors differed from my reconstruction.

Read from the file, not rebuilt from the reviewer's quotes.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

P15 = "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex"
EDITS = []


def edit(old, new, label):
    EDITS.append((old, new, label))


edit(
    r"""dependence, whereas the Level~4 figure includes $\pi$ channels.  At
matched angular content the ordering reverses---Paper~12's own basis
with $|m| \le 1$ reaches 99.1\%---so no coordinate-system advantage
is claimed here.""",
    r"""dependence, whereas the Level~4 figure includes $\pi$ channels.  Nor
is a corrected comparison available in the other direction:\ Paper~12's
own basis with $|m| \le 1$ reaches $99.1\%$, but that figure and the
ones here differ in truncation, cusp treatment and solver class alike,
so no coordinate-system advantage is claimed in either direction.""",
    "D1a: P15 abstract -- reverse-ordering claim removed")

edit(
    r"""basis is $\varphi$-independent and spans only $m_1 = m_2 = 0$.  Run
at matched angular content the ordering reverses.  Paper~12's own
basis with $|m| \le 1$ reaches 99.1\%, against 94.1\% here at
$l_{\max}=4$ and 96.0\% at $l_{\max}=6$ with a cusp correction;\ and a
grid-based prolate spheroidal calculation reaches
99.97\%~\cite{tao_mccurdy_rescigno2010}.""",
    r"""basis is $\varphi$-independent and spans only $m_1 = m_2 = 0$.  Nor
does a corrected comparison run the other way, because no matched pair
exists:\ Paper~12's own basis with $|m| \le 1$ reaches $99.1\%$ at
$(j_{\max},l_{\max}) = (3,3)$, against $94.1\%$ here at $l_{\max}=4$ and
$96.0\%$ at $l_{\max}=6$ with a Schwartz cusp correction and a different
solver class.  A grid-based prolate spheroidal
calculation~\cite{tao_mccurdy_rescigno2010} reaches $99.97\%$, which
settles the narrower point that those coordinates are not cusp-limited
without ordering the two solvers here.""",
    "D1b: P15 SCOPE -- reverse-ordering claim removed")

edit(
    r"""This is the structural advantage of Level~4 coordinates: the cusp is
always a boundary condition, never a coordinate singularity, for
\emph{any} internuclear separation.""",
    r"""This is the structural difference of Level~4 coordinates: the cusp is
always a boundary condition rather than a coordinate surface, for
\emph{any} internuclear separation.  (It is a difference in where the
cusp sits;\ as the scope note in Sec.~\ref{sec:comparison} records, it
does not translate into an accuracy advantage.)""",
    "B7: P15 body -- the source sentence the synthesis copied")

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
