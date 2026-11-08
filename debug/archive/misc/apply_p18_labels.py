"""Fix the two labels I invented instead of reading.

`sec:irreducibility` and `sec:mu_level3` do not exist.  The real labels are
`sec:mu` (the Level-3 mu subsection, L374) and, for the geometric-elevation
obstruction, the paragraph "Irreducibility of the Level 4 piecewise structure"
at L2233 which carries NO label -- so it is referred to by name, not by ref.
`sec:algebraic_curve` (L3198) is the section stating that the pencil structure
is specific to Level 3, which is the better target for the structural point.

Same class of error the C11 gate caught earlier in this sprint: writing a
reference from memory instead of reading it.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

P18 = "papers/group3_foundations/paper_18_exchange_constants.tex"

PAIRS = [
    (r"""is an \emph{algebraic} function over $\mathbb{Q}(\pi,\sqrt2)$---pointwise
diagonalisation there is computational convenience, not necessity
(Sec.~\ref{sec:mu_level3};\ Track~P1).""",
     r"""is an \emph{algebraic} function over $\mathbb{Q}(\pi,\sqrt2)$---pointwise
diagonalisation there is computational convenience, not necessity
(Sec.~\ref{sec:mu} and Sec.~\ref{sec:algebraic_curve};\ Track~P1)."""),
    (r"""both state for $\mu(R)$, and it is withdrawn.""",
     r"""both state for $\mu(R)$, and it is withdrawn."""),
    (r"""other irreducibility results for Level~4 are untouched and should not be
read as withdrawn with it:\ the geometric-elevation obstruction of
Sec.~\ref{sec:irreducibility} (three routes ruled out), and the structural
one just stated (no global $P(\rho,\mu)$, Track~S).""",
     r"""other irreducibility results for Level~4 are untouched and should not be
read as withdrawn with it:\ the geometric-elevation obstruction recorded
under ``Irreducibility of the Level~4 piecewise structure'' below (three
routes ruled out), and the structural one just stated (no global
$P(\rho,\mu)$, Track~S)."""),
]

with io.open(P18, encoding="utf-8") as fh:
    t = fh.read()

n = 0
for old, new in PAIRS:
    if old in t and old != new:
        t = t.replace(old, new, 1)
        n += 1

with io.open(P18, "w", encoding="utf-8") as fh:
    fh.write(t)

print("repaired %d invented reference(s)" % n)
if n < 2:
    print("WARNING: expected 2")
    sys.exit(1)
