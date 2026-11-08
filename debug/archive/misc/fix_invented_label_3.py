"""Third invented reference this session: sec:comparison does not exist in
Paper 15.  Replace the \\ref with a plain pointer rather than guess a label.

Recording the pattern, because three is not an accident: when I write a
cross-reference from memory instead of grepping for the label first, it is
wrong.  The C10 gate catches it, which is why none reached a reader -- but the
cheaper fix is to look the label up before writing it.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

P15 = "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex"

OLD = r"""cusp sits;\ as the scope note in Sec.~\ref{sec:comparison} records, it
does not translate into an accuracy advantage.)"""
NEW = r"""cusp sits;\ as the scope note accompanying the comparison with Paper~12
records, it does not translate into an accuracy advantage.)"""

with io.open(P15, encoding="utf-8") as fh:
    t = fh.read()

if OLD not in t:
    print("FAILED to match")
    sys.exit(1)

with io.open(P15, "w", encoding="utf-8") as fh:
    fh.write(t.replace(OLD, NEW, 1))

print("invented label removed")
