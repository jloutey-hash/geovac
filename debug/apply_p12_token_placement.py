r"""Move the p12-grid-floor-955 token next to the wording it shelters.

C16's MARKER_WINDOW is +-2 lines, deliberately: the +-5 window used to let a
withdrawal note on one claim shelter a LIVE zombie of a different claim a few
lines away (the P7 L127 cross-shelter bug, FULL #6).  My first placement put
both tokens in a header three to five lines above the quoted wording, so the
gate correctly reported the quote as live.

The gate was right and I was wrong;  the fix is to put the token where the
retired wording actually is, not to widen the window.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

P12 = "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex"

OLD = r"""\textbf{[SCOPE 2026-09-14, corrected twice]}
[retracted 2026-09-14: p12-envelope-insensitive]
[retracted 2026-09-14: p12-grid-floor-955]
Two earlier versions of this passage overstated how stable that number is.
The first called the conclusion ``insensitive to the choice'' and then
quoted the figure that refutes it.  The second replaced it with a range,
``$95.5$--$99.1\%$ across the grid'', and called the qualitative conclusion
``robust across the grid''.  Both are wrong, and the second is wrong in the
direction that flatters the result."""

NEW = r"""\textbf{[SCOPE 2026-09-14, corrected twice]}
Two earlier versions of this passage overstated how stable that number is.
The first called the conclusion ``insensitive to the choice''
[retracted 2026-09-14: p12-envelope-insensitive] and then quoted the figure
that refutes it.  The second replaced it with a range, ``$95.5$--$99.1\%$
across the grid'', and called the qualitative conclusion ``robust across the
grid'' [retracted 2026-09-14: p12-grid-floor-955].  Both are wrong, and the
second is wrong in the direction that flatters the result."""

with io.open(P12, encoding="utf-8") as fh:
    t = fh.read()

if OLD not in t:
    print("FAILED to match")
    sys.exit(1)

with io.open(P12, "w", encoding="utf-8") as fh:
    fh.write(t.replace(OLD, NEW, 1))

print("  + tokens moved adjacent to the wording they shelter")
