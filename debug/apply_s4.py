"""S4: the 'every variational point exceeds 98.4%' universal, scoped to what
was actually measured.  On the full (alpha x threshold) grid the weakest
variational value is 95.5%.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

P12 = "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex"

OLD = r"""gap---is insensitive to the choice, since every variational
point exceeds $98.4\%$ and the independent Gaussian route, which is
well conditioned, gives $99.10\%$."""

NEW = r"""gap---is insensitive to the choice:\ every variational
point on the two one-dimensional slices quoted above exceeds $98.4\%$,
the weakest variational value anywhere on the full
$(\alpha, \mathrm{threshold})$ grid is $95.5\%$, and the independent
Gaussian route, which is well conditioned, gives $99.10\%$."""

with io.open(P12, encoding="utf-8") as fh:
    t = fh.read()

if OLD not in t:
    print("FAILED to match")
    sys.exit(1)

with io.open(P12, "w", encoding="utf-8") as fh:
    fh.write(t.replace(OLD, NEW, 1))

print("S4 applied")
