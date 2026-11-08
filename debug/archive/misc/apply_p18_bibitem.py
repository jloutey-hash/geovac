"""Paper 18: add the Paper 12 bibitem the re-pricing now cites.

Title taken verbatim from paper_12_algebraic_vee.tex's own \\title, not
reconstructed -- the C11 gate caught exactly that mistake earlier in this
sprint.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

P18 = "papers/group3_foundations/paper_18_exchange_constants.tex"

ANCHOR = r"""\bibitem{loutey_paper13}"""

NEW = r"""\bibitem{loutey_paper12}
J.~Loutey,
``Algebraic Two-Electron Integrals on the Prolate Spheroidal Lattice,''
GeoVac Paper~12 (2026).

\bibitem{loutey_paper13}"""

with io.open(P18, encoding="utf-8") as fh:
    t = fh.read()

if ANCHOR not in t:
    print("FAILED: anchor not found")
    sys.exit(1)
if r"\bibitem{loutey_paper12}" in t:
    print("already present")
    sys.exit(0)

t = t.replace(ANCHOR, NEW, 1)
with io.open(P18, "w", encoding="utf-8") as fh:
    fh.write(t)
print("Paper 18: loutey_paper12 bibitem added")
