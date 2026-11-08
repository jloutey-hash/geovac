"""C11 defect, self-introduced: the Paper-15 bibitem added to Paper 12 in this
sprint carried a FABRICATED title.

I wrote a plausible-sounding title instead of reading Paper 15's own \\title.
Paper 13 already cites the same paper correctly; that wording is copied here.
Exactly the class C11 exists to catch, caught on my own remediation -- the
"remediated text is not clean text" rule.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

PATH = "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex"

OLD = r"""\bibitem{loutey_paper15}
J.~Loutey,
``Molecule-Frame Hyperspherical Coordinates: the Level-4 Geometry
for Two-Center, Two-Electron Systems,''
GeoVac Paper~15 (2026)."""

NEW = r"""\bibitem{loutey_paper15}
J.~Loutey,
``The Level~4 Natural Geometry: Two-Center Two-Electron Molecules in
Molecule-Frame Hyperspherical Coordinates,''
GeoVac Paper~15 (2026)."""

with io.open(PATH, encoding="utf-8") as fh:
    text = fh.read()

if OLD not in text:
    print("FAILED TO MATCH the fabricated bibitem")
    sys.exit(1)

with io.open(PATH, "w", encoding="utf-8") as fh:
    fh.write(text.replace(OLD, NEW, 1))

print("Paper 12: loutey_paper15 bibitem title corrected to Paper 15's own title")
