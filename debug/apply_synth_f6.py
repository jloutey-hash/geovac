"""Synthesis F6: both summary surfaces give 96.0% as the arc's H2 ceiling and
omit the better prolate-spheroidal number -- the Sec. 9 summary-surface rule
applied in reverse (under-claiming).

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

S = "papers/synthesis/group2_quantum_chemistry_synthesis.tex"
EDITS = []


def edit(old, new, label):
    EDITS.append((old, new, label))


edit(
    r"""for He, $96.0\%$ of $D_e$ for $\mathrm{H}_2$, and---in the""",
    r"""for He, ${\sim}99.1\%$ of $D_e$ for $\mathrm{H}_2$ in prolate
spheroidal coordinates with the azimuthal channels open (against
$96.0\%$ in the molecule-frame hyperspherical treatment, cusp-corrected),
and---in the""",
    "F6a: abstract ceilings list carries the better H2 number")

edit(
    r"""4 & $\mathrm{H}_2$ (2, 2) & molecule-frame hyperspherical & $96.0\%$ of $D_e$ ($l_{\max}=6$, 61 channels) & \cite{loutey_paper15} \\""",
    r"""4 & $\mathrm{H}_2$ (2, 2) & molecule-frame hyperspherical & $96.0\%$ of $D_e$ ($l_{\max}=6$, 61 channels, cusp-corrected) & \cite{loutey_paper15} \\
2$'$ & $\mathrm{H}_2$ (2, 2) & prolate spheroidal, $|m|\le1$ & ${\sim}99.1\%$ of $D_e$ & \cite{loutey_paper12} \\""",
    "F6b: hierarchy table gains the prolate two-electron row")

edit(
    r"""coordinates for $\mathrm{H}_2$ ($96.0\%$ of""",
    r"""coordinates for $\mathrm{H}_2$ (${\sim}99.1\%$ in prolate spheroidal
coordinates with the azimuthal channels open~\cite{loutey_paper12};\
$96.0\%$ of""",
    "F6c: conclusion carries the better H2 number")

with io.open(S, encoding="utf-8") as fh:
    t = fh.read()

applied, failed = [], []
for old, new, label in EDITS:
    if old in t:
        t = t.replace(old, new, 1)
        applied.append(label)
    else:
        failed.append(label)

with io.open(S, "w", encoding="utf-8") as fh:
    fh.write(t)

print("applied %d edits" % len(applied))
for a in applied:
    print("  +", a)
if failed:
    print("UNMATCHED:")
    for f in failed:
        print("  -", f)
    sys.exit(1)
