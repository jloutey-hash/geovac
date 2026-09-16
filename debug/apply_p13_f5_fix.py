"""Two remaining owed items.

(1) Paper 13 L1896-1898 attributes the graph-native FCI convergence floor to
    the cusp.  The failed-approaches ledger records the OPPOSITE, as a
    documented negative: "0.20% floor is small-Z graph-validity-boundary
    artifact (Z_c ~ 1.84), not cusp. Z=10 sign flip confirms."
    (docs/failed_approaches_ledger.md, CUSP-2 row.)  Paper 13 itself cites
    Z_c ~ 1.84 elsewhere.  Since this whole sprint is about a floor wrongly
    attributed to the cusp, the contradiction is not incidental.

(2) Synthesis Table II row 896 certifies prolate spheroidal V_ee as
    "algebraic ... with no quadrature" under a caption that now reads as
    covering the |m| <= 1 headline, which Paper 12 explicitly denies.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

P13 = "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex"
SYN = "papers/synthesis/group2_quantum_chemistry_synthesis.tex"
EDITS = []


def edit(path, old, new, label):
    EDITS.append((path, old, new, label))


edit(
    P13,
    r"""FCI basis invariance was verified (Sprint~3D):
transforming $V_{ee}$ to the graph eigenbasis gives identical energies
(to $< 3 \times 10^{-15}$~Ha), confirming the convergence floor
is from the cusp (embedding exchange constant, Paper~18), not basis
mismatch.""",
    r"""FCI basis invariance was verified (Sprint~3D):
transforming $V_{ee}$ to the graph eigenbasis gives identical energies
(to $< 3 \times 10^{-15}$~Ha), confirming the convergence floor is not
basis mismatch.

\textbf{[WITHDRAWN]} Earlier versions went on to attribute that floor to
the cusp.  The project's own record contradicts it:\ the Schwartz-tail
diagnostic (CUSP-2) found the tail correction \emph{worsens} accuracy and
re-diagnosed the floor as a small-$Z$ graph-validity-boundary artifact
near $Z_c \approx 1.84$---confirmed by a sign flip at $Z = 10$---rather
than as cusp content.  Basis invariance shows what the floor is
\emph{not};\ it does not identify what it is.""",
    "P13: cusp-floor attribution withdrawn (contradicted the ledger)")

edit(
    SYN,
    r"""Prolate spheroidal $V_{ee}$ & algebraic & Neumann expansion, recurrence moments \cite{loutey_paper12} \\""",
    r"""Prolate spheroidal $V_{ee}$ ($\sigma$) & algebraic & Neumann expansion, recurrence moments \cite{loutey_paper12} \\
Prolate spheroidal $V_{ee}$ ($|m|\ge1$) & quadrature & same kernel, ordered-$\xi$ integral by spectral panels \cite{loutey_paper12} \\""",
    "synthesis F5: the algebraic row scoped to sigma; the |m|>=1 cost shown")

by_path = {}
for path, old, new, label in EDITS:
    by_path.setdefault(path, []).append((old, new, label))

applied, failed = [], []
for path, items in by_path.items():
    with io.open(path, encoding="utf-8") as fh:
        t = fh.read()
    for old, new, label in items:
        if old in t:
            t = t.replace(old, new, 1)
            applied.append(label)
        else:
            failed.append(label)
    with io.open(path, "w", encoding="utf-8") as fh:
        fh.write(t)

print("applied %d edits" % len(applied))
for a in applied:
    print("  +", a)
if failed:
    print("UNMATCHED:")
    for f in failed:
        print("  -", f)
    sys.exit(1)
