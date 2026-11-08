"""Dependents of Paper 12's withdrawn cusp diagnosis: Paper 13's motivation
and the group2 synthesis.

The claim was owned by Paper 12, corrected there, and restated in these two
documents in their own words -- the exact shape the Sec. 9 retraction->dependents
rule exists to catch.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

EDITS = []


def edit(path, old, new, label):
    EDITS.append((path, old, new, label))


P13 = "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex"
SYN = "papers/synthesis/group2_quantum_chemistry_synthesis.tex"

# --------------------------------------------------------------- Paper 13
edit(
    P13,
    r"""The motivation for this step comes from Paper~12's
diagnosis~\cite{loutey_paper12}: the prolate spheroidal CI for H$_2$
saturates at 92.4\% of the exact dissociation energy $D_e$, with the
remaining 7.6\% gap attributed to the electron-electron cusp.  The
exact two-electron wavefunction contains non-analytic terms
$r_{12}^{1/2}$ and $r_{12}\ln r_{12}$ near the three-body
coalescence~\cite{Fock1954}, which no polynomial basis in
single-electron coordinates can represent.""",
    r"""The motivation for this step is the structure of the two-electron
wavefunction itself.  Near the three-body coalescence the exact
wavefunction contains non-analytic terms $r_{12}^{1/2}$ and
$r_{12}\ln r_{12}$~\cite{Fock1954}, which no polynomial basis in
single-electron coordinates represents exactly, and for an atom the
coalescence is genuinely three-body---nuclear and interelectronic
singularities meet at a single point.

\emph{Scope note.}  Earlier versions of this paragraph took the
motivation from Paper~12~\cite{loutey_paper12}, whose prolate
spheroidal CI for H$_2$ saturates at 92.4\% of $D_e$, and read that
residual as the cusp demanding a coordinate change.  Paper~12 has
since withdrawn that diagnosis:\ its gap was the absent $m \ne 0$
configurations, and restoring them in the same prolate spheroidal
basis reaches 99.09\%.  Nothing in the present paper depends on the
withdrawn reading---helium's coalescence is three-body whether or
not H$_2$'s residual was misattributed---but the motivation is
restated here so it no longer rests on it.""",
    "P13: motivation no longer rests on Paper 12's withdrawn diagnosis")

# ------------------------------------------------------ synthesis, locus 1
edit(
    SYN,
    r"""eliminated~\cite{loutey_paper12}.  The paper then diagnoses the
remaining $7.6\%$ gap honestly:\ it is a one-electron basis
completeness limit, not an integration error.  The electron--electron
cusp requires the non-analytic terms $r_{12}^{1/2}$ and
$r_{12}\ln r_{12}$ (the Kato/Fock structure~\cite{kato1957}) that no polynomial prolate
spheroidal basis can represent.  This is the framework's embedding-tier
transcendental:\ the cusp is two-dimensional in the angular pair
coordinate and is the canonical example of an exchange constant that no
single algebraic insertion reproduces~\cite{loutey_paper18}.  The
diagnosis identifies the next natural geometry---hyperspherical
coordinates, where the coalescence is a boundary condition rather than
a coordinate singularity.""",
    r"""eliminated~\cite{loutey_paper12}.  The remaining $7.6\%$ gap is a
basis-space limit rather than an integration error---and, as Paper~12
now establishes, a \emph{configuration}-space one.  That basis is
$\varphi$-independent, so it spans only $m_1 = m_2 = 0$;\ but a
$^1\Sigma_g^+$ state constrains only the total $M = m_1 + m_2$, so the
$\pi^2$ and $\delta^2$ configurations---which carry the angular
correlation---are absent by construction.  Restoring them in the same
basis with the same algebraic $V_{ee}$ reaches $99.09\%$ of $D_e$, a
gain of $11.6$~mHa against the $0.34$~mHa available from tripling the
$\sigma$ basis~\cite{loutey_paper12}.

\emph{Withdrawn.}  Earlier versions of this synthesis reported
Paper~12's original diagnosis---that the gap was the electron--electron
cusp demanding non-analytic $r_{12}^{1/2}$ and $r_{12}\ln r_{12}$
terms, and that this identified hyperspherical coordinates as the next
natural geometry.  Paper~12 has withdrawn both.  The cusp remains the
framework's canonical embedding-tier transcendental in Paper~18's
taxonomy~\cite{loutey_paper18}---that classification is about the cusp
itself and does not depend on Paper~12's H$_2$ residual---but it is not
what capped this calculation.""",
    "synthesis: correct the Paper 12 diagnosis")

# ------------------------------------------------------ synthesis, locus 2
edit(
    SYN,
    r"""channels; a complete-basis extrapolation reaches
$\sim\!97\%$)~\cite{loutey_paper15}.  This exceeds the $92.4\%$ that
prolate spheroidal CI with algebraic Neumann integrals achieved
(Paper~12), confirming the cusp-resolution advantage of the
hyperspherical geometry.""",
    r"""channels; a complete-basis extrapolation reaches
$\sim\!97\%$)~\cite{loutey_paper15}.

\emph{Withdrawn.}  Earlier versions read this as exceeding the
$92.4\%$ of prolate spheroidal CI (Paper~12) and confirming a
cusp-resolution advantage for the hyperspherical geometry.  The
comparison was not like-for-like:\ the $96.0\%$ includes $\pi$
channels and the $92.4\%$ does not.  At matched angular content the
ordering reverses---Paper~12's own basis with $|m| \le 1$ reaches
$99.09\%$---so no coordinate-system advantage is claimed.  Both
papers agree on the physics that matters here:\ the $m \ne 0$
channels carry the angular correlation, and neither geometry
gets it without them.""",
    "synthesis: withdraw the cusp-resolution advantage")

by_path = {}
for path, old, new, label in EDITS:
    by_path.setdefault(path, []).append((old, new, label))

failed, applied = [], []
for path, items in by_path.items():
    with io.open(path, encoding="utf-8") as fh:
        text = fh.read()
    for old, new, label in items:
        if old in text:
            text = text.replace(old, new, 1)
            applied.append(label)
        else:
            failed.append(label)
    if not failed:
        with io.open(path, "w", encoding="utf-8") as fh:
            fh.write(text)

if failed:
    print("FAILED TO MATCH:")
    for f in failed:
        print("  -", f)
    sys.exit(1)

print("applied %d edits" % len(applied))
for a in applied:
    print("  +", a)
