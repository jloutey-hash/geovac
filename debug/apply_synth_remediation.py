"""group2 synthesis: the C9 reviewer's findings F1-F7.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

S = "papers/synthesis/group2_quantum_chemistry_synthesis.tex"
EDITS = []


def edit(old, new, label):
    EDITS.append((old, new, label))


# ---- F2: 'Restoring them' drops the |m|<=1 qualifier and the delta limit.
edit(
    r"""correlation---are absent by construction.  Restoring them in the same
basis with the same algebraic $V_{ee}$ reaches $99.1\%$ of $D_e$, a
gain of $11.6$~mHa against the $0.34$~mHa available from tripling the
$\sigma$ basis~\cite{loutey_paper12}.""",
    r"""correlation---are absent by construction.  Restoring them at
$|m| \le 1$ in the same basis, with the same Neumann kernel, reaches
${\sim}99.1\%$ of $D_e$---a gain of $11.6$~mHa against the $0.34$~mHa
available from tripling the $\sigma$ basis~\cite{loutey_paper12}.  The
$|m| = 2$ sector was not obtained there, so the $\delta$ contribution
(${\approx}0.5$~mHa) rests on an independent Gaussian-basis route;\ and
the quoted figure carries a stability envelope, the basis being strongly
linearly dependent at that size.""",
    "F2: |m|<=1 qualifier, delta limit and the envelope restored")

# ---- F1: the Paper 18 independence sentence was a rescue; P18 is now fixed.
edit(
    r"""Paper~12 has withdrawn both.  The cusp remains the
framework's canonical embedding-tier transcendental in Paper~18's
taxonomy~\cite{loutey_paper18}---that classification is about the cusp
itself and does not depend on Paper~12's H$_2$ residual---but it is not
what capped this calculation.""",
    r"""Paper~12 has withdrawn both, and the correction reached
further than this section:\ Paper~18's Level-4
subsection~\cite{loutey_paper18} had built an \emph{irreducibility}
argument on the same saturation---that the product basis is
``structurally blind to three-body coalescence''---and that argument has
been withdrawn with it, taking the qualitative Level-3/Level-4
distinction with it.  What survives is narrower and genuinely
independent:\ the cusp keeps its place as the framework's embedding-tier
transcendental in Paper~18's taxonomy, which is a statement about the
cusp in any product basis and never rested on the H$_2$ residual.  It is
simply not what capped this calculation.""",
    "F1: the rescue sentence replaced by what actually happened")

# ---- F3: dangling anaphora -- 'that geometry' pointed at deleted text.
edit(
    r"""Paper~13~\cite{loutey_paper13} introduces that geometry.  In""",
    r"""Paper~13~\cite{loutey_paper13} introduces hyperspherical coordinates
for helium---on their own grounds, the atomic coalescence being
genuinely three-body, rather than as a consequence of the H$_2$
residual above.  In""",
    "F3: dangling anaphora repaired with the independent grounds")

# ---- F7: 'here' dropped, turning a paper-local abstention into a verdict.
edit(
    r"""ordering reverses---Paper~12's own basis with $|m| \le 1$ reaches
$99.1\%$---so no coordinate-system advantage is claimed.  Both
papers agree on the physics that matters here:\ the $m \ne 0$
channels carry the angular correlation, and neither geometry
gets it without them.""",
    r"""ordering reverses---Paper~12's own basis with $|m| \le 1$ reaches
${\sim}99.1\%$.  Neither paper claims a coordinate-system advantage
over the other, and the two published figures are not a matched pair in
any case:\ the $96.0\%$ is $l_{\max} = 6$ with a Schwartz cusp
correction and a different solver class.  What both papers do agree on
is the physics that matters here:\ the $m \ne 0$ channels carry the
angular correlation, and neither geometry gets it without them.""",
    "F7: paper-local abstention no longer stated as a corpus verdict")

with io.open(S, encoding="utf-8") as fh:
    t = fh.read()

# ---- F4: 96.0% presented as variational at four loci.
n_f4 = 0
for old, new in [
    (r"recovers 96.0\% of $D_e$ for $\mathrm{H}_2$",
     r"recovers 96.0\% of $D_e$ for $\mathrm{H}_2$ (with a Schwartz cusp correction;\ ${\sim}95\%$ pure-variational)"),
    (r"solver that avoids the adiabatic approximation, recovers $96.0\%$ of",
     r"solver that avoids the adiabatic approximation, recovers $96.0\%$ (cusp-corrected;\ ${\sim}95\%$ pure-variational) of"),
]:
    if old in t:
        t = t.replace(old, new, 1)
        n_f4 += 1

applied, failed = [], []
for old, new, label in EDITS:
    if old in t:
        t = t.replace(old, new, 1)
        applied.append(label)
    else:
        failed.append(label)

with io.open(S, "w", encoding="utf-8") as fh:
    fh.write(t)

print("applied %d edits (+%d F4 tier disclosures)" % (len(applied), n_f4))
for a in applied:
    print("  +", a)
if failed:
    print("UNMATCHED (%d):" % len(failed))
    for f in failed:
        print("  -", f)
    sys.exit(1)
