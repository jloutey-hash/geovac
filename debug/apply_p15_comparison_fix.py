"""Paper 15: correct the Paper-12 comparison, which was sigma+pi against
sigma-only and reverses when matched.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
Evidence: debug/sprint_tmr_method_memo.md; Paper 12 Sec. "Restoring the
Azimuthal Channels".
"""
import io
import sys

PATH = "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex"
EDITS = []


def edit(old, new, label):
    EDITS.append((old, new, label))


# ---------------------------------------------------------------- abstract
edit(
    r"""${\sim}95\%$, itself above 92.4\%),
exceeding the 92.4\% achieved by prolate spheroidal CI with algebraic
Neumann integrals (Paper~12).""",
    r"""${\sim}95\%$).
\textbf{[SCOPE]} This figure was previously compared against the
92.4\% of the prolate spheroidal CI with algebraic Neumann integrals
(Paper~12) and the difference read as an advantage of hyperspherical
coordinates.  That comparison was not like-for-like:\ Paper~12's
92.4\% is a $\sigma$-only result, its basis carrying no azimuthal
dependence, whereas the Level~4 figure includes $\pi$ channels.  At
matched angular content the ordering reverses---Paper~12's own basis
with $|m| \le 1$ reaches 99.09\%---so no coordinate-system advantage
is claimed here.""",
    "abstract: withdraw the coordinate-advantage comparison")

# ------------------------------------------------------------- introduction
edit(
    r"""Paper~12~\cite{paper12} demonstrated that prolate spheroidal CI with
algebraic Neumann $V_{ee}$ integrals saturates at 92.4\% of the exact
dissociation energy $D_e$.  The remaining 7.6\% gap was attributed to the
electron--electron cusp, which is a coordinate singularity in
$(\xi_1,\eta_1,\xi_2,\eta_2)$ space and cannot be captured by any
polynomial basis in those coordinates.""",
    r"""Paper~12~\cite{paper12} demonstrated that prolate spheroidal CI with
algebraic Neumann $V_{ee}$ integrals saturates at 92.4\% of the exact
dissociation energy $D_e$ \emph{in a $\sigma$-only basis}.  That
paper originally attributed the remaining 7.6\% gap to the
electron--electron cusp;\ it has since withdrawn that attribution and
identified the gap as the absent $m \ne 0$ configurations, which,
restored in the same basis, take it to 99.09\%.""",
    "intro: restate Paper 12's status")

edit(
    r"""We demonstrate that this approach recovers 96.0\% of $D_e$ at
$l_{\max}=6$ with $\sigma{+}\pi$ channels (61~channels) using a
variational 2D solver with Schwartz cusp correction, exceeding
Paper~12 and confirming the cusp-resolution advantage.""",
    r"""We demonstrate that this approach recovers 96.0\% of $D_e$ at
$l_{\max}=6$ with $\sigma{+}\pi$ channels (61~channels) using a
variational 2D solver with Schwartz cusp correction.  \textbf{[SCOPE]}
Earlier versions read this as exceeding Paper~12 and confirming a
cusp-resolution advantage for hyperspherical coordinates;\ since the
Paper~12 figure it was compared against is $\sigma$-only, that
inference is withdrawn.  The Level~4 result stands on its own terms
as a coupled-channel treatment of the two-center, two-electron
problem.""",
    "intro: withdraw the cusp-resolution advantage claim")

# ------------------------------------------------ the comparison subsection
edit(
    r"""The $\sigma$-only result (87.0\%) falls below Paper~12, but adding
$\pi$~channels pushes the Level~4 result to 94.1\%---a
1.7~percentage-point improvement over Paper~12.  This confirms that
the hyperspherical coordinates' cusp resolution provides a genuine
advantage, but only when the angular basis includes both $\sigma$ and
$\pi$ partial waves.""",
    r"""The $\sigma$-only result (87.0\%) falls below Paper~12, and adding
$\pi$~channels pushes the Level~4 result to 94.1\%.

\textbf{[SCOPE]} Earlier versions read that 1.7-percentage-point
difference as evidence that hyperspherical coordinates resolve the
cusp better than prolate spheroidal ones.  The comparison does not
support it, because the two entries do not contain the same physics:
the 94.1\% includes $\pi$ channels and the 92.4\% does not---Paper~12's
basis is $\varphi$-independent and spans only $m_1 = m_2 = 0$.  Run
at matched angular content the ordering reverses.  Paper~12's own
basis with $|m| \le 1$ reaches 99.09\%, against 94.1\% here at
$l_{\max}=4$ and 96.0\% at $l_{\max}=6$ with a cusp correction;\ and a
grid-based prolate spheroidal calculation reaches
99.97\%~\cite{tao_mccurdy_rescigno2010}.  What the table below
establishes is that $\pi$ channels are worth ${\sim}7$ percentage
points \emph{in this geometry}---which is true, and is the same
lesson Paper~12 learned in its own.  No claim of a coordinate-system
advantage survives.""",
    "comparison subsection: withdraw and state the reversal")

# --------------------------------------------------- pi-channel discussion
edit(
    r"""This improvement is essential: at $l_{\max}=4$, the $\sigma$-only
result (87.0\%) falls below Paper~12's 92.4\%, while the
$\sigma{+}\pi$ result (94.1\%) exceeds it.""",
    r"""This improvement is essential.  It is also, as Paper~12 has since
established in prolate spheroidal coordinates, the dominant missing
piece there too:\ the same $\pi$ configurations are worth 11.6~mHa
in that geometry.  Comparisons between the two papers are therefore
meaningful only at matched $m_{\max}$.""",
    "pi-channel paragraph: matched-comparison caveat")

# ------------------------------------------------------------- conclusion
edit(
    r"""  \item $\pi$-orbital channels ($m \ne 0$) are essential, providing
    ${\sim}7$~percentage points of improvement at each $l_{\max}$.
    Without $\pi$ channels, Level~4 ($\sigma$-only, 87.0\%) falls
    below Paper~12 (92.4\%); with them, it exceeds it (94.1\%).""",
    r"""  \item $\pi$-orbital channels ($m \ne 0$) are essential, providing
    ${\sim}7$~percentage points of improvement at each $l_{\max}$.
    The same holds in prolate spheroidal coordinates, where the
    corresponding configurations are worth 11.6~mHa (Paper~12);\
    comparisons between the two geometries are meaningful only at
    matched $m_{\max}$, and at matched $m_{\max}$ this paper claims
    no advantage over Paper~12.""",
    "conclusion: matched-comparison caveat")

# ------------------------------------------------------------ bibliography
edit(
    r"""\bibitem{paper12}""",
    r"""\bibitem{tao_mccurdy_rescigno2010}
L.~Tao, C.~W. McCurdy, and T.~N. Rescigno,
``Grid-based methods for diatomic quantum scattering problems.
III. Double photoionization of molecular hydrogen in prolate
spheroidal coordinates,''
\emph{Phys. Rev. A} \textbf{82}, 023423 (2010).

\bibitem{paper12}""",
    "bibliography: TMR 2010")

with io.open(PATH, encoding="utf-8") as fh:
    text = fh.read()

applied, failed = [], []
for old, new, label in EDITS:
    if old in text:
        text = text.replace(old, new, 1)
        applied.append(label)
    else:
        failed.append(label)

if failed:
    print("FAILED TO MATCH:")
    for f in failed:
        print("  -", f)
    sys.exit(1)

with io.open(PATH, "w", encoding="utf-8") as fh:
    fh.write(text)

print("applied %d edits" % len(applied))
for a in applied:
    print("  +", a)
