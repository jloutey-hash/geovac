"""Apply the azimuthal-channel correction to Paper 12.

Per memory rule feedback_no_heredoc_backslashes: LaTeX edits go through a
Write-tool script file, never a bash heredoc or inline python -c, because both
halve backslashes.

Evidence: debug/sprint_tmr_method_memo.md (canonical), plus
  debug/prolate_ci_general_m.py   general-m Neumann in Paper 12's own basis
  debug/p12_m_channel_probe.py    independent Gaussian-basis route
  debug/neumann_kernel_check.py   pointwise check of the kernel itself
"""
import io
import sys

PATH = "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex"

EDITS = []


def edit(old, new, label):
    EDITS.append((old, new, label))


# ------------------------------------------------------------------ abstract
edit(
    r"""functions.  The Neumann result plateaus at 92.4\%, confirming that
$V_{ee}$ integration error has been eliminated.  The remaining
7.6\% gap is diagnosed as a one-electron basis completeness limit:
the electron-electron cusp requires non-analytic terms
($r_{12}^{1/2}$, $r_{12}\ln r_{12}$) that no polynomial prolate
spheroidal basis can represent.  This diagnosis identifies the next
natural geometry in the GeoVac program: hyperspherical coordinates
for the three-body coalescence.""",
    r"""functions.  The Neumann result plateaus at 92.4\%, confirming that
$V_{ee}$ integration error has been eliminated.  \textbf{[MEASURED]}
The remaining 7.6\% gap is a basis-space limit, and we identify which
one.  The basis used here is $\varphi$-independent, so it spans only
$m_1 = m_2 = 0$ configurations; but a $^1\Sigma_g^+$ state constrains
the \emph{total} $M = m_1 + m_2$ to zero, not each $m_i$
separately, so the $\pi^2$, $\delta^2$, \ldots\ configurations---every
one of them $^1\Sigma_g^+$---are absent by construction.  Restoring
them in the same basis with the same algebraic $V_{ee}$ gives
$99.09\%$ of $D_e$ at $|m| \le 1$, a gain of $11.6$~mHa that no
growth along the $\sigma$ axis reproduces (tripling the $\sigma$
basis is worth $0.34$~mHa).  \textbf{[WITHDRAWN]} This paper
previously attributed the gap to the electron-electron cusp requiring
non-analytic $r_{12}^{1/2}$ and $r_{12}\ln r_{12}$ terms, and read
that as motivation for hyperspherical coordinates.  Neither holds:
no such term appears in the calculation above, nor in the grid-based
prolate spheroidal treatment of Tao, McCurdy and Rescigno, which
reaches $0.05$~mHa in these same coordinates.""",
    "abstract: diagnosis corrected")

# --------------------------------------------------------------- introduction
edit(
    r"""In the prolate spheroidal
Hylleraas-type CI framework, these integrals are six-dimensional
(three coordinates per electron, with azimuthal symmetry reducing
the effective dimensionality).""",
    r"""In the prolate spheroidal
Hylleraas-type CI framework, these integrals are six-dimensional
(three coordinates per electron).  The $\sigma$-only basis adopted
below reduces the effective dimensionality to four;
Sec.~\ref{sec:gap} shows that this is a \emph{restriction} on the
configuration space, not a symmetry of the $^1\Sigma_g^+$ state, and
measures what it costs.""",
    "intro: azimuthal reduction is a restriction, not a symmetry")

edit(
    r"""Section~\ref{sec:gap} diagnoses the origin of the remaining
7.6\% gap.""",
    r"""Section~\ref{sec:gap} diagnoses the origin of the remaining
7.6\% gap and measures its removal.""",
    "intro: roadmap")

# ------------------------------------------------------- the Neumann kernel
edit(
    r"""\begin{equation}
  \frac{1}{r_{12}} = \frac{2}{R}\sum_{l=0}^{\infty}\sum_{m=0}^{l}
    (2-\delta_{m0})\frac{(l-m)!}{(l+m)!}\,
    P_l^m(\xi_<)\,Q_l^m(\xi_>)\,
    P_l^m(\eta_1)\,P_l^m(\eta_2)\,
    \cos m(\varphi_1 - \varphi_2),
  \label{eq:neumann_full}
\end{equation}
where $\xi_< = \min(\xi_1, \xi_2)$, $\xi_> = \max(\xi_1, \xi_2)$,
$P_l^m$ and $Q_l^m$ are associated Legendre functions of the first
and second kind, and $R$ is the internuclear distance.""",
    r"""\begin{equation}
  \frac{1}{r_{12}} = \frac{2}{R}\sum_{l=0}^{\infty}\sum_{m=0}^{l}
    (2-\delta_{m0})\,(-1)^m\,(2l+1)
    \left[\frac{(l-m)!}{(l+m)!}\right]^{2}
    P_l^m(\xi_<)\,Q_l^m(\xi_>)\,
    P_l^m(\eta_1)\,P_l^m(\eta_2)\,
    \cos m(\varphi_1 - \varphi_2),
  \label{eq:neumann_full}
\end{equation}
where $\xi_< = \min(\xi_1, \xi_2)$, $\xi_> = \max(\xi_1, \xi_2)$,
$P_l^m$ and $Q_l^m$ are associated Legendre functions of the first
and second kind, and $R$ is the internuclear distance.
\textbf{[MEASURED]} Earlier versions of this paper printed
Eq.~\eqref{eq:neumann_full} without the $(-1)^m$ and $(2l+1)$
factors and with the factorial ratio unsquared.  That form does not
reproduce $1/r_{12}$---checked pointwise, it diverges---whereas the
expression above agrees with $1/|\mathbf{r}_1-\mathbf{r}_2|$ to
$2\times10^{-6}$ by $l = 10$ at well-separated $\xi$.  The error was
invisible for a decade of use because every calculation in this
paper used only the $m = 0$ specialisation
Eq.~\eqref{eq:neumann_sigma}, in which all three discrepancies
vanish identically.""",
    "Eq. neumann_full: prefactor corrected + provenance")

edit(
    r"""For $^1\Sigma_g^+$ states of homonuclear diatomics, the
wavefunction has $m = 0$ symmetry.  The azimuthal integration
$\int_0^{2\pi} d\varphi_1 \int_0^{2\pi} d\varphi_2$ projects
onto the $m = 0$ component:""",
    r"""The basis of Eq.~\eqref{eq:basis_function} below carries no
$\varphi$ dependence, so for \emph{that basis} the azimuthal
integration $\int_0^{2\pi} d\varphi_1 \int_0^{2\pi} d\varphi_2$
retains only the $m = 0$ component:""",
    "sec:neumann: restriction stated correctly")

edit(
    r"""Equation~\eqref{eq:neumann_sigma} is the key identity: it
separates the two-electron Coulomb kernel into a sum of products
of functions of individual coordinates.""",
    r"""Equation~\eqref{eq:neumann_sigma} is the key identity for the
$\sigma$-only sector: it separates the two-electron Coulomb kernel
into a sum of products of functions of individual coordinates.
\textbf{[SCOPE]} It is \emph{not} the whole kernel for a
$^1\Sigma_g^+$ state.  Retaining only $m = 0$ is a restriction to
$\sigma$ configurations, not a consequence of $\Sigma$ symmetry: the
$\Sigma$ label fixes $M = m_1 + m_2 = 0$, which $\pi^2$ and
$\delta^2$ configurations satisfy with $m_2 = -m_1 \ne 0$.  The
terms with $m \ne 0$ in Eq.~\eqref{eq:neumann_full} are precisely
what couples those configurations to the $\sigma$ sector.
Sections~\ref{sec:gap} and~\ref{sec:azimuthal} measure the cost of
dropping them.""",
    "sec:neumann: scope note after neumann_sigma")

# ---------------------------------------------------------------- r12 formula
edit(
    r"""\begin{equation}
  r_{12}^2 = \frac{R^2}{4}\bigl[
    (\xi_1^2-1)(1-\eta_1^2) + (\xi_2^2-1)(1-\eta_2^2)
    - 2\sqrt{(\xi_1^2-1)(1-\eta_1^2)(\xi_2^2-1)(1-\eta_2^2)}
    + (\xi_1\eta_1 - \xi_2\eta_2)^2
  \bigr].
  \label{eq:r12_prolate}
\end{equation}""",
    r"""\begin{equation}
  r_{12}^2 = \frac{R^2}{4}\bigl[
    (\xi_1^2-1)(1-\eta_1^2) + (\xi_2^2-1)(1-\eta_2^2)
    - 2\sqrt{(\xi_1^2-1)(1-\eta_1^2)(\xi_2^2-1)(1-\eta_2^2)}
      \,\cos(\varphi_1-\varphi_2)
    + (\xi_1\eta_1 - \xi_2\eta_2)^2
  \bigr].
  \label{eq:r12_prolate}
\end{equation}
(The $\cos(\varphi_1-\varphi_2)$ factor was omitted in earlier
versions of this equation, which therefore held only at
$\varphi_1 = \varphi_2$.)""",
    "Eq. r12_prolate: restore cos(dphi)")

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

print("applied %d edits to %s" % (len(applied), PATH))
for a in applied:
    print("  +", a)
