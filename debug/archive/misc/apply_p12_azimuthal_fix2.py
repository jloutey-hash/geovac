"""Paper 12 azimuthal correction, part 2: sec:gap, new sec:azimuthal,
the hierarchy subsection, and the conclusion.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

PATH = "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex"
EDITS = []


def edit(old, new, label):
    EDITS.append((old, new, label))


# ------------------------------------------------- sec:gap, first subsection
edit(
    r"""The basis functions~\eqref{eq:basis_function} are products
of integer powers of $\xi$ and $\eta$ times a common
exponential $e^{-\alpha\xi}$.  These are smooth, analytic
functions of the single-electron coordinates.  No matter
how many such functions we include, they span only the
polynomial$\times$exponential subspace of $L^2$.

\subsection{The electron-electron cusp}""",
    r"""The basis functions~\eqref{eq:basis_function} are products
of integer powers of $\xi$ and $\eta$ times a common
exponential $e^{-\alpha\xi}$.  These are smooth, analytic
functions of the single-electron coordinates, and---the point that
matters here---they carry no $\varphi$ dependence at all.

\textbf{[MEASURED]} The convergence study itself shows that the gap
is not reachable along the directions the basis does span.  Going
from $N = 27$ to $N = 46$ to $N = 72$ moves the energy
$92.2\% \to 92.4\% \to 92.4\%$:\ $0.34$~mHa for $2.7\times$ the
functions.  Whatever is missing, more $\sigma$ functions do not
supply it.  Sections~\ref{sec:azimuthal_diagnosis}
and~\ref{sec:azimuthal} identify it and remove it.

\subsection{The electron-electron cusp is real, but it is not this
gap}
\label{sec:azimuthal_diagnosis}""",
    "sec:gap: saturation control + retitle")

# ----------------------------------------- close of the cusp subsection
edit(
    r"""$R^{1/2}$ and $R\ln R$ terms are fundamentally non-analytic:
no finite sum of polynomial$\times$exponential functions in
prolate spheroidal coordinates can represent them.""",
    r"""$R^{1/2}$ and $R\ln R$ terms are fundamentally non-analytic:
no finite sum of polynomial$\times$exponential functions in
prolate spheroidal coordinates can represent them \emph{exactly}.

\textbf{[WITHDRAWN]} Earlier versions of this paper concluded from
this that the cusp is what caps the calculation at 92.4\%.  It is
not, and the argument was never sound:\ non-analyticity implies slow
convergence of a partial-wave expansion, not a floor at 7.6\%.  Two
calculations settle it.  Tao, McCurdy and
Rescigno~\cite{tao_mccurdy_rescigno2010}, working in these same
prolate spheroidal coordinates with a polynomial angular basis, no
$r_{12}$ factors and no non-analytic functions, reach
$-1.17442$~Ha---$0.05$~mHa from exact---at angular truncation
$l_{\max} = 6$.  And the calculation reported in
Sec.~\ref{sec:azimuthal} below, run in the basis of
Eq.~\eqref{eq:basis_function} with the same algebraic $V_{ee}$ and
likewise no $r_{12}$ dependence, reaches $99.09\%$.  Neither could
exist if representing the cusp were the binding constraint at this
accuracy.  What the cusp does cost is the slow tail of the
partial-wave series;\ at $|m| \le 1$ and $l \le 3$ that residue is
under one percent of $D_e$.""",
    "sec:gap: withdraw the cusp-as-cause reading")

# ------------------------------------------------- the missing geometry
edit(
    r"""\subsection{The missing geometry}

The 7.6\% gap identifies a \emph{geometric} limitation.  The
prolate spheroidal coordinates $(\xi_1, \eta_1, \xi_2, \eta_2)$
are natural for the \emph{one-electron} nuclear attraction
problem (they separate the two-center Coulomb potential), but
they are not natural for the \emph{two-electron} coalescence
(they do not separate $1/r_{12}$).  The cusp is a property
of the \emph{interelectron} coordinate $r_{12}$, which in
prolate spheroidals is a complicated function of all four
variables:""",
    r"""\subsection{The missing channels}

\textbf{[MEASURED]} The 7.6\% gap is a \emph{configuration-space}
limitation, and an elementary one.  A $^1\Sigma_g^+$ state requires
the \emph{total} azimuthal quantum number $M = m_1 + m_2$ to vanish.
It does not require $m_1 = m_2 = 0$.  The configurations
$\pi^2$ ($m_1 = +1$, $m_2 = -1$), $\delta^2$ ($m_1 = +2$,
$m_2 = -2$), and so on all have $M = 0$ and are all
$^1\Sigma_g^+$;\ they carry the \emph{angular} part of the electron
correlation, the left--right and in--out parts being what a
$\sigma$-only space already describes.  The basis
Eq.~\eqref{eq:basis_function} is $\varphi$-independent and therefore
excludes every one of them, and Eq.~\eqref{eq:neumann_sigma} drops
exactly the kernel terms that would couple them in.  The two
restrictions are consistent with each other, which is why the
calculation converges cleanly to a wrong answer.

Section~\ref{sec:azimuthal} restores the channels and measures the
result:\ $92.4\% \to 99.09\%$ of $D_e$.  What remains after that is
ordinary basis incompleteness---higher $l$, higher $n$, and the
$\delta$ channels, worth a further $0.5$~mHa on an independent
Gaussian-basis estimate.

For completeness, the interelectron coordinate in prolate
spheroidals is a function of all five relative variables:""",
    "sec:gap: missing geometry -> missing channels")

edit(
    r"""The cusp $\partial\Psi/\partial r_{12} \neq 0$ at $r_{12} = 0$
cannot be represented by smooth functions of $\xi_1, \eta_1,
\xi_2, \eta_2$ individually.  It requires either explicit
$r_{12}$ dependence (the traditional approach, but hard to
integrate) or a coordinate system where the coalescence point
has a simple representation.

The natural geometry for the three-body coalescence problem
is the \emph{hyperspherical} coordinate system
$(\rho, \theta, \hat{\Omega})$, where $\rho = (r_1^2 + r_2^2)^{1/2}$
is the hyperradius and the cusp condition becomes a boundary
condition in $\theta$ at $\theta = \pi/4$ (the $r_1 = r_2$
manifold).  In this coordinate system, the Fock
expansion~\eqref{eq:fock_expansion} in $\rho^{1/2}$ and
$\rho\ln\rho$ becomes a local expansion in the natural
variable, amenable to lattice discretization.""",
    r"""Note the $\cos(\varphi_1-\varphi_2)$ dependence:\ the
interelectron distance is not a function of
$\xi_1, \eta_1, \xi_2, \eta_2$ alone.  A $\varphi$-independent basis
cannot resolve it, which is the same restriction stated in
coordinates.

\textbf{[WITHDRAWN]} Earlier versions of this paper closed this
section by identifying hyperspherical coordinates
$(\rho, \theta, \hat{\Omega})$ as the natural geometry the gap
called for, on the grounds that the coalescence cusp has a simple
representation there.  That inference rested on the cusp diagnosis
withdrawn above and does not survive it.  Hyperspherical coordinates
remain the right setting for genuine three-body coalescence
problems, and Paper~13~\cite{loutey_paper13} uses them for helium
on independent grounds;\ what is withdrawn is the claim that
\emph{this} calculation's residual motivates them.""",
    "sec:gap: withdraw the hyperspherical motivation")

# ------------------------------------------- NEW SECTION: the measurement
edit(
    r"""%% ========================================================================
%% VII. IMPLICATIONS FOR THE GEOVAC PROGRAM
%% ========================================================================

\section{Implications for the GeoVac Program}""",
    r"""%% ========================================================================
%% VI-B. RESTORING THE AZIMUTHAL CHANNELS
%% ========================================================================

\section{Restoring the Azimuthal Channels}
\label{sec:azimuthal}

\textbf{[MEASURED]} We extend the basis of
Eq.~\eqref{eq:basis_function} with an azimuthal quantum number,

\begin{equation}
  u_{j l \mu}(\xi,\eta,\varphi) =
    \xi^{\,j}\,\eta^{\,l}\,
    (\xi^2-1)^{\mu/2}(1-\eta^2)^{\mu/2}\,e^{-\alpha\xi},
  \label{eq:basis_mu}
\end{equation}
pairing $m_1 = +\mu$ with $m_2 = -\mu$ through
$\cos\mu(\varphi_1-\varphi_2)$---the $M = 0$, $^1\Sigma_g^+$
combination, which for $\mu = 1$ is the familiar
$\pi_x\pi_x + \pi_y\pi_y$.  Setting $\mu = 0$ recovers
Eq.~\eqref{eq:basis_function} exactly.  The overlap, kinetic and
nuclear-attraction matrices remain \emph{exact}---every integrand is
a polynomial in $\xi$ and $\eta$ times $e^{-2\alpha\xi}$, evaluated
against the moments $A_n$ of Sec.~\ref{sec:auxiliary} and
elementary $\eta$ moments---and all three are diagonal in $\mu$.
$V_{ee}$ uses the full kernel Eq.~\eqref{eq:neumann_full};\ its
$m \neq 0$ terms are what couple different $\mu$.

The $\eta$ integrals remain exact polynomial moments, because the
kernel's $(1-\eta^2)^{m/2}$ and the basis's $(1-\eta^2)^{\mu/2}$
always combine to an integer power.  The sum over $l$ terminates
exactly:\ integrating by parts $m$ times (the boundary terms vanish
since $(1-\eta^2)^{s}$ has a zero of order $s \ge m$ at
$\eta = \pm1$) shows the moment is identically zero for
$l > Q + 2s - m$, with $Q$ the combined $\eta$ power and
$s = (\mu_i + \mu_j + m)/2$.  That rule must be imposed
explicitly rather than left to cancellation:\ the Legendre
derivative coefficients reach $\sim\!10^{10}$, and their
floating-point residue multiplied by the radial integral produced a
$-3\times10^{8}$~Ha artifact before the cutoff was enforced.

\begin{table}[t]
\caption{H$_2$ at $R = 1.4011$~bohr with the azimuthal channels
  restored.  Same basis family, same algebraic $V_{ee}$, $\alpha$
  scanned, canonical orthogonalization throughout.
  $D_e^{\rm exact} = 0.174475$~Ha.}
\label{tab:azimuthal}
\begin{ruledtabular}
\begin{tabular}{crdrd}
  \multicolumn{1}{c}{$(j_{\max}, l_{\max})$}
  & \multicolumn{1}{c}{$N_{\sigma}$}
  & \multicolumn{1}{c}{$D_e\%$ ($\mu = 0$)}
  & \multicolumn{1}{c}{$N_{|m|\le1}$}
  & \multicolumn{1}{c}{$D_e\%$ ($|m| \le 1$)} \\
\hline
$(1,1)$ &   6 & 74.19 &  12 & 81.63 \\
$(2,1)$ &  12 & 75.87 &  24 & 83.47 \\
$(2,2)$ &  27 & 92.25 &  54 & 98.96 \\
$(3,2)$ &  46 & 92.37 &  92 & 99.00 \\
$(3,3)$ &  72 & 92.42 & 144 & 99.09 \\
\end{tabular}
\end{ruledtabular}
\end{table}

Table~\ref{tab:azimuthal} gives the result.  At the largest basis
the azimuthal channels are worth $+11.64$~mHa, taking
$92.42\% \to 99.09\%$ of $D_e$.  The $\mu = 0$ column reproduces the
Neumann column of Table~\ref{tab:convergence} to
$161$, $1.6$, $1.0$, $6.7$ and $58$~$\mu$Ha.

Three controls make the reading unambiguous.  \emph{First}, the
gain is not a basis-count effect:\ tripling the basis along the
$\sigma$ axis is worth $0.34$~mHa (Sec.~\ref{sec:gap}), while
opening the azimuthal axis at fixed $(j_{\max}, l_{\max})$ is worth
$11.6$~mHa.  \emph{Second}, an independent route agrees.  A full CI
in a Cartesian Gaussian basis---different functions, different
integrals, different code---restricted to $m = 0$ orbitals converges
to $92.34\%$, within $0.2$~mHa of this paper's $\sigma$-only value
in a completely unrelated basis;\ releasing $|m| \le 1$ takes it to
$99.10\%$, against $99.09\%$ here.  \emph{Third}, the general-$m$
kernel of Eq.~\eqref{eq:neumann_full} was checked pointwise against
$1/|\mathbf{r}_1-\mathbf{r}_2|$, agreeing to $2\times10^{-6}$.

\textbf{[SCOPE]} Two limits are worth stating plainly.  The
$|m| = 2$ sector was not obtained here:\ the fourth derivatives of
$Q_l$ near $\xi = 1$ lose precision in the present radial
quadrature, so the $\delta$ contribution ($\approx 0.5$~mHa) is
known only from the Gaussian route.  And the $\mu > 0$ radial
integrals are evaluated by a graded-panel spectral quadrature rather
than by the $A_l$, $B_l$, $X_l$ recurrences of
Sec.~\ref{sec:auxiliary};\ the quadrature-free property established
in this paper for $\sigma$ states is therefore \emph{not} yet
extended to $\mu > 0$.  Doing so requires generalising those three
auxiliary tables to associated Legendre functions, which is
well-defined but not attempted here.

\textbf{[MEASURED]} Finally, a caution about the $\sigma$ column
itself.  The overlap matrix of this basis family is strongly
linearly dependent at high powers and a single $\alpha$:\ measured
$\mathrm{cond}(S) = 2.6\times10^{14}$ at $(j,l) = (3,3)$, $\mu = 0$
(the $N = 72$ basis of Table~\ref{tab:convergence}), rising to
$2.0\times10^{16}$ once $\mu = 1$ doubles the basis.  Beyond
$\sim\!10^{16}$ a direct generalised eigensolve is not
trustworthy---at $N = 144$ it returned $-79$~Ha---so
Table~\ref{tab:azimuthal} uses canonical orthogonalization
throughout.  The $-1.161304$~Ha quoted in
Table~\ref{tab:convergence} sits $58\,\mu$Ha below the conditioned
value, i.e.\ the last two of its six decimals reflect linear
dependence rather than physics.  No conclusion in this paper turns
on them.

%% ========================================================================
%% VII. IMPLICATIONS FOR THE GEOVAC PROGRAM
%% ========================================================================

\section{Implications for the GeoVac Program}""",
    "NEW sec:azimuthal with measurement table")

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
