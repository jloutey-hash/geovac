"""DELTA remediation, pass 2: the claims-reviewer findings on Papers 12/13/15.

All of these are defects in THIS sprint's own remediation.  Verified against
primary text by the PM before editing.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

EDITS = []
P12 = "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex"
P13 = "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex"
P15 = "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex"


def edit(path, old, new, label):
    EDITS.append((path, old, new, label))


# ---- F1 LARGE: the withdrawn comparison, sign-flipped, self-contradicting.
edit(
    P12,
    r"""is a statement about what \emph{this paper} licenses, not about
hyperspherical coordinates, whose use for two-electron atoms rests
on separate grounds in Paper~13~\cite{loutey_paper13}.  Where the
two coordinate systems have been compared at matched angular
content, prolate spheroidal is the more accurate for
H$_2$:\ 99.09\% here, against 96.0\% for the molecule-frame
hyperspherical treatment at $l_{\max} = 6$ with a Schwartz cusp
correction~\cite{loutey_paper15}, and 99.97\% for the grid-based
prolate spheroidal calculation of
Ref.~\cite{tao_mccurdy_rescigno2010}.""",
    r"""is a statement about what \emph{this paper} licenses, not about
hyperspherical coordinates, whose use for two-electron atoms rests
on separate grounds in Paper~13~\cite{loutey_paper13}.  What the
present result removes is the \emph{evidence} the fourth entry
rested on, not a verdict in its favour:\ comparisons between the two
geometries are meaningful only at matched $m_{\max}$, and at matched
$m_{\max}$ neither this paper nor Paper~15~\cite{loutey_paper15}
claims an advantage over the other.  (The two published numbers,
99.09\% here and 96.0\% there, are not a matched pair --- the latter
is $l_{\max} = 6$ with a Schwartz cusp correction and a different
solver class --- and neither is the 99.97\% of
Ref.~\cite{tao_mccurdy_rescigno2010}, which is a grid method.  The
honest statement is that the prolate spheroidal geometry is not
disqualified for H$_2$ by a cusp it was said to be unable to
represent.)""",
    "F1: coordinate-system comparison replaced by the matched-m_max standard")

# ---- F2 LARGE: quadrature-freeness now false for the headline, at 4 surfaces.
edit(
    P12,
    r"""integrals are exact within the Neumann truncation order $L$, with
no six-dimensional quadrature and no fitted parameters; the sole
residual numerical integration is a single one-dimensional adaptive
quadrature for the log-singular second-kind Legendre moment $B_l$.""",
    r"""integrals are exact within the Neumann truncation order $L$, with
no six-dimensional quadrature and no fitted parameters; in the
$\sigma$ sector the sole residual numerical integration is a single
one-dimensional adaptive quadrature for the log-singular second-kind
Legendre moment $B_l$.  (The $|m| \le 1$ extension reported below
uses the same kernel but evaluates its ordered-$\xi$ integral by a
two-dimensional spectral quadrature rather than by those
recurrences;\ the quadrature-free property is a $\sigma$-sector
statement throughout this paper.)""",
    "F2a: abstract quadrature claim scoped to sigma")

edit(
    P12,
    r"""Restoring
them in the same basis with the same algebraic $V_{ee}$ gives
$99.09\%$ of $D_e$ at $|m| \le 1$,""",
    r"""Restoring
them in the same basis with the same Neumann kernel gives
$99.09\%$ of $D_e$ at $|m| \le 1$,""",
    "F2b: abstract 'same algebraic V_ee' -> 'same Neumann kernel'")

edit(
    P12,
    r"""recurrence relations on quantum number labels, the sole residual
numerical step being a single one-dimensional adaptive quadrature
for the log-singular second-kind moment $B_l$.""",
    r"""recurrence relations on quantum number labels, the sole residual
numerical step in the $\sigma$ sector being a single
one-dimensional adaptive quadrature for the log-singular
second-kind moment $B_l$ (Sec.~\ref{sec:azimuthal} states what the
$|m| \ge 1$ extension costs).""",
    "F2c: Sec. VII.B scoped")

edit(
    P12,
    r"""six-dimensional quadrature for electron-electron repulsion, leaving
only a single one-dimensional quadrature for the log-singular $B_l$
moment.""",
    r"""six-dimensional quadrature for electron-electron repulsion, leaving
only a single one-dimensional quadrature for the log-singular $B_l$
moment in the $\sigma$ sector.""",
    "F2d: conclusion scoped")

edit(
    P12,
    r"""channels at $|m| \le 1$ in the same basis with the same
    algebraic $V_{ee}$ (Sec.~\ref{sec:azimuthal}), a gain of""",
    r"""channels at $|m| \le 1$ in the same basis with the same
    Neumann kernel (Sec.~\ref{sec:azimuthal}), a gain of""",
    "F2e: conclusion item 6 wording")

# ---- F9 SMALL: Eq. basis_mu has no phi, and 'recovers exactly' is a type error.
edit(
    P12,
    r"""\begin{equation}
  u_{j l \mu}(\xi,\eta,\varphi) =
    \xi^{\,j}\,\eta^{\,l}\,
    (\xi^2-1)^{\mu/2}(1-\eta^2)^{\mu/2}\,e^{-\alpha\xi},
  \label{eq:basis_mu}
\end{equation}""",
    r"""\begin{equation}
  u_{j l \mu}(\xi,\eta) =
    \xi^{\,j}\,\eta^{\,l}\,
    (\xi^2-1)^{\mu/2}(1-\eta^2)^{\mu/2}\,e^{-\alpha\xi},
  \qquad
  \Phi_\mu = u\,u'\cos\mu(\varphi_1-\varphi_2),
  \label{eq:basis_mu}
\end{equation}""",
    "F9a: basis_mu is a one-electron factor; the phi pairing made explicit")

edit(
    P12,
    r"""Setting $\mu = 0$ recovers
Eq.~\eqref{eq:basis_function} exactly.""",
    r"""At $\mu = 0$ the factor $u$ reduces to the one-electron factor of
Eq.~\eqref{eq:basis_function}, and $\Phi_0$ to its symmetrised
product, so the $\mu = 0$ sector is Paper~12's original basis.""",
    "F9b: 'recovers exactly' type mismatch corrected")

# ---- F10 SMALL: kept dimensions + the orthogonalisation threshold.
edit(
    P12,
    r"""$(1,1)$ &   6 & 74.19 &  12 & 81.63 \\
$(2,1)$ &  12 & 75.87 &  24 & 83.47 \\
$(2,2)$ &  27 & 92.25 &  54 & 98.96 \\
$(3,2)$ &  46 & 92.37 &  92 & 99.00 \\
$(3,3)$ &  72 & 92.42 & 144 & 99.09 \\""",
    r"""$(1,1)$ &   6 & 74.19 &  12 & 81.63 \\
$(2,1)$ &  12 & 75.87 &  24 & 83.47 \\
$(2,2)$ &  27 & 92.25 &  54 & 98.96 \\
$(3,2)$ &  46 (46) & 92.37 &  92 (76) & 99.00 \\
$(3,3)$ &  72 (65) & 92.42 & 144 (115) & 99.09 \\""",
    "F10a: kept dimensions added to the table")

edit(
    P12,
    r"""$\alpha$
  scanned, canonical orthogonalization throughout.
  $D_e^{\rm exact} = 0.174475$~Ha.}""",
    r"""$\alpha$
  scanned, canonical orthogonalization throughout (eigenvectors of
  $S$ with eigenvalue below $10^{-11}\lambda_{\max}$ discarded;\ the
  surviving dimension is in parentheses where it differs from $N$).
  $D_e^{\rm exact} = 0.174475$~Ha.}""",
    "F10b: orthogonalisation threshold stated")

# ---- F11 SMALL: duration language + a false scope statement.
edit(
    P12,
    r"""vanish identically.  The error was
invisible for a decade of use because every calculation in this
paper used only the $m = 0$ specialisation
Eq.~\eqref{eq:neumann_sigma}, in which all three discrepancies
vanish identically.""",
    r"""vanish identically.  The error went
undetected because every calculation in \emph{earlier versions} of
this paper used only the $m = 0$ specialisation
Eq.~\eqref{eq:neumann_sigma}, in which all three discrepancies
vanish identically;\ Sec.~\ref{sec:azimuthal} is the first use of
the general form.""",
    "F11: 'a decade of use' (C18 class) removed; scope corrected")

# ---- F8 SMALL: the italicised 'central thesis' is a surviving cusp misattribution.
edit(
    P12,
    r"""The Neumann algebraic approach at 92.2\% with $p = 0$ basis
functions outperforms the $r_{12}$ + numerical approach at
86.8\% with $p \leq 2$ basis functions.  This confirms
the central thesis: \emph{the cusp does not need $r_{12}^p$
basis functions---it needs exact $V_{ee}$ computed from
quantum number algebra.}""",
    r"""The Neumann algebraic approach at 92.2\% with $p = 0$ basis
functions outperforms the $r_{12}$ + numerical approach at
86.8\% with $p \leq 2$ basis functions.  What that comparison
establishes is narrower than it looks:\ \emph{an exact $V_{ee}$ beats
$r_{12}$ basis functions carrying grid error}.  It says nothing
about the cusp either way.  (Earlier versions drew the stronger
conclusion that ``the cusp does not need $r_{12}^p$ basis
functions---it needs exact $V_{ee}$''; that is withdrawn with the
rest of the cusp reading, and Sec.~\ref{sec:azimuthal_diagnosis}
gives the cusp's actual cost here.)""",
    "F8: 'central thesis' cusp misattribution withdrawn")

# ---- F3 LARGE: Paper 13's ABSTRACT still asserts both withdrawn claims.
edit(
    P13,
    r"""Paper~12
identified the electron-electron cusp as the limiting factor in
prolate spheroidal CI for H$_2$ and pointed to hyperspherical
coordinates as the next natural geometry for multi-electron systems.""",
    r"""Hyperspherical
coordinates are the natural setting for a genuine three-body
coalescence, where the nuclear and interelectronic singularities
meet at a point.  (Earlier versions took this motivation from
Paper~12's reading of its H$_2$ residual as an electron-electron
cusp limit;\ Paper~12 has since withdrawn that reading, and nothing
below depends on it.)""",
    "F3: Paper 13 abstract no longer asserts the withdrawn claims")

# ---- F5 SMALL: Paper 15 abstract claims an advantage it disclaims 12 lines later.
edit(
    P15,
    r"""The electron--electron cusp appears at $\alpha=\pi/4$, $\theta_{12}=0$,
independent of molecular geometry---a structural advantage over
prolate spheroidal coordinates, where it is a coordinate singularity.""",
    r"""The electron--electron cusp appears at $\alpha=\pi/4$, $\theta_{12}=0$,
independent of molecular geometry, where in prolate spheroidal
coordinates the same locus is not a coordinate surface.  That is a
structural difference in where the cusp sits, and---see the scope
note below---it does not translate into an accuracy advantage.""",
    "F5: P15 abstract structural clause no longer claims an advantage")

# ---- F6 SMALL: P15 table row + figure caption stage the withdrawn comparison.
edit(
    P15,
    r"""Paper~12 (Neumann $V_{ee}$)                  & 0.1612 & 92.4 \\""",
    r"""Paper~12 (Neumann $V_{ee}$, $\sigma$-only)   & 0.1612 & 92.4 \\
Paper~12 (Neumann $V_{ee}$, $|m|\le1$)       & 0.1729 & 99.1 \\""",
    "F6a: P15 comparison table labels the sigma-only row and adds the matched one")

edit(
    P15,
    r"""The dashed line marks Paper~12's 92.4\%.""",
    r"""The dashed line marks Paper~12's $\sigma$-only 92.4\% (its
  $|m|\le1$ value is 99.1\%, not a like-for-like reference for a
  $\sigma{+}\pi$ curve).""",
    "F6b: P15 figure caption dashed line scoped")

by_path = {}
for path, old, new, label in EDITS:
    by_path.setdefault(path, []).append((old, new, label))

applied, failed = [], []
for path, items in by_path.items():
    with io.open(path, encoding="utf-8") as fh:
        text = fh.read()
    ok = True
    for old, new, label in items:
        if old in text:
            text = text.replace(old, new, 1)
            applied.append(label)
        else:
            failed.append("%s   [%s]" % (label, path))
            ok = False
    if ok:
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
