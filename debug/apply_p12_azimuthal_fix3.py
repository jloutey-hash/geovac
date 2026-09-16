"""Paper 12 azimuthal correction, part 3: hierarchy, conclusion, pointers,
and the two new bibliography entries.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

PATH = "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex"
EDITS = []


def edit(old, new, label):
    EDITS.append((old, new, label))


# --------------------------------------------- convergence-section pointer
edit(
    r"""The diminishing returns beyond $N = 27$ reflect the fact that
additional $\sigma$-type basis functions with higher $\xi$
and $\eta$ powers cannot access the missing physics
(Section~\ref{sec:gap}).""",
    r"""The diminishing returns beyond $N = 27$ reflect the fact that
additional $\sigma$-type basis functions with higher $\xi$
and $\eta$ powers cannot access the missing physics---which is
not in the $\sigma$ sector at all
(Sections~\ref{sec:gap} and~\ref{sec:azimuthal}).""",
    "convergence: pointer to the corrected diagnosis")

# ------------------------------------------------------ hierarchy items 3/4
edit(
    r"""  \item \textbf{Two-center, two-electron} (H$_2$):
    Prolate spheroidal $\times$ prolate spheroidal with
    Neumann-algebraic $V_{ee}$.  One-electron physics is
    exact; electron correlation is resolved to 92.4\%.
    \emph{This paper.}
  \item \textbf{Three-body coalescence}: Hyperspherical
    coordinates $(\rho, \theta, \hat{\Omega})$ for the
    interelectron cusp.  \emph{Identified but not yet
    implemented.}
\end{enumerate}

Each level in this hierarchy addresses a qualitatively different
physical singularity: the nuclear Coulomb singularity ($1/r$),
the two-center potential ($1/r_A + 1/r_B$), the electron-electron
Coulomb singularity ($1/r_{12}$), and the coalescence cusp
(non-analytic $r_{12}$ dependence).  The GeoVac program
identifies each singularity's natural coordinate system and
discretizes the Laplace--Beltrami operator in that system.""",
    r"""  \item \textbf{Two-center, two-electron} (H$_2$):
    Prolate spheroidal $\times$ prolate spheroidal with
    Neumann-algebraic $V_{ee}$.  One-electron physics is
    exact; electron correlation is resolved to 92.4\% in the
    $\sigma$ sector and 99.09\% once the azimuthal channels are
    restored (Sec.~\ref{sec:azimuthal}).  \emph{This paper.}
\end{enumerate}

Each level in this hierarchy addresses a qualitatively different
physical singularity: the nuclear Coulomb singularity ($1/r$),
the two-center potential ($1/r_A + 1/r_B$), and the
electron-electron Coulomb singularity ($1/r_{12}$).  The GeoVac
program identifies each singularity's natural coordinate system and
discretizes the Laplace--Beltrami operator in that system.

\textbf{[WITHDRAWN]} A fourth entry appeared here in earlier
versions---``three-body coalescence: hyperspherical coordinates for
the interelectron cusp, identified but not yet implemented''---
presented as a consequence of this paper's residual.  It is removed,
because the residual was misdiagnosed (Sec.~\ref{sec:gap}).  This
is a statement about what \emph{this paper} licenses, not about
hyperspherical coordinates, whose use for two-electron atoms rests
on separate grounds in Paper~13~\cite{loutey_paper13}.  Where the
two coordinate systems have been compared at matched angular
content, prolate spheroidal is the more accurate for
H$_2$:\ 99.09\% here, against 96.0\% for the molecule-frame
hyperspherical treatment at $l_{\max} = 6$ with a cusp
correction~\cite{loutey_paper15}, and 99.97\% for the grid-based
prolate spheroidal calculation of
Ref.~\cite{tao_mccurdy_rescigno2010}.""",
    "hierarchy: drop the fourth level, add the matched comparison")

# ----------------------------------------------------------- conclusion
edit(
    r"""  \item \textbf{The 7.6\% gap} is precisely diagnosed: it is
    a basis completeness limit, not an integration error.  The
    electron-electron cusp requires non-analytic terms
    ($r_{12}^{1/2}$, $r_{12}\ln r_{12}$) that no polynomial
    prolate spheroidal basis can represent.  This identifies
    hyperspherical coordinates as the next natural geometry
    in the GeoVac hierarchy.
\end{enumerate}""",
    r"""  \item \textbf{The 7.6\% gap} is a basis-space limit, not an
    integration error---and specifically a \emph{configuration}
    limit:\ the $\varphi$-independent basis spans only
    $m_1 = m_2 = 0$, while a $^1\Sigma_g^+$ state constrains only
    $M = m_1 + m_2$.  The $\pi^2$ and $\delta^2$ configurations it
    excludes carry the angular correlation.
  \item \textbf{99.09\% of exact $D_e$} on restoring the azimuthal
    channels at $|m| \le 1$ in the same basis with the same
    algebraic $V_{ee}$ (Sec.~\ref{sec:azimuthal}), a gain of
    $11.6$~mHa, against $0.34$~mHa from tripling the $\sigma$ basis.
    The $\sigma$-only value reproduces independently in a Gaussian
    basis to $0.2$~mHa, and the extended value to $0.01$
    percentage points.
\end{enumerate}

\textbf{[WITHDRAWN]} Earlier versions of this paper concluded that
the residual was the electron-electron cusp demanding non-analytic
$r_{12}^{1/2}$ and $r_{12}\ln r_{12}$ terms, and read that as
identifying hyperspherical coordinates as the next natural geometry.
Both are withdrawn:\ the calculations above contain no such term and
reach 99\%, and the grid-based prolate spheroidal treatment of
Ref.~\cite{tao_mccurdy_rescigno2010} reaches $0.05$~mHa in these
coordinates with a polynomial angular basis.""",
    "conclusion: item 5 corrected, item 6 added, withdrawal stated")

# ------------------------------------------------------------ bibliography
edit(
    r"""\bibitem{Fock1935}
V.~Fock,
``Zur Theorie des Wasserstoffatoms,''""",
    r"""\bibitem{loutey_paper13}
J.~Loutey,
``The Hyperspherical Lattice: Two-Electron Atoms as Coupled
Channel Graphs,''
GeoVac Paper~13 (2026).

\bibitem{loutey_paper15}
J.~Loutey,
``Molecule-Frame Hyperspherical Coordinates: the Level-4 Geometry
for Two-Center, Two-Electron Systems,''
GeoVac Paper~15 (2026).

\bibitem{tao_mccurdy_rescigno2010}
L.~Tao, C.~W. McCurdy, and T.~N. Rescigno,
``Grid-based methods for diatomic quantum scattering problems.
III. Double photoionization of molecular hydrogen in prolate
spheroidal coordinates,''
\emph{Phys. Rev. A} \textbf{82}, 023423 (2010).

\bibitem{Fock1935}
V.~Fock,
``Zur Theorie des Wasserstoffatoms,''""",
    "bibliography: paper 13, paper 15, TMR 2010")

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
