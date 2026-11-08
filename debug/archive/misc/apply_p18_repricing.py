"""Paper 18: re-price the Level-4 exchange-constant subsection.

Its argument rests on Paper 12's withdrawn diagnosis, in the strongest form
anywhere in the corpus, and uses it to carry a TAXONOMY classification.
Verified by the PM: Track M (debug/archive/tracks_misc/track_m_convergence.py)
imports `generate_basis` from `geovac.hylleraas` -- the phi-independent,
sigma-only basis -- and its table reproduces Paper 12's exactly.  So its
l_max scan and its extrapolation ("5% error would require l_max ~ 45") ran
along the ANGULAR axis while the actual limitation was the AZIMUTHAL one.

What is withdrawn: "structurally incomplete", "structurally blind to
three-body coalescence", "cannot represent at any practical basis size",
"irreducible", the 92.5% ceiling as a product-space property, the
"surpassing the prolate CI ceiling" comparison, and the qualitative Level-3
vs Level-4 distinction that rested on all of it.

What survives: the product-space approach is algebraic in the sigma sector;
mu(rho,R) is genuinely transcendental (computed by pointwise diagonalisation)
and keeps its place in the transcendental-content taxonomy; and the r12^2
observation is a real measurement, re-read as evidence about angular
correlation rather than about blindness.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

P18 = "papers/group3_foundations/paper_18_exchange_constants.tex"
EDITS = []


def edit(old, new, label):
    EDITS.append((old, new, label))


# ---- the subsection's opening claim
edit(
    r"""At Level~4 the compact product basis (prolate spheroidal CI) is fully
algebraic but structurally incomplete; the non-compact mol-frame
hyperspherical reparameterization introduces $\mu(\rho,R)$ as the
irreducible cost of accessing three-body coalescence physics that the
compact description cannot represent.""",
    r"""At Level~4 the compact product basis (prolate spheroidal CI) is fully
algebraic in its $\sigma$ sector, and the non-compact mol-frame
hyperspherical reparameterization introduces $\mu(\rho,R)$ as the
transcendental cost of that reparameterization.

\textbf{[WITHDRAWN 2026-09-14]} Earlier versions of this subsection said
more:\ that the product basis is \emph{structurally incomplete}, that
$\mu(\rho,R)$ is \emph{irreducible} because that basis is ``structurally
blind to three-body coalescence,'' and that the reparameterization
``accesses physics the algebraic product space cannot represent at any
practical basis size.''  All of that rested on a saturation which has since
been re-diagnosed, and none of it survives;\ the paragraphs below give the
corrected reading and what it costs this taxonomy.""",
    "P18: opening claim withdrawn")

# ---- Track M re-scoped
edit(
    r"""However, a systematic convergence study (Track~M) demonstrates
that this algebraic approach saturates at $\sim$92.5\% of the
exact dissociation energy $D_e$~\cite{loutey_track_m}.  The
angular basis convergence is rapid---92.3\% at $l_{\max} = 2$
(46 basis functions), 92.5\% at $l_{\max} = 4$ (114 basis
functions)---but successive improvements decay by a factor of
200 from $l = 1 \to 2$ to $l = 2 \to 3$, with the $l = 3 \to 4$
increment contributing only 0.04\%.  Extrapolation in the slow
regime predicts that 5\% error would require $l_{\max} \sim 45$
and 1\% error $l_{\max} \sim 113$, both impractical.  The
effective ceiling is 92.5\% $D_e$.""",
    r"""A systematic convergence study (Track~M) found that this approach
saturates at $\sim$92.5\% of the exact dissociation energy
$D_e$~\cite{loutey_paper12}:\ 92.3\% at $l_{\max} = 2$ (46 basis
functions), 92.5\% at $l_{\max} = 4$ (114), with successive improvements
decaying by a factor of 200 from $l = 1 \to 2$ to $l = 2 \to 3$ and the
$l = 3 \to 4$ increment contributing $0.04\%$.  Extrapolating that decay
predicted $l_{\max} \sim 45$ for $5\%$ error and $l_{\max} \sim 113$ for
$1\%$, both impractical.

\textbf{[MEASURED 2026-09-14]} That saturation is real but it is a
$\sigma$-\emph{sector} saturation, and the extrapolation ran along the
wrong axis.  Track~M scanned $l_{\max}$ in a $\varphi$-independent basis,
which spans only $m_1 = m_2 = 0$;\ a $^1\Sigma_g^+$ state constrains only
the \emph{total} $M = m_1 + m_2$, so every $\pi^2$ and $\delta^2$
configuration---which carry the angular correlation---was absent by
construction.  Opening that axis in the same basis, with the same
algebraic $V_{ee}$, reaches ${\sim}99.1\%$ of $D_e$
(Paper~12~\cite{loutey_paper12}, Sec.~``Restoring the Azimuthal
Channels'').  The $\sigma$ axis does saturate---growing it to $342$
functions still yields $92.40\%$---but it saturates because it is the
wrong axis, not because the product space is incomplete.""",
    "P18: Track M re-scoped as a sigma-sector scan on the wrong axis")

# ---- the 'smoking gun' re-read
edit(
    r"""The smoking gun is explicit correlation: 9~basis functions with
$r_{12}^2$ factors (Hylleraas $p = 2$) achieve 94.7\% $D_e$,
exceeding the 114-function prolate CI plateau.  The residual
7.5\% matches Paper~12's cusp diagnosis (7.6\%): the
electron-electron cusp ($1/r_{12}$ singularity) is a completeness
deficiency of the product basis, not a convergence issue.""",
    r"""The observation that looked like a smoking gun is explicit correlation:\
9~basis functions with $r_{12}^2$ factors (Hylleraas $p = 2$) achieve
$94.7\%$ $D_e$, exceeding the 114-function $\sigma$-only plateau.
\textbf{[MEASURED 2026-09-14]} It is a real measurement and it points
somewhere else than was read.  An $r_{12}$ factor supplies angular
correlation, which is precisely what the absent $m \neq 0$ configurations
were carrying;\ so does opening the azimuthal axis, and more cheaply
still---$54$ functions at $|m| \le 1$ reach $98.95\%$.  Two different
devices supplying the same missing ingredient is evidence about
\emph{what was missing}, not evidence that a basis without either is
structurally deficient.  The numerical coincidence between the residual
$7.5\%$ here and Paper~12's $7.6\%$ is genuine, and both are now
attributed to the same cause:\ the $\sigma$-only restriction.""",
    "P18: the r12 'smoking gun' re-read")

# ---- Paper 15 comparison
edit(
    r"""The mol-frame hyperspherical reparameterization (Paper~15)
resolves this by introducing coordinates $(\rho, R, \alpha)$
where $\rho$ is the electron hyperradius and $\alpha$ the
hyperangle.  The angular eigenvalues $\mu(\rho, R)$ are now a
transcendental function of \textit{two} continuous parameters,
reflecting the simultaneous three-body coalescence in coordinate
and configuration space.  This solver achieves 94.1\% $D_e$,
surpassing the prolate CI ceiling.""",
    r"""The mol-frame hyperspherical reparameterization (Paper~15) introduces
coordinates $(\rho, R, \alpha)$ where $\rho$ is the electron hyperradius
and $\alpha$ the hyperangle.  The angular eigenvalues $\mu(\rho, R)$ are
a transcendental function of \emph{two} continuous parameters.  That
solver achieves $94.1\%$ $D_e$ at $l_{\max} = 4$ with $\sigma$ and $\pi$
channels.  \textbf{[SCOPE]} Earlier versions read this as ``surpassing the
prolate CI ceiling.''  It is not a matched comparison---the $94.1\%$
includes $\pi$ channels and the $92.5\%$ does not---and Paper~15 has
withdrawn the coordinate-system reading it supported.""",
    "P18: Paper 15 comparison scoped")

# ---- the Level 3/4 'qualitative' distinction
edit(
    r"""The distinction between Levels~3 and 4 is qualitative.  At
Level~3, the $S^3 \times S^3$ FCI converges (albeit slowly)---the
exchange constant $\mu(R)$ is a computational convenience that
accelerates convergence but is not fundamentally necessary.  At
Level~4, the prolate product-space CI \textit{saturates}: the
exchange constant $\mu(\rho, R)$ is irreducible because the
product basis is structurally blind to three-body coalescence.
The reparameterization does not merely accelerate convergence;
it accesses physics that the algebraic product space cannot
represent at any practical basis size.""",
    r"""\textbf{[WITHDRAWN 2026-09-14]} Earlier versions closed here with a
\emph{qualitative} distinction between Levels~3 and 4:\ that at Level~3
$\mu(R)$ is a computational convenience while at Level~4 $\mu(\rho,R)$ is
\emph{irreducible}, because the product basis ``saturates'' and is
``structurally blind to three-body coalescence.''  The saturation was the
$\sigma$-only restriction, so that distinction has lost its evidence and
is withdrawn.

\textbf{[OBSERVATION]} What stands in its place is weaker and, on present
evidence, the same at both levels:\ $\mu(\rho,R)$ is a
\emph{reparameterization} constant, as $\mu(R)$ is at Level~3.  Both are
genuinely transcendental---computed by pointwise diagonalisation, with no
known algebraic replacement (Paper~12; the project's algebraic
registry)---and both accelerate convergence in their own coordinates.
Neither is shown to be necessary for the energy.  That places $\mu(\rho,R)$
in this taxonomy by its \emph{transcendental content}, which is unaffected,
rather than by an irreducibility argument, which is not available.

\textbf{[OPEN]} Whether any Level-4 exchange constant is irreducible in the
strong sense is now an open question rather than a settled one.  Answering
it would require a product-space calculation that converges, or provably
fails to converge, with the azimuthal axis open---which the $|m| \le 1$
result begins but does not finish, since $|m| = 2$ was not obtained there.""",
    "P18: Level 3/4 qualitative distinction withdrawn; mu re-priced")

# ---- L1604 'key irreducibility result'
edit(
    r"""single global characteristic polynomial. The product-space CI
    saturation at 92.5\% $D_e$ remains the key irreducibility result.""",
    r"""single global characteristic polynomial. (Earlier versions added
    that ``the product-space CI saturation at 92.5\% $D_e$ remains the
    key irreducibility result'';\ that saturation is now understood as a
    $\sigma$-sector artifact and the irreducibility reading is withdrawn
    --- Sec.~\ref{sec:mu_level4}.)""",
    "P18 L1604: 'key irreducibility result' withdrawn")

# ---- L3021 Class-C statement
edit(
    r"""  \item \textbf{Two-electron energies} computed via FCI on
    $S^3 \times S^3$ are formally Class~S (graph eigenvalues
    of the product-space Hamiltonian), but achieve limited
    accuracy due to the product-basis completeness deficiency
    (Sec.~\ref{sec:mu_level4}).  High-accuracy two-electron
    energies require the Class~C exchange constants.
    The classification thus tracks both the type of observable
    and the accuracy target.""",
    r"""  \item \textbf{Two-electron energies} computed via FCI on
    $S^3 \times S^3$ are formally Class~S (graph eigenvalues
    of the product-space Hamiltonian).  \textbf{[SCOPE 2026-09-14]}
    Earlier versions added that they ``achieve limited accuracy due to
    the product-basis completeness deficiency'' and that ``high-accuracy
    two-electron energies require the Class~C exchange constants.''  The
    completeness deficiency was a $\sigma$-sector restriction
    (Sec.~\ref{sec:mu_level4}), so the \emph{require} does not follow:\ a
    product-space basis with the azimuthal axis open reaches ${\sim}99\%$
    of $D_e$ without any Class~C constant.  What the classification tracks
    is the type of observable;\ whether it also tracks an accuracy target
    is no longer established.""",
    "P18 L3021: Class-C necessity claim scoped")

# ---- the track_m bibitem has no permanent home (Sec. 9 policy)
edit(
    r"""\bibitem{loutey_track_m}
J.~Loutey,
``Track M: Prolate CI $l_{\max}$ Convergence Study,''
GeoVac supplementary material (2026).

""",
    "",
    "P18: track_m bibitem removed (no permanent home; data is in Paper 12)")

with io.open(P18, encoding="utf-8") as fh:
    t = fh.read()

applied, failed = [], []
for old, new, label in EDITS:
    if old in t:
        t = t.replace(old, new, 1)
        applied.append(label)
    else:
        failed.append(label)

with io.open(P18, "w", encoding="utf-8") as fh:
    fh.write(t)

print("applied %d edits" % len(applied))
for a in applied:
    print("  +", a)
if failed:
    print("UNMATCHED (%d):" % len(failed))
    for f in failed:
        print("  -", f)
    sys.exit(1)
