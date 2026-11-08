"""Paper 18, round 2: fix the defects my own re-pricing introduced.

Verified by the PM against primary text before editing:

LARGE-2  I wrote "Both are genuinely transcendental ... with no known
         algebraic replacement (Paper 12; the project's algebraic registry)".
         Paper 18 L404-420 says the OPPOSITE -- mu(R) "is an *algebraic*
         function of R ... a root of a polynomial, not a transcendental
         function" -- and docs/algebraic_registry.md says "algebraic
         (implicit) ... Point-by-point diagonalization is computational
         convenience."  I asserted the negation of my own citation, and cited
         Paper 12, which does not own mu at all.

LARGE-3  "the same at both levels" flattens a real structural distinction:
         the table gives mu(R) "alg. (implicit)" and mu(rho,R) "piecewise",
         and L3193-3199 says the pencil structure is "specific to Level 3"
         while Level 4 eigenvalues "do not satisfy a global characteristic
         polynomial".  The distinction that fell was completeness/necessity;
         the STRUCTURAL one is independent and still live.

LARGE-1  "with the same algebraic V_ee" -- the exact defect I had just fixed
         in Paper 12, reintroduced here, in the paper whose whole function is
         to price departures from algebraicity.

LARGE-4  "require" survived in Claim 4 and the Class-C definition.
LARGE-5  the 94.7% datum lost its only citation when I deleted the track_m
         bibitem; Paper 12 reports 86.8% with 18 functions for the nearest
         comparison, so the number is not merely uncited but contested.
S2       "342 functions still yields 92.40%" came from a reviewer's own
         control run and has no home in the corpus.  Unsourced numbers do not
         go in papers.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

P18 = "papers/group3_foundations/paper_18_exchange_constants.tex"
EDITS = []


def edit(old, new, label):
    EDITS.append((old, new, label))


# ---- LARGE-1 + S2 + the unsourced sigma control
edit(
    r"""construction.  Opening that axis in the same basis, with the same
algebraic $V_{ee}$, reaches ${\sim}99.1\%$ of $D_e$
(Paper~12~\cite{loutey_paper12}, Sec.~``Restoring the Azimuthal
Channels'').  The $\sigma$ axis does saturate---growing it to $342$
functions still yields $92.40\%$---but it saturates because it is the
wrong axis, not because the product space is incomplete.""",
    r"""construction.  Opening that axis in the same basis, with the same
Neumann kernel, reaches ${\sim}99.1\%$ of $D_e$
(Paper~12~\cite{loutey_paper12}, Sec.~``Restoring the Azimuthal
Channels'').  The $\sigma$ axis does saturate---Paper~12's own
convergence table buys $0.34$~mHa for $2.7\times$ the functions---but it
saturates because it is the wrong axis, not because the product space is
incomplete.

\textbf{[SCOPE]} That $99.1\%$ is \emph{not} obtained with the same
\emph{algebraic} $V_{ee}$, and for this paper the distinction is the whole
point.  Paper~12's own scope note records that the $\mu > 0$ radial
integrals leave the $A_l$, $B_l$, $X_l$ recurrences for a two-dimensional
spectral quadrature, so the quadrature-free property holds in the $\sigma$
sector only.  The honest exchange-constant statement is therefore sharper
than the withdrawn one, and points the other way:\ the accuracy formerly
credited to an irreducible \emph{reparameterization} constant is bought,
in the product space, by \emph{leaving the quadrature-free sector}.  That
is a departure priced in transcendental content---exactly the currency
this paper deals in---rather than a structural incapacity.""",
    "LARGE-1 + S2: algebraic-V_ee claim corrected; unsourced 342 removed")

# ---- LARGE-5: the 94.7% datum
edit(
    r"""The observation that looked like a smoking gun is explicit correlation:\
9~basis functions with $r_{12}^2$ factors (Hylleraas $p = 2$) achieve
$94.7\%$ $D_e$, exceeding the 114-function $\sigma$-only plateau.""",
    r"""The observation that looked like a smoking gun is explicit correlation:\
adding $r_{12}$ factors to a small $\sigma$-only basis buys accuracy that
many more $\sigma$ functions do not.  (Earlier versions quoted
``$94.7\%$ with 9 functions'' from an unpublished track note;\ that note
has no permanent home and the figure is not reproduced in Paper~12, whose
nearest comparison reports $86.8\%$ with 18 $r_{12}$-bearing functions
against $92.2\%$ for the algebraic $\sigma$-only route.  The qualitative
point below does not depend on which run is quoted, so the specific
number is withdrawn rather than repaired.)""",
    "LARGE-5: unbacked 94.7% withdrawn rather than left uncited")

# ---- LARGE-2 + LARGE-3: the replacement [OBSERVATION]
edit(
    r"""\textbf{[OBSERVATION]} What stands in its place is weaker and, on present
evidence, the same at both levels:\ $\mu(\rho,R)$ is a
\emph{reparameterization} constant, as $\mu(R)$ is at Level~3.  Both are
genuinely transcendental---computed by pointwise diagonalisation, with no
known algebraic replacement (Paper~12; the project's algebraic
registry)---and both accelerate convergence in their own coordinates.
Neither is shown to be necessary for the energy.  That places $\mu(\rho,R)$
in this taxonomy by its \emph{transcendental content}, which is unaffected,
rather than by an irreducibility argument, which is not available.""",
    r"""\textbf{[OBSERVATION]} What stands in its place is a \emph{structural}
distinction, which is independent of the withdrawn one and survives it
intact.  At Level~3 the angular Hamiltonian is a linear matrix pencil, so
$\mu(R)$ satisfies a global characteristic polynomial $P(R,\mu) = 0$ and
is an \emph{algebraic} function over $\mathbb{Q}(\pi,\sqrt2)$---pointwise
diagonalisation there is computational convenience, not necessity
(Sec.~\ref{sec:mu_level3};\ Track~P1).  At Level~4 the split-region
Legendre expansion breaks the pencil, and $\mu(\rho,R)$ is
piecewise-algebraic with \emph{no} global characteristic polynomial
(Track~S).  That is what separates the two levels in this taxonomy, and it
is why the table below tiers them differently---``alg.\ (implicit)''
against ``piecewise.''

\textbf{[WITHDRAWN]} An earlier version of this paragraph said instead that
both constants are ``genuinely transcendental \ldots\ with no known
algebraic replacement,'' citing the project's algebraic registry.  That is
the \emph{negation} of what the registry and Sec.~\ref{sec:mu_level3} both
state for $\mu(R)$, and it is withdrawn.  What is true of both is narrower:\
neither has been shown \emph{necessary} for the energy, so neither earns
the label irreducible on completeness grounds.""",
    "LARGE-2 + LARGE-3: transcendence claim withdrawn; structural distinction restored")

# ---- the [OPEN] now names what survives
edit(
    r"""\textbf{[OPEN]} Whether any Level-4 exchange constant is irreducible in the
strong sense is now an open question rather than a settled one.  Answering
it would require a product-space calculation that converges, or provably
fails to converge, with the azimuthal axis open---which the $|m| \le 1$
result begins but does not finish, since $|m| = 2$ was not obtained there.""",
    r"""\textbf{[OPEN]} What fell is the \emph{completeness} sense of
irreducibility---that the product space cannot represent the physics.  Two
other irreducibility results for Level~4 are untouched and should not be
read as withdrawn with it:\ the geometric-elevation obstruction of
Sec.~\ref{sec:irreducibility} (three routes ruled out), and the structural
one just stated (no global $P(\rho,\mu)$, Track~S).  What is now open is
narrower:\ whether any Level-4 exchange constant is necessary for the
\emph{energy}.  Answering it would require a product-space calculation
that converges, or provably fails to converge, with the azimuthal axis
open---which the $|m| \le 1$ result begins but does not finish, since
$|m| = 2$ was not obtained there.""",
    "S1: the [OPEN] no longer reads as withdrawing all Level-4 irreducibility")

# ---- LARGE-4: 'require' in the Class-C definition and Claim 4
edit(
    r"""These **require** the higher exchange constants""".replace("**", r"\textbf{").replace(r"\textbf{require\textbf{", "require"),
    "PLACEHOLDER", "unused")
EDITS.pop()

edit(
    r"""These require the higher exchange constants""",
    r"""These are where the higher exchange constants appear""",
    "LARGE-4a: Class-C 'require' softened to 'appear'")

edit(
    r"""Multi-particle correlated observables (Class~C) require the
    embedding, flow, or composition exchange constants of
    Sec.~\ref{sec:taxonomy}.""",
    r"""Multi-particle correlated observables (Class~C) are where the
    embedding, flow, or composition exchange constants of
    Sec.~\ref{sec:taxonomy} appear.  \textbf{[SCOPE 2026-09-14]} An
    earlier wording said they \emph{require} those constants.  That is
    falsified by this paper's own Sec.~\ref{sec:mu_level4}:\ H$_2$'s $D_e$
    is a multi-particle correlated observable and is reached to
    ${\sim}99\%$ in a product space carrying no Class-C constant.  The
    classification tracks which constants \emph{appear}, not which are
    necessary.""",
    "LARGE-4b: Claim 4 necessity wording corrected")

# ---- S6: the table row
edit(
    r"""  4 (prolate CI) & prolate product & (none) & (none) & algebraic \\""",
    r"""  4 (prolate CI) & prolate product & (none) & (none) & algebraic ($\sigma$)\footnotemark[1] \\""",
    "S6: prolate-CI table row scoped to sigma")

# ---- S14: 'No transcendental constant appears anywhere'
edit(
    r"""elements are likewise algebraic.  No transcendental constant
appears anywhere in the computation.""",
    r"""elements are likewise algebraic.  In the $\sigma$ sector no
transcendental constant appears beyond the single one-dimensional
quadrature Paper~12 records for its log-singular $B_l$ moment;\ the
$|m| \ge 1$ extension leaves that sector, as the scope note below states.""",
    "S14: absolute algebraicity claim scoped")

# ---- S4: dangling 'this singularity'
edit(
    r"""The transcorrelated approach (Paper~14, Sec.~IV; Paper~15,
Sec.~VI.J) removes this singularity via a Jastrow similarity""",
    r"""The transcorrelated approach (Paper~14, Sec.~IV; Paper~15,
Sec.~VI.J) removes the $1/r_{12}$ singularity via a Jastrow similarity""",
    "S4: dangling antecedent repaired")

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
