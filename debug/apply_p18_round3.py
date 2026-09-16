"""Paper 18 round 3 -- withdraw, do not replace.

Three rounds running, my edits to this subsection have introduced the next
round's defects, and the shape is the same every time: I withdraw a claim and
then write a REPLACEMENT READING that is itself unbacked.  Round 1 wrote
"structurally blind" out and "genuinely transcendental" in (false, contradicted
by this paper's own Sec. sec:mu).  Round 2 wrote that out and wrote in "a
departure priced in transcendental content" (contradicted by Paper 12, which
calls the mu>0 quadrature an implementation state: "well-defined but not
attempted here") plus "carrying no Class-C constant" (unsupported while the
price question is open) plus "the reading the abstract already states" (the
abstract states the STRONGER observable-indexed reading).

This round removes the replacement narratives instead of writing another one.
Verified against primary text before editing:

LARGE-7  Track M's numbers (114 functions, 92.5%, "factor of 200", the
         l_max ~ 45 / 113 extrapolation) appear NOWHERE in Paper 12 -- I
         redirected the citation to a document that does not contain them when
         I deleted the track_m bibitem.  Repriced to Paper 12's actual ladder.
LARGE-2  abstract L41-44 indexes content to the OBSERVABLE ("of any GeoVac
         observable ... determined entirely by"); the route-relative reading
         denies that indexing.  The appeal is removed and the tension flagged.
LARGE-4  the Class-C DEFINITION still said "These require the higher exchange
         constants" seven lines above the paragraph denying it.
LARGE-3  the conclusion still certified Claim 4 "has held across all tested
         cases".
LARGE-5/6 the two horns.  Primary text resolves toward 6: no price is
         established, so the consolation goes and the counterexample's status
         becomes an open question rather than an assertion.
LARGE-1  a cusp-floor attribution at two loci in THIS paper that Paper 13 has
         [WITHDRAWN] -- my own sweep missed it while editing this file.
SMALL-1  the sec:algebraic_curve bullet says mu(rho,R) IS algebraic at every
         truncation, contradicting the same section's own later bullet.
SMALL-2  "mu(R) is needed to achieve sub-0.1%" is false twice.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

P18 = "papers/group3_foundations/paper_18_exchange_constants.tex"
EDITS = []


def edit(old, new, label):
    EDITS.append((old, new, label))


# ---- LARGE-7: reprice Track M's numbers to what Paper 12 actually contains
edit(
    r"""A systematic convergence study (Track~M) found that this approach
saturates at $\sim$92.5\% of the exact dissociation energy
$D_e$~\cite{loutey_paper12}:\ 92.3\% at $l_{\max} = 2$ (46 basis
functions), 92.5\% at $l_{\max} = 4$ (114), with successive improvements
decaying by a factor of 200 from $l = 1 \to 2$ to $l = 2 \to 3$ and the
$l = 3 \to 4$ increment contributing $0.04\%$.  Extrapolating that decay
predicted $l_{\max} \sim 45$ for $5\%$ error and $l_{\max} \sim 113$ for
$1\%$, both impractical.""",
    r"""Paper~12's convergence study finds that this approach saturates at
$92.4\%$ of the exact dissociation energy $D_e$~\cite{loutey_paper12}:\
$92.37\%$ at 46 basis functions and $92.42\%$ at 72, the last near-tripling
of the basis buying $0.34$~mHa.

\textbf{[WITHDRAWN 2026-09-14]} Earlier versions quoted a finer ladder here
(``92.5\% at $l_{\max}=4$, 114 basis functions'', a decay ``by a factor of
200'', and an extrapolation predicting $l_{\max}\sim45$ for $5\%$ error)
from an unpublished track note.  That note has no permanent home, none of
those figures is reproduced in Paper~12, and the extrapolation is in any
case an extrapolation along the wrong axis---see below.  They are withdrawn
rather than re-sourced.""",
    "LARGE-7: Track M's numbers repriced to Paper 12's actual ladder")

# ---- LARGE-5 + LARGE-6: drop the consolation reading entirely
edit(
    r"""\textbf{[SCOPE]} That $99.1\%$ is \emph{not} obtained with the same
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
    r"""\textbf{[SCOPE]} That $99.1\%$ is \emph{not} obtained with the same
\emph{algebraic} $V_{ee}$:\ Paper~12's own scope note records that the
$\mu > 0$ radial integrals leave the $A_l$, $B_l$, $X_l$ recurrences for a
two-dimensional spectral quadrature, so the quadrature-free property holds
in the $\sigma$ sector only.

\textbf{[OPEN]} Whether that departure carries \emph{transcendental} content
in this paper's sense is not settled here, and the question matters for the
taxonomy.  Paper~12 describes the quadrature as an implementation state
rather than a structural obstruction---generalising the three auxiliary
tables to associated Legendre functions is ``well-defined but not attempted
there''---which would make the $|m|\le1$ route carry no new exchange
constant.  The Level-2 analogue points the other way:\ opening $m \neq 0$
there is exactly how the embedding constant $e^aE_1(a)$ arises
(Sec.~\ref{sec:stieltjes}), and the tier table's ``(none)'' entry for the
prolate product space is a $\sigma$-sector entry with no counterpart row
for $|m| \ge 1$.  An earlier version of this paragraph asserted the first
reading as though it were established;\ it is not, and the taxonomy owes
this row a price or a proof that there is none.""",
    "LARGE-5/6: the consolation reading replaced by the open question it is")

# ---- LARGE-2 + LARGE-4 + LARGE-5: Claim 4 scope, without the false appeal
edit(
    r"""observables (Class~C) are where the embedding, flow, or composition
exchange constants of Sec.~\ref{sec:taxonomy} appear.}

\medskip

\noindent\textbf{[SCOPE 2026-09-14]} An earlier wording said Class-C
observables \emph{require} those constants.  That is falsified by this
paper's own Sec.~\ref{sec:mu_level4}:\ the H$_2$ dissociation energy is a
multi-particle correlated observable, and a prolate product space carrying
no Class-C constant reaches ${\sim}99\%$ of it.  What the classification
tracks is which constants \emph{appear} in a given evaluation route, not
which are unavoidable---the reading the abstract already states.""",
    r"""observables (Class~C) are where the embedding, flow, or composition
exchange constants of Sec.~\ref{sec:taxonomy} appear.}

\medskip

\noindent\textbf{[SCOPE 2026-09-14]} An earlier wording said Class-C
observables \emph{require} those constants, and the abstract still states
the claim in the stronger, observable-indexed form (``the transcendental
content of any GeoVac observable is determined entirely by the type of
projection'').  Sec.~\ref{sec:mu_level4} puts that under strain:\ the
H$_2$ dissociation energy is a multi-particle correlated observable, and a
prolate product space reaches ${\sim}99\%$ of it by a route whose exchange
constant this paper has not yet priced.  Whether that is a counterexample
to the strong form or a route that carries an unpriced constant is the open
question recorded there.  Pending it, the body states the weaker claim---the
classification tracks which constants \emph{appear} in a given evaluation
route---and the two forms are not reconciled.  \textbf{This is a live
tension in the paper's headline claim, not a settled scoping.}""",
    "LARGE-2/5: the false appeal to the abstract removed; the tension stated")

# ---- LARGE-4: the Class-C definition itself
edit(
    r"""    the inter-particle coordinate.  These require the higher exchange
    constants: $e^a E_1(a)$ at Level~2 for $m \neq 0$ states""",
    r"""    the inter-particle coordinate.  The higher exchange constants appear
    here: $e^a E_1(a)$ at Level~2 for $m \neq 0$ states""",
    "LARGE-4: the Class-C definition no longer defines the class by requirement")

# ---- LARGE-3: the conclusion certifies a claim the body records as strained
edit(
    r"""Claim~4 (transcendental content determined by projection type) has
held across all tested cases.""",
    r"""Claim~4 (transcendental content determined by projection type) has held
in its route-relative form across all tested cases;\ its stronger,
observable-indexed form is under strain from the Level-4 azimuthal result
and is not currently reconciled with the body
(Sec.~\ref{sec:observable_classification}).""",
    "LARGE-3: the conclusion no longer certifies Claim 4 unbroken")

# ---- LARGE-1: the cusp-floor zombie, two loci in this paper
edit(
    r"""direct evidence that the cusp is embedding content the graph cannot
absorb, explaining the graph-native CI convergence floor (Paper~13,
Track~DI).""",
    r"""direct evidence that the cusp is embedding content the graph cannot
absorb.  \textbf{[WITHDRAWN]} Earlier versions added that this explains the
graph-native CI convergence floor;\ Paper~13 has since withdrawn that
attribution, re-diagnosing the floor as a small-$Z$ graph-validity-boundary
artifact near $Z_c \approx 1.84$.""",
    "LARGE-1a: cusp-floor attribution withdrawn")

edit(
    r"""FCI basis invariance confirmed: the floor is embedding content
(cusp), not basis mismatch.""",
    r"""FCI basis invariance confirms only what the floor is \emph{not}
(basis mismatch);\ Paper~13 withdraws the cusp attribution.""",
    "LARGE-1b: the second cusp-floor locus")

# ---- SMALL-1: the self-contradicting bullet
edit(
    r"""The adiabatic eigenvalues $\mu(R)$ and $\mu(\rho,R)$ are algebraic
    functions""",
    r"""The adiabatic eigenvalues $\mu(R)$ are algebraic
    functions""",
    "SMALL-1a: mu(rho,R) removed from the algebraic bullet")

# ---- SMALL-2: 'needed to achieve sub-0.1%' is false twice
edit(
    r"""the $\mu(R)$ parameterization is needed to achieve sub-0.1\%
    accuracy.""",
    r"""the $\mu(R)$ parameterization was once thought necessary for sub-0.1\%
    accuracy.  \textbf{[WITHDRAWN]} It is not:\ the adiabatic route that
    carries $\mu(R)$ itself floors at $0.19$--$0.20\%$, and sub-$0.1\%$ is
    reached without it by the 2D variational solver ($0.022\%$ raw).""",
    "SMALL-2: the necessity claim withdrawn")

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
