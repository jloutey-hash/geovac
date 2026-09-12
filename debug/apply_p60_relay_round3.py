"""v5.10.14: the conflict resolved against the relay, and three provenance fixes.

The metric-free/unlinked conflict was posed back to the relayed source with our
projection derivation.  The derivation was accepted and the earlier description
withdrawn:  the -2.90250 calculation was an ORDINARY VARIATIONAL CI, H C = E S C,
and the metric-free cancellation belongs exclusively to the locked isoenergetic
method.  That is a genuine external test of eq:scale_lock -- an attempt to
describe a calculation as simultaneously unlinked and metric-free, posed against
it and retracted.

Three provenance corrections follow from the same round:
 * the variational bound for the FIXED-SCALE metric-free problem is not in the
   canon (whose variational discussion is framed on the scanned k-determinant),
   so it is ours to claim;
 * the identification SW matrix = V_0-weighted overlap is likewise not a
   quotable sentence in the texts -- our formalization, their integrals;
 * molecularly there is NO closed-form counterpart to beta_nu Z = p_kappa/R_nu;
   the betas are eigenvalues of the one-electron molecular Sturmian problem at
   fixed E.  That last one makes the paper's existing "S_SW-orthonormal by
   construction" a derived consequence rather than an asserted fact.
"""
import io

PAP = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
tex = io.open(PAP, encoding="utf-8").read()

# ------------------------------------------- 1. attribution + the retraction
A_OLD = r"""\textbf{[OBSERVATION]} A second consultation of the primary canon, relayed, both
identifies the posing and agrees with the mechanism:\ in that calculation the
scaling parameter was \emph{scanned as a free variational parameter and
minimized}, unlinked from the output eigenvalue --- with the ground state's
difficulty attributed there, as here, to $1s^{2}$ forcing both electrons to
$Q_\nu=p_\kappa/\sqrt2$ and blocking in-out correlation.  That is our
scale-optimized posing, it is what Eq.~\eqref{eq:no_selection} required, and it
is what our ladder independently matches.  We record the attribution as
secondary-source:\ it is consistent with three internal lines of evidence but is
not a quotation, and a primary-source check remains the one thing that could
overturn this section."""

A_NEW = r"""\textbf{[OBSERVATION]} Consultations of the primary canon, relayed,
settle the attribution --- and one of them supplies an external test of
Eq.~\eqref{eq:scale_lock}.  The posing was identified first:\ the scaling
parameter \emph{scanned as a free variational parameter and minimized}, unlinked
from the output eigenvalue, with the ground state's difficulty attributed there,
as here, to $1s^{2}$ forcing both electrons to $Q_\nu=p_\kappa/\sqrt2$.  The same
account also described that calculation as evaluating a \emph{metric-free} matrix
at each step --- which Eq.~\eqref{eq:scale_lock} forbids, since projecting a basis
built at $E_{\rm basis}$ against a target $E_{\rm out}$ leaves the $L^{2}$ overlap
with coefficient $(E_{\rm basis}-E_{\rm out})$, and unlinking the scale is
precisely what makes that coefficient nonzero.  Put back to the source together
with the projection, the derivation was accepted and the description corrected:\
the calculation was an \emph{ordinary variational} configuration interaction,
$HC=E\,SC$, and the metric-free cancellation belongs exclusively to the locked
isoenergetic method.  We record the attribution as secondary-source --- it is not
a quotation, and a primary-source check remains the one thing that could overturn
this section --- but the episode is worth stating for what it is:\ an external
description of a calculation as simultaneously unlinked and metric-free, posed
against Eq.~\eqref{eq:scale_lock} and withdrawn."""
assert A_OLD in tex, "attribution locus not found"
tex = tex.replace(A_OLD, A_NEW, 1)

# --------------------------------------------- 2. the bound is ours to claim
B_OLD = r"""identity makes the bound automatic rather than fortunate---the isoenergetic root
is the lowest root of the fixed-scale generalized problem $H(p_\kappa)C=E\,SC$,
so every point lies above the exact value by construction, not by observation."""
B_NEW = r"""identity makes the bound automatic rather than fortunate---the isoenergetic root
is the lowest root of the fixed-scale generalized problem $H(p_\kappa)C=E\,SC$,
so every point lies above the exact value by construction, not by observation.
\textbf{[OBSERVATION]} The canon does not appear to state this.  Its variational
discussion is framed on the Rayleigh--Ritz principle for the \emph{scanned}
$k$-dependent determinant and does not partition the argument between the locked
and scanned postings;\ a targeted consultation returned no statement that the
fixed-scale metric-free problem bounds from above.  We therefore claim it."""
assert B_OLD in tex, "bound locus not found"
tex = tex.replace(B_OLD, B_NEW, 1)

# ------------------------------- 3. SW identification is ours + molecular beta
S_OLD = r"""That identifies the Shibuya--Wulfman matrix rather than merely naming it:\ for
$V_0=-\sum_A Z_A/|\mathbf{r}-\mathbf{R}_A|$ its off-diagonal entries are exactly
the cross-center nuclear-attraction integrals, which is what
\texttt{geovac/shibuya\_wulfman.py} computes.  The molecular metric is not an
extra structure the method acquires under generalization;\ it \emph{is} $V_0$."""
S_NEW = r"""That identifies the Shibuya--Wulfman matrix rather than merely naming it:\ for
$V_0=-\sum_A Z_A/|\mathbf{r}-\mathbf{R}_A|$ its off-diagonal entries are exactly
the cross-center nuclear-attraction integrals, which is what
\texttt{geovac/shibuya\_wulfman.py} computes.  The molecular metric is not an
extra structure the method acquires under generalization;\ it \emph{is} $V_0$.
\emph{Provenance:} the Shibuya--Wulfman integrals and their momentum-space
evaluation are Avery's;\ the identification of the resulting matrix \emph{as} the
$V_0$-weighted overlap is our formalization --- a targeted consultation could not
locate that sentence in the texts, only the integrals it would describe.

\textbf{[ESTABLISHED, from Avery]} Molecularly there is no closed-form
counterpart to $\beta_\nu Z=p_\kappa/R_\nu$.  Because $V_0$ is multi-center, the
weighting factors are obtained as the exact eigenvalues of the one-electron
molecular Sturmian problem $[-\tfrac12\nabla^{2}-E+\beta v_0]\phi=0$ at fixed
$E$, with a many-electron configuration's $\beta_\nu$ following from the
isoenergetic condition on that spectrum.  This is what makes the molecular
orbitals of Sec.~\ref{sec:molecular} $S_{\rm SW}$-orthonormal \emph{by
construction} rather than by arrangement:\ potential-weighted orthonormality is
automatic for eigenvectors of a common Sturmian problem with distinct $\beta$,
and $S_{\rm SW}$ is the potential doing the weighting.  Two facts the paper had
been carrying separately are one fact."""
assert S_OLD in tex, "SW locus not found"
tex = tex.replace(S_OLD, S_NEW, 1)

io.open(PAP, "w", encoding="utf-8").write(tex)
print("paper: retraction recorded; bound claimed; SW provenance + molecular beta added")
