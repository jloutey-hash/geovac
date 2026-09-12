"""v5.10.13: the general-V0 form, the attribution, and the named obstruction.

Three edits, all following a second relayed consultation of Avery's canon whose
load-bearing leg was independently re-derived and numerically verified here:

 1. The -2.90250 withdrawal gains an attribution.  The scaling parameter was
    scanned as a free variational parameter and minimised, not locked to the
    eigenvalue -- which is what eq:no_selection already required, and what our
    own scale-optimised ladder already matched.
 2. eq:general_v0 added:  V C = V_0 B C, every L2 overlap cancelling for ANY
    local V_0, orthonormal configurations or not.  Its atomic specialisation
    reproduces eq:secular bit-identically.  Molecularly it identifies what the
    Shibuya-Wulfman matrix IS -- V_0's cross-center elements, which is literally
    what geovac/shibuya_wulfman.py computes.
 3. The [OPEN] V_0-shape question gains its named obstruction:  losing 1/r
    breaks the Fock/hyperspherical mapping and with it the closed-form
    multi-center and inter-electron integrals.  A wall with a mechanism, not
    unexplored ground.
"""
import io

PAP = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
tex = io.open(PAP, encoding="utf-8").read()

# ------------------------------------------------------------- 1. attribution
W_OLD = r"""matches instead our scale-optimized ladder ($1.64$~mHa at $K=100$).  Which posing
produced it --- or which basis, the genuine shared-scale Coulomb--Sturmians being
the other candidate --- could not be confirmed against the primary source, and we
flag it as the single result that would overturn this section were it to prove
otherwise."""

W_NEW = r"""matches instead our scale-optimized ladder ($1.64$~mHa at $K=100$).
\textbf{[OBSERVATION]} A second consultation of the primary canon, relayed, both
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
assert W_OLD in tex, "withdrawal locus not found"
tex = tex.replace(W_OLD, W_NEW, 1)

# ---------------------------------------------------------- 2. eq:general_v0
M_OLD = r"""with $S$ the Shibuya--Wulfman matrix~\cite{averyphd}, whose condition number"""
M_NEW = r"""with $S$ the Shibuya--Wulfman matrix~\cite{averyphd}, whose condition number"""
assert M_OLD in tex

ANCHOR = r"""The clean atomic result does not fully transfer to molecules.  The multi-center
one-electron isoenergetic problem is again \emph{generalized},"""
INSERT = r"""\textbf{[INTERNAL THEOREM]} \emph{The general-$V_0$ form, and what the metric
is.}  Before specializing, the general statement costs two lines and is worth
having, because it says what the molecular metric actually is.  Projecting
$\langle\Phi_\mu|H-E|\Phi_\nu\rangle=0$ and substituting the Sturmian equation on
the ket, $(-\tfrac12\sum_j\nabla_j^{2}-E)|\Phi_\nu\rangle=-\beta_\nu
V_0|\Phi_\nu\rangle$, gives
\begin{equation}
  \sum_\nu\Big(\langle\Phi_\mu|V|\Phi_\nu\rangle-\beta_\nu
    \langle\Phi_\mu|V_0|\Phi_\nu\rangle\Big)C_\nu=0,
  \qquad\text{that is}\qquad
  \mathbf{V}\,\mathbf{C}=\mathbf{V}_0\,\mathbf{B}\,\mathbf{C},
  \label{eq:general_v0}
\end{equation}
$\mathbf{B}=\mathrm{diag}(\beta_\nu)$.  Every $L^{2}$ overlap cancels
\emph{identically} --- for any local $V_0$, and whether or not the configurations
are orthonormal.  \textbf{[MEASURED]} Atomically
$\langle\Phi_\mu|V_0|\Phi_\nu\rangle=-Z R_\nu\delta_{\mu\nu}$ and
$\beta_\nu Z=p_\kappa/R_\nu$, so $\mathbf{V}_0\mathbf{B}=-p_\kappa\mathbb{1}$ and
Eq.~\eqref{eq:general_v0} collapses to Eq.~\eqref{eq:secular}:\ verified to
$6\times10^{-10}$, with the reconstruction $Z\,\mathrm{diag}\,R_\nu-G$
bit-identical to the assembled $M$.  \textbf{[OBSERVATION]} Molecularly
$\mathbf{V}_0$ is \emph{not} diagonal, so the matrix on the right-hand side is
the weighting potential itself --- not the identity, and not the $L^{2}$ overlap.
That identifies the Shibuya--Wulfman matrix rather than merely naming it:\ for
$V_0=-\sum_A Z_A/|\mathbf{r}-\mathbf{R}_A|$ its off-diagonal entries are exactly
the cross-center nuclear-attraction integrals, which is what
\texttt{geovac/shibuya\_wulfman.py} computes.  The molecular metric is not an
extra structure the method acquires under generalization;\ it \emph{is} $V_0$.

The clean atomic result does not fully transfer to molecules.  The multi-center
one-electron isoenergetic problem is again \emph{generalized},"""
assert ANCHOR in tex, "molecular section anchor not found"
tex = tex.replace(ANCHOR, INSERT, 1)

# -------------------------------------------------------- 3. the [OPEN] wall
O_OLD = r"""$\beta_\nu Z_w = p_\kappa/R_\nu$ eliminates $Z_w$ identically, so the lever must
change the shape of $V_0$ and not merely its strength."""
O_NEW = r"""$\beta_\nu Z_w = p_\kappa/R_\nu$ eliminates $Z_w$ identically, so the lever must
change the shape of $V_0$ and not merely its strength.  That lever has a named
price, and it is not small:\ a screened or mean-field $V_0$ improves the physical
quality of each configuration, but losing the $1/r$ form breaks the mapping to
hyperspherical harmonics under the Fock projection, and with it the closed-form
evaluation of the multi-center and inter-electron integrals --- the property this
entire encoding rests on.  The open question is therefore sharper than ``does a
better $V_0$ exist'':\ it is whether any $V_0$ of different radial shape
\emph{preserves the Fock mapping}.  If none does, the trade is a wall with a
mechanism rather than unexplored ground."""
assert O_OLD in tex, "[OPEN] locus not found"
tex = tex.replace(O_OLD, O_NEW, 1)

io.open(PAP, "w", encoding="utf-8").write(tex)
print("paper: attribution + eq:general_v0 + [OPEN] obstruction applied")
