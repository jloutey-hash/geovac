r"""Round-3 DELTA claims remediation: Papers 12/13/15 + synthesis + field guide
+ the claim-test matrix.

Every OLD string below was read back out of the file with sed before this script
was written.  Every finding was re-verified against primary text before being
accepted -- including the two that turned out to be MY OWN round-2 prose
(M-5 and M-6), which is where the largest defects keep landing.

PRINCIPLE, THIRD ROUND RUNNING: withdraw, do not replace.  Rounds 1 and 2 each
withdrew a claim and then wrote a replacement reading that became the next
round's defect.  Nothing below asserts a new mechanism.  Where two figures
conflict and I cannot resolve them from the corpus, the note says
[UNRECONCILED] and names both readings rather than picking one.

Findings NOT applied here, and why:
  L-4 and L-3 are genuine two-way numeric conflicts that need the solver to
  settle.  They are flagged in-paper as UNRECONCILED rather than guessed, and
  raised to the PI in the session summary.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

P12 = "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex"
P13 = "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex"
P15 = "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex"
SY = "papers/synthesis/group2_quantum_chemistry_synthesis.tex"
FG = "papers/synthesis/geovac_field_guide.tex"
MX = "docs/claim_test_matrix.md"

EDITS = []


def edit(path, old, new, label):
    EDITS.append((path, old, new, label))


# =====================================================================
# PAPER 12
# =====================================================================

# M-1: the two most-read surfaces quote the top of the envelope, bare.
edit(P12,
     r"""$99.1\%$ of $D_e$ at $|m| \le 1$, a gain of $11.6$~mHa that no
growth along the $\sigma$ axis reproduces (tripling the $\sigma$
basis is worth $0.34$~mHa).""",
     r"""$99.1\%$ of $D_e$ at $|m| \le 1$ (stability envelope
$99.0$--$99.1\%$;\ Sec.~\ref{sec:azimuthal}), a gain of $11.6$~mHa that
no growth along the $\sigma$ axis reproduces (near-tripling the
$\sigma$ basis is worth $0.34$~mHa).""",
     "M-1a/N-3a: abstract gains the envelope; 'tripling' -> 'near-tripling'")

edit(P12,
     r"""  \item \textbf{99.1\% of exact $D_e$} on restoring the azimuthal
    channels at $|m| \le 1$ in the same basis with the same
    Neumann kernel (Sec.~\ref{sec:azimuthal}), a gain of
    $11.6$~mHa, against $0.34$~mHa from tripling the $\sigma$ basis.""",
     r"""  \item \textbf{99.1\% of exact $D_e$} (stability envelope
    $99.0$--$99.1\%$) on restoring the azimuthal
    channels at $|m| \le 1$ in the same basis with the same
    Neumann kernel (Sec.~\ref{sec:azimuthal}), a gain of
    $11.6$~mHa, against $0.34$~mHa from near-tripling the $\sigma$
    basis.""",
     "M-1b/N-3b: conclusion gains the envelope; 'tripling' corrected")

# M-2: the stated residual budget is a third of the actual residual.
edit(P12,
     r"""ordinary basis incompleteness---higher $l$, higher $n$, and the
$\delta$ channels, worth a further $0.5$~mHa on an independent
Gaussian-basis estimate.""",
     r"""ordinary basis incompleteness---higher $l$, higher $n$, and the
$\delta$ channels.  That remainder is ${\approx}1.6$~mHa in total, of
which the $\delta$ channels account for ${\approx}0.5$~mHa on an
independent Gaussian-basis estimate;\ the balance is not itemised
here.""",
     "M-2: residual budget no longer reads as 0.5 mHa from exact")

# M-5: a false universal inside my own round-2 envelope paragraph.
edit(P12,
     r"""points fail---including $\alpha = 1.00$, which is the natural default and
the value every smaller basis in this paper uses;\ $\alpha = 1.15$ and""",
     r"""points fail---including $\alpha = 1.00$;\ $\alpha = 1.15$ and""",
     "M-5: false universal removed (alpha is optimised per basis size, L601)")

# M-6: 'insensitive to the choice', refuted by the number offered for it.
edit(P12,
     r"""section---that the azimuthal channels close essentially all of the
$7.6\%$ gap---is insensitive to the choice:\ every variational
point on the two one-dimensional slices quoted above exceeds $98.4\%$,
the weakest variational value anywhere on the full
$(\alpha, \mathrm{threshold})$ grid is $95.5\%$, and the independent
Gaussian route, which is well conditioned, gives $99.10\%$.""",
     r"""section---that the azimuthal channels, and not growth along the
$\sigma$ axis, are what close the gap---is robust across the grid.  The
\emph{size} of the closure is not.

\textbf{[SCOPE 2026-09-14]} An earlier version of this sentence called
the conclusion ``insensitive to the choice'' and then quoted the figure
that refutes it.  Every variational point on the two one-dimensional
slices above exceeds $98.4\%$, but the weakest variational value
anywhere on the full $(\alpha, \mathrm{threshold})$ grid is $95.5\%$,
which closes about $41\%$ of the $7.6$-point gap rather than
``essentially all'' of it.  The defensible statement is that the
qualitative effect is robust while the quantitative value ranges over
$95.5$--$99.1\%$ across the grid, the independent and well-conditioned
Gaussian route giving $99.10\%$.""",
     "M-6: the self-refuting 'insensitive' replaced by what was measured")

# M-7: the quadrature scoping missed two in-paper surfaces.
edit(P12,
     r"""Eq.~\eqref{eq:basis_function} with the same algebraic $V_{ee}$ and""",
     r"""Eq.~\eqref{eq:basis_function} with the same Neumann kernel and""",
     "M-7a: L804 'same algebraic V_ee' -> 'same Neumann kernel'")

edit(P12,
     r"""restored.  Same basis family, same algebraic $V_{ee}$, $\alpha$""",
     r"""restored.  Same basis family, same Neumann kernel, $\alpha$""",
     "M-7b: tab:azimuthal caption, same correction")

# =====================================================================
# PAPER 13
# =====================================================================

# N-8: a spectral bound stated as a measured accuracy, attributed to the graph.
edit(P13,
     r"""loutey_paper7}.  The discrete graph Laplacian reproduces
    hydrogenic energies to $< 0.1\%$~\cite{loutey_paper0,
    loutey_paper1}.""",
     r"""loutey_paper7}.  \textbf{[SCOPE 2026-09-14]} Earlier versions said
    here that the discrete graph Laplacian ``reproduces hydrogenic
    energies to $<0.1\%$''.  That states as a measured accuracy what is
    a property of the construction:\ $E_0 = \kappa \lambda_{\max}$
    holds by construction, so any such figure is a bound on the
    truncation deficit of the spectrum, not an accuracy against
    experiment~\cite{loutey_paper0, loutey_paper1}.""",
     "N-8: the 'reproduces to <0.1%' paraphrase withdrawn (bound, not accuracy)")

# L-2: the corrected cusp framing never reached Paper 13.  Three loci.
edit(P13,
     r"""    \emph{This paper.}  Hyperspherical coordinates $(R, \alpha,
    \theta_{12})$ place the electron-electron cusp at a boundary
    condition rather than a coordinate singularity.""",
     r"""    \emph{This paper.}  Hyperspherical coordinates $(R, \alpha,
    \theta_{12})$ place the electron-electron cusp on a coalescence
    manifold, where it becomes a boundary condition on the angular
    problem at fixed $R$.""",
     "L-2a: intro locus no longer contrasts against 'a coordinate singularity'")

edit(P13,
     r"""a manifold in the five-dimensional hyperangular space, not a point
singularity---the crucial structural advantage of these coordinates.""",
     r"""a manifold in the five-dimensional hyperangular space, not a point
singularity---a structural feature of these coordinates.

\textbf{[SCOPE 2026-09-14]} Earlier versions called this a ``crucial
structural advantage''.  It is a structural \emph{difference}:\ Paper~15
records that it does not translate into an accuracy advantage, and
Paper~12 reaches $99.1\%$ of H$_2$'s $D_e$ in coordinates where the same
locus is not a coordinate surface at all.""",
     "L-2b: 'crucial structural advantage' -> difference, with the scope note")

edit(P13,
     r"""coalescence manifold, rather than imposing conditions on the
hyperradial equation.  This is the key advantage over single-electron
coordinate systems: the cusp becomes a \emph{boundary condition on
the angular problem} at fixed $R$, not a singularity in the radial
coordinate.""",
     r"""coalescence manifold, rather than imposing conditions on the
hyperradial equation:\ the cusp becomes a \emph{boundary condition on
the angular problem} at fixed $R$, not a singularity in the radial
coordinate.  This is a difference from single-electron coordinate
systems, not a demonstrated advantage over them.""",
     "L-2c: 'key advantage over' -> difference from")

# M-3: six uncaveated 0.05% loci survive the abstract fix.
edit(P13,
     r"""yields $E = -2.9052$~Ha for helium (0.05\% error), substantially
improving on all previous GeoVac results for this system.""",
     r"""yields $E = -2.9052$~Ha for helium, nominally $0.05\%$ from the
exact value.  That figure is \emph{non-variational}---the energy lies
below the exact one---and is superseded as an accuracy by the properly
variational two-dimensional treatment at $0.022\%$.""",
     "M-3a: L130, the 'substantially improving on all previous' universal")

edit(P13,
     r"""is surprisingly the most accurate at 0.05\% error.  Higher""",
     r"""is surprisingly the closest to exact, at a nominal $0.05\%$---a
non-variational figure, so ``closest'' here is not ``best
converged''.  Higher""",
     "M-3b: L361, 'the most accurate' universal")

edit(P13,
     r"""energy ($N_\alpha = 200$, $N_R = 3000$).  The $l_{\max} = 0$
result is optimal because higher-$l$ channels amplify
finite-difference errors near the boundaries.}""",
     r"""energy ($N_\alpha = 200$, $N_R = 3000$).  The $l_{\max} = 0$
result is numerically closest to exact because higher-$l$ channels
amplify finite-difference errors near the boundaries;\ it is
non-variational, so it is not the best-converged entry.}""",
     "M-3c: tab:lmax caption, 'optimal'")

edit(P13,
     r"""Table~\ref{tab:comparison} compares this result with other GeoVac
methods for helium.  The hyperspherical solver is both more accurate
and faster than every previous approach, while using a dramatically
smaller matrix.""",
     r"""Table~\ref{tab:comparison} compares this result with other GeoVac
methods for helium.  The hyperspherical solver is faster than every
previous approach and uses a dramatically smaller matrix;\ its nominal
$0.05\%$ is non-variational and so is not comparable as an accuracy to
the variational entries beside it.""",
     "M-3d: L541, 'more accurate than every previous approach'")

edit(P13,
     r"""demonstrate that the adiabatic hyperspherical method achieves 0.05\%
accuracy for helium.""",
     r"""demonstrate that the adiabatic hyperspherical method reaches a
nominal, non-variational $0.05\%$ for helium.""",
     "M-3e: L951, 'achieves 0.05% accuracy'")

# M-4: the 0.05% is given two incompatible causes.
edit(P13,
     r"""This sparsity is ultimately responsible for the 0.05\% accuracy
of the single-channel ($l = 0$) adiabatic approximation: the
ground channel is only weakly coupled to excited channels.""",
     r"""This sparsity is why the single-channel ($l = 0$) adiabatic
approximation works as well as it does:\ the ground channel is only
weakly coupled to excited channels.

\textbf{[WITHDRAWN 2026-09-14]} Earlier versions made the sparsity
responsible for the $0.05\%$ figure specifically.  This paper's own
closing section attributes that number instead to fortuitous
cancellation between the adiabatic approximation error and the
finite-difference discretization error.  The two readings are not
compatible---one makes the number a structural property of the
framework, the other an artifact---and the sparsity attribution is
withdrawn.""",
     "M-4: the sparsity-causes-0.05% attribution withdrawn")

# M-8: an unverifiable provenance detail inside the round-2 [WITHDRAWN].
edit(P13,
     r"""near $Z_c \approx 1.84$---confirmed by a sign flip at $Z = 10$---rather
than as cusp content.""",
     r"""near $Z_c \approx 1.84$---rather than as cusp content.""",
     "M-8: the unsourced 'sign flip at Z = 10' clause dropped")

# =====================================================================
# PAPER 15
# =====================================================================

# L-1: 'completely decoupled' contradicts the preceding paragraph, the
# measurement, and Paper 12.
edit(P15,
     r"""\paragraph{$\sigma$--$\pi$ decoupling.}
Because the nuclear coupling is diagonal in $m$, and the electron-electron
coupling also conserves the total $M = m_1 + m_2$, the $\sigma$ ($m = 0$)
and $\pi$ ($|m| = 1$) sectors are completely decoupled in the angular
eigenvalue problem.  $\pi$ channels contribute a constant
${\sim}6.6$~percentage-point additive offset to $D_e$ regardless of
$\sigma$ $l_{\max}$ (Table~\ref{tab:extended_convergence}):""",
     r"""\paragraph{The $\pi$ offset.}
\textbf{[WITHDRAWN 2026-09-14]} Earlier versions of this paragraph
stated that the $\sigma$ ($m = 0$) and $\pi$ ($|m| = 1$) sectors are
``completely decoupled in the angular eigenvalue problem'', on the
grounds that the nuclear coupling is diagonal in $m$ and the $e$--$e$
coupling conserves $M = m_1 + m_2$.  That does not follow, and the
paragraph immediately above says as much:\ $(m_1, m_2) = (0,0)$ and
$(+1,-1)$ both carry $M = 0$, so conserving $M$ is precisely what places
them in the \emph{same} block, and the $e$--$e$ multipole expansion
couples them.  The measurement points the same way---decoupled blocks
would make the ground state a minimum over blocks, so adding $\pi$
channels could not lower it by a roughly constant amount, which is what
is observed.  What survives is the empirical regularity, stated as one:\
$\pi$ channels contribute a near-constant additive offset to $D_e$
across $\sigma$ $l_{\max}$ (Table~\ref{tab:extended_convergence}):""",
     "L-1a: the sigma-pi decoupling claim withdrawn; the regularity kept")

# N-5: 'constant' and the 7.2 pp outlier at the headline truncation.
edit(P15,
     r"""This constant offset means $\pi$ convergence is independent of $\sigma$
convergence; the $\pi$ contribution is set by the frozen $l_{\max} = 2$
$\pi$ channels throughout the extended study.""",
     r"""The offset is near-constant over $l_{\max} = 2$--$5$ and rises to
$7.2$~pp at $l_{\max} = 6$, the truncation that generates this paper's
headline figure, so ``constant'' is an approximation rather than an
exact regularity.  The $\pi$ contribution is in any case set by the
frozen $l_{\max} = 2$ $\pi$ channels throughout the extended study,
which is a property of the methodology and not a demonstrated
independence of $\pi$ convergence from $\sigma$ convergence.""",
     "N-5: 'constant' scoped; the independence inference withdrawn")

edit(P15,
     r"""$\pi$ channels contribute a constant ${\sim}6.6$~pp offset
independent of $\sigma$ $l_{\max}$, consistent with complete
$\sigma$--$\pi$ decoupling in the angular eigenvalue problem.""",
     r"""$\pi$ channels contribute a near-constant ${\sim}6.6$~pp offset
across $\sigma$ $l_{\max}$, rising to $7.2$~pp at $l_{\max} = 6$.
\textbf{[WITHDRAWN 2026-09-14]} Earlier versions read this as
consistent with complete $\sigma$--$\pi$ decoupling;\ the two sectors
are coupled through the $e$--$e$ multipole expansion, and the
regularity is reported here as an observation without that
explanation.""",
     "L-1b: the conclusion's decoupling restatement withdrawn")

# M-9: an ordering clause survived the round-2 removal.
edit(P15,
     r"""The $\sigma$-only result (87.0\%) falls below Paper~12, and adding
$\pi$~channels pushes the Level~4 result to 94.1\%.""",
     r"""The $\sigma$-only result (87.0\%) falls below Paper~12's
$\sigma$-only $92.4\%$---at a different truncation and in a different
solver class---and adding $\pi$~channels pushes the Level~4 result to
94.1\%.""",
     "M-9: the surviving ordering clause scoped to its truncation/solver")

# M-10: the 93.6-vs-94.1 explanation names one of two differences.
edit(P15,
     r"""  $\sigma{+}\pi$ uses $m_{\max} = 1$ with $\pi$ channels frozen at
  their $l_{\max} = 2$ values --- which is why this column's
  $l_{\max} = 4$ entry ($93.6\%$) sits below the $94.1\%$ of
  Table~\ref{tab:comparison}, where the $\pi$ channels are solved
  at full $l_{\max}$;\ the two are different calculations, not a
  discrepancy.""",
     r"""  $\sigma{+}\pi$ uses $m_{\max} = 1$ with $\pi$ channels frozen at
  their $l_{\max} = 2$ values.  This column's $l_{\max} = 4$ entry
  ($93.6\%$) and the $94.1\%$ of Table~\ref{tab:comparison} differ in
  \emph{two} respects, not one:\ the $\pi$ channels are frozen here and
  solved at full $l_{\max}$ there, and this table is the 2D+cusp solver
  while that one is the plain variational 2D solver.  The two effects
  partly cancel, so the difference between the entries measures
  neither alone;\ the two are different calculations, not a
  discrepancy.""",
     "M-10: the caption now names both differences")

# L-3: the delta row's channel count cannot fit its stated baseline.
edit(P15,
     r"""At $l_{\max} = 4$, adding $\delta$ channels ($m_{\max} = 2$,
8~additional channels beyond the 29 $\sigma{+}\pi$ channels)
yields a modest $+0.65$~percentage-point gain, from $87.0\%$ to
$87.6\%$ relative to $\sigma$-only at the same $l_{\max}$ (the figure
quoted in Table~\ref{tab:extended_convergence}'s footnote).""",
     r"""At $l_{\max} = 4$, adding $\delta$ channels ($m_{\max} = 2$) yields a
modest $+0.65$~percentage-point gain, from $87.0\%$ to $87.6\%$
\emph{relative to $\sigma$-only} at the same $l_{\max}$ (the figure
quoted in Table~\ref{tab:extended_convergence}'s footnote).

\textbf{[UNRECONCILED 2026-09-14]} The channel count attached to this
row does not fit that baseline.  Table~\ref{tab:extended_convergence}
lists $N_{\rm ch} = 37$, and an earlier version of this sentence read
that as ``8 additional channels beyond the 29 $\sigma{+}\pi$
channels''---but a $\sigma{+}\pi{+}\delta$ calculation at
$l_{\max} = 4$ cannot return $87.6\%$ when its own $\sigma{+}\pi$
subset returns $93.6\%$.  Either the row is $\sigma{+}\delta$, which
would be $13 + 8 = 21$ channels rather than 37, or the $87.6\%$ belongs
to a different calculation.  The two readings are not resolved here,
and the row should not be relied on until they are.""",
     "L-3a: the delta row flagged UNRECONCILED, both readings named")

edit(P15,
     r"""  \item \textbf{Higher $m$ channels.}  $\delta$ channels ($m_{\max} = 2$)
    add $+0.65$~pp at $l_{\max} = 4$ but at $8.7\times$ cost.
    CBS extrapolation including $\delta$ would push the limit above 97\%.""",
     r"""  \item \textbf{Higher $m$ channels.}  $\delta$ channels ($m_{\max} = 2$)
    add $+0.65$~pp at $l_{\max} = 4$ but at $8.7\times$ cost.
    \textbf{[WITHDRAWN 2026-09-14]} Earlier versions carried that gain
    into this budget and concluded that CBS extrapolation including
    $\delta$ ``would push the limit above 97\%''.  The $+0.65$~pp was
    measured against $\sigma$-only, not against the $\sigma{+}\pi$
    baseline this budget uses, so the transfer is invalid;\ and the row
    it comes from is itself unreconciled (see the $\delta$-channel cost
    paragraph above).  No $\delta$ contribution to the $\sigma{+}\pi$
    gap is claimed here.""",
     "L-3b: the invalid transfer of the delta gain into the sigma+pi budget")

# N-6: an unbacked necessity claim, now also false across geometries.
edit(P15,
     r"""Reaching $> 99\%$ would require higher $m$ channels ($m_{\max} \ge 2$)
at extended $l_{\max}$, where the cost scaling becomes impractical.""",
     r"""Reaching $> 99\%$ along this route would need higher $m$ channels
($m_{\max} \ge 2$) at extended $l_{\max}$, where the cost scaling
becomes impractical.  That is an extrapolation from the CBS estimate
rather than a proven necessity, and it is not a necessity in general:\
Paper~12 reaches $99.1\%$ at $|m| \le 1$ in a different geometry.""",
     "N-6: the 'would require' necessity scoped to this route")

# L-4: two values for the same quantity; the abstract follows from neither.
edit(P15,
     r"""At $l_{\max} = 4$, the correction differential (between $R_{\rm eq}$ and
dissociation) is 0.39~mHa.  The cusp correction is a property of the
partial-wave basis truncation, not of the hyperradial solution method.""",
     r"""At $l_{\max} = 4$, the correction differential (between $R_{\rm eq}$ and
dissociation) is 0.39~mHa.  The cusp correction is a property of the
partial-wave basis truncation, not of the hyperradial solution method.

\textbf{[UNRECONCILED 2026-09-14]} That measurement and the
${\sim}1.7$~mHa / ${\sim}1.0$~percentage-point figure quoted above do
not sit together under this paper's own $1/(l_{\max} + 1/2)^4$
scaling:\ propagating $0.39$~mHa from $l_{\max} = 4$ gives
${\approx}0.09$~mHa, or ${\approx}0.05$~pp, at $l_{\max} = 6$---not
$1.0$~pp.  The abstract's pure-variational ``${\sim}95\%$'' is obtained
by subtracting the $1.0$~pp figure from $96.0\%$;\ under the scaling
law it would instead be ${\approx}96.0\%$.  One of the two figures is
wrong.  Which one is not settled here, and the pure-variational value
should be read as uncertain over ${\sim}95$--$96\%$ until it is.""",
     "L-4a: the two cusp-correction magnitudes flagged UNRECONCILED")

edit(P15,
     r"""(61~channels, CBS extrapolation ${\sim}97\%$;\ the pure-variational value is
${\sim}95\%$).""",
     r"""(61~channels, CBS extrapolation ${\sim}97\%$;\ the pure-variational
value is ${\sim}95$--$96\%$, this paper's two estimates of the cusp
contribution being unreconciled).""",
     "L-4b: the abstract's pure-variational figure widened to the honest range")

# =====================================================================
# SYNTHESIS + FIELD GUIDE + CLAIM MATRIX
# =====================================================================

# M-12: the third scope limit is missing from the synthesis.
edit(SY,
     r"""$V_{ee}$ matrix is exact within the Neumann truncation order with no
six-dimensional quadrature and no fitted parameters~\cite{loutey_paper12}.""",
     r"""$V_{ee}$ matrix is exact within the Neumann truncation order with no
six-dimensional quadrature and no fitted parameters~\cite{loutey_paper12}.
That property is established in the $\sigma$ sector;\ the azimuthal
extension described below does not inherit it.""",
     "M-12a: the quadrature-free claim scoped to sigma at its first surface")

edit(SY,
     r"""$|m| = 2$ sector was not obtained there, so the $\delta$ contribution
(${\approx}0.5$~mHa) rests on an independent Gaussian-basis route;\ and
the quoted figure carries a stability envelope, the basis being strongly
linearly dependent at that size.""",
     r"""$|m| = 2$ sector was not obtained there, so the $\delta$ contribution
(${\approx}0.5$~mHa) rests on an independent Gaussian-basis route;\ the
quoted figure carries a stability envelope, the basis being strongly
linearly dependent at that size;\ and the quadrature-free property does
not extend to it---the $\mu > 0$ radial integrals are evaluated by
spectral quadrature rather than by the closed-form recurrences, so the
headline figure is not itself a quadrature-free result.""",
     "M-12b: the third scope limit added to the azimuthal paragraph")

edit(SY,
     r"""available from tripling the $\sigma$ basis~\cite{loutey_paper12}.""",
     r"""available from near-tripling the $\sigma$ basis~\cite{loutey_paper12}.""",
     "N-3c: synthesis 'tripling' corrected")

# N-10: an uncaveated non-variational extrapolation in a summary table.
edit(FG,
     r"""3 & He (1-center, 2e) & Hyperspherical & $0.004\%$ (cusp) & 13 \\""",
     r"""3 & He (1-center, 2e) & Hyperspherical & $0.004\%$ (cusp-corrected;\ raw variational $0.022\%$) & 13 \\""",
     "N-10: field-guide He row carries the raw variational figure")

# M-11: the matrix carries the exact statement the paper corrected.
edit(MX,
     r"""new 2026-09-14. The m=0 tests could never see this: all three discrepancies vanish at m=0 |""",
     r"""new 2026-09-14. The m=0 tests could never see the (-1)^m or the squaring of the factorial ratio, both of which are trivial at m=0. The (2l+1) does NOT vanish at m=0; eq:neumann_sigma carries it because that independent derivation supplies it |""",
     "M-11: the matrix's false 'all three vanish at m=0' corrected")

# ---------------------------------------------------------------------
by_path = {}
for path, old, new, label in EDITS:
    by_path.setdefault(path, []).append((old, new, label))

applied, failed = [], []
for path, items in by_path.items():
    with io.open(path, encoding="utf-8") as fh:
        t = fh.read()
    for old, new, label in items:
        if old in t:
            t = t.replace(old, new, 1)
            applied.append(label)
        else:
            failed.append("%s   [%s]" % (label, path))
    with io.open(path, "w", encoding="utf-8") as fh:
        fh.write(t)

print("applied %d of %d edits" % (len(applied), len(EDITS)))
for a in applied:
    print("  +", a)
if failed:
    print("")
    print("UNMATCHED (%d):" % len(failed))
    for f in failed:
        print("  -", f)
    sys.exit(1)
