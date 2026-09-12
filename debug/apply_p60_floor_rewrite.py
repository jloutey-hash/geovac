"""Paper 60: retire the angular attribution of the accuracy floor.

Three loci carried the claim that the floor is angular truncation and that the
family "cannot reach chemical accuracy at any K":  the abstract, the Sec.4 floor
paragraph, and the conclusion.  The span-vs-posing diagnostic (2026-09-08)
refutes both halves -- the floor is the scale lock, and it is a GROUND-STATE
pathology that 2^1S does not share at identical encoding cost.

Numbers cited here are measured, not fitted, wherever a fitted floor would be
window-sensitive:  the excited-state ladder value at the largest computed basis
(K=202) is used rather than its extrapolated floor.
"""
import io

REG = "debug/qa/numeric_registry.py"
PAP = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"

# ------------------------------------------------------------------ registry
NEW = '''    "p60_exc_gap_k202": dict(
        value=1.786, convention="constant: mHa above the exact He 2^1S energy "
                                "(-2.145974046 Ha) reached by the METRIC-FREE "
                                "isoenergetic posing, full s+p+d+f, K=202. "
                                "MEASURED ladder value, not an extrapolated floor",
        q=None,
        provenance="MEASURED 2026-09-08, debug/p60_excited_ladder.py. Companion "
                   "ground-state value at the SAME K and the same ||M||_1=192.9 is "
                   "7.158 mHa. Ladder: 1.972 (K=74), 1.899 (100), 1.849 (130), "
                   "1.813 (164), 1.786 (202). A free-floor fit over that window "
                   "returns 1.65 mHa but is window-sensitive, so the measured "
                   "endpoint is what is cited in the paper.",
        aliases={1.972: "K=74", 1.813: "K=164"}),
'''
src = io.open(REG, encoding="utf-8").read()
anchor = '    "p60_cond_S_converged": dict('
assert anchor in src and "p60_exc_gap_k202" not in src
src = src.replace(anchor, NEW + anchor, 1)
# annotate the floor entry with its newly-identified mechanism
old_prov = '"MEASURED 2026-09-07 on the converged ladder K = 74..514, fit "'
assert old_prov in src
src = src.replace(
    old_prov,
    '"MECHANISM (2026-09-08): the floor is the SCALE LOCK, not angular "\n'
    '                   "truncation -- freeing lambda over the identical span "\n'
    '                   "reaches 1.28 mHa at K=130. And it is GROUND-STATE "\n'
    '                   "specific: 2^1S sits at 1.786 mHa at K=202. "\n'
    '                   "MEASURED 2026-09-07 on the converged ladder K = 74..514, fit "', 1)
io.open(REG, "w", encoding="utf-8").write(src)
print("registry: p60_exc_gap_k202 added; p60_energy_floor provenance annotated")

tex = io.open(PAP, encoding="utf-8").read()

# ------------------------------------------------------------------ abstract
A_OLD = r"""The cheap rule is also the one that stops converging---it saturates
$4.0\times$ above chemical accuracy at any basis size---so the two facts are
one fact.  The result is therefore a measured trade-off between cost growth
and attainable accuracy, not an unconditional efficiency claim.  The
development follows."""

A_NEW = r"""The cheap rule also carries an accuracy floor, whose mechanism we
locate exactly:\ the metric-free form holds \emph{iff} the basis scale is locked
to the eigenvalue, $\lambda=p_\kappa=\sqrt{-2E}$, and that lock is not the
variational optimum.  Solving the \emph{same} span variationally with the scale
freed reaches chemical accuracy, at the price of the entire encoding advantage
($1$-norm growth from $K^{0.72}$ to $K^{2.75}$).  The floor is moreover a
\emph{ground-state} pathology rather than a property of the method:\ because
$\|M\|_1$ does not depend on which root is extracted, every $^{1}S$ state costs
the same to encode, and at $K=202$ the ground state sits $4.49\times$ above
chemical accuracy while $2\,^{1}S$ sits $\gvq{p60_exc_gap_k202}{1.12}\times$
above it, the posing cost falling threefold to fourfold per rung up the $^{1}S$
ladder.  The result is therefore a measured trade-off between cost growth,
attainable accuracy and \emph{which state is targeted}, not an unconditional
efficiency claim.  The development follows."""

assert A_OLD in tex, "abstract locus not found"
tex = tex.replace(A_OLD, A_NEW, 1)

# ------------------------------------------------------------------ Sec. 4
S_OLD = r"""of the deficit.  \textbf{The floor is $4.0\times$ chemical accuracy
($1.594$~mHa), so this family cannot reach chemical accuracy at any $K$.}  The
angular detail that would close the gap is exactly what $l_{\max}$ is capped
against --- which is to say the cheap growth and the saturation are the same
fact seen twice.  \textbf{[OPEN]} Whether the floor can be lowered while the
sublinear growth is preserved --- for instance by a \emph{fixed}, non-growing
increment of high-$l$ functions --- is open, and is the natural next question
this section raises.  As it stands the fixed-$l_{\max}$ family is best read as a
\emph{fixed-accuracy tier}:\ sublinear cost at a known and unimprovable
$\sim\!6.4$~mHa, rather than an approximation that converges."""

S_NEW = r"""of the deficit.  The floor is $4.0\times$ chemical accuracy
($1.594$~mHa).

\textbf{[INTERNAL THEOREM]} \emph{Where the floor comes from.}  It is neither
angular truncation nor a deficiency of the span.  Write the configurations at a
free global scale $\lambda$, which the Goscinskian construction fixes to
$p_\kappa$.  Because each $\Phi_\nu$ is hydrogenic at $Q_\nu=\lambda/R_\nu$, the
one-body Coulomb metric is exactly diagonal,
\begin{equation}
  W_{\mu\nu}\equiv\Big\langle\Phi_\mu\Big|\sum_j r_j^{-1}\Big|\Phi_\nu\Big\rangle
    = R_\nu\,\delta_{\mu\nu},
  \label{eq:W_diagonal}
\end{equation}
verified entrywise to $5\times10^{-11}$.  Hence the kinetic matrix at unit scale
is $T=\mathbb{1}-S/2$, and the variational problem over the same span is
$H(\lambda)C=E\,SC$ with
$H(\lambda)=\lambda^{2}(\mathbb{1}-S/2)+\lambda(-Z\,\mathrm{diag}\,R_\nu+G)$ and
$G_{\mu\nu}=\langle\Phi_\mu|r_{12}^{-1}|\Phi_\nu\rangle$.  Putting
$E=-\lambda^{2}/2$ cancels every $S$ term \emph{identically} and leaves
$(Z\,\mathrm{diag}\,R_\nu-G)C=\lambda C$, which is Eq.~\eqref{eq:secular}.  The
metric-free posing is therefore the variational problem of the same span,
evaluated at the one scale where the $L^{2}$ metric drops out:
\begin{equation}
  \text{metric-free}\iff E=-\lambda^{2}/2\iff\lambda=p_\kappa=\sqrt{-2E}.
  \label{eq:scale_lock}
\end{equation}
\textbf{[MEASURED]} That lock is not the variational optimum.  Freeing $\lambda$
over the identical span reaches $0.15$~mHa of the exact $s$-limit at $K=136$
($s$-only, where $-2.879029$~Ha is known independently) against $4.43$~mHa
locked, and $\gvq{p60_span_deficit_spdf}{1.28}$~mHa of the exact energy at
$K=130$ ($spdf$) against $7.46$~mHa locked.  The span is not the limitation; the
lock is.  Two corollaries.  The variational bound is automatic rather than
fortunate --- $E_{\rm iso}$ is the lowest root of $H(p_\kappa)C=E\,SC$, so no
point can fall below the exact value.  And freeing the scale forfeits the
encoding advantage entirely:\ $\|M\|_1\sim K^{0.72}$ becomes
$\|H(\lambda^{\ast})\|_1\sim K^{1.95}$ and
$\|S^{-1/2}HS^{-1/2}\|_1\sim K^{2.75}$, with $\mathrm{cond}(S)\sim K^{0.94}$ --- a
factor $5\times10^{3}$ at $K=136$.  \emph{Metric-free posing, sublinear $1$-norm
and accuracy floor are one fact, not three.}

\textbf{[ESTABLISHED, from Avery]} \emph{Why the ground state in particular.}
Helium's ground state needs \emph{in-out} radial correlation, supplied
variationally by split-shell $1s\,1s'$ functions carrying two independent
exponents.  A Goscinskian $1s^{2}$ configuration places both electrons at $n=1$
and hence at one exponent $Q_\nu/n_j$ --- for \emph{any} weighting potential,
since a separable $V_0$ still assigns a single $\beta_\nu$ per configuration.
Excited configurations obtain two scales free from $n_a\neq n_b$.
\textbf{[MEASURED]} The state-dependence this predicts is large, and it is free:
$\|M\|_1$ is a property of $M$, not of which root is extracted, so every $^{1}S$
state costs the same to block-encode and only the delivered accuracy differs.
At $K=202$ ($\|M\|_1=192.9$) the ground-state error is $7.16$~mHa ($4.49\times$
chemical accuracy) against $\gvq{p60_exc_gap_k202}{1.786}$~mHa ($1.12\times$)
for $2\,^{1}S$.  Reference-free --- the posing cost
$E_{\rm iso}-\min_\lambda E_{\rm var}$ requires no known limit --- the $^{1}S$
ladder at $K=105$ reads $\gvq{p60_posing_cost_ground}{4.21}$,
$\gvq{p60_posing_cost_exc}{0.98}$, $0.32$ and $0.13$~mHa across the first four
roots, a threefold-to-fourfold reduction per rung.
\textbf{[OBSERVATION]} The shrinking Rydberg gaps do not spoil this at chemical
accuracy:\ resolving the $k$-th root needs
$\epsilon_p=\epsilon_E/p_\kappa\approx8\times10^{-4}$ against gaps of $0.336$,
$0.041$, $0.014$ and $0.006$ --- margins of $254$, $27$, $8.6$ and $3.9$ --- so
the query count $\|M\|_1/\epsilon_p\approx2\times10^{5}$ is state-independent
through $4\,^{1}S$.  The unmeasured cost of an interior root is the overlap of a
cheap trial state, not spectral resolution.

\textbf{[OPEN]} Two questions replace the one this section used to raise.
Whether a weighting potential of different \emph{radial shape} can reposition
the exponent distribution enough to lower the ground-state floor while
preserving Eq.~\eqref{eq:W_diagonal} --- a constant rescaling cannot, since
$\beta_\nu Z_w = p_\kappa/R_\nu$ eliminates $Z_w$ identically, so the lever must
change the shape of $V_0$ and not merely its strength.  And whether the
excited-state operating point survives the state-preparation cost that
Eq.~\eqref{eq:scale_lock} says nothing about.  As it stands the family is best
read as a \emph{fixed-accuracy tier whose accuracy depends on the target
state}:\ sublinear cost at ${\sim}6.4$~mHa for the ground state and
${\sim}1.8$~mHa for $2\,^{1}S$, bought at the same price."""

assert S_OLD in tex, "Sec.4 locus not found"
tex = tex.replace(S_OLD, S_NEW, 1)

# ------------------------------------------------------------------ conclusion
C_OLD = r"""growth, and is paid for by a $6.4$~mHa accuracy floor that no basis size in
that family crosses.  Read as a fixed-accuracy tier it is a real and cheap
operating point;\ read as a route to chemical accuracy it is not one, and
whether the floor can be lowered without losing the growth is the open
question we leave."""

C_NEW = r"""growth, and is paid for by an accuracy floor whose mechanism is
Eq.~\eqref{eq:scale_lock}:\ the metric-free form exists only at the scale
$\lambda=p_\kappa$, and that scale is not the variational optimum.  The floor is
a \emph{ground-state} pathology --- $6.4$~mHa there, but ${\sim}1.8$~mHa for
$2\,^{1}S$ at identical encoding cost, falling threefold to fourfold per rung up
the $^{1}S$ ladder.  Read as a ground-state route to chemical accuracy it is not
one;\ read as a cheap fixed-accuracy tier, and especially as an
\emph{excited-state} method --- which is where the generalized-Sturmian
literature has always been strongest --- it is a real operating point.  The open
question we leave is whether a weighting potential of different radial shape can
lower the ground-state floor without breaking Eq.~\eqref{eq:W_diagonal}."""

assert C_OLD in tex, "conclusion locus not found"
tex = tex.replace(C_OLD, C_NEW, 1)

io.open(PAP, "w", encoding="utf-8").write(tex)
print("paper: abstract + Sec.4 floor paragraph + conclusion rewritten")
