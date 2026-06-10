"""Phase C: rewrite Paper 38 to the unconditional theorem.

Replaces the named-gaps section with the translation-seminorm
framework (kernel proposition, injectivity lemma, lifted-state lemma,
unconditional theorem, remarks), and updates abstract / intro theorem
/ L5 conditionality accordingly. All replacements uniqueness-asserted.
"""
import io

P = r'papers/group1_operator_algebras/paper_38_su2_propinquity_convergence.tex'
s = io.open(P, encoding='utf-8').read()
fails = []


def rep(old, new, tag):
    global s
    if s.count(old) == 1:
        s = s.replace(old, new)
    else:
        fails.append((tag, s.count(old)))


# ---------------------------------------------------------------- 1.
NEW_SECTION = r'''\subsection{The unconditional theorem:\ compression and lifted
states under the translation seminorm}\label{sec:named_gaps}

Both gaps that conditioned earlier formulations of the main theorem
are closed by one reframing:\ metrize the truncated state space by
the \emph{left-translation Lipschitz seminorm} rather than by the
truthful-Dirac commutator seminorm. On the continuum the two agree
exactly (Lemma~\ref{lem:continuum_lip});\ at finite cutoff the
translation seminorm is non-degenerate while the Dirac-commutator
seminorm is not (Remark~\ref{rem:dirac_degeneracy}). All proofs in
this subsection are self-contained;\ the kernel quantities
$\KS_{\nmax}$, $Z_{\nmax}$, $\gamma_{\nmax}$ are those of
\S~\ref{sec:L2} below (forward references).

\begin{definition}[translation seminorm]
\label{def:translation_seminorm}
Let $\lambda_g$ denote left translation on $\sthree = \SU(2)$, acting
on the spinor space $L^2(\sthree, \Sigma) \cong L^2(\SU(2)) \otimes
\C^2$ (left trivialization) by $U_g = \lambda_g \otimes 1$, and
$\rho_g(T) = U_g T U_g^*$. Since $\DCH$ is built from the
left-invariant frame, $[U_g, \DCH] = 0$;\ hence $U_g$ preserves every
truncation $\Hilb_{\nmax}$ and conjugation preserves the operator
system:\ $\rho_g(P_{\nmax} M_f P_{\nmax}) = P_{\nmax}
M_{\lambda_g f} P_{\nmax}$. For $T \in \Op_{\nmax}$ and
$f \in C(\sthree)$ define
\begin{equation}
L_{\nmax}(T) := \sup_{g \ne e}
\frac{\opnorm{\rho_g(T) - T}}{d(e, g)},
\qquad
L(f) := \sup_{g \ne e}
\frac{\norm{\lambda_g f - f}_{L^\infty}}{d(e, g)},
\end{equation}
with $d$ the round geodesic distance (bi-invariant).
\end{definition}

\begin{lemma}[continuum identification]\label{lem:continuum_lip}
For $f \in C^\infty(\sthree)$,
$L(f) = \mathrm{Lip}(f) = \norm{\nabla f}_{L^\infty}
= \opnorm{[\DCH, M_f]}$.
\end{lemma}

\begin{proof}
$\norm{\lambda_g f - f}_\infty = \sup_x |f(g^{-1}x) - f(x)| \le
\mathrm{Lip}(f)\, d(g^{-1}x, x) = \mathrm{Lip}(f)\, d(e, g)$ by
bi-invariance, so $L(f) \le \mathrm{Lip}(f)$;\ conversely left
translations act transitively and realize every short geodesic
($d(g^{-1}x, x) = d(e, g)$ along the orbit), so the supremum
recovers $\mathrm{Lip}(f)$. The identities $\mathrm{Lip}(f) =
\norm{\nabla f}_\infty = \opnorm{[\DCH, M_f]}$ on the continuum are
standard (Connes' distance formula on a closed spin manifold).
\end{proof}

\begin{lemma}[per-band injectivity]\label{lem:band_injectivity}
For $1 \le N \le 2\nmax - 1$ let $E_N \subset C^\infty(\sthree)$ be
the degree-$(N{-}1)$ harmonic band ($\dim E_N = N^2$), an irreducible
$\mathrm{SO}(4)$-module of type $V_{(N-1)/2} \otimes V_{(N-1)/2}$.
The compression $\iota_N : E_N \to B(\Hilb_{\nmax})$,
$f \mapsto P_{\nmax} M_f P_{\nmax}$, is injective.
\end{lemma}

\begin{proof}
$\iota_N$ is $\mathrm{SO}(4)$-equivariant for the conjugation action
(rotations implemented on $\Hilb_{\nmax}$ commute with $P_{\nmax}$).
Since $E_N$ is irreducible, $\iota_N$ is either zero or injective
(Schur). It is nonzero:\ every multiplier matrix $\Mat_{NLM}$ with
$N \le 2\nmax - 1$ is nonzero, by non-vanishing of the stretched
$3$-$Y$ matrix element within the top kept shell (Avery--Wen--Avery
closed form~\cite{avery_wen_avery1986});\ verified
symbolically/numerically for every band at $\nmax \le 5$, where the
per-band rank additionally equals the full $N^2$
(\texttt{debug/p38\_g1g2\_band\_diagnostics.py}, frozen in
\texttt{tests/test\_p38\_action\_seminorm.py}).
\end{proof}

\begin{proposition}[kernel condition on the truthful substrate]
\label{prop:kernel_condition}
$L_{\nmax}(T) = 0$ if and only if $T \in \C 1$.
\end{proposition}

\begin{proof}
$L_{\nmax}(T) = 0$ means $\rho_g(T) = T$ for all $g$. Write $T =
P_{\nmax} M_f P_{\nmax}$ with $f$ band-limited to $N \le 2\nmax - 1$
(every element of $\Op_{\nmax}$ has this form). Then $P_{\nmax}
M_{\lambda_g f - f} P_{\nmax} = 0$;\ left translation preserves each
harmonic band, so $\lambda_g f - f$ is band-limited and
Lemma~\ref{lem:band_injectivity} gives $\lambda_g f = f$ for all
$g$;\ transitivity forces $f$ constant. The converse is immediate.
Finite dimensionality plus the kernel condition make
$(S(\Op_{\nmax}), \mathrm{MK}_{L_{\nmax}})$ a compact metric space of
finite diameter, and Lemma~\ref{lem:continuum_lip} provides the same
on the continuum side --- the standard Monge--Kantorovich setting on
both ends.
\end{proof}

\begin{lemma}[spinor window inclusion and the lifted state]
\label{lem:lifted_state}
Let $J = (\nmax - 1)/2$, let
$h = \sum_{j \le J} \sqrt{2j+1}\,\chi_j$ be the central spectral
Fej\'er amplitude of Definition~\ref{def:central_fejer} (the scalar
window then has exactly $\nmax$ levels, and $|h|^2 / Z_{\nmax}$ is
the Fej\'er kernel $\KS_{\nmax}$), fix a unit spinor $\chi \in
\C^2$, and set $\xi := (h \otimes \chi)/\sqrt{Z_{\nmax}}$. Then:
\begin{enumerate}\setlength{\itemsep}{2pt}
\item[(a)] (\emph{Exact window inclusion.}) $\xi \in \Hilb_{\nmax}$.
Under $\mathrm{SU}(2)_L \times \mathrm{SU}(2)_R$,
\begin{equation*}
L^2(\SU(2)) \otimes \C^2 \;=\; \bigoplus_{j} V_j \otimes
\bigl(V_{j+1/2} \oplus V_{j-1/2}\bigr),
\end{equation*}
and these summands are precisely the Camporesi--Higuchi
eigenspaces~\cite{camporesi_higuchi1996}:\ $V_j \otimes V_{j+1/2}$
is the positive shell $n = 2j$ and $V_j \otimes V_{j-1/2}$ the
negative shell $n = 2j-1$, with dimensions $(n+1)(n+2)$ matching
shell by shell. Hence $h \otimes \chi$ lies in the union of the
\emph{complete} shells $n \le 2J = \nmax - 1$ --- exact containment,
no spillover.
\item[(b)] (\emph{Lifted state.}) For every $f \in C(\sthree)$,
$\langle \xi, (P_{\nmax} M_f P_{\nmax})\,\xi\rangle =
\int_{\sthree} f\, \KS_{\nmax}\, d\mu_{\mathrm{Haar}}$:\ the induced
measure is exactly the Fej\'er measure, with first moment
$\gamma_{\nmax}$ (Lemma~\ref{lem:L2}(d)).
\item[(c)] (\emph{Dual map.}) $\upsilon : \Op_{\nmax} \to
C(\sthree)$, $\upsilon(T)(g) := \langle U_g \xi,\, T\, U_g
\xi\rangle$, is unital and positive, satisfies the tautological
contraction $L(\upsilon(T)) \le L_{\nmax}(T)$, and acts on compressed
multipliers as Fej\'er smoothing:\
$\upsilon(P_{\nmax} M_f P_{\nmax}) = \KS_{\nmax} * f$.
\end{enumerate}
\end{lemma}

\begin{proof}
(a) is the displayed decomposition ($\C^2$ in the right tensor slot;\
$V_j \otimes V_{1/2} = V_{j+1/2} \oplus V_{j-1/2}$), with the
eigenvalue and dimension identifications of
Camporesi--Higuchi~\cite{camporesi_higuchi1996}. (b):\ $\xi \in
\mathrm{ran}\, P_{\nmax}$ absorbs the projections, so the pairing is
$\int f\, |h|^2\, |\chi|^2 / Z_{\nmax} = \int f\, \KS_{\nmax}$.
(c):\ positivity and unitality are immediate from the vector-state
form. The contraction:\ for $q = g' g^{-1}$,
$|\upsilon(T)(g) - \upsilon(T)(g')| = |\langle \xi,
(\rho_{g^{-1}}(T) - \rho_{g'^{-1}}(T))\, \xi\rangle| \le
\opnorm{\rho_q(T) - T} \le L_{\nmax}(T)\, d(e, q) =
L_{\nmax}(T)\, d(g, g')$ by bi-invariance. The smoothing identity:\
$\upsilon(P M_f P)(g) = \int f(x)\, \KS_{\nmax}(g^{-1}x)\, dx =
(\KS_{\nmax} * f)(g)$, using centrality and inversion-invariance of
$\KS_{\nmax}$;\ verified to machine precision ($\le 3 \times
10^{-16}$ across bands and random group elements,
\texttt{debug/p38\_g1g2\_scalar\_prototype.py}).
\end{proof}

\begin{theorem}[Main theorem, unconditional]
\label{thm:main_unconditional}
For every $\nmax \ge 1$,
\begin{equation}\label{eq:unconditional_bound}
d_{\mathrm{GH}}\Bigl( \bigl(S(\Op_{\nmax}),
\mathrm{MK}_{L_{\nmax}}\bigr),\;
\bigl(S(C(\sthree)), \mathrm{MK}_{L}\bigr) \Bigr)
\;\le\; \gamma_{\nmax}
\;=\; \frac{4}{\pi}\,\frac{\log \nmax}{\nmax} + O(1/\nmax).
\end{equation}
\end{theorem}

\begin{proof}
Let $S : C(\sthree) \to \Op_{\nmax}$, $S(f) = P_{\nmax} M_f
P_{\nmax}$ (unital completely positive) and let $\upsilon$ be the
dual map of Lemma~\ref{lem:lifted_state}. Four facts:
\begin{enumerate}\setlength{\itemsep}{1pt}
\item[(i)] $L_{\nmax}(S(f)) \le L(f)$:\ $\rho_g(S(f)) - S(f) =
S(\lambda_g f - f)$ and compression contracts operator norms.
\item[(ii)] $L(\upsilon(T)) \le L_{\nmax}(T)$
(Lemma~\ref{lem:lifted_state}(c)).
\item[(iii)] $\norm{\upsilon(S(f)) - f}_\infty =
\norm{\KS_{\nmax} * f - f}_\infty \le \gamma_{\nmax}\, L(f)$ (the
good-kernel estimate in the proof of Lemma~\ref{lem:L4}(c), at the
function level).
\item[(iv)] $\opnorm{S(\upsilon(T)) - T} \le \gamma_{\nmax}\,
L_{\nmax}(T)$:\ on $T = P M_f P$,
\begin{equation*}
S(\upsilon(T)) = P\, M_{\KS_{\nmax} * f}\, P
= \int_{\SU(2)} \rho_g(T)\, \KS_{\nmax}(g)\, dg =: \Phi(T)
\end{equation*}
--- the operator-side conjugation average against the same Fej\'er
measure --- and $\opnorm{\Phi(T) - T} \le \int \opnorm{\rho_g(T) -
T}\, \KS_{\nmax}(g)\, dg \le L_{\nmax}(T) \int d(e, g)\,
\KS_{\nmax}(g)\, dg = \gamma_{\nmax}\, L_{\nmax}(T)$.
\end{enumerate}
By (i)--(ii) the state-space pullbacks $S^* : S(\Op_{\nmax}) \to
S(C(\sthree))$ and $\upsilon^* : S(C(\sthree)) \to S(\Op_{\nmax})$
are $1$-Lipschitz for the Monge--Kantorovich metrics;\ by
(iii)--(iv) they are $\gamma_{\nmax}$-almost mutually inverse, e.g.\
$d(S^*\upsilon^*\omega, \omega) = \sup_{L(f) \le 1}
|\omega(\upsilon(S(f))) - \omega(f)| \le \gamma_{\nmax}$. The
correspondence $\{(\omega, \upsilon^*\omega)\} \cup \{(S^*\phi,
\phi)\}$ then has distortion $\le 2\gamma_{\nmax}$ (the standard
four-pair approximate-isometry computation), whence
$d_{\mathrm{GH}} \le \gamma_{\nmax}$. The rate is
Lemma~\ref{lem:L2}(d.ii)--(d.iii). \qedhere
\end{proof}

\begin{remark}[no inverse estimate;\ the Berezin pair as refinement]
\label{rem:no_inverse}
Both almost-inverse defects in the proof are Fej\'er smoothings at
the same moment $\gamma_{\nmax}$ --- one on functions, one as the
conjugation average $\Phi$ --- so no partial-inverse or transference
estimate is required anywhere:\ the dual-reach gap of earlier
formulations is dissolved rather than closed. The Berezin map
$B_{\nmax}$ of Definition~\ref{def:berezin} (Lemmas~L4--L5) remains
as an explicit forward-reconstruction refinement;\ it is no longer
load-bearing for the distance bound. The scalar $C(G)$ case of this
compression/lifted-state pattern, for general compact metric groups
with the same single-moment bound, is due to
Gaudillot-Estrada--van~Suijlekom~\cite{gaudillot_vs2023};\ the
contribution here is the spinor-window transport (the exact-fit
inclusion of Lemma~\ref{lem:lifted_state}(a)), per-band
multiplicity-one injectivity on the chirality-doubled system, the
closed-form rate moment, and the $4/\pi$ asymptotic constant. All
assembly steps above are proved self-contained.
\end{remark}

\begin{remark}[degeneracy of the Dirac-commutator seminorm]
\label{rem:dirac_degeneracy}
The truthful chirality-diagonal commutator seminorm
$\opnorm{[\DCH, \cdot]}$ on $\Op_{\nmax}$ has kernel strictly larger
than the scalars at every finite cutoff ($10/14$ multipliers at
$\nmax = 2$, $26/55$ at $\nmax = 3$, bit-exact;\ every top-band
multiplier compresses shell-diagonally), so it does not metrize the
truncated state space. An engineered off-diagonal modification of
$\DCH$ restores its kernel condition ($1/14$, $1/55$) but is no
longer needed for the main theorem. Frozen falsifiers:\
\texttt{tests/test\_p45\_kplus\_degeneracy.py},
\texttt{tests/test\_p38\_action\_seminorm.py}.
\end{remark}

\begin{remark}[history]\label{rem:history38}
An earlier draft of this paper (May 2026;\ repository history and the
corresponding Zenodo deposit) stated the main theorem in the
Latr\'emoli\`ere propinquity, with a Lemma~L2 that conflated the
kernel's Plancherel mass distribution with its multiplier symbol and
a dual-reach argument resting on that conflation. A verification
pass (2026-06-09;\ repository \texttt{CHANGELOG}, v3.106.0) produced
an intermediate conditional form;\ the present unconditional form
(2026-06-10) follows from the translation-seminorm reframing of this
subsection. The companion Lorentzian paper proves the related
$K^{+}$-compression degeneracy theorem.
\end{remark}'''

OLD_SECTION = r'''\subsection{Named gaps and history}\label{sec:named_gaps}

The main theorem is conditional on two named gaps.
\begin{enumerate}
\item[(G1)] \emph{Dual-direction reach.} The estimate
$\mathrm{reach}_P \le \gamma_{\nmax}$ requires a bounded partial
inverse of the Berezin map on the relevant window, which no forward
bound supplies (an inverse is governed by the multiplier-symbol
\emph{minimum};\ on the full support window the inverse norm grows
like $O(\nmax^{2})$). The correct mechanism is the
antiderivative transference of
Leimbach--van~Suijlekom~\cite{leimbach_vs2024} (Lemmas~3.4--3.5
there), whose conclusion for scalar sectors on compact metric groups
is covered by
Gaudillot-Estrada--van~Suijlekom~\cite{gaudillot_vs2023};\ its
adaptation to the chirality-doubled spinor substrate is open.
\item[(G2)] \emph{Kernel condition and the Dirac substrate.} For the
truthful chirality-diagonal Camporesi--Higuchi Dirac, the Lipschitz
seminorm $\opnorm{[\DCH, \cdot]}$ on the truncated multiplier system
has kernel far exceeding the scalars (verified bit-exactly:\ $10/14$
multipliers at $\nmax = 2$, $26/55$ at $\nmax = 3$), so the truncated
object is not a compact quantum metric space and its state-space
metric is degenerate. The kernel condition holds for an engineered off-diagonal
modification of $\DCH$ with E1-selection-pattern couplings
(cf.~\cite{paper32}~\S~III;\ kernel $= \C 1$ exactly:\ $1/14$,
$1/55$), at the cost that the modified operator is a
hand-parameterized stand-in for the genuine spinor-bundle
off-diagonal structure;\ no proved truthful-to-off-diagonal
comparison is currently available, and L3's constant must be
re-established on that substrate. Frozen falsifier:\
\texttt{tests/test\_p45\_kplus\_degeneracy.py} in the repository.
\end{enumerate}

\begin{remark}[History and retraction]\label{rem:history38}
An earlier draft of this paper (May 2026;\ archived in the repository
history and the corresponding Zenodo deposit) stated the main theorem
unconditionally in the Latr\'emoli\`ere propinquity, with a
Lemma~L2 that conflated the kernel's Plancherel mass distribution
with its multiplier symbol and a dual-reach argument resting on that
conflation. A verification pass (2026-06-09;\ audit record in the
repository \texttt{CHANGELOG}, v3.106.0) produced the present
corrected, conditional form;\ the companion Lorentzian paper proves
the related $K^{+}$-compression degeneracy theorem and explains why
the Lorentzian extension of the present construction failed.
\end{remark}'''

rep(OLD_SECTION, NEW_SECTION, 'named-gaps-section')

# ---------------------------------------------------------------- 2.
rep(r'''We prove — subject to two named gaps (\S~\ref{sec:named_gaps}) —
that the Connes--van~Suijlekom truncations of the round
$\sthree = \SU(2)$ Camporesi--Higuchi spectral triple converge to the
continuum spectral triple in van~Suijlekom's state-space
Gromov--Hausdorff distance~\cite{vs2021_jgp}, with explicit asymptotic
rate
\(
\Lambda(\Tcal_{\nmax}, \Tcal_{\sthree}) = O\!\bigl(\log \nmax / \nmax\bigr)
\)
and asymptotic constant $4/\pi$. The two gaps:\ the dual-direction
reach estimate (antiderivative transference, adapted to the spinor
substrate), and the kernel condition
$L(a) = 0 \Leftrightarrow a \in \C 1$, which fails for the truthful
chirality-diagonal Dirac and is verified numerically only for the
engineered off-diagonal Dirac.
The proof proceeds via five lemmas:''',
    r'''We prove, unconditionally, that the Connes--van~Suijlekom
truncations of the round $\sthree = \SU(2)$ Camporesi--Higuchi
spectral triple converge to the continuum spectral triple in
van~Suijlekom's state-space Gromov--Hausdorff
distance~\cite{vs2021_jgp}, with explicit rate
\(
d_{\mathrm{GH}} \le \gamma_{\nmax}
= (4/\pi + o(1)) \log \nmax / \nmax .
\)
The truncated state spaces are metrized by the left-translation
Lipschitz seminorm, which equals the Dirac/metric Lipschitz constant
on the continuum and whose kernel at every finite cutoff is exactly
the scalars (the truthful Dirac-commutator seminorm itself
degenerates at finite cutoff;\ a remark quantifies this). The proof
is a compression/lifted-state argument in the pattern of the
compact-metric-group theorem of
Gaudillot-Estrada--van~Suijlekom~\cite{gaudillot_vs2023}, with
self-contained proofs:\ the dual map evaluates conjugated operators
against the central spectral Fej\'er vector state, which embeds
\emph{exactly} into the spinor window via the
$\mathrm{SU}(2)_L \times \mathrm{SU}(2)_R$ block decomposition
$V_j \otimes (V_{j \pm 1/2})$ of the Camporesi--Higuchi eigenspaces,
and both almost-inverse defects collapse to Fej\'er smoothing — one
on functions, one as operator conjugation-averaging — at the common
moment $\gamma_{\nmax}$. The contribution relative to the scalar
theory:\ the spinor-window transport, per-band multiplicity-one
injectivity on the chirality-doubled multiplier system, the
closed-form rate moment, and the $4/\pi$ asymptotic constant.
Supporting structure is organized in five lemmas:''',
    'abstract')

# ---------------------------------------------------------------- 3.
rep(r'''\begin{theorem}[Main theorem, conditional;\ cf.~\S\ref{sec:main_theorem}]
\label{thm:main_intro}
Subject to the two named gaps of \S~\ref{sec:named_gaps}
(dual-direction reach;\ kernel condition on the off-diagonal Dirac
substrate), the truncated triples $\Tcal_{\nmax}$ converge to
$\Tcal_{\sthree}$ in van~Suijlekom's state-space Gromov--Hausdorff
distance~\cite{vs2021_jgp}, with explicit rate''',
    r'''\begin{theorem}[Main theorem;\ cf.~\S\ref{sec:named_gaps} and
\S\ref{sec:main_theorem}]
\label{thm:main_intro}
Unconditionally, the truncated triples $\Tcal_{\nmax}$ converge to
$\Tcal_{\sthree}$ in van~Suijlekom's state-space Gromov--Hausdorff
distance~\cite{vs2021_jgp}, under the translation-seminorm
metrization of \S~\ref{sec:named_gaps}, with explicit rate''',
    'thm-intro')

# ---------------------------------------------------------------- 4.
rep(r'''The present paper treats that case, subject to the two named gaps
of \S~\ref{sec:named_gaps}.''',
    r'''The present paper treats that case;\ the theorem is
unconditional under the translation-seminorm metrization of
\S~\ref{sec:named_gaps}.''', 'treats-case')

# ---------------------------------------------------------------- 5.
rep(r'''reach$_{B}$, reach$_{P}$, height$_{B}$, height$_{P}$ for the four
approximation estimates the pair supplies. For the truncated-side
Monge--Kantorovich metric to be non-degenerate the kernel condition
of \S~\ref{sec:named_gaps}(G2) must hold, which conditions the
statement on the off-diagonal Dirac substrate.''',
    r'''reach$_{B}$, reach$_{P}$, height$_{B}$, height$_{P}$ for the four
approximation estimates the pair supplies. The truncated-side
Monge--Kantorovich metric is non-degenerate unconditionally under the
translation-seminorm metrization
(Proposition~\ref{prop:kernel_condition});\ the distance bound itself
is Theorem~\ref{thm:main_unconditional}, and the pair below supplies
the explicit forward-reconstruction refinement.''', 'L5-preamble')

# ---------------------------------------------------------------- 6.
rep(r'''\begin{lemma}[L5:\ distance bound from the approximation pair,
conditional]\label{lem:L5}
Subject to the two named gaps of \S~\ref{sec:named_gaps}
(the dual-direction reach estimate G1;\ the kernel condition G2,
conditioning the statement on the off-diagonal Dirac substrate), the
pair''',
    r'''\begin{lemma}[L5:\ Berezin-pair refinement of the forward
direction]\label{lem:L5}
Within the framework of \S~\ref{sec:named_gaps} (kernel condition
Proposition~\ref{prop:kernel_condition};\ distance bound
Theorem~\ref{thm:main_unconditional}), the pair''', 'lem-L5')

# ---------------------------------------------------------------- 7.
rep(r'''\textbf{(Named gap G1.)} The dual estimate
$\mathrm{reach}_P \le \gamma_{\nmax}$ cannot be inferred from any
forward cb-norm bound:\ a partial inverse is governed by the symbol
\emph{minimum}, and on the full support window the inverse norm grows
like $O(\nmax^{2})$. The correct route is the
antiderivative-transference mechanism of
Leimbach--van~Suijlekom~\cite{leimbach_vs2024} (Lemmas~3.4--3.5),
whose conclusion for scalar sectors on compact metric groups is
covered by Gaudillot-Estrada--van~Suijlekom~\cite{gaudillot_vs2023};\
its adaptation to the chirality-doubled spinor substrate is the open
named gap on which this lemma is conditional. If the adapted constant
exceeds $1$, the right-hand side of (\ref{eq:L5_bound}) degrades by
that $O(1)$ factor;\ the rate class $O(\log\nmax/\nmax)$ is
unaffected.''',
    r'''\textbf{(Dissolved.)} The dual estimate is no longer needed:\
Theorem~\ref{thm:main_unconditional} obtains the distance bound via
the lifted-state map, whose almost-inverse defect is the conjugation
average $\Phi$ at the same moment $\gamma_{\nmax}$
(Remark~\ref{rem:no_inverse}) — no partial inverse of $B_{\nmax}$,
and no transference estimate, enters anywhere. (For the record:\ a
partial inverse of $B_{\nmax}$ on the full support window would be
governed by the symbol \emph{minimum} and grow like $O(\nmax^2)$,
which is why no forward bound could supply it.) The pair $(B_{\nmax},
P_{\nmax})$ retains its role as the explicit forward-reconstruction
refinement.''', 'reachP')

# ---------------------------------------------------------------- 8.
rep(r'''\begin{theorem}[Main theorem;\ rate-controlled, conditional on the
two named gaps of \S~\ref{sec:named_gaps}]\label{thm:main}
Subject to the dual-direction reach estimate (G1) and the kernel
condition on the off-diagonal Dirac substrate (G2),
the truncated Camporesi--Higuchi spectral triples $\Tcal_{\nmax}$
converge to the round-$\sthree$ Camporesi--Higuchi spectral triple
$\Tcal_{\sthree}$ in van~Suijlekom's state-space Gromov--Hausdorff
distance~\cite{vs2021_jgp}:''',
    r'''\begin{theorem}[Main theorem;\ rate-controlled restatement]
\label{thm:main}
Unconditionally (Theorem~\ref{thm:main_unconditional}), the
truncated Camporesi--Higuchi spectral triples $\Tcal_{\nmax}$
converge to the round-$\sthree$ Camporesi--Higuchi spectral triple
$\Tcal_{\sthree}$ in van~Suijlekom's state-space Gromov--Hausdorff
distance~\cite{vs2021_jgp}, under the translation-seminorm
metrization of \S~\ref{sec:named_gaps}:''',
    'thm-main')

# ---------------------------------------------------------------- 9.
rep(r'''\begin{proposition}[Limit identification, conditional]
\label{prop:limit_identification}
Under the same conditionality as Theorem~\ref{thm:main}, the
state-space GH limit of the sequence $(\Tcal_{\nmax})_{\nmax \ge 1}$
is''',
    r'''\begin{proposition}[Limit identification]
\label{prop:limit_identification}
The state-space GH limit of the sequence
$(\Tcal_{\nmax})_{\nmax \ge 1}$ is''', 'limit-id')

io.open(P, 'w', encoding='utf-8').write(s)
print('FAILS:', fails if fails else 'none')
