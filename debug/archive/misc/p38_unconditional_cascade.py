"""Cascade the Paper 38 unconditional theorem to P45, P32, claims
register, README, field guide, and the N1 outreach note."""
import io

fails = []


def rep(path, old, new, tag):
    s = io.open(path, encoding='utf-8').read()
    if s.count(old) == 1:
        io.open(path, 'w', encoding='utf-8').write(s.replace(old, new))
    else:
        fails.append((tag, s.count(old)))


P45 = r'papers/group1_operator_algebras/paper_45_lorentzian_propinquity.tex'
P32 = r'papers/group1_operator_algebras/paper_32_spectral_triple.tex'
FG = r'papers/synthesis/geovac_field_guide.tex'
N1 = r'docs/outreach/note_n1_su2_truncations.tex'
CR = r'docs/claims_register.md'
RM = r'README.md'

# ----------------------------------------------------------- Paper 45
rep(P45, r'''and a conditional spatial statement
(Proposition~\ref{prop:spatial_conditional}):\ the spatial truncations
converge in van~Suijlekom's state-space Gromov--Hausdorff
distance~\cite{vs2021_jgp,gaudillot_vs2023} at the companion-paper
rate $O(\log\nmax/\nmax)$, conditional on the Paper~38 erratum repairs,
with the temporal factor a metrically
invisible spectator carried $\JL$-equivariantly.''',
    r'''and a spatial convergence statement
(Proposition~\ref{prop:spatial_conditional}):\ the spatial truncations
converge in van~Suijlekom's state-space Gromov--Hausdorff
distance~\cite{vs2021_jgp,gaudillot_vs2023} at the companion-paper
rate $O(\log\nmax/\nmax)$ — unconditional, via the companion paper's
translation-seminorm theorem — with the temporal factor a metrically
invisible spectator carried $\JL$-equivariantly.''', 'p45-abstract')

rep(P45, r'''(\S~\ref{sec:withdrawn_assembly}), and the conditional spatial
statement (Proposition~\ref{prop:spatial_conditional}).''',
    r'''(\S~\ref{sec:withdrawn_assembly}), and the spatial convergence
statement (Proposition~\ref{prop:spatial_conditional}), unconditional
via the companion paper's translation-seminorm theorem.''',
    'p45-list')

rep(P45, r'''annihilation theorem, dissects the failure of the compression
assembly, and states the conditional spatial convergence statement.''',
    r'''annihilation theorem, dissects the failure of the compression
assembly, and states the spatial convergence statement.''',
    'p45-outline')

rep(P45, r'\subsection{Conditional spatial statement}\label{sec:main_theorem}',
    r'\subsection{The spatial convergence statement}\label{sec:main_theorem}',
    'p45-subsec')

rep(P45, r'''\begin{proposition}[Spatial sector, conditional]
\label{prop:spatial_conditional}
Fix $\Nt \ge 1$, $T > 0$. Conditional on the companion $\sthree$
paper's spatial Gromov--Hausdorff convergence statement
(van~Suijlekom state-space distance~\cite{vs2021_jgp};\ its two named
gaps:\ the dual-direction reach estimate via antiderivative
transference~\cite{leimbach_vs2024,gaudillot_vs2023}, and the kernel
condition on the engineered off-diagonal Dirac substrate), the
spatial truncations underlying''',
    r'''\begin{proposition}[Spatial sector]
\label{prop:spatial_conditional}
Fix $\Nt \ge 1$, $T > 0$. By the companion $\sthree$ paper's
(unconditional) spatial Gromov--Hausdorff convergence theorem
(van~Suijlekom state-space distance~\cite{vs2021_jgp}, with the
truncated state space metrized by the left-translation Lipschitz
seminorm;\ compression/lifted-state argument in the pattern
of~\cite{gaudillot_vs2023}), the spatial truncations underlying''',
    'p45-prop')

rep(P45, r'''all estimates are the
Paper~38 estimates verbatim, subject to that paper's two named gaps
(which is the conditionality).''',
    r'''all estimates are the
Paper~38 estimates verbatim, now unconditional under that paper's
translation-seminorm metrization.''', 'p45-proof')

rep(P45, r'''is nonetheless a well-defined closed-form expression;\ its surviving
meaning is the conditional spatial statement of
Proposition~\ref{prop:spatial_conditional}, to which we now turn.''',
    r'''is nonetheless a well-defined closed-form expression;\ its surviving
meaning is the spatial convergence statement of
Proposition~\ref{prop:spatial_conditional}, to which we now turn.''',
    'p45-eq-tail')

rep(P45, r'''The ``$\Lprop$ bound'' column is the
closed-form expression $\Cthreejoint\cdot\gammajoint$, whose surviving
meaning is the conditional spatial rate of
Proposition~\ref{prop:spatial_conditional};''',
    r'''The ``$\Lprop$ bound'' column is the
closed-form expression $\Cthreejoint\cdot\gammajoint$, whose surviving
meaning is the spatial rate of
Proposition~\ref{prop:spatial_conditional};''', 'p45-caption')

# ----------------------------------------------------------- Paper 32
rep(P32, r'''but \emph{fails} for the truthful chirality-diagonal CH Dirac used in
the statement of the theorem ($10/14$, $26/55$); the
truthful$\leftrightarrow$offdiag bridge is a named gap.  The spatial
convergence content is expected to survive in the van~Suijlekom
state-space framework.
\end{remark}''',
    r'''but \emph{fails} for the truthful chirality-diagonal CH Dirac
($10/14$, $26/55$).  Items (ii)--(iv) are resolved (2026-06-10) by
the translation-seminorm reframing of the companion $S^3$ paper:\ the
truncated state space is metrized by the left-translation Lipschitz
seminorm (equal to the Lipschitz constant on the continuum), whose
kernel is exactly the scalars on the \emph{truthful} substrate
(Schur multiplicity-one plus per-band injectivity of band-limited
compression);\ the dual direction is supplied by an exact-fit spinor
lifted state (the Fej\'er vector $h \otimes \chi$, lying in the
window by the $V_j \otimes V_{j \pm 1/2}$ block decomposition of the
Camporesi--Higuchi eigenspaces), and both almost-inverse defects are
Fej\'er smoothings at the common moment $\gamma_{n_{\max}}$ --- no
partial inverse or transference estimate remains.  The convergence
theorem is unconditional in van~Suijlekom's state-space framework at
rate $(4/\pi + o(1))\log n_{\max}/n_{\max}$.
\end{remark}''', 'p32-caveat')

rep(P32, r'''This closure claim is qualified by the
2026-06-09 status caveat (Remark~\ref{rem:gh_status_caveat}):\ the
proof chain is under repair, with the spatial content expected to
survive in the van~Suijlekom state-space framework.''',
    r'''The proof chain was corrected in June 2026
(Remark~\ref{rem:gh_status_caveat}):\ the theorem stands
unconditionally in van~Suijlekom's state-space framework under the
translation-seminorm metrization.''', 'p32-wh1')

# ----------------------------------------------------- claims register
rep(CR, r'''| 8 | State-space GH convergence of truncated Camporesi–Higuchi (Dirac) triples on S³, rate (4/π)log n/n | Paper 38 | CONDITIONAL | Gaps: G1 dual-direction reach (transference not adapted to spinor substrate); G2 kernel condition (holds for engineered off-diagonal Dirac 1/14, 1/55; fails truthful CH 10/14, 26/55). Scalar-sector convergence for compact groups is **prior art** (Gaudillot-Estrada–van Suijlekom, arXiv:2310.14733) |''',
    r'''| 8 | State-space GH convergence of truncated Camporesi–Higuchi (Dirac) triples on S³, rate (4/π)log n/n, translation-seminorm metrization | Paper 38 | INTERNAL THEOREM (unconditional) | Former gaps closed 2026-06-10: kernel condition proved on the truthful substrate (Schur multiplicity-one + per-band injectivity; `tests/test_p38_action_seminorm.py`); dual direction via exact-fit spinor lifted state — both almost-inverse defects are Fejér smoothings at γ_n, no inverse estimate remains. Scalar-sector convergence for compact groups is **prior art** (Gaudillot-Estrada–van Suijlekom, arXiv:2310.14733); contribution = spinor transport + per-band injectivity + explicit 4/π rate |''',
    'cr-row8')

# --------------------------------------------------------------- README
rep(RM, r'''- **Paper 38** State-space Gromov–Hausdorff convergence of truncated Dirac triples on SU(2), explicit 4/π rate — **conditional on two named gaps** (dual-direction reach; kernel condition on the off-diagonal Dirac substrate)''',
    r'''- **Paper 38** State-space Gromov–Hausdorff convergence of truncated Dirac triples on SU(2), explicit 4/π rate — **unconditional** (2026-06-10: translation-seminorm metrization; kernel condition proved on the truthful substrate; dual direction via an exact-fit spinor lifted state)''',
    'rm-bullet')

rep(RM, r'''| **38** | **SU(2) state-space GH convergence** | **Five-lemma proof, 4/π rate; conditional on two named gaps (see claims register)** |''',
    r'''| **38** | **SU(2) state-space GH convergence** | **Unconditional; 4/π rate; compression/lifted-state proof, translation seminorm (see claims register)** |''',
    'rm-table')

# ---------------------------------------------------------- field guide
rep(FG, r'''with the constant $4/\pi$ universal across compact connected Lie
groups~\cite{loutey_paper40}.  A June 2026 adversarial verification
pass left the theorem standing \emph{conditionally}, with two named
gaps (a dual-direction reach estimate, and a kernel condition tied to
the choice of Dirac substrate;\ see Paper~38~\S~1 ``Named gaps and
history'' and the claims register,
\texttt{docs/claims\_register.md}).''',
    r'''with the constant $4/\pi$ universal across compact connected Lie
groups~\cite{loutey_paper40}.  A June 2026 adversarial verification
pass first reduced the theorem to a conditional form with two named
gaps, and the gaps were then closed by a translation-seminorm
reframing (the kernel condition is proved on the truthful substrate;\
the dual direction uses an exact-fit spinor lifted state, and both
almost-inverse defects are Fej\'er smoothings at the same moment):\
the theorem is unconditional in van~Suijlekom's state-space
framework.  See Paper~38~\S~1 and the claims register,
\texttt{docs/claims\_register.md}.''', 'fg-passage')

rep(FG, r'''    Math.OA arc.  Headline:\ conditional SU(2) state-space
    Gromov--Hausdorff convergence with explicit $4/\pi$ rate
    (Paper~38, two named gaps);''',
    r'''    Math.OA arc.  Headline:\ unconditional SU(2) state-space
    Gromov--Hausdorff convergence with explicit $4/\pi$ rate
    (Paper~38, translation-seminorm metrization);''', 'fg-group1')

rep(FG, r'''    Papers~38, 45, and 32.  Paper~38 is the conditional SU(2)
    state-space Gromov--Hausdorff convergence theorem (two named
    gaps);''',
    r'''    Papers~38, 45, and 32.  Paper~38 is the unconditional SU(2)
    state-space Gromov--Hausdorff convergence theorem
    (translation-seminorm metrization, compression/lifted-state
    proof);''', 'fg-readers')

# ------------------------------------------------------------ N1 note
rep(N1, r'''specialists. Verification status is flagged sentence by sentence;\
two named gaps are exactly that --- gaps.''',
    r'''specialists. Verification status is flagged sentence by
sentence;\ proofs are self-contained in the full write-ups, and
anything verified only numerically is flagged as such.''', 'n1-open')

rep(N1, r'''\section{Result A:\ conditional convergence (two named gaps)}

\begin{theorem}[conditional]
Subject to gaps G1--G2 below, the truncations of the chirality-doubled
Dirac triple above converge to $\mathcal{T}_{S^3}$ in van~Suijlekom's
state-space Gromov--Hausdorff distance, with
$d \le C_3\,\gamma_n = \tfrac{4\log n}{\pi n} + O(1/n)$.
\end{theorem}

\noindent\textbf{G1 (dual-direction reach).} The reverse
approximation direction needs a bounded partial inverse of $B_n$ on
the relevant window;\ no forward bound supplies it (the inverse is
governed by the symbol \emph{minimum};\ on the full support window it
grows like $O(n^2)$). The natural route is the antiderivative
transference of Leimbach--van~Suijlekom \cite{LvS2024} (their
Lemmas~3.4--3.5), whose conclusion for scalar sectors on compact
metric groups is covered by \cite{GEvS2023};\ its adaptation to the
spinor-multiplier substrate is open.

\smallskip
\noindent\textbf{G2 (kernel condition).} For the truthful
chirality-diagonal Dirac, the seminorm $\opnorm{[\DCH, \cdot]}$ on
the truncated multiplier system has kernel far exceeding the scalars
(bit-exactly:\ $10/14$ basis multipliers at cutoff $n = 2$, $26/55$
at $n = 3$), so the truncated object is not a compact quantum metric
space. An engineered off-diagonal modification of $\DCH$ restores
the kernel condition exactly ($1/14$, $1/55$ --- identity only), but
no comparison theorem between the truthful and modified operators is
currently available, and the Lipschitz constant must be
re-established on the modified substrate.''',
    r'''\section{Result A:\ unconditional convergence}

\begin{theorem}
The truncations of the chirality-doubled Dirac triple above converge
to $\mathcal{T}_{S^3}$ in van~Suijlekom's state-space
Gromov--Hausdorff distance, with the truncated state spaces metrized
by the left-translation Lipschitz seminorm
$L_n(T) = \sup_{g \ne e} \norm{U_g T U_g^* - T} / d(e, g)$ (equal to
the metric Lipschitz constant on the continuum):
$d_{\mathrm{GH}} \le \gamma_n = \tfrac{4 \log n}{\pi n} + O(1/n)$.
\end{theorem}

\noindent\emph{Proof shape} (self-contained in the full write-up).
The kernel condition $L_n(T) = 0 \Leftrightarrow T \in \C 1$ holds on
the truthful substrate:\ left translation preserves harmonic bands,
band-limited compression is injective band-by-band (Schur
multiplicity-one;\ the nonvanishing matrix element is the stretched
$3$-$Y$ coefficient), and transitivity does the rest. The dual map
evaluates conjugated operators against the Fej\'er vector state
$\xi = (h \otimes \chi)/\sqrt{Z}$, which lies \emph{exactly} in the
Dirac window:\ under $\mathrm{SU}(2)_L \times \mathrm{SU}(2)_R$ the
spinor space is $\bigoplus_j V_j \otimes (V_{j+1/2} \oplus
V_{j-1/2})$ and these summands are precisely the Dirac eigenspaces.
Both almost-inverse defects collapse to Fej\'er smoothing — one on
functions, one as the operator conjugation-average — at the common
moment $\gamma_n$, so no partial-inverse or transference estimate
appears anywhere. The truthful Dirac-\emph{commutator} seminorm, by
contrast, degenerates at finite cutoff (kernel $10/14$ at $n = 2$,
$26/55$ at $n = 3$, bit-exact) — it is the wrong finite-cutoff
metrization, and this degeneracy is quantified in the write-up. The
scalar $C(G)$ pattern is \cite{GEvS2023};\ the contribution here is
the spinor-window transport, the per-band injectivity, the
closed-form rate moment, and the $4/\pi$ constant.''', 'n1-resultA')

rep(N1, r'''\begin{question}
For G1:\ is the \cite{LvS2024}/\cite{GEvS2023} transference the right
route on a spinor-multiplier substrate, and is there a standard
perturbation argument for G2's truthful-vs-modified Dirac comparison
(the modification is bounded but not small)?
\end{question}''',
    r'''\begin{question}
For Result A:\ is the finite-cutoff degeneracy of the
Dirac-commutator seminorm — and its repair by the translation
seminorm — known or recorded? And is the exact-fit spinor-window
embedding of the Fej\'er vector state (via
$V_j \otimes V_{j \pm 1/2}$) the construction you would have used,
or is there a slicker route to the dual map on bundle-valued
truncations?
\end{question}''', 'n1-q2')

print('FAILS:', fails if fails else 'none')
