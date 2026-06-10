"""Apply Phase-2d standalone-audit rewords to Papers 38 and 45.

Run from repo root. Uniqueness-asserted replacements; misses reported.
"""
import io

fails = []


def rep(path, old, new, tag):
    s = io.open(path, encoding='utf-8').read()
    if s.count(old) == 1:
        io.open(path, 'w', encoding='utf-8').write(s.replace(old, new))
    else:
        fails.append((tag, s.count(old)))


P38 = r'papers/group1_operator_algebras/paper_38_su2_propinquity_convergence.tex'
P45 = r'papers/group1_operator_algebras/paper_45_lorentzian_propinquity.tex'

# --- Paper 38 ---

rep(P38, r'''The kernel condition holds for the engineered
off-diagonal Dirac of the R3.5 construction (kernel $= \C 1$
exactly:\ $1/14$, $1/55$), at the cost that the off-diagonal operator
is a hand-parameterized stand-in for the genuine spinor-bundle
off-diagonal structure;''',
    r'''The kernel condition holds for an engineered off-diagonal
modification of $\DCH$ with E1-selection-pattern couplings
(cf.~\cite{paper32}~\S~III;\ kernel $= \C 1$ exactly:\ $1/14$,
$1/55$), at the cost that the modified operator is a
hand-parameterized stand-in for the genuine spinor-bundle
off-diagonal structure;''', 'p38-G2')

rep(P38, r'''(repository:\ \texttt{debug/p45\_adversarial\_L2\_memo.md})''',
    r'''(supporting computation in the project repository,
\texttt{debug/p45\_adversarial\_L2\_memo.md})''', 'p38-L2proof')

rep(P38, r'''The closed form of $b$ remains open, and
the constant joins a small documented list of irreducible-but-natural
GeoVac framework constants alongside $K = \pi(B+F-\Delta)$
of~\cite{paper2}.''',
    r'''The closed form of $b$ remains open, and the constant joins a
small documented list of natural-but-unidentified constants in this
programme, alongside the numerical coincidence formula studied
in~\cite{paper2}.''', 'p38-constants')

rep(P38, r'''\textbf{(ii) Lorentzian extension.} The Camporesi--Higuchi spectral
triple is Riemannian, with KO-dimension $3 \pmod 8$. The Lorentzian
analog --- via Wick rotation of $\sthree$ to $H^3$ or to a static
patch of $\mathrm{de Sitter}_3$ --- carries different KO-dimension and
sign convention, and may admit a different propinquity theory entirely.
We are not aware of a Lorentzian propinquity framework in the
literature; the topic is open.

\emph{Update (2026-05-16):} a partial Riemannian-side
operator-system-level result is now in place via the
\emph{modular Hamiltonian closure} (Sprint~L1, see GeoVac framework
documentation Paper~32~\S~VIII.F): the modular Hamiltonian $K$ on the
truncated Camporesi--Higuchi spectral triple $\mathcal{T}_{\nmax}$
satisfies the period closure $\sigma_{2\pi}(O) = O$ bit-exactly at
$\nmax \in \{2, 3, 4, 5\}$ for the hemispheric-wedge-restricted state
and all four physical witnesses (Bisognano--Wichmann, Hartle--Hawking,
Sewell, Unruh). This lifts the four-witness Wick-rotation theorem from
``structural correspondence'' to ``literal identification at the
operator-system level (Riemannian)''. The closure is finite-cutoff
exact (not a propinquity limit-statement), because the geometric BW
realization $K = J_{\mathrm{polar}}$ has integer spectrum on the
full-Dirac basis, so $\mathrm{e}^{\mathrm{i}\cdot 2\pi \cdot n} = 1$
holds for all integer $n$. The genuine Lorentzian extension to
signature $(3, 1)$ via the Bizi--Brouder--Besnard~2018 Krein-space
spectral triple, with corresponding Lorentzian propinquity, remains
the open question of this subsection; Sprint L1's Riemannian closure
is the prerequisite that confirms the operator-system identification
is sound before lifting signature.

\emph{Update (2026-05-16, same day, evening):} Sprint L1-tighten
extends the Sprint L1 closure from the geometric BW-$\alpha$ realization
($K = J_{\mathrm{polar}}$) to the load-bearing Tomita--Takesaki
BW-$\gamma$ construction: $K_{\mathrm{TT}} = -\log \Delta$ extracted
from the polar decomposition $S = J_{\mathrm{TT}} \Delta^{1/2}$ on the
GNS Hilbert--Schmidt space $H_{\mathrm{GNS}} = M_{\dim_W}(\mathbb{C})$.
Verdict at every tested $\nmax \in \{2, 3, 4, 5\}$ and every physical
witness:\ \emph{UNIFIED\_STRONG} closure. The two constructions are
operator-action-level conjugates of each other ($\sigma_t^{\mathrm{TT}}
= \sigma_{-t}^{\alpha}$ bit-exact at general $t$, identical at the
period $t = 2\pi$), and the state-dependent $J_{\mathrm{TT}}(X) =
\rho^{1/2} X^* \rho^{-1/2}$ verifies $J_{\mathrm{TT}}^2 = +I$ at the
full Tomita level. The L1 ``BW-$\alpha$-only'' caveat is removed: the
four-witness Wick-rotation theorem is closed at signature $(3, 0)$ via
both the geometric and the Tomita-Takesaki readings. See GeoVac
framework documentation Paper~32~\S~VIII.F (Sprint L1-tighten extension
paragraph) and
\texttt{debug/l1\_tighten\_tomita\_results\_memo.md}.''',
    r'''\textbf{(ii) Lorentzian extension.} The Camporesi--Higuchi spectral
triple is Riemannian, with KO-dimension $3 \pmod 8$. On the
Riemannian side, the modular-Hamiltonian structure of the truncations
is closed at the operator-system level:\ the modular Hamiltonian on
$\mathcal{T}_{\nmax}$ satisfies the period closure
$\sigma_{2\pi}(O) = O$ bit-exactly at $\nmax \in \{2, 3, 4, 5\}$ for
the hemispheric-wedge-restricted state and all four physical
witnesses (Bisognano--Wichmann, Hartle--Hawking, Sewell, Unruh), in
both the geometric realization ($K = J_{\mathrm{polar}}$, integer
spectrum, so $e^{i \cdot 2\pi n} = 1$ exactly) and the
Tomita--Takesaki realization ($K_{\mathrm{TT}} = -\log\Delta$ on the
GNS Hilbert--Schmidt space), which are flow-conjugate,
$\sigma_t^{\mathrm{TT}} = \sigma_{-t}^{\alpha}$;\ see the four-witness
companion paper~\cite{paper42} and \cite{paper32}~\S~VIII.F. The
genuinely Lorentzian extension at signature $(3, 1)$ via Krein-space
spectral triples is \emph{open}, and is now known to be obstructed
along its most natural route:\ compressing a Krein-self-adjoint
product Dirac to the Krein-positive subspace annihilates the spatial
Dirac, so the compressed Lipschitz seminorm vanishes identically on
the natural operator system and no quantum-metric structure survives.
The companion Lorentzian paper~\cite{paper45} proves this degeneracy
theorem and states the repaired target (a Toeplitz temporal
multiplier algebra in the Connes--van~Suijlekom circle pattern, and a
Krein-compatible Lipschitz seminorm that does not compress the
Dirac).''', 'p38-lorentzian')

rep(P38, r'''\bibitem[Loutey, internal Paper~24]{paper24}''',
    r'''\bibitem[Loutey, internal Paper~32]{paper32}
J.~Loutey,
\newblock \emph{The GeoVac spectral triple:\ construction, Connes
axiom audit, and the $\pi$-source case-exhaustion theorem},
\newblock GeoVac internal preprint (2026), DOI via Zenodo.

\bibitem[Loutey, internal Paper~42]{paper42}
J.~Loutey,
\newblock \emph{Tomita--Takesaki modular structure on truncated
$\SU(2)$ spectral triples:\ four-witness Wick-rotation identification
at finite cutoff},
\newblock GeoVac internal preprint (2026), DOI via Zenodo.

\bibitem[Loutey, internal Paper~45]{paper45}
J.~Loutey,
\newblock \emph{Lorentzian propinquity on truncated
$\SU(2) \times U(1)_T$ Krein spectral triples:\ a degeneracy theorem},
\newblock GeoVac internal preprint (2026), DOI via Zenodo.

\bibitem[Loutey, internal Paper~24]{paper24}''', 'p38-bibitems')

# --- Paper 45 ---

rep(P45, r'''This paper closes the convergence-theorem leg of the
Lorentzian-extension program for truncated GeoVac spectral triples,
following''',
    r'''This paper concerns the convergence-theorem leg of the
Lorentzian-extension programme for truncated spectral triples of the
Geometric Vacuum (GeoVac) framework~\cite{paper32} --- and proves
that the leg's most natural device fails structurally --- following''',
    'p45-intro')

rep(P45, r'''both with the corrected Lemma~L2 convention
of \S~\ref{sec:L2}''',
    r'''both with the mass-distribution/multiplier-symbol
distinction of Lemma~\ref{lem:L2}''', 'p45-l2conv')

rep(P45, r'''The driver
\texttt{debug/l3b\_2\_sub\_sprint\_D\_compute.py} evaluates the panel''',
    r'''The driver
\texttt{debug/l3b\_2\_sub\_sprint\_D\_compute.py} (project repository)
evaluates the panel''', 'p45-driver')

rep(P45, r'''(reproduced from \texttt{debug/data/l3b\_2\_sub\_sprint\_D.json};\ see''',
    r'''(data and driver in the project repository,
\texttt{debug/data/l3b\_2\_sub\_sprint\_D.json};\ see''', 'p45-caption')

rep(P45, r'''reproduced from \texttt{debug/data/l3b\_2\_sub\_sprint\_D.json}.
The corresponding computational driver is''',
    r'''reproduced from the project repository
(\texttt{debug/data/l3b\_2\_sub\_sprint\_D.json});\ the corresponding
computational driver is''', 'p45-appB')

print('FAILS:', fails if fails else 'none')
