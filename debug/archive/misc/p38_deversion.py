"""One-shot de-versioning edits for Paper 38.

Run from repo root: python debug/p38_deversion.py
Each replacement asserts uniqueness; misses reported, not fatal.
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


# 1. Title
rep(r'''\title{State-space Gromov--Hausdorff convergence of \\
truncated Camporesi--Higuchi spectral triples on $\sthree$\\
\large (v2.0 erratum of ``Latr\'emoli\`ere propinquity convergence of
truncated Camporesi--Higuchi spectral triples on $\sthree$'')}''',
    r'''\title{State-space Gromov--Hausdorff convergence of \\
truncated Camporesi--Higuchi spectral triples on $\sthree$}''', 'title')

# 2. Date
rep(r'\date{May 2026 (v1 preprint);\ v2.0 erratum:\ June 9, 2026}',
    r'\date{June 2026 (preprint)}', 'date')

# 3. Abstract parenthetical
rep(r'''We prove — subject to two named gaps recorded in the v2.0 erratum
(2026-06-09) — that the Connes--van~Suijlekom truncations of the
round $\sthree = \SU(2)$ Camporesi--Higuchi spectral triple converge
to the continuum spectral triple in van~Suijlekom's state-space
Gromov--Hausdorff distance~\cite{vs2021_jgp}, with explicit asymptotic
rate
\(
\Lambda(\Tcal_{\nmax}, \Tcal_{\sthree}) = O\!\bigl(\log \nmax / \nmax\bigr)
\)
and asymptotic constant $4/\pi$. (Version~1 stated the result in the
Latr\'emoli\`ere propinquity;\ the erratum withdraws that framework
attribution, corrects the v1 Lemma~L2(b)--(c) — the value
$2/(\nmax+1)$ is the kernel's Plancherel mass-distribution maximum,
not a cb-norm;\ the smoothing map is unital completely positive with
cb-norm $1$ — and records the two gaps:\ the dual-direction reach
estimate, and the kernel condition $L(a) = 0 \Leftrightarrow a \in
\C 1$, which fails for the truthful chirality-diagonal Dirac and is
verified numerically only for the engineered off-diagonal Dirac.)
The proof proceeds via five lemmas:''',
    r'''We prove — subject to two named gaps (\S~\ref{sec:named_gaps}) —
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
The proof proceeds via five lemmas:''', 'abstract')

# 4. Erratum section -> named-gaps pointer; delete section, content moves to sec:named_gaps
rep(r'''% =====================================================================
\section*{Erratum and corrections (v2.0, 2026-06-09)}
% =====================================================================

An adversarial verification pass (corpus hardening sprint, 2026-06-09;\
canonical record \texttt{debug/sprint\_p45\_hardening\_phase1\_memo.md}
and \texttt{debug/p45\_adversarial\_L2\_memo.md} in the project
repository) found that version~1 of this paper contained one wrong
lemma, one wrong framework attribution, and two unproved steps. The
main convergence statement survives in corrected, conditional form.

\begin{enumerate}
\item \textbf{Lemma~L2(b)--(c) corrected.} Version~1 asserted that the
Fej\'er smoothing map acts with symbol $(2j+1)/Z_{\nmax}$ and has
cb-norm $2/(\nmax+1)$, citing Bo\.zejko--Fendler. The quantity
$(2j+1)/Z$ is the kernel's Plancherel \emph{mass distribution};\ the
true multiplier symbol $\sigma(j')$ equals $1$ at the trivial
representation (forced by the normalization asserted in L2(a)),
satisfies $0 \le \sigma \le 1$, and is supported on $j' \le
2j_{\max}$ — at $\nmax = 2$ the values are $(1, \sqrt2/3, 2/9)$,
refuting the v1 claim $(1/3, 2/3, 0)$ bit-exactly. The smoothing map
is unital completely positive with cb-norm exactly $1$;\ the
Bo\.zejko--Fendler citation (a discrete-groups theorem) is withdrawn.
The convergence-controlling moment $\gamma_{\nmax}$ (L2(d):\ closed
form, $4/\pi$ asymptote, uniform bound) is a property of the kernel
itself and is unaffected.
\item \textbf{Framework attribution corrected.} The distance actually
bounded by the $(B_{\nmax}, P_{\nmax})$ approximation pair is
van~Suijlekom's state-space Gromov--Hausdorff
distance~\cite{vs2021_jgp} (the device of
Leimbach--van~Suijlekom~\cite{leimbach_vs2024}), not Latr\'emoli\`ere's
propinquity:\ Latr\'emoli\`ere's framework defines the distance via
tunnels carried by quantum metric spaces with quantum isometries and
imposes hypotheses (Leibniz Lip-norms on C$^{*}$-algebras;\ kernel
condition) that are not established for the truncated operator
system. The title, theorem statements, and L5 are restated
accordingly.
\item \textbf{Named gap 1:\ dual-direction reach.} The v1
reach$_{P}$ argument inferred a bounded partial inverse of
$B_{\nmax}$ from the v1 ``cb-norm'' $2/(\nmax+1)$;\ a forward bound
cannot control a partial inverse, whose norm is governed by the
symbol \emph{minimum} on the relevant window. The corrected route is
the antiderivative-transference mechanism of
Leimbach--van~Suijlekom~\cite{leimbach_vs2024} (Lemmas~3.4--3.5
there), whose conclusion for scalar sectors on compact metric groups
is covered by Gaudillot-Estrada--van~Suijlekom~\cite{gaudillot_vs2023};\
its adaptation to the chirality-doubled spinor substrate is open.
\item \textbf{Named gap 2:\ kernel condition and the Dirac choice.}
For the truthful chirality-diagonal Camporesi--Higuchi Dirac, the
Lipschitz seminorm $\opnorm{[\DCH, \cdot]}$ on the truncated
multiplier system has kernel far exceeding the scalars (verified
bit-exactly:\ $10/14$ multipliers at $\nmax = 2$, $26/55$ at
$\nmax = 3$), so the truncated object is not a compact quantum metric
space and its state-space metric is degenerate. The kernel condition
holds for the engineered off-diagonal Dirac of the project's R3.5
construction (kernel $= \C 1$ exactly:\ $1/14$, $1/55$), at the cost
that the off-diagonal operator is a hand-parameterized stand-in for
the genuine spinor-bundle off-diagonal structure, and L3's constant
must be re-established for it. The v1 text bridged this by an
unproved ``bounded perturbation robustness'' assertion, now
withdrawn;\ the theorem below is stated as conditional on this gap.
\item \textbf{Novelty repositioning.} State-space GH convergence of
scalar spectral truncations for all compact metric groups (including
$\SU(2)$) is published in
Gaudillot-Estrada--van~Suijlekom~\cite{gaudillot_vs2023}. The
residual contribution of this paper is the spinor substrate, the
closed-form rate moment $\gamma_{\nmax}$, and the $4/\pi$ constant —
not the existence of convergence for $\SU(2)$ scalars.
\end{enumerate}

\noindent Downstream:\ Papers~39/40 inherit item~1 (Paper~40 carries
its own erratum remark);\ Papers~45--49 are descoped separately
(Paper~45 v2, whose Theorem~3.2 [$K^{+}$ annihilation] also explains
why the Lorentzian extension of the present construction failed). The
frozen falsifier for item~4 is
\texttt{tests/test\_p45\_kplus\_degeneracy.py}.

% =====================================================================
\section{Introduction}\label{sec:intro}''',
    r'''% =====================================================================
\section{Introduction}\label{sec:intro}''', 'erratum-section')

# 5. thm:main_intro
rep(r'''\begin{theorem}[Main theorem, conditional;\ cf.~\S\ref{sec:main_theorem}
and the Erratum]
\label{thm:main_intro}
Subject to the two named gaps of the v2.0 erratum (dual-direction
reach;\ kernel condition on the off-diagonal Dirac substrate), the
truncated triples $\Tcal_{\nmax}$ converge to $\Tcal_{\sthree}$ in
van~Suijlekom's state-space Gromov--Hausdorff
distance~\cite{vs2021_jgp}, with explicit rate''',
    r'''\begin{theorem}[Main theorem, conditional;\ cf.~\S\ref{sec:main_theorem}]
\label{thm:main_intro}
Subject to the two named gaps of \S~\ref{sec:named_gaps}
(dual-direction reach;\ kernel condition on the off-diagonal Dirac
substrate), the truncated triples $\Tcal_{\nmax}$ converge to
$\Tcal_{\sthree}$ in van~Suijlekom's state-space Gromov--Hausdorff
distance~\cite{vs2021_jgp}, with explicit rate''', 'thm-intro')

# 6. Intro CvS/GE-vS paragraph
rep(r'''explicit convergence to the continuum in a quantitative quantum
Gromov--Hausdorff distance as an open question, deferring it to
``elsewhere.'' The convergence itself, at the state-space level, is
established in van~Suijlekom~\cite{vs2021_jgp} for the circle and
sphere, and — published before v1 of the present paper, identified in
the 2026-06-09 audit — for \emph{all compact metric groups} by
Gaudillot-Estrada--van~Suijlekom~\cite{gaudillot_vs2023};\ see Erratum
item~5 for the resulting novelty repositioning.''',
    r'''explicit convergence to the continuum in a quantitative quantum
Gromov--Hausdorff distance as an open question, deferring it to
``elsewhere.'' The convergence itself, at the state-space level, is
established in van~Suijlekom~\cite{vs2021_jgp} for the circle and
sphere, and for \emph{all compact metric groups} by
Gaudillot-Estrada--van~Suijlekom~\cite{gaudillot_vs2023};\ the
residual contribution of the present paper is therefore the spinor
substrate, the closed-form rate moment, and the explicit $4/\pi$
constant — not the existence of convergence for $\SU(2)$ scalars.''',
    'intro-gevs')

# 7. "treats that case"
rep(r'''The present paper treats that case, subject to the erratum's two
named gaps.''',
    r'''The present paper treats that case, subject to the two named gaps
of \S~\ref{sec:named_gaps}.''', 'treats-case')

# 8. Insert named-gaps subsection before Setup
rep(r'''% =====================================================================
\section{Setup}\label{sec:setup}''',
    r'''\subsection{Named gaps and history}\label{sec:named_gaps}

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
metric is degenerate. The kernel condition holds for the engineered
off-diagonal Dirac of the R3.5 construction (kernel $= \C 1$
exactly:\ $1/14$, $1/55$), at the cost that the off-diagonal operator
is a hand-parameterized stand-in for the genuine spinor-bundle
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
\end{remark}

% =====================================================================
\section{Setup}\label{sec:setup}''', 'named-gaps')

# 9. L2(b) tag
rep(r'''\item[(b)] (Mass distribution and multiplier symbol;\ \emph{corrected
2026-06-09}.)''',
    r'''\item[(b)] (Mass distribution and multiplier symbol.)''', 'L2b-tag')

# 10. L2(b) trailing parenthetical
rep(r''' (The v1
statement assigned the mass distribution as the symbol;\ it is
refuted at the trivial representation by (a) itself.)''', '', 'L2b-paren')

# 11. L2(c) tag + parenthetical
rep(r'''\item[(c)] (Cb-norm;\ \emph{corrected 2026-06-09}.) $T_{\KS_{\nmax}}$
is unital completely positive (convolution against a probability
density), hence $\Vert T_{\KS_{\nmax}}\Vert_{\mathrm{cb}} = 1$. (The
v1 value $2/(\nmax+1)$ is the mass-distribution maximum of (b), not
the cb-norm of any map used in this paper.)''',
    r'''\item[(c)] (Cb-norm.) $T_{\KS_{\nmax}}$
is unital completely positive (convolution against a probability
density), hence $\Vert T_{\KS_{\nmax}}\Vert_{\mathrm{cb}} = 1$. (The
mass-distribution maximum $2/(\nmax+1)$ of (b) is \emph{not} the
cb-norm of any map used in this paper.)''', 'L2c')

# 12. L2 proof BF sentence
rep(r'''on the truncated side the induced central multiplier
acts block-scalar on $\bigoplus_{\pi}\Mat_{d_{\pi}}(\C)$ with cb-norm
$\sup_{\pi}|\sigma(\pi)| = 1$ — elementary, requiring neither
amenability nor Bo{\.z}ejko--Fendler~\cite{bozejko_fendler1991}
(whose discrete-groups theorem the v1 proof cited inapplicably).''',
    r'''on the truncated side the induced central multiplier
acts block-scalar on $\bigoplus_{\pi}\Mat_{d_{\pi}}(\C)$ with cb-norm
$\sup_{\pi}|\sigma(\pi)| = 1$ — an elementary statement requiring no
amenability input.''', 'L2-proof-bf')

# 13. L2 proof verification sentence
rep(r'''unitality at $j' = 0$ is
$c_{0} = \sum_{j}(2j+1) = Z$, and the $\nmax = 2$ values were
verified bit-exactly against independent Weyl-integration quadrature
(2026-06-09 audit).''',
    r'''unitality at $j' = 0$ is
$c_{0} = \sum_{j}(2j+1) = Z$, and the $\nmax = 2$ values are verified
bit-exactly against independent Weyl-integration quadrature
(repository:\ \texttt{debug/p45\_adversarial\_L2\_memo.md}).''',
    'L2-proof-verif')

# 14. def:berezin
rep(r'''\begin{definition}[Berezin map, convolution form;\ \emph{corrected
2026-06-09}]\label{def:berezin}''',
    r'''\begin{definition}[Berezin map, convolution form]\label{def:berezin}''',
    'berezin-tag')

rep(r''' (The v1 definition used the
mass-distribution weights $N/Z_{\nmax}$ on $N \le \nmax$;\ with those
weights the map is not unital — $B(\mathbf{1}) = Z^{-1} I$ — and the
approximate-identity property fails maximally at $f = \mathbf{1}$.
The convolution form is the definition all proofs below actually
use;\ note its image reaches multiplier labels up to the achievable
envelope $N \le 2\nmax - 1$, not $N \le \nmax$ as v1 implied.)''',
    r''' (A sum form weighted by the mass distribution $N/Z_{\nmax}$ in
place of $\sigma_{\nmax}(N)$ would not be unital —
$B(\mathbf{1}) = Z^{-1} I$ — and would fail the approximate-identity
property maximally at $f = \mathbf{1}$. Note the image reaches
multiplier labels up to the achievable envelope $N \le 2\nmax - 1$.)''',
    'berezin-paren')

# 15. L5 preamble
rep(r'''\textbf{(v2 framework correction.)} Version~1 framed this assembly in
Latr\'emoli\`ere's quantum Gromov--Hausdorff propinquity for metric
spectral triples~\cite{latremoliere2018, latremoliere_metric_st_2017},
treating $(B_{\nmax}, P_{\nmax})$ as a ``direct tunnel.'' In
Latr\'emoli\`ere's actual framework the infimum runs over tunnels
carried by quantum metric spaces with quantum isometries, under
hypotheses (Leibniz Lip-norms on C$^{*}$-algebras;\ kernel condition)
not established for the truncated operator system;\ a pair of UCP
maps is not a tunnel there. The distance the pair certifies is
van~Suijlekom's state-space Gromov--Hausdorff
distance~\cite{vs2021_jgp}, via two-directional approximation data in
the Leimbach--van~Suijlekom pattern~\cite{leimbach_vs2024}. We write''',
    r'''A remark on framework before assembling. In Latr\'emoli\`ere's
quantum Gromov--Hausdorff propinquity for metric spectral
triples~\cite{latremoliere2018, latremoliere_metric_st_2017}, the
distance is an infimum over tunnels carried by quantum metric spaces
with quantum isometries, under hypotheses (Leibniz Lip-norms on
C$^{*}$-algebras;\ kernel condition) not established for the
truncated operator system;\ a pair of UCP maps is not a tunnel there.
The distance the pair $(B_{\nmax}, P_{\nmax})$ certifies is
van~Suijlekom's state-space Gromov--Hausdorff
distance~\cite{vs2021_jgp}, via two-directional approximation data in
the Leimbach--van~Suijlekom pattern~\cite{leimbach_vs2024}. We write''',
    'L5-preamble')

rep(r'''reach$_{B}$, reach$_{P}$, height$_{B}$, height$_{P}$ for the four
approximation estimates the pair supplies. For the truncated-side
Monge--Kantorovich metric to be non-degenerate the kernel condition
of Erratum item~4 must hold, which conditions the statement on the
off-diagonal Dirac substrate.''',
    r'''reach$_{B}$, reach$_{P}$, height$_{B}$, height$_{P}$ for the four
approximation estimates the pair supplies. For the truncated-side
Monge--Kantorovich metric to be non-degenerate the kernel condition
of \S~\ref{sec:named_gaps}(G2) must hold, which conditions the
statement on the off-diagonal Dirac substrate.''', 'L5-kernel-ref')

# 16. lem:L5 statement
rep(r'''\begin{lemma}[L5:\ distance bound from the approximation pair,
conditional]\label{lem:L5}
Subject to the two named gaps of the v2.0 erratum (the dual-direction
reach estimate of item~3;\ the kernel condition of item~4, conditioning
the statement on the off-diagonal Dirac substrate), the pair''',
    r'''\begin{lemma}[L5:\ distance bound from the approximation pair,
conditional]\label{lem:L5}
Subject to the two named gaps of \S~\ref{sec:named_gaps}
(the dual-direction reach estimate G1;\ the kernel condition G2,
conditioning the statement on the off-diagonal Dirac substrate), the
pair''', 'lem-L5')

# 17. reach_P paragraph
rep(r'''\textbf{(v2:\ named gap.)} The v1 text asserted the dual estimate
$\mathrm{reach}_P \le \gamma_{\nmax}$ ``via the L2(c) cb-norm
equality'' — an invalid inference (Erratum item~3:\ a forward cb-norm
bound cannot control a partial inverse, whose norm is set by the
symbol minimum;\ on the full support window the inverse norm grows
like $O(\nmax^{2})$). The corrected route is the
antiderivative-transference mechanism of
Leimbach--van~Suijlekom~\cite{leimbach_vs2024} (Lemmas~3.4--3.5),
whose conclusion for scalar sectors on compact metric groups is
covered by Gaudillot-Estrada--van~Suijlekom~\cite{gaudillot_vs2023};\
its adaptation to the chirality-doubled spinor substrate is the open
named gap on which this lemma is conditional. If the adapted constant
exceeds $1$, the right-hand side of (\ref{eq:L5_bound}) degrades by
that $O(1)$ factor;\ the rate class $O(\log\nmax/\nmax)$ is
unaffected.''',
    r'''\textbf{(Named gap G1.)} The dual estimate
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
unaffected.''', 'reachP')

# 18. thm:main
rep(r'''\begin{theorem}[Main theorem;\ rate-controlled, conditional on the
erratum's two named gaps]\label{thm:main}
Subject to the dual-direction reach estimate (Erratum item~3) and the
kernel condition on the off-diagonal Dirac substrate (Erratum item~4),
the truncated Camporesi--Higuchi spectral triples $\Tcal_{\nmax}$''',
    r'''\begin{theorem}[Main theorem;\ rate-controlled, conditional on the
two named gaps of \S~\ref{sec:named_gaps}]\label{thm:main}
Subject to the dual-direction reach estimate (G1) and the kernel
condition on the off-diagonal Dirac substrate (G2),
the truncated Camporesi--Higuchi spectral triples $\Tcal_{\nmax}$''',
    'thm-main')

# 19. why_su2 note
rep(r'''Definition~\ref{def:central_fejer}. (\emph{v2 note:}\ for central
multipliers the cb-norm computation is in fact elementary —
block-scalar action — on any compact group;\ the v1 appeal to
Bo{\.z}ejko--Fendler~\cite{bozejko_fendler1991} is withdrawn per the
Erratum, and the genuine non-abelian difficulty is the
Clebsch--Gordan structure of the multiplier symbol, not the cb-norm.)''',
    r'''Definition~\ref{def:central_fejer}. (For central multipliers the
cb-norm computation is elementary — block-scalar action — on any
compact group;\ the genuine non-abelian difficulty is the
Clebsch--Gordan structure of the multiplier symbol, not the cb-norm.)''',
    'why-su2')

# 20. Remove now-orphan Bozejko-Fendler bibitem
rep(r'''\bibitem[Bo{\.z}ejko \& Fendler(1991)]{bozejko_fendler1991}
M.~Bo{\.z}ejko and G.~Fendler,
\newblock \emph{Herz--Schur multipliers and uniformly bounded
representations of discrete groups},
\newblock Arch.~Math. \textbf{57} (1991), 290--298.

''', '', 'bf-bibitem')

io.open(P, 'w', encoding='utf-8').write(s)
print('FAILS:', fails if fails else 'none')
