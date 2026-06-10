"""Descope corrections to the Group 1 synthesis (Phase 2e follow-through)."""
import io

P = r'papers/synthesis/group1_operator_algebras_synthesis.tex'
s = io.open(P, encoding='utf-8').read()
fails = []


def rep(old, new, tag):
    global s
    if s.count(old) == 1:
        s = s.replace(old, new)
    else:
        fails.append((tag, s.count(old)))


# 1. Status note after abstract
rep(r'''\end{abstract}''',
    r'''\end{abstract}

\bigskip\noindent\textbf{Status note (June 2026).}  An adversarial
verification pass retracted the Lorentzian convergence claims of this
arc:\ the K$^{+}$-compression device underlying Papers~45--46 and the
bridge constructions of Papers~48--49 is degenerate (the compressed
Lipschitz seminorm vanishes identically;\ Paper~45's annihilation
theorem), and Paper~38's Riemannian convergence theorem is now stated
conditionally with two named gaps.  This synthesis is corrected at
the highest-visibility points;\ where narrative detail below still
reflects the pre-retraction state, the Status notes in Papers~38 and
42--49 and the claims register
(\texttt{docs/claims\_register.md}) are authoritative.''', 'status-note')

# 2. Intro "Fourth, ..." passage
rep(r'''Fourth, a $K^{+}$-restricted
weak-form Lorentzian propinquity convergence theorem holds on
truncated $\SU(2)\times\Uone_{t}$ Krein spectral triples (Paper~45) --
to the author's knowledge the first such convergence result in the
math.OA literature, closing a question deferred by Connes and van
Suijlekom \cite{connes_vs2021}.  The convergence has since been
strengthened in two directions:\ Paper~46 closes the strong-form
Lorentzian propinquity convergence (without $K^{+}$ restriction) on
the natural and enlarged operator-system substrates, and Paper~47
closes the de-compactification limit $T \to \infty$ at the
norm-resolvent / spectral level via a two-rate hybrid composite.''',
    r'''Fourth, the natural route to a Lorentzian convergence theorem on
truncated $\SU(2)\times\Uone_{t}$ Krein spectral triples ---
compressing to the Krein-positive subspace --- is closed by a
degeneracy theorem (Paper~45):\ the compression annihilates the
spatial Dirac and the restricted Lipschitz seminorm vanishes
identically.  An earlier convergence claim through this device,
together with the strong-form and bridge extensions built on it
(Papers~46, 48, 49), was retracted in June 2026;\ the question
deferred by Connes and van Suijlekom \cite{connes_vs2021} remains
open, now with a sharpened specification.  Paper~47's
de-compactification limit $T \to \infty$ survives at the
norm-resolvent / spectral level.''', 'intro-fourth')

# 3. First-claim paragraph
rep(r'''This is, to the author's knowledge, the first
Lorentzian propinquity convergence theorem in the math.OA literature
(Paper~45~\cite{paper45}, Theorem~5.1).  Connes and van~Suijlekom
\cite{connes_vs2021} deferred this question three times to
``elsewhere''; the present arc closes it at the K$^{+}$-weak-form
level on truncated Krein spectral triples.''',
    r'''The corresponding convergence claim (formerly Theorem~5.1 of
Paper~45) was retracted in June 2026:\ the quantity $\Lprop$ is
degenerate --- the K$^{+}$ compression annihilates the spatial Dirac
and the restricted Lipschitz seminorm vanishes identically
(Paper~45's annihilation theorem~\cite{paper45}).  The question
Connes and van~Suijlekom \cite{connes_vs2021} deferred to
``elsewhere'' remains open, now with a sharpened specification.''',
    'first-claim')

# 4. "first Lorentzian instance"
rep(r'''The Paper~38
$\SU(2)$ result and the Paper~40 generalisation are the first
non-abelian non-flat instances at quantitative rate; the Paper~45
result is the first Lorentzian instance.''',
    r'''The Paper~38
$\SU(2)$ result and the Paper~40 generalisation are non-abelian
non-flat instances at quantitative rate (conditional;\ scalar-sector
convergence for compact groups is due to Gaudillot-Estrada and
van~Suijlekom, arXiv:2310.14733);\ on the Lorentzian side the
Paper~45 contribution is a degeneracy theorem, not a convergence
instance.''', 'first-instance')

# 5. Theorem block replacement
rep(r'''\begin{theorem}[Paper~45, Theorem~5.1]
\label{thm:p45_main}
For the truncated $\SU(2) \times \Uone_{T}$ Krein spectral triples
$\Tcal^{L}_{\nmax,\Nt,T}$, with $\DL$ as in
\eqref{eq:lorentzian_dirac} and $\Krein_{\nmax,\Nt} =
\HGV^{\nmax} \otimes \C^{\Nt}$,
\begin{equation}
\Lprop(\Tcal^{L}_{\nmax,\Nt,T}, \Tcal^{L}_{\Manifold})
\;\le\;
C_{3}^{\mathrm{joint}} \cdot
\gamma^{\mathrm{joint}}_{\nmax,\Nt,T}
\;\longrightarrow\; 0
\quad \text{as } (\nmax, \Nt) \to (\infty, \infty),
\label{eq:p45_main}
\end{equation}
with $C_{3}^{\mathrm{joint}} \le 1$ asymptotically tight (inherited
from Paper~38 L3 via the joint Lichnerowicz structural identity) and
\begin{equation}
\gamma^{\mathrm{joint}}_{\nmax,\Nt,T}
\;=\;
O\!\bigl(\log \nmax / \nmax \;+\; T / \Nt\bigr).
\label{eq:p45_rate}
\end{equation}
\end{theorem}''',
    r'''\begin{theorem}[Paper~45, annihilation theorem]
\label{thm:p45_main}
The device of Definition~\ref{def:k_plus_prop} is degenerate.  For
$\DL$ as in \eqref{eq:lorentzian_dirac}, Krein-self-adjointness of
the anti-Hermitian spatial block forces
$\{\JL, \DGV \otimes I\} = 0$, hence
$P^{+}(\DGV \otimes I)P^{+} = 0$:\ the compression annihilates the
spatial Dirac.  The surviving temporal block is momentum-diagonal and
commutes with every momentum-diagonal temporal multiplier, so the
restricted Lipschitz seminorm vanishes identically on the operator
system and $\Lprop$ of \eqref{eq:k_plus_prop_def} is not a quantum
Gromov--Hausdorff-type distance.  (Bit-exact:\ restricted-seminorm
kernel $42/42$ and $275/275$ multipliers at
$(\nmax, \Nt) = (2,3), (3,5)$;\ frozen falsifier in the repository
test suite.)
\end{theorem}

An earlier revision of the arc claimed a convergence theorem through
this device (formerly Theorem~5.1 of Paper~45);\ the claim was
retracted in June 2026.  The rate expression
\begin{equation}
C_{3}^{\mathrm{joint}} \cdot \gamma^{\mathrm{joint}}_{\nmax,\Nt,T}
\;=\;
O\!\bigl(\log \nmax / \nmax \;+\; T / \Nt\bigr)
\label{eq:p45_rate}
\end{equation}
survives as a well-defined closed form whose meaning is a
\emph{conditional spatial} statement (Paper~38, two named gaps), the
temporal factor being a metrically invisible spectator.''',
    'p45-theorem')

# 6. Simplifications intro
rep(r'''The proof of Theorem~\ref{thm:p45_main} transfers the Paper~38
five-lemma chain to the joint $\SU(2) \times \Uone$ setting.  Three
structural simplifications arise that are not visible in the
Riemannian setting, and are recorded as substantive new content of
the Paper~45 proof:''',
    r'''The original transfer of the Paper~38 five-lemma chain to the
joint $\SU(2) \times \Uone$ setting recorded three structural
``simplifications.''  Under the degeneracy analysis each is
reinterpreted:\ they are \emph{symptoms} of the temporal factor's
metric invisibility, not conveniences:''', 'simplifications-intro')

# 7. cross-term paragraph: append diagnosis sentence
rep(r'''The temporal direction contributes nothing to the joint
commutator, and $C_{3}^{\mathrm{joint}} = C_{3}^{\SU(2)}$ verbatim.''',
    r'''The temporal direction contributes nothing to the joint
commutator, and $C_{3}^{\mathrm{joint}} = C_{3}^{\SU(2)}$ verbatim.
Read correctly, \eqref{eq:vanishing_cross_term} is the
\emph{diagnosis}:\ the momentum-diagonal temporal algebra is
invisible to the Dirac commutator.''', 'cross-term')

# 8. cb-norm paragraph
rep(r'''\paragraph{$\Nt$-independent joint cb-norm.}
By the Bo\.zejko--Fendler central-multiplier equality
\cite{bozejko_fendler1991} on the amenable group product
$\SU(2) \times \Uone$, the joint cb-norm of the central spectral
Fej\'er multiplier is
\begin{equation}
\cbnorm{S_{K^{\mathrm{joint}}}}
\;=\;
\frac{2}{\nmax + 1},
\label{eq:joint_cb_norm}
\end{equation}
\emph{independent} of $\Nt$.  This reflects the trivial
central-multiplier structure of the compact abelian factor $\Uone$,
and means that the $\Nt$-dependence in
\eqref{eq:p45_rate} lives entirely in the joint mass-concentration
moment (the $T/\Nt$ Fej\'er-on-the-circle rate), not in the
cb-norm height bound.''',
    r'''\paragraph{Corrected cb-norm convention.}  The smoothing map of
the joint central spectral Fej\'er kernel is unital completely
positive, with cb-norm exactly $1$;\ the quantity $2/(\nmax+1)$
formerly quoted as the ``joint cb-norm'' is the maximum of the SU(2)
Plancherel \emph{mass distribution}, a different invariant (the
Bo\.zejko--Fendler appeal is withdrawn;\ the correct statement is
elementary).  The $\Nt$-dependence of the rate
\eqref{eq:p45_rate} lives entirely in the joint mass-concentration
moment (the $T/\Nt$ Fej\'er-on-the-circle rate).''', 'cb-norm')

# 9. Riemannian-limit paragraph: append caveat
rep(r'''agrees with the Riemannian Paper~38 bound bit-exactly, providing the
load-bearing falsifier check for the joint construction.''',
    r'''agrees with the Riemannian Paper~38 bound bit-exactly.  Under the
degeneracy analysis this recovery is an automatic formula identity
(the panel quantity only ever contained the spatial rate formula),
not an operator-level falsifier;\ the falsifier with teeth is the
restricted-seminorm check of Theorem~\ref{thm:p45_main}.''',
    'riemannian-limit')

# 10. Honest scope opening
rep(r'''Theorem~\ref{thm:p45_main} is, to the author's knowledge, the first
published Lorentzian propinquity convergence theorem in the math.OA
literature.  Several aspects of its scope require explicit comment.''',
    r'''Theorem~\ref{thm:p45_main} closes the natural device;\ the
Lorentzian convergence question itself remains open.  Several aspects
of scope require explicit comment.''', 'honest-scope')

# 11. Q1 item correction
rep(r'''The K$^{+}$-restriction reduces the indefinite Krein product to a
Hilbert inner product on $\Kplus$ and lets the standard Latr\'emoli\`ere
machinery apply verbatim.''',
    r'''The K$^{+}$-restriction reduces the indefinite Krein product to a
Hilbert inner product on $\Kplus$ --- but annihilates the spatial
Dirac (Theorem~\ref{thm:p45_main}), so no quantum-metric machinery
applies on the compression.''', 'q1-item')

# 12-13. P48/P49 subsection flags
rep(r'''\label{sec:closures_paper48}

Paper~48~\cite{paper48} constructs''',
    r'''\label{sec:closures_paper48}

\emph{Status (June 2026):\ the metric-level content of this
subsection is descoped with the Paper~45 retraction;\ the bridge
design survives as a conditional proposal.  See the Status notes in
Papers~48--49 and \texttt{docs/claims\_register.md}.}

Paper~48~\cite{paper48} constructs''', 'p48-flag')

rep(r'''\label{sec:closures_paper49}

Paper~49~\cite{paper49} closes''',
    r'''\label{sec:closures_paper49}

\emph{Status (June 2026):\ the $\Lambda$-inheritance and
K$^{+}$-framed claims of this subsection are descoped with the
Paper~45 retraction;\ the cocycle-deficit algebra and the OSLPLS
category design survive.  See Paper~49's Status note and
\texttt{docs/claims\_register.md}.}

Paper~49~\cite{paper49} closes''', 'p49-flag')

io.open(P, 'w', encoding='utf-8').write(s)
print('FAILS:', fails if fails else 'none')
