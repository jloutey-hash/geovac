"""One-shot de-versioning edits for Paper 45, part 2 (section 6 onward).

Run from repo root: python debug/p45_deversion_part2.py
Each replacement asserts uniqueness; misses are reported, not fatal.
"""
import io

P = r'papers/group1_operator_algebras/paper_45_lorentzian_propinquity.tex'
s = io.open(P, encoding='utf-8').read()
fails = []


def rep(old, new, tag):
    global s
    if s.count(old) == 1:
        s = s.replace(old, new)
    else:
        fails.append((tag, s.count(old)))


rep(r'''\begin{remark}[Erratum record for v1 Lemma~L2]\label{rem:cb_Nt_indep}
Version~1 asserted $\cbnorm{S_{\Kjoint}} = 2/(\nmax+1)$, citing a
Bo\.zejko--Fendler central-multiplier
equality~\cite{bozejko_fendler1991} (a theorem about Herz--Schur
multipliers and uniformly bounded representations of \emph{discrete}
groups, inapplicable here) and a ``Pisier~Ch.~8 Thm.~8.10'' that does
not state the claimed equality. The value $2/(\nmax+1)$ is the
maximum of the SU(2) Plancherel \emph{mass distribution}
(Lemma~\ref{lem:L2}(c) without the $\Uone$ factor) — a correct
supremum of a quantity that is not the cb-norm of any map used in the
proofs;\ the v1 ``$\Nt$-independence of the cb-norm'' arose only from
mixing the mass convention on the SU(2) factor with the symbol
convention on the $\Uone$ factor. The correct cb-norm statement,
(\ref{eq:L2_main}), is elementary and requires neither amenability
nor Bo\.zejko--Fendler. Two downstream consequences:\ (i) the v1
``contractivity'' item L4(b) is corrected to $\cbnorm{\Bjoint} \le 1$;\
(ii) the v1 reach$_{P}$ argument, which inferred a bounded partial
inverse of $\Bjoint$ from the smallness of the v1 ``cb-norm,'' is
invalid — a forward bound cannot control a partial inverse, whose
norm is governed by the symbol \emph{minimum} on the relevant window
(see \S~\ref{sec:withdrawn_assembly}).
\end{remark}''',
    r'''\begin{remark}[The mass maximum is not a cb-norm]\label{rem:cb_Nt_indep}
The mass-distribution maximum $2/(\nmax+1)$ of Lemma~\ref{lem:L2}(c)
is not the cb-norm of any map appearing in this paper;\ the smoothing
map is unital completely positive with cb-norm $1$
(Lemma~\ref{lem:L2}(b)), an elementary statement requiring no
amenability input. Conflating the two objects has two consequences
worth flagging:\ (i) it would assign the UCP map $\Bjoint$ a
``cb-norm'' below $1$, which is impossible;\ (ii) it would suggest a
bounded partial inverse of $\Bjoint$ ``at the same scale'' — invalid,
since an inverse is governed by the symbol \emph{minimum} on the
relevant window (see \S~\ref{sec:withdrawn_assembly}).
\end{remark}''', 'cb-remark')

rep(r'''\begin{remark}[v2 correction:\ the v1 sum form is not equivalent]
\label{rem:berezin_inequivalence}
Version~1 defined $\Bjoint$ by a sum form,
$\sum_{N,L,M,q} m^{\SU(2)}(N)\, m \text{-or-} \sigma^{\Uone}(q)\,
c_{N,L,M,q}\,(M^{\spat}_{N,L,M} \otimes M^{\temp}_{q})$ with
$M^{\temp}_{q}$ a momentum-\emph{diagonal} matrix, and asserted
equivalence with (\ref{eq:joint_berezin}). The two forms are
inequivalent on the temporal factor:\ under the convolution form, the
image of the Fourier mode $e^{iqt}$ is
$\sigma^{\Uone}(q)$ times the \emph{truncated shift} (Toeplitz)
operator on $\C^{\Nt}$, which is not momentum-diagonal for $q \ne 0$;\
the v1 sum form instead assigned it a momentum-diagonal image, which
is not the compression of any multiplication operator. Two
consequences. (i) The weights:\ the v1 sum form used the SU(2)
\emph{mass distribution} $m^{\SU(2)}(N)$ where the convolution form
acts by the multiplier symbol $\sigma^{\SU(2)}(N)$;\ with the mass
weights the sum form is not unital
($\Bjoint(\mathbf{1}) = Z^{-1} I \ne I$), maximally violating the
approximate-identity property at $f = \mathbf{1}$. (ii) The target:\
the convolution form (\ref{eq:joint_berezin}) — the only coherent
definition, and the one all proofs below actually use — has temporal
image in the Toeplitz compression algebra, which is \emph{not
contained} in the operator system $\Op^{L}$ of
\S~\ref{sec:op_system}, whose temporal factor is momentum-diagonal.
This target mismatch is part of the degeneracy diagnosis of
\S~\ref{sec:main}:\ the v1 operator system is not the algebra the
Berezin map naturally reconstructs, and the momentum-diagonal choice
is precisely what makes the temporal factor invisible to the Dirac
commutator (Lemma~\ref{lem:L3}). Rebuilding the construction with the
Toeplitz temporal algebra is the reopened target of
\S~\ref{sec:open}.
\end{remark}''',
    r'''\begin{remark}[Sum-form pitfalls:\ weights and targets]
\label{rem:berezin_inequivalence}
Two natural-looking sum-form alternatives to
(\ref{eq:joint_berezin}) fail, and both failures are instructive.
(i) \emph{Weights:}\ a sum form weighted by the mass distribution
$m^{\SU(2)}(N)$ in place of the multiplier symbol
$\sigma^{\SU(2)}(N)$ is not unital
($\Bjoint(\mathbf{1}) = Z^{-1} I \ne I$) and maximally violates the
approximate-identity property at $f = \mathbf{1}$.
(ii) \emph{Targets:}\ under the convolution form, the image of the
Fourier mode $e^{iqt}$ is $\sigma^{\Uone}(q)$ times the
\emph{truncated shift} (Toeplitz) operator on $\C^{\Nt}$ — not
momentum-diagonal for $q \ne 0$ — so $\Bjoint$'s temporal image is
\emph{not contained} in the operator system $\Op^{L}$ of
\S~\ref{sec:op_system}, whose temporal factor is momentum-diagonal;\
and a sum form assigning $e^{iqt}$ a momentum-diagonal image is not
the compression of any multiplication operator. This target mismatch
is part of the degeneracy diagnosis of \S~\ref{sec:main}:\ the
momentum-diagonal operator system is not the algebra the Berezin map
naturally reconstructs, and that choice is precisely what makes the
temporal factor invisible to the Dirac commutator
(Lemma~\ref{lem:L3}). Rebuilding the construction with the Toeplitz
temporal algebra is the reopened target of \S~\ref{sec:open}.
\end{remark}''', 'berezin-remark')

rep(r''' \emph{(v1 claimed $\cbnorm{\Bjoint} \le
2/(\nmax+1)$;\ corrected per Lemma~\ref{lem:L2} and
Remark~\ref{rem:cb_Nt_indep} — a UCP map cannot have cb-norm below
$1$.)}''', '', 'L4b-paren')

rep(r'''\emph{(v1 asserted this for all $f$;\ for temporally-varying $f$ the
convolution-form image is Toeplitz on the temporal factor, where
$[\gamma^{0} \otimes \partial_{t}, \cdot]$ does not vanish and the
required estimate is precisely the Connes--van~Suijlekom
$\Sone$-type Toeplitz commutator bound~\cite{connes_vs2021} — part of
the reopened target, \S~\ref{sec:open}.)}''',
    r'''\emph{(For temporally-varying $f$ the convolution-form image is
Toeplitz on the temporal factor, where
$[\gamma^{0} \otimes \partial_{t}, \cdot]$ does not vanish;\ the
required estimate is precisely the Connes--van~Suijlekom
$\Sone$-type Toeplitz commutator bound~\cite{connes_vs2021} — part of
the reopened target, \S~\ref{sec:open}.)}''', 'L4d-paren')

rep(r'\section{The annihilation theorem and the withdrawn v1 assembly}\label{sec:main}',
    r'\section{The annihilation theorem and the failure of the compression assembly}\label{sec:main}',
    'sec6-header')

rep(r'''level, so the $\Kplus$-compression of the pair is well-defined.
Version~1 asserted that the compressed pair is ``a valid
Latr\'emoli\`ere tunneling pair'' between the compressed triples of
Definition~\ref{def:weak_form_propinquity}. That assertion fails twice
over:\ the cited framework contains no such device (Erratum item~2),
and — more fundamentally — the compressed triples are not quantum
metric spaces at all, as we now prove.''',
    r'''level, so the $\Kplus$-compression of the pair is well-defined. One
might hope the compressed pair to be a tunneling pair between the
compressed triples of Definition~\ref{def:weak_form_propinquity} in
some Riemannian framework. It is not, twice over:\ no applicable
framework admits a UCP-pair tunnel on these objects, and — more
fundamentally — the compressed triples are not quantum metric spaces
at all, as we now prove.''', 'postdef')

rep(r'''\subsection{The withdrawn v1 assembly}\label{sec:withdrawn_assembly}

For the record, and so that citations of the v1 quantities resolve to
an honest account, we summarize what version~1 asserted and the status
of each piece. Version~1 claimed the bound''',
    r'''\subsection{Anatomy of the failed assembly}\label{sec:withdrawn_assembly}

We record where each piece of the natural assembly stands — both
because the surviving estimates are genuine analytic facts about the
kernels, and because the failure points specify the repair
(\S~\ref{sec:open}). The assembly one would attempt is the bound''',
    'anatomy-open')

rep(r'''for the ``direct UCP tunneling pair'' $(\Bjoint, \Pjoint)$, attributing
(\ref{eq:propinquity_bound}) and the four constituent definitions to
``Thm~5.5'' and ``Def.~3.4 / Def.~3.5'' of
Latr\'emoli\`ere~\cite{latremoliere_metric_st_2017}. \emph{No such
numbered results exist in that paper}, whose framework (tunnels
carried by quantum metric spaces with quantum isometries) contains no
UCP-pair device;\ the nearest genuine framework is van~Suijlekom's
state-space Gromov--Hausdorff bound via approximation maps in both
directions~\cite{vs2021_jgp,leimbach_vs2024}. Independently of the
attribution, Theorem~\ref{thm:kplus_annihilation} shows the left-hand
side of (\ref{eq:propinquity_bound}) is not a distance, so the v1
``Proposition (reach and height bounds)''\label{prop:reach_height}
and the v1 main theorem are withdrawn. Piece by piece, measured
against the \emph{classical} joint-gradient unit ball (an auxiliary
third seminorm, not the Lipschitz ball of either triple — itself a
structural defect of the v1 bookkeeping):''',
    r'''for the UCP pair $(\Bjoint, \Pjoint)$, with the four constituents
measured against the \emph{classical} joint-gradient unit ball (an
auxiliary third seminorm — not the Lipschitz ball of either triple,
itself a structural defect of any such bookkeeping). No
four-constituent UCP-pair tunnel exists in Latr\'emoli\`ere's
framework~\cite{latremoliere_metric_st_2017} (tunnels there are
carried by quantum metric spaces with quantum isometries);\ the
nearest genuine framework is van~Suijlekom's state-space
Gromov--Hausdorff bound via approximation maps in both
directions~\cite{vs2021_jgp,leimbach_vs2024}. Independently of the
framework question, Theorem~\ref{thm:kplus_annihilation} shows the
left-hand side of (\ref{eq:propinquity_bound}) is not a distance, so
no assembly through Definition~\ref{def:weak_form_propinquity} can
succeed.\label{prop:reach_height} Piece by piece:''', 'anatomy-body')

rep(r'''\item \emph{reach$_{P}$ (dual roundtrip).} Withdrawn. The v1 argument
inferred a bounded partial inverse of $\Bjoint$ from the v1
``cb-norm'' $2/(\nmax+1)$;\ a forward bound cannot control a partial
inverse, whose norm is governed by the multiplier-symbol
\emph{minimum} on the relevant window (at full support the inverse
norm grows like $O(\nmax^{2})$). A corrected dual-direction estimate
is the antiderivative-transference mechanism of
Leimbach--van~Suijlekom~\cite{leimbach_vs2024} (Lemmas~3.4--3.5
there);\ for scalar sectors on compact metric groups the conclusion is
covered by Gaudillot-Estrada--van~Suijlekom~\cite{gaudillot_vs2023}.
Adapting it to the chirality-doubled spinor substrate is a named gap
inherited by Paper~38's erratum.''',
    r'''\item \emph{reach$_{P}$ (dual roundtrip).} No bounded partial
inverse of $\Bjoint$ can be inferred from any forward bound:\ an
inverse is governed by the multiplier-symbol \emph{minimum} on the
relevant window (at full support the inverse norm grows like
$O(\nmax^{2})$). The correct dual-direction mechanism is the
antiderivative transference of
Leimbach--van~Suijlekom~\cite{leimbach_vs2024} (Lemmas~3.4--3.5
there);\ for scalar sectors on compact metric groups the conclusion
is covered by Gaudillot-Estrada--van~Suijlekom~\cite{gaudillot_vs2023}.
Adapting it to the chirality-doubled spinor substrate is a named gap
shared with the companion $\sthree$ paper.''', 'reachP')

rep(r'''\item \emph{height$_{B}$ (Lipschitz distortion).} Withdrawn — false as
stated. Counterexample:\ for purely temporal $f = \mathbf{1} \otimes
f_{t}$ with $\norm{\partial_{t} f_{t}}_{\infty} = 1$, the classical
side gives $\norm{f}_{\mathrm{Lip}} = 1$ while
$\opnorm{[\DL, \Bjoint(f)]}$ vanishes on the momentum-diagonal
reading (and the claimed bound $\Cthreejoint\gammajoint \to 0$), so
the distortion is $\ge 1$ and does not tend to $0$. The spatial
restriction of the statement survives via Lemma~\ref{lem:L4}(d).''',
    r'''\item \emph{height$_{B}$ (Lipschitz distortion).} Fails:\ for
purely temporal $f = \mathbf{1} \otimes f_{t}$ with
$\norm{\partial_{t} f_{t}}_{\infty} = 1$, the classical side gives
$\norm{f}_{\mathrm{Lip}} = 1$ while $\opnorm{[\DL, \Bjoint(f)]}$
vanishes on the momentum-diagonal reading, so the distortion is
$\ge 1$ and cannot tend to $0$. The spatial restriction survives via
Lemma~\ref{lem:L4}(d).''', 'heightB')

rep(r'''\begin{theorem}[v1 main theorem — \textbf{withdrawn}]\label{thm:main}
\emph{(For the record;\ see Erratum.)} Version~1 asserted:\ for all
$\nmax \ge 1$, $\Nt \ge 1$, $T > 0$,
\begin{equation}\label{eq:main_bound}
   \Lprop\!\bigl(\Tcal^{L}_{\nmax,\Nt,T},\;
                 \Tcal^{L}_{\Manifold}\bigr)
   \;\le\; \Cthreejoint(\nmax, \Nt) \cdot \gammajoint_{\nmax, \Nt, T}
   \;\xrightarrow[(\nmax, \Nt)\to(\infty, \infty)]{}\; 0.
\end{equation}
This statement is withdrawn:\ its left-hand side is not a quantum
Gromov--Hausdorff-type distance
(Theorem~\ref{thm:kplus_annihilation}(iii)), its proof imported
nonexistent results (Erratum item~2), and one of its four constituent
bounds is false (height$_{B}$ above). The right-hand side is a
well-defined closed-form rate expression whose surviving meaning is
given by Proposition~\ref{prop:spatial_conditional}.
\end{theorem}''',
    r'''The assembly therefore produces no distance bound. The rate
expression it would have produced,
\begin{equation}\label{eq:main_bound}
   \Cthreejoint(\nmax, \Nt) \cdot \gammajoint_{\nmax, \Nt, T}
   \;=\; O\!\Bigl(\tfrac{\log\nmax}{\nmax} + \tfrac{T}{\Nt}\Bigr)
   \;\xrightarrow[(\nmax, \Nt)\to(\infty, \infty)]{}\; 0,
\end{equation}
is nonetheless a well-defined closed-form expression;\ its surviving
meaning is the conditional spatial statement of
Proposition~\ref{prop:spatial_conditional}, to which we now turn.''',
    'thm-main-delete')

rep(r'''The per-harmonic functional
form, the asymptotic limit $1^{-}$, and the rate expressions of the
withdrawn v1 statement (Theorem~\ref{thm:main}) are insensitive to
this distinction.''',
    r'''The per-harmonic functional
form, the asymptotic limit $1^{-}$, and the rate expressions of
\S~\ref{sec:withdrawn_assembly} (eq.~\ref{eq:main_bound}) are
insensitive to this distinction.''', 'envelope-ref')

rep(r'''\begin{remark}[No limit identification]\label{rem:limit_identification}
Version~1 included a ``limit identification'' remark asserting, via
the nonexistent ``Latr\'emoli\`ere Thm~5.5,'' that
$\Lambdab$-convergence identifies the limit on the
Wasserstein--Kantorovich state space of
$\Tcal^{+}_{\Manifold} = (\Pplus C^{\infty}(\Manifold)\Pplus,
\Kplus_{\infty}, \Pplus D^{L}_{\infty} \Pplus)$. That remark is
withdrawn on both grounds:\ the citation is empty, and the candidate
limit object is itself degenerate — the compressed continuum Dirac
loses its spatial part by the same anticommutation mechanism as
Theorem~\ref{thm:kplus_annihilation}(i), so its Lipschitz seminorm
kernel contains all time-independent functions and its
Monge--Kantorovich structure does not separate states.
\end{remark}''',
    r'''\begin{remark}[The compressed continuum object is also degenerate]
\label{rem:limit_identification}
There is no limit identification to state:\ the candidate limit
object $\Tcal^{+}_{\Manifold} = (\Pplus C^{\infty}(\Manifold)\Pplus,
\Kplus_{\infty}, \Pplus D^{L}_{\infty} \Pplus)$ is itself degenerate
— the compressed continuum Dirac loses its spatial part by the same
anticommutation mechanism as
Theorem~\ref{thm:kplus_annihilation}(i), so its Lipschitz seminorm
kernel contains all time-independent functions and its
Monge--Kantorovich structure does not separate states.
\end{remark}''', 'limit-id')

rep(r'''\begin{remark}[The v1 ``bookkeeping'' framing, inverted]
\label{rem:bookkeeping}
Version~1 described the assembly as ``the bookkeeping layer of the
proof shape,'' with ``no new analytical input.'' The 2026-06-09 audit
inverted this assessment:\ the assembly layer is precisely where the
construction failed — a degenerate seminorm, a nonexistent imported
theorem, and a false constituent bound — while the kernel-level
analytic facts (Lemma~\ref{lem:L4}(a)--(c),(e)) survive. The lesson
recorded for the project's verification protocol:\ ``bookkeeping''
steps that import external theorems are exactly the steps that
require line-level hypothesis checking.
\end{remark}''',
    r'''\begin{remark}[Assembly layers are not mere bookkeeping]
\label{rem:bookkeeping}
In proof shapes of this kind it is tempting to regard the assembly
layer as bookkeeping, with the analytical content concentrated in the
kernel-level lemmas. The present construction inverts that
assessment:\ the kernel-level facts
(Lemma~\ref{lem:L4}(a)--(c),(e)) are robust, while the assembly layer
is where the construction fails — a degenerate seminorm, no
admissible tunnel device, and a false constituent bound. Assembly
steps that import external theorems are exactly the steps that
require line-level hypothesis checking.
\end{remark}''', 'bookkeeping')

rep(r'''Version~1 composed (\ref{eq:main_bound}) with the norm-resolvent
convergence $\Tcal^{L}_{\Manifold} \to
\Tcal^{L}_{\sthree\times\R_{t}}$ of Paper~47~\cite{paper47}~\S~4 into
a ``two-rate hybrid convergence'' to the non-compact carrier. The
outer (norm-resolvent) arrow is operator-level, does not use the
Lipschitz seminorm, and survives this erratum with rate
$O(e^{-|\mathrm{Im}\,z|\,T/2})$ on compact-temporal-support
subspaces. The inner (metric) arrow is withdrawn together with
(\ref{eq:main_bound});\ the metric-level statement on both the
compact and non-compact carriers reverts to open (Paper~47 carries
its own erratum notice).''',
    r'''A two-rate hybrid composite to the non-compact carrier
$\sthree\times\R_{t}$ has two arrows. The outer (norm-resolvent)
arrow of Paper~47~\cite{paper47}~\S~4 is operator-level, does not use
the Lipschitz seminorm, and is unaffected by the degeneracy, with
rate $O(e^{-|\mathrm{Im}\,z|\,T/2})$ on compact-temporal-support
subspaces. The inner (metric) arrow would require a Q1-repaired
construction;\ the metric-level statement on both the compact and
non-compact carriers is open (Paper~47 carries its own Status note).''',
    'two-rate')

rep(r'''\textbf{(v2 reinterpretation.)} The driver
\texttt{debug/l3b\_2\_sub\_sprint\_D\_compute.py} evaluates the panel''',
    r'''The driver
\texttt{debug/l3b\_2\_sub\_sprint\_D\_compute.py} evaluates the panel''',
    'panel-intro')

rep(r'''at $T = 2\pi$. The 2026-06-09 audit established that the tabulated
``$\Lprop$ bound'' values are \emph{evaluations of the closed-form
rate expressions}''',
    r'''at $T = 2\pi$. The tabulated
``$\Lprop$ bound'' values are \emph{evaluations of the closed-form
rate expressions}''', 'panel-audit')

rep(r'''(iv) the v1 constituent
inequalities measured against the auxiliary classical-gradient ball.''',
    r'''(iv) the assembly constituent
inequalities of \S~\ref{sec:withdrawn_assembly} against the auxiliary
classical-gradient ball.''', 'panel-iv')

rep(r'''\caption{Rate-formula panel (v2 reinterpretation;\ v1 presented this
as the ``numerical verification panel for the main theorem'').''',
    r'''\caption{Rate-formula panel.''', 'caption-head')

rep(r'''mass-distribution
maximum}, not a cb-norm (Erratum item~3;\ the cb-norm of the smoothing
map is $1$ at every cell).''',
    r'''mass-distribution
maximum}, not a cb-norm (the cb-norm of the smoothing map is $1$ at
every cell;\ Lemma~\ref{lem:L2}).''', 'caption-cb')

rep(r'''\emph{(v1 stated this as equality of propinquities;\ neither side is
an evaluated metric — see \S~\ref{sec:panel}.)}''',
    r'''\emph{(Neither side is an evaluated metric — see
\S~\ref{sec:panel}.)}''', 'riemlimit-paren')

rep(r'''\textbf{(v2 correction on falsifier strength.)} Version~1 presented
this reduction as ``the load-bearing falsifier of the construction.''
The 2026-06-09 audit showed it is far weaker than advertised:\ because
the panel quantity only ever contained the spatial rate formula
(\S~\ref{sec:panel}), the $\Nt = 1$ reduction is an automatic formula
identity — it is \emph{insensitive} to the degeneracy of the metric
construction, which it failed to detect. A falsifier with actual teeth
for this construction is the operator-level seminorm check of
Corollary~\ref{cor:annihilation_numerics}, which detects the failure
immediately. The corresponding $\Nt = 1$ reductions in Paper~42
(Tomita--Takesaki construction) and Paper~43 (Lorentzian Dirac) are
genuinely operator-level and are not affected by this correction.''',
    r'''\textbf{On falsifier strength.} Because the panel quantity only
ever contains the spatial rate formula (\S~\ref{sec:panel}), this
$\Nt = 1$ reduction is an automatic formula identity — it is
\emph{insensitive} to the degeneracy of the metric construction,
which it cannot detect. A falsifier with actual teeth for this
construction is the operator-level seminorm check of
Corollary~\ref{cor:annihilation_numerics}, which detects the failure
immediately. The corresponding $\Nt = 1$ reductions in Paper~42
(Tomita--Takesaki construction) and Paper~43 (Lorentzian Dirac) are
genuinely operator-level.''', 'falsifier-strength')

rep(r'''The rate-formula ratio $1.3223 / 2.0746 = 0.6374$ is monotone
decreasing across the panel and matches the single-factor Paper~38 and
Paper~39~\cite{paper39} ratios bit-identically. Version~1 read this
match as confirmation that ``the SU(2) factor dominates the joint
rate.'' The v2 reading is sharper and less flattering:\ the match is
bit-identical \emph{because the tabulated quantity is the spatial
formula} — the temporal factor never enters the maximum at these
cells, and (by Theorem~\ref{thm:kplus_annihilation}) never enters the
operator-level seminorm at all.''',
    r'''The rate-formula ratio $1.3223 / 2.0746 = 0.6374$ is monotone
decreasing across the panel and matches the single-factor Paper~38
and Paper~39~\cite{paper39} ratios bit-identically — \emph{because
the tabulated quantity is the spatial formula}:\ the temporal factor
never enters the maximum at these cells, and (by
Theorem~\ref{thm:kplus_annihilation}) never enters the operator-level
seminorm at all.''', 'asymp-consistency')

rep(r'''Version~1 presented this paper as the Lorentzian \emph{convergence}
leg of that arc, with papers (42, 43, 45) jointly constituting ``the
Lorentzian-extension closure'' of the Riemannian quartet (38, 39, 40).
After the erratum, the honest statement is:\ Papers~42/43 close
modular-structure statements at the operator level (unaffected here),
while the convergence leg is \emph{not closed} — the present paper now
proves that the v1 device for it cannot work
(Theorem~\ref{thm:kplus_annihilation}) and reopens the problem
(\S~\ref{sec:open}).''',
    r'''Papers~42/43 close modular-structure statements at the operator
level;\ those are unaffected by the degeneracy. The convergence leg
of the arc, by contrast, is \emph{not closed}:\ the present paper
proves that the natural device for it cannot work
(Theorem~\ref{thm:kplus_annihilation}) and states the open problem
(\S~\ref{sec:open}).''', 'four-witness')

rep(r'''The present paper works on
truncated Krein spectral triples — and, after the erratum, proves a
\emph{negative} structural result about the obvious operator-algebraic
convergence device on them, rather than a convergence theorem.''',
    r'''The present paper works on
truncated Krein spectral triples — and proves a \emph{negative}
structural result about the obvious operator-algebraic convergence
device on them, rather than a convergence theorem.''', 'ms-remark')

rep(r'''The 2026-06-09 erratum converts
v1's ``weak-form vs strong-form'' split into a single reopened
problem, now with a sharper specification supplied by the failure
analysis.''',
    r'''The degeneracy theorem collapses the former ``weak-form vs
strong-form'' split into a single open problem, with a sharper
specification supplied by the failure analysis.''', 'Q1')

rep(r'''On the v1 substrate the
``temporal contribution to the Lipschitz seminorm'' is identically
zero for \emph{all} multipliers (not only pure tensors), so questions
about joint sharpness are moot until Q1 supplies a construction in
which the temporal direction carries metric weight.''',
    r'''On the momentum-diagonal substrate the temporal contribution to
the Lipschitz seminorm is identically zero for \emph{all} multipliers
(not only pure tensors), so questions about joint sharpness are moot
until Q1 supplies a construction in which the temporal direction
carries metric weight.''', 'Q4')

rep(r'''The v2 content obeys the structural-skeleton scope discipline of the
GeoVac framework in its negative form:''',
    r'''The results here obey the structural-skeleton scope discipline of
the GeoVac framework in its negative form:''', 'scope')

rep(r'''\textbf{(v2 rewrite.)} The v1 version of this appendix asserted
$\int_{\SU(2)} \KSU\,\chi_{j'}\,dg = (2j'+1)/Z$ for $j' \le j_{\max}$
— a false identity (already at $\nmax = 2$, $j' = \tfrac12$:\ the
true value is $2\sqrt2/3 \ne 2/3$) that conflated the kernel's
Plancherel mass distribution with its character coefficients. We
record the corrected objects.''',
    r'''We record the kernel's two Plancherel-side objects. A caution
worth stating explicitly:\ the identity
$\int_{\SU(2)} \KSU\,\chi_{j'}\,dg = (2j'+1)/Z$ — which would make
the character coefficients equal the mass distribution — is
\emph{false} (already at $\nmax = 2$, $j' = \tfrac12$ the left side
is $2\sqrt2/3 \ne 2/3$);\ the two objects are related but distinct,
as follows.''', 'appA-intro')

rep(r'\begin{lemma}[Mass distribution and multiplier symbol, corrected]',
    r'\begin{lemma}[Mass distribution and multiplier symbol]', 'appA-title')

rep(r'''The $\nmax = 2$
values are direct, and were verified bit-exactly against independent
Weyl-integration quadrature in the 2026-06-09 audit
(\texttt{debug/p45\_adversarial\_L2\_memo.md}).''',
    r'''The $\nmax = 2$
values are direct, and are verified bit-exactly against independent
Weyl-integration quadrature in the repository
(\texttt{debug/p45\_adversarial\_L2\_memo.md}).''', 'appA-proof')

rep(r'''(15+ test functions) verifies the formula-level content:\ the
tunneling-pair construction, the v1 constituent inequalities measured
against the auxiliary classical-gradient ball
(\S~\ref{sec:withdrawn_assembly}), the rate-formula panel, the''',
    r'''(15+ test functions) verifies the formula-level content:\ the
tunneling-pair construction, the assembly constituent inequalities of
\S~\ref{sec:withdrawn_assembly} (classical-gradient ball), the
rate-formula panel, the''', 'appB')

rep(r'''Fix $\Nt \ge 1$, $T > 0$. Conditional on Paper~38's spatial
Gromov--Hausdorff convergence statement as repaired under its own
erratum (van~Suijlekom state-space distance~\cite{vs2021_jgp};\
corrected Lemma~L2;\ the dual-direction reach estimate supplied by the
antiderivative-transference mechanism~\cite{leimbach_vs2024,gaudillot_vs2023};\
kernel condition holding on the engineered off-diagonal Dirac
substrate), the spatial truncations underlying''',
    r'''Fix $\Nt \ge 1$, $T > 0$. Conditional on the companion $\sthree$
paper's spatial Gromov--Hausdorff convergence statement
(van~Suijlekom state-space distance~\cite{vs2021_jgp};\ its two named
gaps:\ the dual-direction reach estimate via antiderivative
transference~\cite{leimbach_vs2024,gaudillot_vs2023}, and the kernel
condition on the engineered off-diagonal Dirac substrate), the
spatial truncations underlying''', 'prop-conditional')

rep(r'''all estimates are the Paper~38
estimates verbatim, subject to that paper's erratum repairs
(which is the conditionality).''',
    r'''all estimates are the Paper~38
estimates verbatim, subject to that paper's two named gaps (which is
the conditionality).''', 'prop-proof')

io.open(P, 'w', encoding='utf-8').write(s)
print('FAILS:', fails if fails else 'none')
