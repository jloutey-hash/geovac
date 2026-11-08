# Paper 30 Wilson merge: material cut from Paper 30

Verbatim passages removed from `papers/group5_qed_gauge/paper_30_su2_wilson.tex`
when Papers 25 and 41 were merged into it (2026-11-07, PI-approved). Line numbers
refer to the pre-merge file (which already carried the PM's same-day cubic-term fix).
Papers 25 and 41 are archived whole and are not excerpted here. Several passages
below were rewritten rather than dropped; the new wording is in the paper.

## Lines 26-86: Former abstract (five action principles, S^3 = SU(2) framing, ST-SU3 sibling paragraph)

```latex
\begin{abstract}
Paper~25~\cite{paper25} placed the GeoVac Hopf graph inside the
Wilson--Hodge lattice-gauge vocabulary with a $U(1)$ connection on
edges, and explicitly left as open (\S VII.A) the question of whether
a non-abelian gauge structure is natural on the same graph.  This paper
answers that question for the compact simply-connected non-abelian Lie
group of rank~1, $\mathrm{SU}(2)$.  The coincidence $S^3 =
\mathrm{SU}(2)$ (as manifolds) makes this case structurally special:
the underlying manifold of the GeoVac graph \emph{is} the gauge group.
Taking link variables $U_e \in \mathrm{SU}(2)$, plaquettes as the
primitive closed non-backtracking walks of the Ihara zeta
setup~\cite{paper29}, and the Wilson action
$S_W[U] = \beta \sum_P \bigl(1 - \tfrac{1}{2}\mathrm{Re}\,\mathrm{tr}\,
U_P\bigr)$, we obtain a gauge-invariant finite-dimensional lattice
gauge theory with Haar-normalizable partition function.
Three structural results follow.
(1)~\emph{$U(1)$ reduction to Paper 25.}  The diagonal (maximal-torus)
sector of $\mathrm{SU}(2)$ link variables reproduces Paper~25's
$U(1)$ Wilson--Hodge action to machine precision.  Paper~25 is the
Cartan-subalgebra projection of the present theory; the abelian
construction it documents is a structural sub-sector of a genuinely
non-abelian one.
(2)~\emph{The weak-coupling kinetic term completes the Hodge-1
Laplacian.}  Linearizing $S_W$ around the trivial vacuum $U_e = I$
produces, at second order in the gauge field, the plaquette (curl)
quadratic form $(\beta/8) \sum_a \langle A^a, B_2 B_2^{\top} A^a\rangle$
on $\mathfrak{su}(2)$-valued $1$-cochains, with per-plaquette
per-component coefficient $1/8 = 1/(4 N_c)$.  This curl form is
\emph{not} a scalar multiple of Paper~25's edge Laplacian
$L_1 = B^{\top} B$ (the two have complementary supports: the plaquette
columns of $B_2$ span $\ker L_1$ exactly); it is the co-exact part that completes
$L_1$ to the full discrete Hodge-1 Laplacian
$B^{\top} B + B_2 B_2^{\top}$ of the 2-complex with the Ihara
plaquettes attached as faces, exactly the face structure
Paper~25~\S II anticipates.
(3)~\emph{Plaquette enumeration and Monte Carlo Wilson-loop data.}
On the $S^3$ Coulomb graph at $n_{\max} = 2, 3, 4$, the first-Betti
counts are $\beta_1 = 0, 2, 8$ and the primitive 4-cycle counts are
$0, 2, 8$; at $n_{\max} = 3$ the two plaquettes form a cycle basis, so
the plaquette variables are independent Haar elements and the Wilson
expectation is exactly $\langle W \rangle = I_2(\beta)/I_1(\beta)$
($0.124, 0.240, 0.433, 0.719$ at $\beta = 0.5, 1, 2, 5$), which long
Metropolis chains reproduce.
The scope is a compact-graph non-abelian lattice gauge construction,
not a Yang--Mills mass-gap claim on $\mathbb{R}^4$; the paper closes
with a synthesis observation that the GeoVac framework now carries
\emph{five} action principles---maximum entropy on nodes, Rayleigh--Ritz
on node-wavefunctions, Wilson action on gauge edges, Ihara path-integral
on walks, and spectral action on Dirac spectrum---suggesting---though we do not claim---a framework-wide
organizing principle of least action.
An $\mathrm{SU}(3)$ sibling on the Bargmann--Segal $S^5$ graph
(Paper~24) has been constructed in Sprint~ST-SU3
(\texttt{debug/st\_su3\_wilson\_memo.md}, May~2026), confirming
Proposition~1's analog (Cartan reduction to $U(1)\times U(1)$) and
Proposition~2's analog (kinetic coefficient $1/12 = 1/(4 N_c)$) at machine
precision; matter coupling on $S^5$ encounters the Sprint~5 Track~S5
Clebsch--Gordan-intertwiner obstruction, so that of Paper~24's
Coulomb/HO asymmetry layers, Wilson gauge alone is universal while
Wilson gauge with natural matter coupling is Coulomb-specific (the
spectrum-computing and calibration-$\pi$ layers are withdrawn).
\end{abstract}
```

## Lines 112-118: Introduction: S^3 = SU(2) framing paragraph

```latex
The $S^3$ case is categorically different.  As manifolds,
$S^3 = \mathrm{SU}(2)$: the underlying manifold of the GeoVac graph
\emph{is} the gauge group.  Taking $\mathrm{SU}(2)$-valued link
variables is therefore not a choice made by matching representation
content across shells; it is the natural object carried by the
manifold itself.  This is the cleanest geometric setting in which a
non-abelian Wilson lattice gauge theory can be built on a GeoVac graph.
```

## Lines 120-142: Introduction: 'What is new' list ((i) Ihara-walk plaquettes, (ii) Cartan = Paper 25, (iii) Hodge completion)

```latex
This paper reports the construction.  The mathematics is
standard~\cite{wilson1974,creutz1983}\footnote{Concurrent work in
the Bonn/DESY Hamiltonian-lattice program develops SU(2) Hamiltonian
lattice gauge theory using discrete partitionings of the gauge group,
with motivation toward quantum simulation: see
Garofalo et al.~\cite{garofalo2023_su2_s3} for the SU(2) Hamiltonian
formulation based on $S_3$ partitionings (here $S_3$ denotes the
symmetric group, distinct from the manifold $S^3$ used here), and
Jakobs et al.~\cite{jakobs2025_su2_dynamics} for dynamical
continuum-limit studies.}: link variables
$U_e \in \mathrm{SU}(2)$, plaquette holonomies $U_P = \prod_{e \in P}
U_e$, Wilson action
$S_W = \beta \sum_P (1 - \tfrac{1}{2} \mathrm{Re}\,\mathrm{tr}\,U_P)$,
node-local gauge transformations $U_e \mapsto g_{s(e)}^{\phantom{\dagger}}
U_e\, g_{t(e)}^{\dagger}$, and Haar measure on each link.  What is new
is (i)~the choice of plaquettes as the primitive closed non-backtracking
walks of Paper~29~\cite{paper29}'s Ihara zeta machinery, rather than the
elementary four-cycles of a hypercubic lattice, (ii)~the confirmation
that the maximal-torus reduction is \emph{exactly} Paper~25's $U(1)$
construction, and (iii)~the identification of the weak-coupling kinetic term as the
plaquette (curl) quadratic form that completes Paper~25's edge
Laplacian $L_1$ to the full discrete Hodge-1 Laplacian of the
plaquette-filled 2-complex.
```

## Lines 144-160: Former Scope paragraph

```latex
\paragraph*{Scope.}  This is an observation paper.  What it establishes
is a self-consistent non-abelian lattice gauge structure on the finite
compact GeoVac Hopf graph, together with three structural results.
What it does \emph{not} establish is a continuum-limit Yang--Mills mass
gap on $\mathbb{R}^4$, confinement (at the sizes tested the plaquettes
form a cycle basis, so Wilson loops factorize into independent
plaquette factors and the area/perimeter-law diagnostic carries no
dynamical information), or any derivation of $\alpha$ (Paper~2 is an Observation).
The rhetoric follows the CLAUDE.md~\S1.5 rule: the mathematical
construction is sharp; the physical interpretation is modest.  The
construction is a non-abelian instance of the Marcolli--van Suijlekom
gauge-network framework~\cite{marcolli_vs2014}, which defines finite
spectral triples on graphs whose spectral action is a Wilson lattice
gauge theory; Perez-Sanchez~\cite{perez_sanchez2025}
establishes that the continuum limit in that lineage is Yang--Mills
without Higgs, consistent with the scope statements of
Section~\ref{sec:scope}.
```

## Lines 162-176: Former Organization paragraph

```latex
\paragraph*{Organization.}  Section~\ref{sec:recap} recaps the
Paper~25 Wilson--Hodge vocabulary in the form needed.
Section~\ref{sec:su2} states the $\mathrm{SU}(2)$ link variables,
plaquette structure, Wilson action, and gauge transformation.
Section~\ref{sec:characters} gives the character expansion of the
partition function on $\mathrm{SU}(2)$.  Section~\ref{sec:results}
reports the three structural results: $U(1)$-reduction exactness,
the kinetic-term / Hodge-completion result, and plaquette/Wilson-loop
data.
Section~\ref{sec:action_principles} places the construction inside a
framework-wide synthesis: the GeoVac project now exhibits five
independent action principles, each organizing a different sector of
its data.  Section~\ref{sec:scope} states explicitly what is and is
not claimed (Yang--Mills adjacency).  Section~\ref{sec:open} lists
open questions, and Section~\ref{sec:conclusion} concludes.
```

## Lines 211-249: Plaquettes as primitive closed non-backtracking walks; former Table 1 (N_4, N_6, N_8 columns)

```latex
The orienting convention for plaquettes differs here from the ``elementary
4-cycle'' convention of a hypercubic lattice.  We take plaquettes to be
the \emph{primitive closed non-backtracking walks} of the graph, i.e.\
the Ihara cycles of Paper~29~\cite{paper29}.  A primitive closed
non-backtracking walk on $L$ edges is a sequence $e_1 e_2 \cdots e_L$
with $\mathrm{head}(e_i) = \mathrm{tail}(e_{i+1})$, no $e_{i+1}$ equal
to the reverse of $e_i$ (including across the closure $e_L \to e_1$),
and not a cyclic rotation of a shorter closed walk.  This is the
natural Wilson generalization on a non-hypercubic graph: primitive
closed non-backtracking walks are exactly the irreducible ``boundary
of a disk'' objects whose holonomies enter the Wilson action.

Table~\ref{tab:plaquettes} collects the plaquette counts on the
$S^3$ Coulomb graph at $n_{\max} = 2, 3, 4$.  At $n_{\max} = 2$
the graph is a forest ($\beta_1 = 0$), so there are no plaquettes
and the $\mathrm{SU}(2)$ Wilson theory is trivial ($Z(\beta) = 1$).
At $n_{\max} = 3$ the theory is non-trivial for the first time: two
primitive 4-cycles in the $p$-block, one primitive 6-cycle, one
primitive 8-cycle.  At $n_{\max} = 4$ the graph has
$\beta_1 = 8$ and a richer plaquette structure.

\begin{table}[t]
\centering
\caption{Plaquette enumeration on the $S^3$ Coulomb graph.  $V, E, c$
are vertex, edge, and component counts; $\beta_1 = E - V + c$ is the
first Betti number; $N_L$ is the number of primitive closed
non-backtracking walks of length $L$ (one representative per
unoriented class).}
\label{tab:plaquettes}
\begin{tabular}{rrrrrrrr}
\toprule
$n_{\max}$ & $V$ & $E$ & $c$ & $\beta_1$ & $N_4$ & $N_6$ & $N_8$ \\
\midrule
2 & 5 & 3 & 2 & 0 & 0 & 0 & 0 \\
3 & 14 & 13 & 3 & 2 & 2 & 1 & 1 \\
4 & 30 & 34 & 4 & 8 & 8 & 7 & 18 \\
\bottomrule
\end{tabular}
\end{table}
```

## Lines 260-265: Shared-edge paragraph after Eq. (p_block_plaquettes)

```latex
These share the undirected edge $\{(2, 1, 0), (3, 1, 0)\}$, but each
also contains an edge the other lacks, and with $\beta_1 = 2$ they form
a cycle basis.  Gauge-fixing along a spanning tree therefore makes the
two plaquette holonomies independent Haar-distributed variables: the
shared edge introduces no correlation, and the plaquette-independent
character expansion of Section~\ref{sec:characters} is exact here.
```

## Lines 398-431: Character-coefficient table at n_max = 3 (Table 2) and its paragraph

```latex
\paragraph*{Coefficient table at $n_{\max} = 3$.}  For the two primitive
4-cycles of Eq.~(\ref{eq:p_block_plaquettes}), the exact
partition function (the plaquettes are independent)
$Z_0(\beta) = e^{-\beta p}\bigl[d_0(\beta)\bigr]^p$ with $p = 2$
gives, together with the first excited coefficient $d_1(\beta) =
(4/\beta) I_2(\beta)$,

\begin{table}[h]
\centering
\caption{Character coefficients and the exact partition function
on the $n_{\max} = 3$ Hopf graph (the two plaquettes are independent).
Numerical data:
\texttt{debug/data/su2\_wilson\_gauge\_results.json}.}
\label{tab:char_coeffs}
\begin{tabular}{rrrrr}
\toprule
$\beta$ & $Z_0(\beta)$ & $d_0(\beta)$ & $d_1(\beta)$ & $d_1/d_0$ \\
\midrule
 0.1 & $8.21\times 10^{-1}$ & $1.0013$ & $0.0500$ & $0.050$ \\
 0.5 & $3.91\times 10^{-1}$ & $1.0316$ & $0.2552$ & $0.247$ \\
 1.0 & $1.73\times 10^{-1}$ & $1.1303$ & $0.5430$ & $0.480$ \\
 2.0 & $4.63\times 10^{-2}$ & $1.5906$ & $1.3779$ & $0.866$ \\
 5.0 & $4.30\times 10^{-3}$ & $9.7343$ & $14.0045$ & $1.439$ \\
10.0 & $5.88\times 10^{-4}$ & $534.20$ & $912.61$  & $1.708$ \\
\bottomrule
\end{tabular}
\end{table}

The entries of Table~\ref{tab:char_coeffs} set the scale at which
higher representations matter in a general partition function.  For
the single-plaquette Wilson loop on this graph they do not enter:
because the two plaquettes of Eq.~(\ref{eq:p_block_plaquettes}) are independent,
$\langle W \rangle = I_2(\beta)/I_1(\beta)$ exactly at every $\beta$
(Section~\ref{sec:results_MC}).
```

## Lines 481-501: U(1)-reduction verification paragraph and Observation 'Paper 25 is the Cartan-subalgebra sector'

```latex
The proposition is verified numerically in
\texttt{tests/test\_su2\_wilson\_gauge.py}, class
\texttt{TestU1Reduction}: for a triangle plaquette with random
diagonal phases $\{\theta_e\} = \{0.3, 0.7, -0.5\}$ and $\beta = 1.7$,
the two action evaluations agree to $|S_W - S_{U(1)}| < 10^{-12}$
(machine precision, 170 lines of test code).  A fuller
random-Haar-diagonal check at $\beta = 2$ on $n_{\max} = 3$ gives
$|S_W - S_{U(1)}| = 0$ exactly.

\begin{observation}[Paper 25 is the Cartan-subalgebra sector]
\label{obs:cartan_subalgebra}
Paper~25's abelian $U(1)$ Wilson--Hodge construction is structurally
the Cartan-subalgebra projection of the present non-abelian
$\mathrm{SU}(2)$ theory on the same graph.  The question ``is
Paper~25 abelian because the underlying geometry forces it to be, or
because abelian is one choice out of several?'' has a definite
answer: Paper~25 is the abelian sub-sector of a genuinely
non-abelian lattice gauge theory on the same graph.  This elevates
Paper~25's $U(1)$ construction from a freestanding abelian model to
the commutative-reduction sector of a natural non-abelian one.
\end{observation}
```

## Lines 503-517: Remark: higher-rank generalization (ST-SU3 Cartan details)

```latex
\begin{remark}[higher-rank generalization]
\label{rem:higher_rank}
On a rank-$r$ gauge group with maximal torus $\mathrm{T} = U(1)^r$,
the analog of Proposition~\ref{prop:u1_reduction} produces $r$
independent (but coupled, via a $\det = 1$ constraint when present)
$U(1)$ sectors.  The $\mathrm{SU}(3)$ case ($r = 2$) is computed in
Sprint~ST-SU3 (May~2026,
\texttt{debug/st\_su3\_wilson\_memo.md}~\S~3): restricting Wilson
$\mathrm{SU}(3)$ links on the Bargmann--Segal $S^5$ graph to the
maximal torus $T \cong U(1) \times U(1)$ produces a coupled
$U(1) \times U(1)$ Wilson action with character
$(1/3)\bigl(\cos a_P + \cos b_P + \cos(a_P + b_P)\bigr)$, verified
at machine precision over random Cartan trials at $N_\text{max} = 2,
3$ (max $|\Delta S| < 1.5 \times 10^{-11}$ at $7\,765$ plaquettes).
\end{remark}
```

## Lines 538-548: Opening of Result 2 (Paper 25 sections II-III framing)

```latex
Paper~25~\S II--III introduces the edge Laplacian $L_1 = B^{\top} B$
as the ``discrete Hodge-1 Laplacian'' whose kernel dimension equals
$\beta_1$ and whose nonzero spectrum coincides with the nonzero
spectrum of $L_0$ (by singular-value decomposition of $B$), and notes
that on a bare 1-complex only the \emph{exact} part
$d_0 d_0^*$ of the manifold Hodge-1 operator
$\Delta_1 = d_0 d_0^* + d_1^* d_1$ is available --- the co-exact part
requires faces.  The present section identifies the \emph{kinetic
term} of the $\mathrm{SU}(2)$ Wilson action at weak coupling and
shows it is exactly that missing co-exact (curl) part, built on the
Ihara plaquettes as faces.
```

## Lines 583-613: Original d_1 = B^T paragraph and Hodge discussion around Eq. (wilson_quadratic)

```latex
The $O(A^3)$ term does not vanish: the commutator part of $\log U_P$
feeds the cubic vertex $\propto \epsilon^{abc}$
(Observation~\ref{obs:L1_structural}).  The
sum $\sum_{e \in P} A_e^{a}$ is exactly the signed sum of the gauge
field over the plaquette boundary, which is the action of the
discrete exterior derivative $d_1 = B^{\top}$ (acting on 1-cochains;
in fact the $d_1$ that maps 1-cochains to 2-cochains, but on a
triangle-free 1-complex we can read it as a signed sum over
plaquettes).  Therefore
\begin{equation}
S_W[U] = \frac{\beta}{8} \sum_P \Bigl(\sum_{e \in P} A_e^{a}\Bigr)^2
+ O(A^3)
= \frac{\beta}{8}\, A^{a\top} M A^{a} + O(A^3),
\label{eq:wilson_quadratic}
\end{equation}
where $M = B_2 B_2^{\top}$ is the plaquette-incidence bilinear form
on the edge space, built from the signed edge--face incidence matrix
$B_2$ whose columns are the (signed) indicator vectors of the Ihara
plaquettes.  In Hodge terms, $M = d_1^* d_1$ is the \emph{co-exact
(curl) part} of the Hodge-1 Laplacian of the 2-complex obtained by
attaching the plaquettes as faces.  It is \emph{not} proportional to
Paper~25's $L_1 = B^{\top} B = d_0 d_0^*$: the two quadratic forms
have complementary supports.  On the $n_{\max} = 3$ graph
($E = 13$, $\beta_1 = 2$) one computes $\operatorname{rank} M = 2 =
\beta_1$, $\operatorname{rank} L_1 = 11 = E - \beta_1$, and
$L_1 B_2 = 0$ exactly --- the plaquette columns span $\ker L_1$ (the
cycle space), which is precisely where the divergence form $L_1$
vanishes.  Their sum $L_1 + M$ has full rank $E$ and is the complete
discrete Hodge-1 Laplacian of the plaquette-filled 2-complex
(verified in \texttt{tests/test\_su2\_wilson\_gauge.py}, class
\texttt{TestWeakCouplingKinetic}).  The structural statement is:
```

## Lines 639-666: Provenance paragraph (withdrawn proportionality to L_1) and Observation 'completes the Hodge pair'

```latex
\textbf{Provenance: MEASURED.}  The $1/8$ coefficient and the
quadratic form are verified numerically against exact
matrix-exponential holonomies, both on a triangle plaquette and on
the full $n_{\max} = 3$ graph (relative deviation below
$5 \times 10^{-3}$ at $\epsilon = 0.01$), in
\texttt{tests/test\_su2\_wilson\_gauge.py}::\texttt{TestWeakCoupling%
Kinetic} --- the $\mathrm{SU}(2)$ mirror of the $\mathrm{SU}(3)$
$1/12$ test in \texttt{tests/test\_su3\_wilson\_s5.py}.  An earlier
version of this proposition asserted that the kinetic form
``coincides with $L_1 = B^{\top} B$ up to a positive scalar
multiple''; that assertion was wrong (the two forms have
complementary supports, as the rank computation above shows) and is
withdrawn --- the correct statement is the completion statement,
which matches the face-structure reading already recorded in
Paper~25~\S II.

\begin{observation}[The Wilson kinetic sector completes the Hodge pair]
\label{obs:L1_structural}
At weak coupling the non-abelian Wilson theory contributes the
co-exact half of the Hodge-1 pair: the free (Gaussian) sector is the
curl form $B_2 B_2^{\top}$ on the cycle space, while Paper~25's
$L_1 = B^{\top} B$ is the exact half, whose role is combinatorial
(divergence / gradient sector) rather than plaquette-kinetic.  The
non-abelian self-couplings (cubic and quartic $A$-vertices) appear at
order $A^3$ and $A^4$ in the Wilson action expansion and involve the
structure constants of $\mathfrak{su}(2)$, namely $\epsilon^{abc}$;
they are parametrically suppressed at weak coupling.
\end{observation}
```

## Lines 668-693: Remark: universal SU(N) coefficient (original wording)

```latex
\begin{remark}[universal $\mathrm{SU}(N)$ kinetic coefficient]
\label{rem:sun_coefficient}
With generators $T^a = \lambda^a/2$ normalized to
$\mathrm{tr}(T^a T^b) = \tfrac{1}{2} \delta^{ab}$, the second-order
expansion of the Wilson action
$S_W[U] = \beta \sum_P (1 - (1/N_c)\,\mathrm{Re}\,\mathrm{tr}\,U_P)$
on a general gauge group $\mathrm{SU}(N_c)$ produces a per-plaquette
per-component coefficient of $1/(4 N_c)$.  The $\mathrm{SU}(2)$ value
$1/8$ derived above is the $N_c = 2$ specialization, verified
numerically in \texttt{tests/test\_su2\_wilson\_gauge.py}::%
\texttt{TestWeakCouplingKinetic::test\_quadratic\_expansion\_numerically}.
The $\mathrm{SU}(3)$ value $1/12$ is verified independently in
Sprint~ST-SU3 (\texttt{debug/st\_su3\_wilson\_memo.md}~\S~4)
symbolically and numerically on the Bargmann--Segal $S^5$ graph
($120$ degrees of freedom at $N_\text{max} = 2$, relative error
$\sim O(\epsilon^2)$;
frozen test \texttt{tests/test\_su3\_wilson\_s5.py}::%
\texttt{TestWeakCouplingKinetic}).
The mechanism: more colors $\Rightarrow$ smaller per-component
coefficient because the $(1/N_c)\,\mathrm{Re}\,\mathrm{tr}$
normalization in the Wilson action divides by $N_c$ while the
trace-of-pair normalization holds $\mathrm{tr}(T^a T^a) = 1/2$ fixed.
The decoupling of components at quadratic order is a universal
$\mathrm{SU}(N_c)$ feature; non-abelian self-interactions enter at
$O(A^3)$ via the structure constants $f^{abc}$ in every case.
\end{remark}
```

## Lines 695-772: Result 3: Monte Carlo Wilson-loop table (unconverged MC column), confinement paragraph, Paper 29 Ihara-count consistency

```latex
\subsection{Result 3: Plaquette counts and Wilson-loop data}
\label{sec:results_MC}

Monte Carlo with Metropolis updates on the full
$\mathrm{SU}(2)^{|E|}$ configuration space provides direct access to
the Wilson-loop expectation value $\langle W_C \rangle = \bigl\langle
\tfrac{1}{2}\mathrm{Re}\,\mathrm{tr}\,U_C \bigr\rangle$.  We run at
$n_{\max} = 3$ (13 forward links), 2000 samples after 500 thermalization
steps, seed = 12, for a Wilson loop coinciding with one of the two
primitive 4-cycles in the $p$-block.

\begin{table}[h]
\centering
\caption{Wilson-loop expectation on the $n_{\max} = 3$ Coulomb graph,
loop = first $p$-block 4-cycle
[$(2,1,-1) \to (2,1,0) \to (3,1,0) \to (3,1,-1)$].  MC: Metropolis
sampling, 2000 samples, stderr reported; the seeded chain is
bit-deterministic but not converged (its deviations from the exact
column exceed the quoted stderr, which ignores autocorrelation).
Exact: $\langle W \rangle = I_2(\beta) / I_1(\beta)$, from the
independence of the two plaquettes.}
\label{tab:wilson_loop}
\begin{tabular}{rrrr}
\toprule
$\beta$ & $\langle W \rangle_{\mathrm{MC}}$ (2000 samples) & $\pm\,$stderr & $\langle W \rangle_{\text{exact}}$ \\
\midrule
0.5 & $0.126$ & $0.011$ & $0.124$ \\
1.0 & $0.158$ & $0.011$ & $0.240$ \\
2.0 & $0.332$ & $0.009$ & $0.433$ \\
5.0 & $0.681$ & $0.006$ & $0.719$ \\
\bottomrule
\end{tabular}
\end{table}

The exact column follows because $P_A$ and $P_B$ of
Eq.~(\ref{eq:p_block_plaquettes}) form a cycle basis ($\beta_1 = 2$):
after gauge-fixing on a spanning tree each plaquette holonomy is an
independent Haar variable, and
$\langle \tfrac12 \mathrm{tr}\,U \rangle$ under the weight
$e^{(\beta/2)\mathrm{tr}\,U}$ is $I_2(\beta)/I_1(\beta)$.  The short
chain's deviations ($34\%$ at $\beta = 1$, $23\%$ at $\beta = 2$) are
sampling error, not shared-edge correlation: three independent chains
of 20\,000 samples give $0.232$--$0.249$ at $\beta = 1$ and
$0.429$--$0.452$ at $\beta = 2$, scattered about the exact values.
$\langle W \rangle$ is monotonic in $\beta$, rising from $0$ at strong
coupling to $1$ at weak coupling.

\paragraph*{Confinement diagnostic on a factorizing graph.}
Confinement in Wilson's sense is the statement that $\langle W_C
\rangle$ obeys an \emph{area law} $\langle W_C \rangle \sim
e^{-\sigma \cdot \mathrm{Area}(C)}$ at strong coupling, crossing over
to a \emph{perimeter law} at weak coupling.  The diagnostic requires
a family of Wilson loops indexed by a geometric area.  On the
$n_{\max} = 3$ Coulomb graph the 6-cycle bounding $P_A \cup P_B$ gives
a loop of area $2$, but because the plaquettes are independent its
expectation is exactly $(I_2/I_1)^2$: an area law that holds at every
$\beta$ by factorization, with no dynamical content and no
area--perimeter crossover.  \emph{This is a structural feature of the
small graph, not a statement about the field theory}.  It persists at
every $n_{\max}$: each $l$-block of the Coulomb graph is a planar
$(n_{\max}-l)\times(2l+1)$ grid whose unit squares form a cycle basis
($\beta_1 = \sum_l 2l(n_{\max}-l-1)$, giving $2$ and $8$ at
$n_{\max} = 3, 4$), so with the 4-plaquette action the theory is the
exactly solvable two-dimensional lattice gauge theory and every
contractible loop of area $A$ has
$\langle W \rangle = (I_2/I_1)^{A}$ [INTERNAL THEOREM].  To see
confinement in the conventional sense, one needs plaquettes that do
not factorize, which this label graph does not supply.

\paragraph*{Consistency with Paper 29's Ihara cycle counts.}  The
primitive 4-cycle count at $n_{\max} = 3$ ($N_4 = 2$) and the
higher-length counts ($N_6 = 1$, $N_8 = 1$) match the M\"obius-inverted
Ihara primitive-walk counts of Paper~29~\cite{paper29} up to the
orientation factor $N_{\text{M\"obius}}(L) = 2 \cdot N_{\text{unoriented}}(L)$
(each undirected cycle gives clockwise and counterclockwise primitive
walks; we keep one representative per unoriented class).  Paper~29's
plaquette-enumeration machinery is reused here without modification;
the Ihara cycles of one paper are the Wilson plaquettes of the other.
```

## Lines 774-856: Section: Principle of Least Action, framework-wide synthesis (five action principles, packing-Wilson duality)

```latex
% ===================================================================
\section{Principle of Least Action: A Framework-Wide Synthesis}
\label{sec:action_principles}
% ===================================================================

The Wilson action of Eq.~(\ref{eq:wilson_action}) organizes the
gauge-sector dynamics of the GeoVac graph as a variational principle
on $\mathrm{SU}(2)^{|E|}$.  This is the fifth independent action
principle in the GeoVac framework, and it is useful to place it in
context.

Over the course of the project, variational or action structures have
emerged on five distinct data sectors of the Hopf graph.  Each is
tightly localized to one kind of structure and does not substitute
for any of the others.  We collect them in
Table~\ref{tab:action_principles}.

\begin{table}[t]
\centering
\small
\caption{Five action principles in the GeoVac framework.  Each
organizes a different sector of the graph data.  The Wilson-action
entry is the contribution of this paper.}
\label{tab:action_principles}
\begin{tabular}{>{\raggedright\arraybackslash}p{3.2cm}>{\raggedright\arraybackslash}p{3.1cm}>{\raggedright\arraybackslash}p{3.0cm}>{\raggedright\arraybackslash}p{3.2cm}}
\toprule
Principle & Data sector & Extremum principle & Reference \\
\midrule
Maximum entropy on packing & Nodes (candidate shell centers) & Maximize entropy under isotropy (Axiom 2) & Paper~0~\cite{paper0} \\
\addlinespace
Rayleigh--Ritz on graph Hamiltonian & Node wavefunctions (0-cochains) & Minimize $\langle \psi | H | \psi \rangle / \langle \psi | \psi \rangle$ & Papers~1, 7~\cite{paper7} \\
\addlinespace
Wilson lattice gauge action & Edge link variables (1-cochains, SU(2)-valued) & Minimize $S_W[U] = \beta \sum_P (1 - \tfrac{1}{2}\mathrm{Re}\,\mathrm{tr}\,U_P)$ & Paper~25~\cite{paper25} $U(1)$; present paper SU(2) \\
\addlinespace
Ihara zeta path integral & Primitive closed non-backtracking walks & $\log \zeta_G(s)^{-1} = \sum_{[C]} s^{\ell(C)}$ as weighted walk sum & Paper~29~\cite{paper29} \\
\addlinespace
Spectral action on Dirac spectrum & Spinor-bundle spectrum ($|\lambda_n|, g_n$) & $\mathrm{Tr}\,f\bigl(D / \Lambda\bigr)$ over modes & Paper~28~\cite{paper28} \\
\bottomrule
\end{tabular}
\end{table}

\paragraph*{What each sector is about.}
The node sector carries the scalar matter content (the wavefunction
amplitudes) and hosts two action principles: Paper~0's packing
axioms as a max-entropy statement (``maximize entropy under
isotropy''), which \emph{builds} the graph, and the Rayleigh--Ritz
principle on the built graph's Hamiltonian, which computes the
graph's spectrum.  The edge sector carries gauge content (Paper~25's
$U(1)$ phases, this paper's SU(2) link variables) and hosts the
Wilson action.  The walk sector aggregates edge data into
Ihara-cycle weights and hosts a discrete path integral
(log-partition function = walk sum).  The Dirac spectrum, finally,
is the specific data set on which the spectral action of
noncommutative-geometry type is computed; Paper~28's QED spectral
sums are literally this.

\paragraph*{The packing--Wilson duality.}
Two of the five principles stand out as structural duals.  Paper~0's
Axiom 2 (geometric indifference) is a max-entropy principle on node
data: among all candidate discrete packings compatible with Axiom 1
(binary distinguishability), the entropy-maximizing one is selected.
The Wilson action is a min-action principle on edge data: among all
SU(2) link configurations, those minimizing plaquette curvature are
selected in the continuum ($\beta \to \infty$) limit.  The two
principles live on complementary cochain degrees (0-form vs 1-form)
and extremize in opposite directions (max vs min), but they are
structurally compatible: Paper~0's max-entropy on nodes fixes
\emph{which} graph the Wilson action lives on, and the Wilson action
fixes \emph{how} gauge dynamics runs on it.  Nothing in the
mathematics forbids a joint variational principle on the full
(node, edge) data; we note this for future work but do not pursue it
here.

\paragraph*{What this section does and does not claim.}
The observation is descriptive: the GeoVac framework, as it currently
stands, presents five independent variational structures, each
localized to a different sector of the graph data.  We do not claim
that these five can be reduced to a single grand principle, nor that
one is more fundamental than the others.  Each is a tool matched to a
data type.  The utility of naming them together is that it makes the
project's full toolkit visible in one place: the framework is not
reducible to ``diagonalize a graph Laplacian,'' because its data has
more kinds of structure than nodes alone.
```

## Lines 859-920: Section: Yang-Mills adjacency, what is and is not claimed

```latex
\section{Yang--Mills Adjacency: What Is and Is Not Claimed}
\label{sec:scope}
% ===================================================================

The present paper inherits the language of Wilson lattice gauge theory.
The rhetoric rule of CLAUDE.md~\S1.5 requires that we be explicit about
scope.  We enumerate what is done and what is not.

\paragraph*{What this paper establishes.}  A gauge-invariant, finite,
Haar-normalizable non-abelian lattice gauge theory on the finite
compact GeoVac Hopf graph at any $n_{\max}$, with (i) link variables
in $\mathrm{SU}(2)$, (ii) plaquettes as primitive closed
non-backtracking walks, (iii) Wilson action and partition function
well-defined without renormalization, (iv) exact $U(1)$ reduction to
Paper~25 in the maximal-torus limit, (v) the weak-coupling kinetic
term identified as the co-exact (curl) completion of
$L_1 = B^{\top} B$, with universal coefficient $1/8 = 1/(4N_c)$.
All 36 tests in \texttt{tests/test\_su2\_wilson\_gauge.py} pass,
including the numerical $U(1)$-reduction test, end-to-end gauge
invariance at machine precision (asserted $< 10^{-13}$), the kinetic
quadratic-form tests, and the exact Wilson-loop values $I_2/I_1$.

\paragraph*{What this paper does NOT establish.}
\begin{itemize}
\item \emph{Yang--Mills mass gap on $\mathbb{R}^4$.}  The Clay
  Millennium Problem requires a continuum-field-theory construction
  of pure Yang--Mills on 4D Minkowski or Euclidean space, proof of a
  positive mass gap, and proof that no interacting scale-invariant
  continuum limit exists.  Our construction lives on a fixed finite
  compact graph; there is no continuum limit of $\mathbb{R}^4$ inside
  it.
\item \emph{Confinement in the conventional sense.}  As shown in
  Section~\ref{sec:results_MC}, the plaquettes form a cycle basis of
  each planar $l$-block at every $n_{\max}$, so Wilson loops obey
  $\langle W \rangle = (I_2/I_1)^{\mathrm{Area}}$ exactly: an area law
  by factorization, which says nothing about confinement.
\item \emph{Lorentzian signature or dynamics.}  This is a Euclidean
  lattice gauge theory on a spatial slice: no time direction, no
  Wick rotation, no $S$-matrix, no scattering.  Dynamics in the
  usual sense is absent.
\item \emph{Renormalization or continuum limit.}  The bare Wilson
  action is the whole story on a finite graph.  No continuum SU(2)
  Yang--Mills limit on $S^3$ is expected: the label graph's edge
  Laplacian does not converge to the $S^3$ Laplacian.
\item \emph{No derivation of $\alpha$.}  Paper~2's
  combination rule $K = \pi(B + F - \Delta)$ is not refined by the
  present construction.  Paper~25~\S IV read the three ingredients as
  Casimir-matter trace / Fock-Dirichlet gauge zeta / Dirac-boundary
  count; the non-abelian SU(2) Wilson structure does not produce a
  new $K$-like object or a non-abelian ``fine-structure constant''
  coupling, and we make no claim in that direction.
\end{itemize}

\paragraph*{Relationship to the Clay problem.}
The mathematical vocabulary is shared (lattice gauge theory,
non-abelian link variables, plaquette action, Wilson loops), but the
setting differs categorically.  The Clay problem is about continuum
field theory on a non-compact 4D spacetime; the present construction
is about discrete gauge theory on a compact 3D-graph spatial slice.
We use Wilson's formalism because it is the correct tool for
non-abelian lattice gauge theories on finite graphs; we do not
imply that the construction transports to $\mathbb{R}^4$.
```

## Lines 923-1011: Open questions: matter coupling, larger graphs, SU(2) Ihara zeta, SU(2) analog of K, transverse photon sector, continuum limit

```latex
\section{Open Questions}
\label{sec:open}
% ===================================================================

\subsection{Matter coupling: Dirac fermions on the Hopf graph}

The GeoVac framework's Tier-2 relativistic pipeline~\cite{paper14}
builds spin-ful composed qubit Hamiltonians in the
$(\kappa, m_j)$ Dirac basis on $S^3$.  These are matter sectors on
the same graph whose gauge sector we have built here.  Coupling the
two---Dirac fermion matter on nodes, $\mathrm{SU}(2)$ Wilson gauge on
edges---would be the first matter$+$gauge lattice theory in the
GeoVac project.  The construction is standard in principle (covariant
hopping amplitudes $\bar\psi_i U_{ij} \psi_j$), but on a
triangle-free simplicial 1-complex with Dirac spinors it would
require the Szmytkowski angular labels of Paper~14's Tier~2 to be
reconciled with the Wilson link-variable orientation convention.  A
clean demonstrator would be: does the spinor spectrum feel the
non-trivial $\mathrm{SU}(2)$ topology (i.e., is the Dirac operator
sensitive to the $\beta_1 = 2$ harmonic 1-forms at $n_{\max} = 3$ when
the gauge field is non-trivially topological)?

\subsection{Larger graphs and the confinement diagnostic}

With the 4-plaquette action this question is closed at every
$n_{\max}$: the plaquettes form a cycle basis of each planar $l$-block,
so Wilson loops factorize and obey $(I_2/I_1)^{\mathrm{Area}}$ exactly
(Section~\ref{sec:results_MC}), with no area--perimeter crossover.  A
confinement diagnostic with dynamical content needs an action whose
plaquettes do not factorize; the label graph, whose $l$-blocks are
disconnected, does not supply one.

\subsection{SU(2) analog of the Ihara zeta}

Paper~29 built the (scalar-valued) Ihara zeta of the GeoVac Hopf
graphs from primitive closed non-backtracking walk counts.  A natural
non-abelian generalization would weight each primitive walk $[C]$ by
a character $\chi_j(U_{[C]})$ of its $\mathrm{SU}(2)$ holonomy,
producing a character-valued zeta function.  Whether the resulting
object has a well-defined Ihara--Bass-type determinantal expression,
or whether it satisfies a graph-Riemann Hypothesis analog, is
speculative.  We record the question without pursuing it.

\subsection{SU(2) analog of Paper 2's $K = \pi(B + F - \Delta)$}

Paper~25~\S VII.A suggested that a non-abelian analog of the Paper~2
combination rule might exist on the $S^3$ graph.  The SU(2) Wilson
structure of the present paper does not by itself produce one; the
Wilson action is coupling-insensitive in a way that Paper~2's $K$ is
not.  A more promising route (if any) would couple the SU(2) Wilson
gauge to Dirac-spinor matter (first open question above) and examine
the combined spectral-gauge invariants.  We do not claim such a route
will succeed; we record it as a natural follow-up within the
Paper-2 observation-tier program.

\subsection{Plaquettes as the transverse photon sector}

The plaquettes enumerated in this paper define a boundary operator
$d_1\colon C^1 \to C^2$ whose transpose gives the co-exact Laplacian
$L_1^{\mathrm{co}} = d_1 d_1^\top$ on edges.  This is the
\emph{transverse} (physical photon) sector of the graph Hodge
decomposition, complementing the exact (longitudinal) sector
$L_1^{\mathrm{ex}} = B^\top B$ used in Paper~25 and the graph-native
QED of Paper~28.

Paper~28's pendant-edge theorem gives
$\Sigma_{\mathrm{scalar}}(\mathrm{GS}) = 2(n_{\max}-1)/n_{\max}$ for
the scalar propagator $G = (B^\top B)^+$, while the transverse
propagator $G_T = (d_1 d_1^\top)^+$ gives $\Sigma_T(\mathrm{GS}) = 0$.
That zero is not a ground-state selection rule: the $l = 0$ block is a
path, every $s$-edge is a bridge lying in no closed walk, and $G_T$
vanishes on all of them, so $\Sigma_T(n\mathrm{S}) = 0$ for every $n$.
Nor is there a continuum zero for it to recover, since the physical
vector vertex couples the ground state (Paper~28, Theorem on the
self-energy, as re-scoped there).  The plaquettes therefore supply the
transverse sector of the graph Hodge decomposition, but no selection
rule of continuum QED is recovered through them.

\subsection{Continuum limit}

The continuum limit of lattice gauge theory on refined Alon--Boppana
sequences of graphs converging to a target manifold is a technical
question that has been addressed in specific cases (hypercubic
lattices converging to $\mathbb{R}^4$, spherical lattices converging
to $S^3$).  The GeoVac label graph is not such a sequence: its edge
Laplacian splits into $n_{\max}$ disconnected $l$-blocks and does not
converge to the $S^3$ Laplacian (the documented graph-to-$S^3$ wall),
so no continuum Yang--Mills limit on $S^3$ is expected from this
construction.
```

## Lines 1013-1089: Open question: SU(3) Wilson sibling on S^5 (ST-SU3; Coulomb/HO layer reading)

```latex
\subsection{SU(3) Wilson sibling on $S^5$}
\label{sec:su3_sibling}

Paper~25~\S~VII.A also flagged as open whether a non-abelian extension
of the GeoVac Hopf $U(1)$ Wilson--Hodge structure transfers to the
Bargmann--Segal $S^5$ graph (Paper~24).  Sprint~5 Track~S5
(\texttt{debug/s5\_gauge\_structure\_memo.md}) returned a clean
negative for the \emph{adjacency-preserving} SU(3) action on the
$(N, 0)$ representation tower: shell-shifting transitions are
Clebsch--Gordan intertwiners between distinct SU(3) irreps, not group
elements, so an SU(3) action on the $(N, 0)$ tower fails the Wilson
``fixed group on every link'' requirement.

Sprint~ST-SU3 (May~2026,
\texttt{debug/st\_su3\_wilson\_memo.md}, with module
\texttt{geovac/su3\_wilson\_s5.py} and tests
\texttt{tests/test\_su3\_wilson\_s5.py} 49/49 pass, two slow-marked) tested a
structurally different question: \emph{Wilson lattice gauge theory}
with link variables $U_e \in \mathrm{SU}(3)$ as external 3$\times$3
matrices in the fundamental representation, independent of the shell
index $N$.  The construction parallels the present paper's
$\mathrm{SU}(2)$ on $S^3$.  Three structural results mirror those of
the present paper:

\begin{itemize}
\setlength\itemsep{0pt}
\item \textbf{Cartan reduction.}  Restricting Wilson SU(3) links to
the maximal torus $T \cong U(1)\times U(1)$ produces a coupled
$U(1)\times U(1)$ Wilson action with character
$(1/3)\bigl(\cos a_P + \cos b_P + \cos(a_P+b_P)\bigr)$ (the third
phase forced by $\det = 1$).  Verified at machine precision over
random Cartan trials at $N_\text{max} = 2, 3$
(max $|\Delta S| < 1.5\times 10^{-11}$ at $7\,765$ plaquettes).
The rank-2 generalization of Proposition~\ref{prop:u1_reduction} per
Remark~\ref{rem:higher_rank}.
\item \textbf{Kinetic coefficient.}  The weak-coupling expansion gives
$1/12 = 1/(4 \cdot 3)$ per plaquette per generator component, the
$N_c = 3$ specialization of Remark~\ref{rem:sun_coefficient}'s
universal $1/(4 N_c)$ formula.  Eight $\mathfrak{su}(3)$ components
are kinetic-decoupled at quadratic order; non-abelian self-couplings
enter at $O(A^3)$ via the structure constants $f^{abc}$.
\item \textbf{Plaquette tabulation.}  $V, E, \beta_1$ at
$N_\text{max} = 2, 3, 4$ are $(10, 15, 6)$, $(20, 42, 23)$,
$(35, 90, 56)$, with primitive 4-cycle counts $N_4 = 15, 88, 272$
--- roughly an order of magnitude denser than the $S^3$ graph at
matched cutoff.  Wilson loops at $N_\text{max} = 2$ with local SU(3)
updates give $\langle W \rangle = 0.05, 0.10, 0.26, 0.65, 0.89$ at
$\beta = 0.5, 1, 2, 5, 10$, systematically lower than the present
paper's SU(2) values at small $\beta$ (consistent with the standard
$N_c$-color randomization of the trace).
\end{itemize}

The Sprint~5 Track~S5 Clebsch--Gordan obstruction is \emph{bypassed}
at the gauge level (Wilson links live in the fundamental SU(3) rep,
not in the $(N, 0)$ shell rep) but \emph{rediscovered} at the
matter-coupling level: coupling Wilson SU(3) to natural shell matter
(in the $(N, 0)$ irreps of growing dimension) requires a CG
intertwiner from $\textbf{3}\otimes (N, 0)$ to $(N+1, 0)$ that is
exactly the Sprint~5 obstruction.  Wilson SU(3) on Bargmann is
therefore a self-consistent \emph{pure-gauge} theory; the natural
matter sector cannot be coupled in the strict Wilson sense.

The structural reading: $S^3 = \mathrm{SU}(2)$ makes the
present paper's case special in a way that does not transfer to the
$\mathrm{SU}(3)$ on $S^5$ analogue.  On $S^3$, matter (spin-$1/2$
fermions, Paper~14 Tier~2) lives in the same representation as the
link variables, so the gauge--matter coupling is natural (an argument,
not yet a coupled construction; Paper~25); on $S^5$,
the matter representations $(N, 0)$ have dimensions that grow with
$N$, so a uniform-link-rep coupling is impossible.  This sharpens
Paper~24's Coulomb/HO asymmetry from three layers
(spectrum-computing $L_0$ and calibration $\pi$, both since
withdrawn, and Wilson gauge with matter coupling) plus the universal Wilson--Hodge vocabulary, to
\emph{four layers}: Wilson gauge alone is universal; only Wilson
gauge \emph{with natural matter coupling} is Coulomb-specific.
Sprint~ST-SU3 contributes the third layer on the gauge-with-matter
side.  The reading of Paper~24 \S~V is updated accordingly.
```

## Lines 1092-1140: Conclusion (four-roles S^3 framing, five action principles)

```latex
\section{Conclusion}
\label{sec:conclusion}
% ===================================================================

Sprint 4 Track RH-Q constructed SU(2) Wilson lattice gauge theory on
the GeoVac Hopf graph (module \texttt{geovac/su2\_wilson\_gauge.py},
36 tests passing, 2026-07-04 count).  The present paper formalizes that construction
and places it in context relative to Paper~25's abelian $U(1)$
predecessor.

Three structural results summarize the content:
(1)~the maximal-torus reduction of the $\mathrm{SU}(2)$ Wilson action
recovers Paper~25's $U(1)$ Wilson--Hodge action exactly, so
Paper~25 is the Cartan-subalgebra sector of a non-abelian theory;
(2)~the weak-coupling kinetic term of the non-abelian Wilson action
is the co-exact (curl) quadratic form on the Ihara plaquettes ---
coefficient $1/8 = 1/(4N_c)$ per plaquette per
$\mathfrak{su}(2)$ component --- which completes the edge Laplacian
$L_1 = B^{\top} B$ that Paper~25 introduced to the full discrete
Hodge-1 Laplacian, with non-abelian self-interactions entering at
cubic order in the gauge field;
(3)~plaquette enumeration on the $S^3$ Coulomb graph matches
Paper~29's Ihara cycle counts, and because the plaquettes form a cycle
basis of each planar $l$-block, Wilson loops factorize exactly,
$\langle W \rangle = (I_2/I_1)^{\mathrm{Area}}$, which long Monte Carlo
chains reproduce.

We close with a framework observation.  The $S^3 = \mathrm{SU}(2)$
coincidence is structurally overdetermined: $S^3$ appears in the
framework as Fock's projection of the Coulomb problem (Paper~7), as
the underlying manifold of the $(n, l, m_l)$ Hopf graph (Paper~14),
as the spin-1/2 carrier of the Dirac sector (Paper~14 Tier-2), and
now as the gauge group of the non-abelian Wilson construction.  Four
independent structural roles on the same 3-manifold, each needing
only the simplest compact simply-connected non-abelian Lie group of
rank~1.  We do not read this as evidence that $S^3$ is fundamental in
any ontological sense (CLAUDE.md~\S1.5 rhetoric); we do observe that
$S^3$ is, in a specific technical sense, the minimal object that
supports all four roles simultaneously, and we record the convergence
as a fact.

The project now contains five action principles---Paper~0's
max-entropy on nodes, Rayleigh--Ritz on the graph Hamiltonian,
Paper~25's and this paper's Wilson action on edges, Paper~29's Ihara
zeta path integral on walks, and Paper~28's spectral action on the
Dirac spectrum---each matched to a distinct sector of the graph
data.  None reduces to any of the others.  The Wilson action has,
with this paper, a genuine non-abelian sibling to its $U(1)$
predecessor; the rest of the table is left as is.
```

## Lines 1161-1163: Bibitem fegan1983 (dropped)

```latex
\bibitem{fegan1983}
H.~D.~Fegan, ``The heat equation on a compact Lie group,''
\textit{Trans.~Amer.~Math.~Soc.}\ \textbf{246}, 339--357 (1978).
```

## Lines 1238-1241: Bibitem rh_q_memo (debug/ citation, dropped)

```latex
\bibitem{rh_q_memo}
GeoVac Project, ``Track RH-Q: SU(2) Wilson Lattice Gauge Theory on
the $S^3 = \mathrm{SU}(2)$ GeoVac Hopf Graph,''
\texttt{debug/su2\_wilson\_gauge\_memo.md}, April 2026.
```
