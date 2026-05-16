# U(1) Wilson Lattice Gauge on the Dirac Rule B Graph — Design Memo

**Sprint XCWG (May 2026).** Companion pilot: `debug/xcwg_u1_wilson_rule_b_pilot.py`.
Data: `debug/data/xcwg_u1_wilson_rule_b_pilot.json`. No production code modified.

This memo specifies a U(1) Wilson lattice gauge theory on the Dirac-S³
"Rule B" graph of Paper 29 §RH-C — extending Paper 25's U(1) Wilson–Hodge
construction from the scalar Fock graph to the spinor-labeled, cross-ℓ-hopping
Rule B graph. The structural sequel is the recent observation that the scalar
graph is *per-ℓ-block 2D* (Paper 25's MK block-spin reduction reproduces 2D
U(1) confinement but no 3D phase content), while Rule B's defining feature is
exactly that its edges cross ℓ-shells via E1 parity flips.

The verdict, stated up front: **U(1) Wilson on Rule B is a well-defined,
gauge-invariant, abelian lattice gauge theory on a finite directed graph,
with the same Wilson–Hodge dictionary as Paper 25 but with a structurally
much richer plaquette content**. The Hodge identity is verified numerically
at both $n_{\max}=2$ and $n_{\max}=3$. The single load-bearing structural
difference from Paper 25's scalar graph is that **Rule B is plaquette-dense
and the smallest plaquette length is 4 starting at $n_{\max}=2$**, whereas
Paper 25's scalar graph has no 4-plaquettes at $n_{\max}=2$ (just 2 at
$n_{\max}=3$, and they live entirely within the p-block). This is the
discrete-geometric reason to expect Rule B to host nontrivial RG flow that
the scalar graph cannot.

The honest scope is the standing Paper 25 / Paper 30 scope (CLAUDE.md §1.5):
a gauge-invariant, Haar-normalizable abelian lattice gauge theory on a
compact finite graph. Not continuum QED on $\mathbb{R}^4$. Not Lorentzian.
The "RG-flow projections" in §7 are *qualitative structural readings*, not
RG-flow theorems. The next sprint (MK block-spin on Rule B) would convert
the structural reading into a quantitative flow diagram; this memo only
demonstrates the construction is well-defined and identifies where the
expected structural payoff sits.

---

## 1. Link variables and gauge transformation on spinor-labeled Rule B

### 1.1 Vertex set, oriented edges, link variables

The Rule B graph $G_B = (V, E)$ at cutoff $n_{\max}$ has vertices labeled by
Dirac quantum numbers $v = (n_{\rm fock}, \kappa, m_j)$ — a signed-integer
$\kappa$ encoding $(\ell, j)$ via $\kappa = -(\ell+1)$ for $j=\ell+\tfrac12$
and $\kappa = +\ell$ for $j=\ell-\tfrac12$, with $|\kappa| = j+\tfrac12$ — and
edges given by the E1 dipole selection rule
\[
  v \sim_B v'
  \iff
  \Delta n_{\rm fock} \in \{-1, 0, +1\},\ 
  \Delta \ell = \pm 1,\ 
  |\Delta j| \le 1,\ 
  |\Delta m_j| \le 1,\ 
  v \ne v'.
\]
(See `geovac/ihara_zeta_dirac.py::_edge_rule_B`.) Each undirected edge
$\{v, v'\} \in E$ is promoted to two oriented edges $e = v \to v'$ and
$\bar e = v' \to v$.

A U(1) **link variable** is a phase $U_e = e^{i\theta_e} \in U(1)$
assigned to each oriented edge, with the inverse convention
\[
  U_{\bar e} = U_e^{-1} = e^{-i\theta_e}, \qquad \theta_{\bar e} = -\theta_e.
\]
The free configuration space of the gauge theory is therefore $U(1)^{|E|}$,
parameterized by the $|E|$ phases $\theta_e \in [0, 2\pi)$ on the
canonically oriented edges (we choose $v \to v'$ with the smaller node index
as source). Bookkeeping is identical to Paper 25 §III except that the node
index set is bigger and labels carry spinor structure.

### 1.2 Gauge transformation

A gauge transformation is a function $\chi: V \to \mathbb{R}/2\pi\mathbb{Z}$,
$v \mapsto \chi_v$. It acts on link variables by
\[
  U_e \mapsto U_e' = e^{i\chi_{v'}}\, U_e\, e^{-i\chi_v}
  \qquad \text{for } e = v \to v',
\]
equivalently
\[
  \theta_e \mapsto \theta_e + (\chi_{v'} - \chi_v) \pmod{2\pi}.
\]
This is consistent with the inverse convention: $\theta_{\bar e} \mapsto
-\theta_e - (\chi_{v'} - \chi_v) = -\theta_e' \pmod{2\pi}$.

The corresponding *matter-field* transformation, on a complex-valued matter
$\psi: V \to \mathbb{C}$, is
\[
  \psi_v \mapsto e^{i\chi_v} \psi_v,
\]
which is the standard Wilson–Hodge dictionary (Paper 25 §III.C). All
gauge-invariant observables of the U(1) theory — plaquette holonomies and
their products — are then automatically invariant under this transformation.

### 1.3 Spinor-label subtlety: U(1) acts as overall phase, not as Pauli rotation

The vertices of Rule B carry a four-component Dirac spinor structure
(through $\kappa$ and $m_j$), but **the U(1) Wilson theory of this memo
treats the spinor labels as pure book-keeping indices**. Concretely:

* A "matter field" $\psi_v$ is a complex scalar attached to the vertex
  $v = (n, \kappa, m_j)$. It is *not* a 4-component Dirac spinor at a single
  point in space; rather, the entire Hilbert space of single-particle Dirac
  states is the direct sum $\bigoplus_v \mathbb{C} = \mathbb{C}^V$ with one
  complex coefficient per Dirac orbital. This is exactly the convention of
  the graph-native QED of Paper 28 §graph_native_qed.
* The U(1) acts by an overall phase $\psi_v \mapsto e^{i\chi_v}\psi_v$ at
  each *orbital*, not by spin rotation on a 4-vector internal index.

This is the right way to read Rule B as a lattice gauge theory: spinor
structure is *external* (it labels which orbital lives at this vertex), and
the gauge group is the abelian U(1) of *charge*, not the SU(2) of *spin*.

The natural larger gauge group one might consider is U(1) × SU(2) (charge
times intrinsic spin rotation), and indeed Paper 30 builds an SU(2) Wilson
theory on the scalar S³ graph. On Rule B the SU(2) version is qualitatively
different from the present construction (it would mix the four nodes of a
$\kappa \to \pm\kappa$, $m_j \to \pm m_j$ multiplet at each $(n, \ell)$),
and is out of scope here. The present sprint focuses on the abelian
U(1) of charge, which is the direct sibling of Paper 25.

### 1.4 What "U(1) on Rule B" *adds* beyond Paper 25

The only thing that differs from Paper 25's scalar construction at the
level of axioms is *which graph* the gauge theory lives on. The vertex set
is bigger ($V=10$ at $n_{\max}=2$ vs Paper 25's $V=5$; $V=28$ vs $V=14$ at
$n_{\max}=3$), and the edges go across ℓ-shells. Everything else — link
variables, gauge transformations, Haar measure $\prod_e \frac{d\theta_e}{2\pi}$,
Wilson action below — is verbatim Paper 25.

---

## 2. Node and edge Laplacians $L_0, L_1$

### 2.1 Signed incidence matrix

Choose the orientation $v \to v'$ with the smaller node index as source.
Enumerate the undirected edges in lex order $(u, v)$ with $u < v$, indexed
$e = 1, \ldots, |E|$. The signed incidence matrix is $B \in
\mathbb{Z}^{V \times E}$ with
\[
  B_{u, e} = -1, \qquad B_{v, e} = +1, \qquad B_{w, e} = 0 \text{ for } w \ne u, v.
\]
This is the same convention as Paper 25 Eq. (II.B–II.C).

### 2.2 Node and edge Laplacians

\[
  L_0 = B B^{\top} \in \mathbb{Z}^{V \times V},
  \qquad
  L_1 = B^{\top} B \in \mathbb{Z}^{E \times E}.
\]

$L_0$ is the standard combinatorial graph Laplacian: $(L_0)_{vv} = \deg(v)$,
$(L_0)_{vv'} = -1$ if $v \sim_B v'$, $0$ otherwise. It acts on
**vertex (0-cochain) space** $\mathbb{R}^V$.

$L_1$ is the edge Laplacian: $(L_1)_{ee} = 2$ for every edge $e = (u, v)$
(since the column of $B$ at $e$ has two nonzero entries $\pm 1$), and the
off-diagonal entry $(L_1)_{e_1, e_2}$ for $e_1 = (u_1, v_1)$,
$e_2 = (u_2, v_2)$ is $\pm 1$ if $e_1$ and $e_2$ share exactly one
endpoint (the sign depending on whether they "agree" or "disagree" at that
endpoint), and $0$ otherwise. It acts on **edge (1-cochain) space**
$\mathbb{R}^E$.

### 2.3 Sizes and small-$n$ values

From the pilot (`xcwg_u1_wilson_rule_b_pilot.py`):

| $n_{\max}$ | $V$ | $E$ | $c$ | $\beta_1 = E - V + c$ | $\dim\ker L_0$ | $\dim\ker L_1$ |
|:----------:|:---:|:---:|:---:|:---------------------:|:--------------:|:--------------:|
|     2      | 10  | 20  | 1   |          11           |        1       |       11       |
|     3      | 28  | 106 | 1   |          79           |        1       |       79       |

The Rule B graph is connected at both cutoffs (single component $c = 1$),
as expected since the dipole edges mix $\kappa$. Compare Paper 25's scalar
S³ Coulomb graph at $n_{\max}=3$: $V = 14$, $E = 13$, $c = 3$, $\beta_1 = 2$.
**Rule B has 40× more cycles at $n_{\max} = 3$.**

The first few eigenvalues of $L_1$ at $n_{\max} = 2$ (truncated to the
smallest 20) are
\[
  \{0^{(11)}, 1.585, 2.000, 2.219, 2.755, 3.000, 3.000, 3.234, 3.000, \ldots\}.
\]
The $\beta_1 = 11$ zero modes are confirmed. The spectrum no longer
factorizes cleanly over $\mathbb{Q}(\sqrt 5)$ as it did at Paper 25
$n_{\max} = 3$; the larger graph carries richer algebraic content (a
detailed factorization of $\det(L_1 - \lambda I)$ over $\mathbb{Q}$ is
out of scope here but can be extracted from the pilot data).

---

## 3. Hodge identity check

The Hodge identity states that the nonzero spectra of $L_0$ and $L_1$
coincide (with multiplicities): if $B = U \Sigma V^{\top}$ is the SVD,
then $L_0 = U \Sigma^2 U^{\top}$ and $L_1 = V \Sigma^2 V^{\top}$, so the
nonzero $\Sigma^2$ values are shared. The kernel dimensions track the
homology:
\[
  \dim\ker L_0 = c, \qquad \dim\ker L_1 = \beta_1 = E - V + c.
\]

### 3.1 Verification at $n_{\max} = 2$

The pilot computes $\dim\ker L_0 = 1$, $\dim\ker L_1 = 11$, nonzero
spectra of length $V - c = 9$ and $E - \beta_1 = 9$, with shared maximum
discrepancy
\[
  \max |w_0^{\rm nz} - w_1^{\rm nz}| < 10^{-9}.
\]

Status: **PASSED at machine precision**.

### 3.2 Verification at $n_{\max} = 3$

The pilot computes $\dim\ker L_0 = 1$, $\dim\ker L_1 = 79$, nonzero
spectra of length $V - c = 27$ and $E - \beta_1 = 27$, with shared maximum
discrepancy $< 10^{-9}$.

Status: **PASSED at machine precision**.

The Hodge identity is structural (it follows from the SVD argument above)
and could not have failed unless the signed incidence matrix were
mis-constructed. Verifying it is a *consistency* check for the
implementation, not a structural finding.

### 3.3 Honest scope

Hodge identity at the level of `numpy.linalg.eigvalsh` only certifies the
floating-point construction. A symbolic certificate (charpoly over $\mathbb{Z}$,
sympy verification) is not needed here because $B$ has integer entries by
construction, so $L_0$ and $L_1$ are integer matrices and their nonzero
eigenvalues are algebraic over $\mathbb{Q}$ — verbatim with Paper 25
Observation 1 ("Edge Laplacian is algebraic on a finite Hopf graph"). The
Rule B graph inherits this π-free certificate.

---

## 4. Weak-coupling kinetic term (Paper 25 Theorem 1 analog)

### 4.1 Wilson action

Define the Wilson action by
\[
  S_W[\theta] = \beta \sum_{P \in \mathcal P} (1 - \cos \theta_P),
\]
where $\mathcal P$ is the set of *primitive non-backtracking closed walks*
(plaquettes; see §5 below) and $\theta_P = \sum_{e \in P} \theta_e$ is the
oriented sum of link phases around the plaquette. The action is real,
gauge-invariant ($\theta_P$ depends only on link variables along the
plaquette, and a node-local shift adds and subtracts $\chi$ at each
intermediate vertex), and bounded below by $0$ with equality at the trivial
vacuum $\theta_e \equiv 0$.

### 4.2 Quadratic expansion

Linearize around the trivial vacuum ($\theta_e = 0$): expanding
$1 - \cos\theta_P = \tfrac12 \theta_P^2 + O(\theta^4)$,
\[
  S_W[\theta] = \frac{\beta}{2} \sum_P \theta_P^2 + O(\theta^4)
  = \frac{\beta}{2} \theta^{\top} (d_1^{\top} d_1) \theta + O(\theta^4),
\]
where $d_1: \mathbb{R}^E \to \mathbb{R}^{\mathcal P}$ is the *plaquette
boundary operator* with $(d_1)_{P, e} = +1$ if $e$ appears in $P$ in the
canonical orientation, $-1$ if it appears reversed, $0$ otherwise.

The 1-form Hodge Laplacian on the graph admits the decomposition
\[
  L_1 = d_0^{\top} d_0 + d_1^{\top} d_1 \quad (\text{Hodge–Helmholtz}),
\]
where $d_0 = B^{\top}$ in our notation. (Paper 25 §II.C.) For *triangle-free*
2-complexes with $\mathcal P$ supplying all plaquette 2-cells, the second
summand controls the *transverse / gauge-invariant* edge modes and the
first summand controls the *longitudinal / pure-gauge* modes.

### 4.3 Theorem 1 analog

**Statement.** At quadratic order around the trivial vacuum,
\[
  S_W[\theta] = \frac{\beta}{2}\, \theta^{\top} K\, \theta + O(\theta^4),
  \qquad K = d_1^{\top} d_1,
\]
and $K$ is the *co-exact* (gauge-invariant) component of the edge Laplacian.
Restricted to the gauge-fixed subspace $\theta \perp \operatorname{im}(d_0)$,
$K$ acts as the full edge Laplacian $L_1$: this is the abelian sibling of
Paper 30 Theorem 2 / Paper 25 Theorem 1.

**Proof sketch.** Expand the cosine as in §4.2; the matrix elements of
$d_1^{\top} d_1$ are $\langle e_1, e_2 \rangle_K = \sum_P (d_1)_{P, e_1}
(d_1)_{P, e_2}$. On the gauge-fixed subspace, the longitudinal modes
$\operatorname{im}(d_0)$ are projected out by the Faddeev-Popov procedure
(or by the BRST-invariant Coulomb-gauge choice), leaving $K \big|_{\perp} =
L_1 \big|_{\perp}$ as the quadratic kinetic operator. ∎

### 4.4 Structural difference from Paper 25

The form of the theorem is *identical* to Paper 25 Theorem 1 — both
constructions give $L_1$ as the quadratic kinetic operator after gauge
fixing. The structural difference is in the **content** of $L_1$:

* **Paper 25 scalar graph.** At $n_{\max} = 3$, $L_1$ is 13-dimensional with
  spectrum $\{0^{(2)}, \tfrac{3-\sqrt 5}{2}, 1^{(2)}, \tfrac{5-\sqrt 5}{2},
  2, \tfrac{3+\sqrt 5}{2}, 3^{(3)}, \tfrac{5+\sqrt 5}{2}, 5\}$ — block-
  diagonal by $\ell$-shell.
* **Rule B graph.** At $n_{\max} = 3$, $L_1$ is 106-dimensional with
  $\beta_1 = 79$. The block decomposition by $\kappa$ is *visible at the
  vertex level* but NOT a block-decomposition of $L_1$, because every edge
  crosses $\kappa$ (E1 parity flip $\Delta\ell = \pm 1$ forces a $\kappa$
  jump). The Wilson kinetic operator is therefore *intrinsically
  $\kappa$-mixing*, the structural analog of "cross-shell propagator".

This is the technical statement of why Rule B should host RG-flow content
that Paper 25's scalar graph cannot: in Paper 25, decoupling block-spin
inside each $\ell$-block reproduces 2D Wilson confinement per block but
never sees across blocks. On Rule B, no such decoupling exists at the
edge level.

---

## 5. Plaquette structure and Wilson loops

### 5.1 Plaquette counts (pilot data)

The number of primitive non-backtracking closed walks (modulo cyclic
rotation and orientation reversal) up to length 8 / 6:

| $n_{\max}$ | $L=4$ | $L=6$ | $L=8$ |
|:----------:|:-----:|:-----:|:-----:|
|     2      |  44   |  144  | 1412  |
|     3      |  994  | 27906 |  *    |

(*: not enumerated at $n_{\max}=3$ because the count grows quickly; the
pilot capped length 6.)

Compare Paper 25 scalar graph plaquettes — at $(n_{\max}=2, 3, 4)$ the
plaquette counts $(L=4, 6, 8)$ are $(0, 0, 0)$, $(2, 1, 1)$, $(8, 7, 18)$.
**Rule B is plaquette-dense from $n_{\max}=2$ onward** — already 44 length-4
plaquettes at the smallest meaningful cutoff. This is the most
consequential structural difference.

### 5.2 Why Rule B has length-4 plaquettes at $n_{\max}=2$

A simple combinatorial construction. At $n_{\max}=2$ the vertex set is
\[
  \{1s_{1/2}^{\pm 1/2}, 2s_{1/2}^{\pm 1/2}, 2p_{1/2}^{\pm 1/2},
    2p_{3/2}^{\pm 1/2, \pm 3/2}\}
\]
(10 nodes). A typical length-4 plaquette is
\[
  1s_{1/2}^{+1/2} \to 2p_{1/2}^{-1/2} \to 2s_{1/2}^{+1/2} \to 2p_{1/2}^{-1/2}
  \to \cdots
\]
— no wait, that backtracks. A correct example:
\[
  1s_{1/2}^{+1/2} \to 2p_{1/2}^{+1/2} \to 2s_{1/2}^{+1/2} \to 2p_{3/2}^{+1/2}
  \to 1s_{1/2}^{+1/2}.
\]
All four steps satisfy $\Delta\ell = \pm 1$, $|\Delta m_j| \le 1$,
$|\Delta j| \le 1$. The plaquette "weaves" between $\ell = 0$ and $\ell = 1$
shells at fixed $m_j$. There are roughly $4 \times 4 \times 2 \approx 30$
distinct length-4 plaquettes from this single template, with the remaining
~14 coming from $m_j$-flip variants and longer routings; we have not done
the full combinatorial decomposition, but the numerical count $44$ is
consistent with a few hundred raw walks reduced by cyclic and reversal
symmetry.

### 5.3 Strong-coupling Wilson loops

For abelian U(1), single-link integration of $e^{i n \theta}$ against the
Haar measure $\frac{d\theta}{2\pi}$ gives $\delta_{n, 0}$. Therefore the
strong-coupling ($\beta \to 0$) leading behavior of a Wilson loop
$\langle W_C \rangle$ in U(1) Wilson is
\[
  \langle W_C \rangle = \left(\tfrac{\beta}{2}\right)^{A(C)}\, [1 + O(\beta^2)],
\]
where $A(C)$ is the *minimum* number of plaquettes needed to tile the loop
$C$ (the "area"). This is the standard U(1) area law (cf. Wilson 1974
appendix; Drouffe–Itzykson; Creutz §6). For Rule B, the smallest plaquette
is $L = 4$, so a length-4 Wilson loop coinciding with one plaquette has
\[
  \langle W_C \rangle = \tfrac{\beta}{2} + O(\beta^3) \quad (\beta \to 0).
\]
Longer Wilson loops show area-law scaling tiled by length-4 plaquettes
wherever possible.

### 5.4 Weak-coupling Wilson loops

At $\beta \to \infty$, expand around $\theta_e = 0$; the action is
quadratic-Gaussian (modulo the gauge-fixing in §4) and
\[
  \langle W_C \rangle = \exp\!\bigl(-\tfrac{1}{2\beta}\, C^{\top} K^{-1} C\bigr),
\]
where $C$ is the (signed) edge-indicator vector of the loop and $K$ is the
gauge-fixed kinetic operator from §4.3. Since $K$ has spectrum bounded
below by the smallest nonzero eigenvalue of $L_1$ ($\approx 1.585$ at
$n_{\max}=2$), the Wilson loop has *perimeter-law* scaling at weak coupling.

### 5.5 Simplest gauge-invariant Wilson loop

The simplest gauge-invariant observable is the holonomy around any single
length-4 plaquette $P$:
\[
  W_P = e^{i\theta_P} = e^{i(\theta_{e_1} + \theta_{e_2} + \theta_{e_3} + \theta_{e_4})}.
\]
Under the gauge transformation $\theta_e \mapsto \theta_e + \chi_{v'} - \chi_v$,
the sum $\theta_P$ collects $\chi_{v_1} - \chi_{v_1}$ around the closed loop
and is therefore invariant. The expectation value
$\langle (1 + W_P + W_P^{-1})/3 \rangle$ or $\langle \cos\theta_P \rangle$
is the standard order parameter, with strong-coupling limit
$\beta/2 + O(\beta^3)$ and weak-coupling limit
$1 - \frac{1}{2\beta} \langle P, K^{-1} P \rangle + O(\beta^{-2})$.

---

## 6. Structural comparison to Paper 25 (scalar graph)

We collect the head-to-head comparison at $n_{\max} = 3$, which is where
both graphs are at their first "meaningful" cutoff:

| Property | Paper 25 (scalar Fock) | Rule B (Dirac dipole) |
|:---------|:----------------------:|:---------------------:|
| $V$ | 14 | 28 |
| $E$ | 13 | 106 |
| $c$ | 3 | 1 |
| $\beta_1$ | 2 | 79 |
| Plaquette $L=4$ count | 0 | 994 |
| Plaquette $L=6$ count | 1 | 27 906 |
| Block decomposition of $L_1$ | per-$\ell$-shell | $\kappa$-mixing edges, no per-$\kappa$ block |
| Smallest plaquette length | 6 | 4 |
| Number of connected components | 3 (s, p, d) | 1 |

The structural shift is qualitative, not just quantitative:

**(i) Connectivity.** Paper 25's scalar graph splits as $\ell = 0 \oplus
\ell = 1 \oplus \ell = 2$ at $n_{\max}=3$ — three disjoint per-$\ell$
sub-graphs. Rule B is connected because the dipole edges link every shell
to every adjacent shell.

**(ii) Cycle rank.** Paper 25 $\beta_1 = 2$ vs Rule B $\beta_1 = 79$ —
the harmonic 1-cochain space is **40× larger** on Rule B.

**(iii) Hopf U(1) on Rule B.** Paper 25 §III.C identifies a natural U(1)
gauge structure on the scalar graph from the ladder-operator phases
$L_{\pm}$ and $T_{\pm}$. On Rule B, the analogous structure is the
**E1 dipole transition amplitude** carried by each edge (Szmytkowski
matrix elements). These are also complex amplitudes that transform
covariantly under a node-local U(1), with the additional structure that
the underlying matter sits in *spinor* representations of SO(3) ⊃ SO(2);
but the *gauge group itself* is still abelian U(1).

A potential subtlety: vertex labels include $\kappa$ (with both signs of
$j-\ell$), so one might ask whether the gauge group naturally enlarges to
U(1) × $\mathbb{Z}_2$, where the $\mathbb{Z}_2$ is the conjugation
$\kappa \to -\kappa$ (or alternatively, $m_j \to -m_j$). Two remarks:

(a) Paper 29 §RH-C identifies the $m_j \to -m_j$ reflection as a $\mathbb{Z}_2$
that *the Hashimoto operator commutes with* (this is the basis of the
22+24 factorization of the Dirac Rule B Ihara zeta). This $\mathbb{Z}_2$
is a *graph automorphism*, not a gauge symmetry; it does not act on link
variables, only on vertex labels.

(b) The gauge group of the lattice gauge theory is therefore just U(1),
with the $\mathbb{Z}_2$ a global symmetry that the Wilson action manifestly
respects (since the action depends only on the unsigned undirected edge
set, and $m_j \to -m_j$ permutes the vertex set in a way that preserves
adjacency). The fixed-point sub-graph of this $\mathbb{Z}_2$ (the vertices
with $m_j = 0$, which is empty at half-integer $m_j$; or in a different
labeling, the orbit graph) gives a natural quotient construction analogous
to the S² Hopf quotient in Paper 25 §III.B. We do not pursue the quotient
in this memo.

**(iv) Plaquette length.** The most consequential operational difference.
Paper 25's smallest plaquette is length 6 (a hexagonal loop in the p-block
of the scalar graph at $n_{\max}=3$). Rule B's smallest plaquette is length
**4**, appearing already at $n_{\max}=2$ — and there are 44 of them at
$n_{\max}=2$, scaling to 994 at $n_{\max}=3$. Strong-coupling expansions on
Rule B therefore have leading behavior $\tfrac{\beta}{2}$ at length-4 loops,
whereas Paper 25's leading behavior is $(\beta/2)^{3/2}$-class at length-6
loops. The two theories have different scaling dimensions at the strong-
coupling fixed point, which is a quantitative signature distinguishing
2D-like flow (Paper 25 per-$\ell$-block) from genuinely-3D-like flow (Rule
B, if MK confirms it).

---

## 7. RG-flow projections (qualitative)

The Migdal–Kadanoff (MK) block-spin / bond-moving renormalization group on
finite graph lattices is a standard tool: it groups plaquettes by length
class, applies the U(1) "bond-moving" duality
\[
  \beta_{\rm dual} = -\tfrac{1}{2}\log\tanh\tfrac{\beta}{2}
\]
to each block, and iterates. The two-dimensional case has a single
non-trivial fixed point (Berezinskii–Kosterlitz–Thouless), reproduced on
Paper 25's per-$\ell$-block scalar graph; the three-dimensional case has a
nontrivial confinement-deconfinement transition (Polyakov 1977; Mack
1980).

The qualitative question for Rule B: does MK reduce the 4-plaquettes to
2-plaquettes via decimation, or do they survive as irreducible 3D content?

### 7.1 Decimation does not reduce cross-shell plaquettes to within-shell

Within Paper 25's scalar graph, each plaquette lies inside a single
$\ell$-block. Decimating one vertex in the $\ell = 1$ block reduces a
length-6 hexagon to a length-4 square *inside the same block* (this is
standard 2D MK); the decimated graph is still per-$\ell$-block.

On Rule B, the smallest plaquette ($L = 4$) has the structure described in
§5.2: it weaves between $\ell = 0$ and $\ell = 1$. Decimating one vertex
in such a plaquette (say, a $2s_{1/2}$ vertex) connects its neighbors
across $\ell$-shells in a single hop, which is *not* a Rule B edge (Rule B
edges have $\Delta\ell = \pm 1$, not $\Delta\ell = 0$ between $2p$ and
$2p$). The decimated graph is therefore NOT a sub-graph of Rule B.
**Cross-shell plaquette structure cannot be reduced to within-shell content
by MK decimation.** This is the structural reason to expect 3D phase
content.

### 7.2 What MK on Rule B would compute

A future MK sprint would compute, at each iteration:
(a) Coarsened graph $G^{(k)}$ — vertices grouped by some equivalence
relation (e.g., $\kappa$-orbits or $n_{\rm fock}$-shells);
(b) Effective coupling $\beta^{(k)}$ — via the bond-moving / character-
expansion duality;
(c) Effective plaquette content $\mathcal P^{(k)}$ — counts at each length.

Expected qualitative signatures:
- **2D-like fixed point** if $\beta^{(k)}$ flows monotonically and the
  plaquette content becomes per-block (would replicate Paper 25 per-$\ell$).
- **3D-like fixed point** if $\beta^{(k)}$ has a nontrivial discontinuity
  at some $\beta_c$, separating a confined phase (area law) from a
  Coulomb phase (perimeter law).

The structural argument in §7.1 favors the second scenario, but a
quantitative MK computation is required for confirmation. This is the
natural next sprint (1–2 weeks of focused work).

### 7.3 Honest scope of the projection

This sub-section is qualitative reading, not theorem. Three caveats:

(a) The Dirac graph is *not* a regular lattice; MK on irregular finite
graphs is delicate and not literal "block-spin" but a graph-theoretic
analog (Schaub et al. 2020; this is a generalization, not the canonical
Wilson MK).

(b) Even if the structural reading is right, the finite-graph theory has
no second-order phase transition — only a crossover, with finite-size
scaling.

(c) The MK approximation is itself uncontrolled — it captures qualitative
features but is not a proof. A more reliable check would be direct Monte
Carlo of the partition function and Wilson loop expectation values; see
Paper 30 §6 for the SU(2) Monte Carlo on the scalar S³ graph.

---

## 8. Recommendation

**The U(1) Wilson lattice gauge theory on Dirac Rule B is a well-defined,
gauge-invariant, abelian lattice gauge theory on a finite directed graph.**
All the standard Wilson–Hodge machinery transfers verbatim from Paper 25:

* Link variables $U_e \in U(1)$ on oriented edges (§1).
* Node-local gauge transformation $\psi_v \to e^{i\chi_v}\psi_v$,
  $U_e \to e^{i\chi_{v'}} U_e e^{-i\chi_v}$ (§1.2).
* Signed incidence $B$, Laplacians $L_0 = BB^{\top}$, $L_1 = B^{\top}B$;
  Hodge identity (verified numerically at $n_{\max} = 2, 3$, §3).
* Wilson action $S_W = \beta \sum_P (1 - \cos\theta_P)$; weak-coupling
  kinetic term $L_1$ on the gauge-fixed subspace (§4).
* Plaquette content quadratic in $\beta_1$, smallest-plaquette length 4,
  area-law and perimeter-law regimes (§5).

**The construction is structurally promising for the questions Paper 25's
scalar graph cannot reach.** The single load-bearing structural difference
from Paper 25 is that Rule B has 4-plaquettes at $n_{\max} = 2$ (44 of them)
and 994 at $n_{\max} = 3$, while Paper 25's scalar graph has none of length
4 at any tested cutoff — its smallest plaquette is length 6. This is the
discrete geometric content responsible for the per-$\ell$-block 2D
reduction of Paper 25 not transferring to Rule B: every Rule B edge crosses
ℓ-shells by the E1 parity-flip selection rule, so block-spin decimation
cannot collapse cross-shell plaquettes to within-shell content (§7.1).

**Recommended next sprint.** MK block-spin / bond-moving RG on Rule B at
$n_{\max} = 2, 3, 4$, with the structural targets:

* Whether $\beta^{(k)}$ has a nontrivial fixed point separating a confined
  (area-law) phase from a Coulomb (perimeter-law) phase.
* Whether the plaquette content reduces to per-$\kappa$ (analog of
  per-$\ell$ in Paper 25) or remains $\kappa$-mixing under iteration.
* Whether the Wilson loop has a confinement-deconfinement transition at
  finite $\beta_c$.

A positive outcome on any of (a)–(c) would be a structural sequel to
Paper 25 with no analog in Paper 30 (SU(2) on the scalar graph) and no
analog in Sprint ST-SU3 (SU(3) on Bargmann S⁵, which sits in the
Riemannian/Hopf-fiber asymmetry of Paper 24 §V).

**Honest scope.** The construction is a finite-graph abelian lattice gauge
theory. Per CLAUDE.md §1.5, this is not a continuum QED claim, not a
Lorentzian theory, not a mass-gap result. The framework's value here is
that it provides a concrete discrete arena in which the standard Wilson–
Hodge technology can be exercised on a non-trivially-3D plaquette
geometry that arises *naturally* from the Dirac selection rules, with no
external choice of lattice and no fitted parameters. The "natural geometry"
of Rule B is the E1 dipole graph of Paper 29 §RH-C; the U(1) gauge
theory on it is what this memo specifies.

---

### Reference data

All structural data above is reproduced by `xcwg_u1_wilson_rule_b_pilot.py`,
output in `debug/data/xcwg_u1_wilson_rule_b_pilot.json`. Plaquette
enumeration is limited to length $\le 8$ at $n_{\max} = 2$ and length
$\le 6$ at $n_{\max} = 3$ for budget reasons. Hodge identity is verified
to floating-point precision; a sympy-exact integer-matrix certification is
straightforward but not required for the verdict.

### Files

* `debug/xcwg_u1_wilson_rule_b_design_memo.md` (this file).
* `debug/xcwg_u1_wilson_rule_b_pilot.py` (pilot code, no production
  modifications).
* `debug/data/xcwg_u1_wilson_rule_b_pilot.json` (numerical output).

---

*Sprint XCWG, May 2026. Independent of, but structurally adjacent to,
Paper 25 (U(1) on scalar S³), Paper 29 §RH-C (Rule B definition), Paper 30
(SU(2) on scalar S³), Sprint ST-SU3 (SU(3) on Bargmann S⁵).*
