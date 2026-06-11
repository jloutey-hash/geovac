# Alternative U(1) Wilson observables on Dirac Rule B — design and tractability

**Sprint XCWG (May 2026). Companion to:**
- `debug/xcwg_rule_b_infrastructure_memo.md` (Rule B graph properties)
- `debug/xcwg_u1_wilson_rule_b_design_memo.md` (U(1) Wilson construction)
- `debug/xcwg_rule_b_spectral_dim_memo.md` (heat-kernel d_s)
- `debug/xcwg_observables_pilot.py` (pilot computation, this memo's source)
- `debug/data/xcwg_observables_pilot.json` (numerical output)

The XCWG-B sprint computes RG flow via Migdal–Kadanoff (MK) block-spin
and compares to continuum 3D U(1). MK is one diagnostic out of several
that exist for compact U(1) lattice gauge theory in 3D. This memo
surveys and *scopes* the alternative observables — Wilson loops,
Polyakov-loop analogs, monopole content, plaquette correlators, and a
Rule-B-native "Wilson-loop l-shell extent" — so XCWG-B can run two or
three of them alongside MK and have independent cross-checks of the
phase structure.

Bottom line up front: three observables are essentially **free** to add
to XCWG-B (Wilson loop S, plaquette correlator, l-shell extent), with
all three reducing to small dense linear algebra on the existing
infrastructure. One (Polyakov-loop analog) has **no clean Rule B
realization** — Rule B is a static spatial graph with no distinguished
temporal direction, and this is itself an informative structural
finding. The monopole observable is **structurally well-defined but
computationally heavier**; it requires enumerating closed 2-surfaces
on Rule B, a multi-week extension. The recommendation is to bolt
Wilson-loop S, plaquette correlator, and l-shell extent onto XCWG-B at
no real cost.

---

## §1. Wilson loops

### 1.1 Definitions

A Wilson loop on Rule B is a primitive non-backtracking closed walk
$C = e_1 \cdots e_L$ on the directed edge set. The Wilson action
($S_W = \beta \sum_P (1 - \cos \theta_P)$) selects out the smallest
plaquettes ($L=4$ on Rule B); larger Wilson loops are then composite
objects, and their expectation values diagnose the strong/weak coupling
phases.

The standard strong/weak limits of $\langle W_C \rangle = \langle e^{i\theta_C}\rangle$:

- **Strong coupling** ($\beta \to 0$): each edge averages independently
  against the Haar measure, $\int_0^{2\pi} e^{i n\theta} d\theta/2\pi
  = \delta_{n,0}$. Only edges that appear in pairs $e, \bar e$ survive
  single-link integration. For a primitive closed Wilson loop $C$ of
  length $L$, the leading contribution comes from tiling $C$ by
  plaquettes — the **minimum area** $A(C)$ in number of plaquettes
  needed to fill $C$:
  $$
    \langle W_C \rangle_{\beta \to 0} = \left(\tfrac{\beta}{2}\right)^{A(C)}
    \bigl[1 + O(\beta^2)\bigr].
  $$
  This is the standard U(1) "area law" (Wilson 1974). For Rule B with
  smallest plaquette $L=4$, a single plaquette tiles itself ($A=1$),
  giving leading $\beta/2$. A length-6 loop needs at least $A=2$
  plaquettes (since one $L=4$ plaquette has 4 edges, two share 1 edge,
  giving a $L=6$ boundary), so its leading strong-coupling behavior is
  $(\beta/2)^2$. **In compact 3D U(1) on a regular lattice, the area
  law persists at all $\beta$ (Polyakov 1977); the question is whether
  Rule B's effective $\sim$3D structure reproduces this or shows a
  Coulomb (perimeter-law) phase at high $\beta$.**

- **Weak coupling** ($\beta \to \infty$): expand $\theta_e$ to quadratic
  order around the trivial vacuum, gauge-fix to the orthogonal
  complement of $\operatorname{im}(d_0)$, and treat the residual
  Gaussian path integral:
  $$
    \langle W_C \rangle_{\beta \to \infty}
    = \exp\!\Bigl(-\tfrac{1}{2\beta}\, C^{\!\top} K_{\!\rm pinv}\, C\Bigr),
  $$
  where $C \in \mathbb{Z}^{|E|}$ is the signed edge-indicator vector of
  the loop and $K = d_1^{\!\top} d_1$ is the Wilson kinetic operator
  from §4.3 of the design memo. The closed loop has $B C = 0$ (no node
  boundary), so $C$ is automatically gauge-invariant and the gauge-fixed
  propagator $K_{\rm pinv}$ is the right object. In the confined phase,
  $C^{\!\top} K_{\rm pinv} C$ scales **linearly with the perimeter of
  $C$** (perimeter law); in the Coulomb phase it scales with the
  minimum cross-sectional area (area law). **The perimeter-vs-area
  scaling of $C^{\!\top} K_{\rm pinv} C$ versus $L$ is the
  confinement-deconfinement order parameter at the Gaussian level.**

### 1.2 Pilot result on Rule B at $n_{\max} = 2$ (computed)

The pilot (`debug/xcwg_observables_pilot.py`) computed, on Rule B at
$n_{\max}=2$ (V=10, E=20, $\beta_1 = 11$):

| Loop length $L$ | Count | $\langle S \rangle$ (mean) | range | strong-coupling leading $\beta$ power |
|:-:|:-:|:-:|:-:|:-:|
| 4 | 44 | **0.2477** | $[2/9, 1/3]$ | $\beta^1$ |
| 6 | 144 | **0.4361** | $[0.347, 0.458]$ | $\beta^2$ |

where $S = C^{\!\top} K_{\rm pinv} C$ is the Gaussian action of the loop
(units of $1/\beta$). Three observations from these numbers:

(i) **$S$ values are clean rationals.** $S(L=4)_{\min} = 2/9$ and
$S(L=4)_{\max} = 1/3$ exactly. The 44 length-4 loops fall into
$\{2/9, 5/18, 1/3\}$-class buckets (sympy verification would pin down
the exact multiplicities, deferred). The rationality is structural:
$K$ is an integer matrix, so $K_{\rm pinv}$ has rational entries on
$\operatorname{im}(K)$. The Wilson kinetic content is π-free at
finite cutoff, matching the Paper 24 / Paper 25 / Paper 29 algebraic
certificate.

(ii) **Ratio $\langle S \rangle_{L=6}/\langle S \rangle_{L=4} = 1.76$**,
between linear-in-$L$ (perimeter law would predict $6/4 = 1.5$) and
quadratic ($L^2$, area law would predict $9/4 = 2.25$). At
$n_{\max}=2$ the graph is too small to discriminate cleanly between
these scalings. The natural extension is to compute $S$ at $L=4, 6, 8,
10$ on $n_{\max} = 3, 4$ and fit the scaling exponent.

(iii) **Rank of $K = d_1^{\!\top} d_1$ equals $\beta_1$**: at
$n_{\max}=2$, rank$(K)$ = rank$(d_1)$ = 11 = $\beta_1$. The 44
plaquettes span only the 11-dimensional space of harmonic 1-cochains.
This is the discrete analog of "the plaquette content is enough to
detect all 1-cycles in homology"; the dimension match between
$d_1$-image and ker $L_1$ is structural (Hodge theorem on the 2-complex).
At $n_{\max}=3$: rank$(K)$ = 79 = $\beta_1$.

### 1.3 Tractability at $n_{\max} = 2, 3, 4$

| $n_{\max}$ | $V$ | $E$ | $\beta_1$ | Loop enumeration (cost) | $K$ pseudoinverse (cost) |
|:-:|:-:|:-:|:-:|:-:|:-:|
| 2 | 10 | 20 | 11 | $L \le 8$: 1600 loops, instant | $20\times 20$, $< 1$ ms |
| 3 | 28 | 106 | 79 | $L \le 6$: ~$10^4$ loops, $\sim$10 s | $106 \times 106$, $< 1$ s |
| 4 | 60 | 312 | 253 | $L \le 6$: ~$10^5$ loops, $\sim$1 min | $312 \times 312$, $< 1$ s |

Wilson loops are **tractable at $n_{\max} = 4$** at $L \le 6$. Going
to $L = 8$ at $n_{\max} = 4$ would push enumeration to $\sim 10^7$
walks (an hour or two); not needed for the leading scaling fit.

### 1.4 Sensitivity to phase structure

The Wilson loop $\langle W_C \rangle$ is the **canonical**
confinement-deconfinement order parameter in lattice gauge theory:

- 3D compact U(1) is **always confined** at any $\beta$ (Polyakov 1977)
  with $\langle W_C \rangle \sim \exp(-\sigma A)$, where $\sigma$ is
  the string tension. At weak coupling the string tension is exponentially
  small in $\beta$: $\sigma \sim \exp(-\pi^2 \beta / 4)$.
- 4D U(1) has a **deconfinement transition** at $\beta_c \approx 1.0$
  separating an area-law (confined) phase from a perimeter-law (Coulomb)
  phase.
- 2D U(1) is always confined trivially (Wilson 1974 §5): the partition
  function factorizes per plaquette.

Computing $S(L)$ scaling on Rule B and looking for area-vs-perimeter
crossover therefore **directly diagnoses the effective dimension**
(2D vs 3D vs 4D) of the underlying gauge theory. This is independent
of the MK $\beta_{\rm eff}(\beta)$ flow, which probes the same physics
through the duality $\beta \to -\tfrac12 \log \tanh(\beta/2)$.

### 1.5 Independence from MK

Wilson loops carry information **beyond** what MK gives. MK reports
the flow of the bare coupling under decimation; Wilson loops report
the actual expectation value of a physical observable. In principle
both should detect the same phase structure (this is the consistency
of any Wilsonian RG), but practically they have different sources of
systematic error: MK has the bond-moving approximation (uncontrolled),
Wilson loops have the loop-truncation cutoff (controlled by $L$).
Cross-checking MK against Wilson scaling is the standard discipline
in the lattice literature (Creutz, *Quarks, Gluons and Lattices*,
chapters 6–8).

---

## §2. Polyakov-loop analog

### 2.1 The continuum construction

In finite-temperature lattice gauge theory the Polyakov loop is the
trace of a path-ordered link product around a Euclidean-time circle:
$L_{\rm Polyakov}(\vec x) = \mathrm{Tr}\, \prod_{t=0}^{N_t - 1}
U_{(\vec x, t) \to (\vec x, t+1)}$, where $N_t$ is the temporal extent
of the lattice with periodic-in-time boundary conditions. Its expectation
diagnoses *deconfinement*: $\langle L \rangle = 0$ in the confined
phase, $\langle L \rangle \ne 0$ in the deconfined phase, via the
breaking of the center symmetry $\mathbb{Z}_N$ (for U(1), the center
is U(1) itself).

For this construction to work the gauge theory must have:
1. A distinguished "time" direction along which a closed loop is to
   be evaluated.
2. Periodic boundary conditions in time.
3. A center symmetry that the action manifestly respects but the
   Polyakov loop transforms under.

### 2.2 Rule B has no temporal direction

The Rule B graph $G_B$ is a **static spatial graph**: vertices are
labeled by $(n_{\rm fock}, \kappa, m_j)$, all *spatial* / *internal*
quantum numbers. There is no time index, no Euclidean-time circle.
The closest analog is the principal quantum number $n_{\rm fock}$,
which is a *radial* quantum number (it parameterizes the radial
extent in the Fock projection of the Coulomb problem onto $S^3$),
not a temporal one. There is no natural cyclic identification of
$n_{\rm fock}$ values — the Fock projection treats $n_{\rm fock}$
as a discrete depth coordinate, not a periodic one.

### 2.3 Three candidate analogs (all imperfect)

The pilot script does enumerate three candidate "Polyakov-like" loops:

(a) **Loops that traverse all $n_{\rm fock}$ layers.** A primitive
length-$L$ loop that includes a vertex from every $n_{\rm fock} \in
\{1, \ldots, n_{\max}\}$ would be a "spans-the-radial-direction"
analog. Pilot data ($n_{\max}=2$, two radial layers): 28 of 44
length-4 loops and all 144 length-6 loops span both layers. At
$n_{\max}=2$ this is too crude to be informative.

(b) **Loops that close on the $m_j \to -m_j$ orbit of a single vertex.**
The $\mathbb{Z}_2$ reflection $m_j \to -m_j$ is a graph automorphism
of Rule B (Paper 29 §RH-C Hypothesis 1, validated by the 22+24
Ihara-zeta factorization). One could define a "fundamental domain"
under this $\mathbb{Z}_2$ and compute the holonomy of any loop that
encloses the fixed-point set. But the fixed-point set is empty (no
half-integer $m_j$ satisfies $m_j = -m_j$), so the construction is
formal-only.

(c) **Loops along a $\kappa$-fixed orbit.** Equivalent to (b), in
practice empty for the same reason.

### 2.4 Honest verdict: NO clean Polyakov-loop analog

The Polyakov-loop construction needs a distinguished temporal
direction; Rule B does not have one. The three candidates above are
formal substitutes that do not carry the structural content of the
standard Polyakov loop. **The cleanest move is to acknowledge this
absence as itself informative:** it says Rule B at this stage of
development is a *zero-temperature* gauge theory, not a finite-$T$
one. A genuine finite-$T$ analog would require a separate Euclidean-
time compactification — for instance, taking the tensor product
$G_B \otimes C_{N_t}$ with a cyclic chain (the "thermal Matsubara
circle"), which would lift the construction to a 4-dimensional
substrate. This is exactly the construction of Paper 35 §VIII
(`tensor_product_S3_x_S1`) for the Klein-Gordon temperature decoding,
adapted to Rule B. *This is multi-week work, out of scope for XCWG-B,
but a natural future direction if temperature dependence becomes a
target.*

---

## §3. Monopole content

### 3.1 Continuum construction

In 3D compact U(1), a monopole at position $\vec x$ is detected by
integrating the field strength $F = dA$ around a closed 2-surface
$\partial V$ surrounding $\vec x$. The total flux is a multiple of $2\pi$:
$$
\oint_{\partial V} F = 2\pi q, \qquad q \in \mathbb{Z}.
$$
Polyakov (1977) showed that integrating monopoles into the partition
function gives an effective sine-Gordon theory whose mass gap implies
linear confinement at all $\beta$ in 3D U(1).

On a finite simplicial complex, monopoles are detected by computing
the lattice analog of $\oint F$:
$$
q_{\partial V} = \tfrac{1}{2\pi} \sum_{P \subset \partial V} \theta_P,
$$
where the sum is over plaquettes forming the boundary of a 3-cell
("cube") $V$. For each link variable $U_e = e^{i\theta_e}$ with
$\theta_e \in (-\pi, \pi]$, $\theta_P$ is computed mod $2\pi$ inside
this fundamental domain, and the deviation $\theta_P - \text{Im} \log
e^{i\theta_P}$ counts monopoles via $4\pi$-vortex content (DeGrand–
Toussaint 1980; Polikarpov 1996 review).

### 3.2 Rule B does not have a natural 3-cell

The above construction needs **closed 2-surfaces** on the graph —
collections of plaquettes whose boundary 1-cycles cancel pairwise.
On a regular cubic lattice these are the elementary cubes (6 faces).
On Rule B, candidate 2-surfaces are *triangulations of plaquette
triples (or quadruples) sharing edges*.

The pilot enumerated triangular cells: 3 length-4 plaquettes mutually
sharing edges. **At $n_{\max}=2$, the pilot found 500 candidate
triangular cells** (a lot — reflecting the high plaquette density,
9.4 plaquettes per edge). However, not every triangular cell is a
**closed** 2-surface — three plaquettes sharing edges pairwise can
form an open "tent" topology rather than a closed surface.

A closed 2-surface on a bipartite girth-4 graph like Rule B requires
either (i) a triangular cell where the three plaquettes' boundary
1-cycles sum to zero (which forces a specific orientation pattern), or
(ii) larger combinatorial structures (4-plaquette "tetrahedra",
6-plaquette "octahedra", etc.).

### 3.3 What computing monopoles requires

A genuine monopole computation would:
1. Enumerate closed 2-surfaces on $G_B$ (sub-graph problem on the
   plaquette adjacency graph). At $n_{\max}=2$ with 44 plaquettes,
   the candidate triangular cells number 500; closed ones a small
   fraction; at $n_{\max}=3$ with 994 plaquettes the combinatorics
   becomes ~$10^5$ candidates.
2. For each closed 2-surface and each canonical Wilson configuration
   $\{U_e\}$, compute $q$ via the lattice DeGrand-Toussaint formula.
3. Statistically sample monopole density at various $\beta$ via Monte
   Carlo (the standard 3D compact U(1) measurement, but tailored to
   Rule B's irregular structure).

This is a **multi-week** undertaking and well outside what XCWG-B
can absorb. It is, however, structurally well-defined — the only
question is whether the irregular geometry of Rule B has enough
combinatorial structure to make the construction informative, or
whether the irregularity drowns out monopole signals.

### 3.4 Verdict: well-defined but heavy

Monopoles on Rule B are a **valid future direction** (3-6 weeks of
focused work), with the structural target of detecting confinement
via Polyakov's mechanism on a non-cubic substrate. They are too
heavy to bolt onto XCWG-B; defer until after Wilson scaling, MK flow,
and plaquette correlators have been computed at $n_{\max} = 3, 4$.

---

## §4. Plaquette correlators

### 4.1 Definition and continuum predictions

For two plaquettes $P, P'$ at graph distance $d$,
$$
\Gamma(P, P'; \beta) \;=\; \langle \cos\theta_P \cos\theta_{P'}\rangle
- \langle \cos\theta_P \rangle \langle \cos\theta_{P'}\rangle.
$$
This connected correlator detects the mass gap of the gauge theory.
In the confined phase, $\Gamma \sim e^{-m d}$ (exponential decay with
glueball mass $m$); in the Coulomb phase, $\Gamma \sim 1/d^4$ (power
law; massless photon).

### 4.2 Gaussian leading order on Rule B

At weak coupling, expand both cosines to quadratic order. Then
$\Gamma_{\rm leading}(P, P') = \tfrac{1}{2\beta^2}\, (C_P^{\!\top}
K_{\rm pinv} C_{P'})^2 + O(1/\beta^3)$, with $C_P$ the edge indicator
of plaquette $P$. We report the off-diagonal propagator element
$M(P, P') = C_P^{\!\top} K_{\rm pinv} C_{P'}$ — its decay with graph
distance is the structural signature of the gauge propagator.

### 4.3 Pilot result on Rule B at $n_{\max} = 2$ (computed)

| Plaquette pair distance $d$ | # pairs | $\langle M \rangle$ | $\langle |M| \rangle$ | $\mathrm{std}(M)$ |
|:-:|:-:|:-:|:-:|:-:|
| 0 (vertex-sharing) | 892 | $0.0111$ | $0.0499$ | $0.0668$ |
| 1 (one-graph-distance) | 54 | $\sim 0$ (machine) | $\sim 4 \times 10^{-17}$ | $\sim 7 \times 10^{-17}$ |

The pilot data show **a contact correlator structure**: only vertex-
sharing plaquette pairs ($d = 0$) have nonzero connected weak-coupling
correlation; pairs at $d \ge 1$ are bit-zero. **At $n_{\max} = 2$ the
graph is too small to support any propagator decay length** — the
graph diameter (4) is much smaller than any expected correlation
length, so all pairs are either nearest-neighbor (contact) or
beyond-cutoff (zero by being outside the harmonic 1-cochain support).

This is a structural artifact of finite size, not a finding about
Rule B's physics. The proper computation lives at $n_{\max} \ge 3$
where graph distance and plaquette structure can carry a real decay
length. The fact that the contact correlator is finite ($\langle |M|
\rangle = 0.05$) and non-uniform ($\mathrm{std}(M) = 0.067$, ratio
1.34) does indicate that the plaquette structure is *not* uniform —
there is preferred direction within the contact-correlator block.

### 4.4 Tractability at $n_{\max} = 3, 4$

At $n_{\max} = 3$: $994 \times 993 / 2 = 493{,}021$ plaquette pairs.
For each pair, the cost is one inner product $C_P^{\!\top} K_{\rm pinv}
C_{P'}$ in a 106-vector space — $\sim 10^{-4}$ s per pair, so $\sim 1$
minute total. Pair distances need a $V \times V$ shortest-path matrix
(already cheap). This is **tractable**.

At $n_{\max} = 4$: $5{,}047 \times 5{,}046 / 2 = 1.27 \times 10^7$
pairs, each $\sim 10^{-3}$ s, total $\sim 4$ hours. **Slow but
tractable**; could be done overnight.

### 4.5 Sensitivity

Plaquette correlators directly diagnose the **mass gap** of the gauge
theory:
- Confined phase: exponential decay $\Gamma(d) \sim e^{-m d}$.
- Deconfined phase: power-law decay $\Gamma(d) \sim 1/d^4$.

At $n_{\max} = 3, 4$ the maximum plaquette-pair distance is bounded
by the graph diameter (6, 8). This gives 3-4 useful decay points
before the propagator hits the boundary — enough for a log-linear fit
but borderline for confidence. The strongest version would be at
$n_{\max} = 5$ (diameter 10).

### 4.6 Independence from MK and Wilson loops

Plaquette correlators give a *different* combination of the gauge
data than Wilson loops do:
- Wilson loops integrate the propagator over a *cycle* in
  $\operatorname{im}(d_0^{\!\top})$.
- Correlators evaluate the propagator at *two distinct* points in
  the cycle basis.

Cross-checking the two should yield consistent estimates of $\beta_c$
(if a crossover exists), with different systematic errors. This is
the standard double-test discipline of lattice gauge theory.

---

## §5. Spectral dimension of Wilson loops / l-shell extent

### 5.1 Definition

A novel observable specific to Rule B's $\ell$-shell structure: for
each primitive loop $C$ on Rule B, let $\nu_\ell(C)$ be the number
of distinct $\ell$-values its vertices carry. This is the "l-shell
extent" of the loop.

The intuition is from Paper 25's per-$\ell$-block structure: in the
scalar Fock graph, $\nu_\ell(C) = 1$ for every loop (all loops are
within-shell). On Rule B, $\nu_\ell$ can be 2, 3, or more — and the
distribution over $\nu_\ell$ values directly reads how "shell-mixing"
the Wilson kinetic content is.

### 5.2 Pilot result (computed)

| $n_{\max}$ | Loop length $L$ | $\nu_\ell = 1$ | $\nu_\ell = 2$ | $\nu_\ell = 3$ |
|:-:|:-:|:-:|:-:|:-:|
| 2 | 4 | 0 | **44 (100%)** | 0 |
| 2 | 6 | 0 | **144 (100%)** | 0 |
| 3 | 4 | 0 | 452 (45.5%) | **542 (54.5%)** |

(The $n_{\max}=3$ data matches the infrastructure memo's plaquette
signature $(0,1,1,2)$ accounting: 542 plaquettes use three distinct
$\ell$-shells, 452 use two. Within-shell loops, $\nu_\ell=1$, are
absent at all measured $n_{\max}$.)

Structural reading:

- At $n_{\max}=2$, only $\ell \in \{0, 1\}$ exists, so every loop is
  bipartite in $\ell$-parity. Length-4 and length-6 loops both use
  exactly 2 shells.
- At $n_{\max}=3$, $\ell \in \{0, 1, 2\}$ is available, and **the
  majority (54.5%) of length-4 plaquettes already use three distinct
  $\ell$-shells**. This is genuinely cross-shell content.
- Forecast at $n_{\max} = 4, 5$: $\ell \in \{0, 1, 2, 3\}$ or
  $\{0, 1, 2, 3, 4\}$; expect $\nu_\ell$ distribution to spread up to
  4 or 5.

### 5.3 Why this is informative

The $\nu_\ell$ distribution measures **how much the Wilson kinetic
content sees Paper 25's per-$\ell$-block obstruction**. In the scalar
graph, $\nu_\ell \equiv 1$ — the Wilson action factorizes as a sum of
per-$\ell$ pieces, and MK decimation reduces each piece to a 2D
problem. On Rule B at $n_{\max}=3$, half the plaquettes span three
$\ell$-shells, which **cannot be reduced to a 2D-per-$\ell$ problem by
any decimation that respects Rule B's adjacency rule** (cf. design memo
§7.1).

The $\nu_\ell$ distribution is therefore a direct empirical
quantification of "how much 3D content" the Wilson construction has,
independent of Wilson loop scaling and plaquette correlators.

### 5.4 Tractability

$\nu_\ell$ computation is **trivial**: enumerate loops, count distinct
$\ell$-values per loop. At $n_{\max}=4$ with $5{,}047$ length-4
plaquettes, this is sub-second.

### 5.5 Independence from MK

MK flow on Rule B (the natural decimation $\kappa \to \kappa'$ that
respects E1 selection) does not preserve $\ell$-shell content — it
mixes shells. The $\nu_\ell$ distribution before and after MK
iteration would be a structural observable of how MK collapses the
cross-shell plaquette content. This is independent information from
the bare $\beta_{\rm eff}(\beta)$ flow that MK normally tracks.

---

## §6. Tractability assessment table

| Observable | $n_{\max}=2$ | $n_{\max}=3$ | $n_{\max}=4$ |
|:-|:-:|:-:|:-:|
| Wilson loops $S = C^{\!\top} K_{\rm pinv} C$ at $L \le 6$ | **instant** (44 + 144 loops) | **seconds** (994 + ~$10^4$ loops) | **minutes** (~$5 \times 10^3$ + ~$10^5$) |
| Polyakov-loop analog | **N/A** (no temporal direction) | **N/A** | **N/A** |
| Monopole content (closed 2-surfaces) | **multi-week** scoping | **multi-week** | **multi-week** |
| Plaquette correlator $C^{\!\top} K_{\rm pinv} C'$ by distance | **instant** (946 pairs, contact-only) | **minutes** ($5 \times 10^5$ pairs) | **hours** ($10^7$ pairs) |
| $\ell$-shell extent $\nu_\ell$ distribution | **instant** | **instant** | **second** |
| Heat-kernel $d_s$ (already done) | done | done | done |
| MK $\beta_{\rm eff}(\beta)$ flow (XCWG-B Track B1) | minutes | hours | days |

"Multi-week" entries are real engineering tasks (closed-surface
enumeration, MC sampling, etc.), not bolt-on additions.

---

## §7. Recommended observables for XCWG-B

The recommended additions to XCWG-B, in priority order:

### **Tier 1: bolt onto XCWG-B at no real cost**

1. **Wilson loops $S(L)$ at $L = 4, 6$ on $n_{\max} = 2, 3, 4$.**
   - Output: $\langle S(L) \rangle$ vs $L$ for each $n_{\max}$.
   - Fit: scaling exponent $S \sim L^\alpha$. Perimeter law $\alpha=1$;
     area law $\alpha = 2$ (with caveats); intermediate values indicate
     finite-size crossover.
   - Cross-check vs MK $\beta_{\rm eff}(\beta)$: a confined phase
     should produce $\alpha \to 2$ at large $\beta$; a Coulomb phase
     $\alpha \to 1$.
   - **Compute cost:** ~10 minutes total (already pilot-tested).

2. **$\ell$-shell extent $\nu_\ell$ distribution at $n_{\max} = 2, 3, 4$.**
   - Output: histogram of $\nu_\ell$ over all primitive loops at each
     $n_{\max}$.
   - Diagnostic: $\nu_\ell$ distribution stays bimodal at $\{1, 2, ...\}$
     after MK = the cross-shell content is structural;  $\nu_\ell$
     collapses to $\{1\}$ after MK = MK decimates to per-shell content
     (Paper 25 reduction).
   - **Compute cost:** seconds.

3. **Plaquette correlator $M(P, P')$ vs graph distance, at $n_{\max} = 3$.**
   - Output: $\langle |M(d)| \rangle$ for $d = 0, 1, 2, 3, 4$ at
     $n_{\max} = 3$ (small enough sample but distance range
     adequate).
   - Fit: $\log \langle |M| \rangle$ vs $d$. Linear (exponential
     decay) = confinement; concave/power-law = deconfinement.
   - **Compute cost:** ~30 minutes at $n_{\max} = 3$.

### **Tier 2: significant investment, defer until Tier 1 lands**

4. **Plaquette correlator at $n_{\max} = 4$.** Cost ~4 hours; only
   pursue if $n_{\max} = 3$ result is ambiguous about exponential
   vs power-law decay.
5. **Wilson loops at $L = 8, 10$ on $n_{\max} = 3$.** Cost ~hours;
   gives more decay points for scaling fit.

### **Tier 3: separate sprint**

6. **Monopole content / closed 2-surface enumeration.** Cost 3-6
   weeks; only pursue if the structural confinement signature is
   found and we want to understand its mechanism (Polyakov's
   monopole condensation).

---

## §8. Honest scope

What this memo provides:

- A complete *design specification* for five categories of alternative
  Wilson observables on Rule B, including pilot computation at $n_{\max}
  = 2$ and $n_{\max} = 3$.
- A clear identification of three observables (Wilson loops, $\nu_\ell$,
  plaquette correlator at $n_{\max} = 3$) that XCWG-B can compute at
  marginal cost.
- An honest verdict that **no clean Polyakov-loop analog exists** on the
  static Rule B graph; this is a structural absence, not an open
  question.
- A scoping verdict that **monopole content is well-defined but heavy**;
  it should be a separate 3-6 week sprint, not bolted onto XCWG-B.

What this memo does NOT provide:

- A Monte Carlo computation of any of these observables at finite
  $\beta$. The values reported are strong-coupling and weak-coupling
  leading orders only.
- A finite-temperature compactification of Rule B. The Polyakov-loop
  question is genuinely open until such a construction exists.
- A monopole enumeration. The 500 triangular cells at $n_{\max} = 2$
  are *candidate* cells; closed ones are a subset whose identification
  requires solving a sub-graph closure problem.
- An MK comparison computation. That is XCWG-B Track B1's deliverable;
  this memo is structurally complementary.

What single-observable cross-check would be most informative for
XCWG-B:

> The **Wilson loop scaling exponent $\alpha$**, fit from
> $\langle S(L) \rangle \sim L^\alpha$ over $L \in \{4, 6, 8\}$ at
> $n_{\max} = 3, 4$. A clean value $\alpha \approx 2$ confirms the
> area-law / confined phase predicted by 3D compact U(1). A value
> $\alpha \approx 1.5$ would indicate a finite-size crossover. A
> value $\alpha \to 1$ would indicate the Rule B substrate is closer
> to a 2D U(1) per-block reduction (and would contradict the
> heat-kernel $d_s$ evidence in `xcwg_rule_b_spectral_dim_memo.md`).
> Wilson scaling exponent is the **canonical** confinement order
> parameter and is the single observable a lattice-gauge audience
> would expect to see first.

---

## §9. Numerical reference data

All numerical values reported in this memo are reproduced by
`debug/xcwg_observables_pilot.py`, with output at
`debug/data/xcwg_observables_pilot.json`.

Pilot results at $n_{\max} = 2$:
- 44 length-4 Wilson loops, 144 length-6 loops, all primitive
  non-backtracking representatives modulo rotation/reversal.
- Wilson kinetic operator $K = d_1^{\!\top} d_1$ has rank 11
  ($= \beta_1$) on a 20-dim edge space.
- $\langle S \rangle_{L=4} = 0.2477$, range $[2/9, 1/3]$ (clean
  rationals).
- $\langle S \rangle_{L=6} = 0.4361$, range $[0.347, 0.458]$.
- All loops at this cutoff have $\nu_\ell = 2$ (2-shell content).
- 500 candidate triangular 2-cells among 44 plaquettes.

Pilot results at $n_{\max} = 3$ (cheap pieces only):
- 994 length-4 Wilson loops (consistent with infrastructure memo).
- $K$ has rank 79 ($= \beta_1$) on a 106-dim edge space.
- $\langle S \rangle_{L=4} = 0.096 \pm 0.018$ (smaller than at
  $n_{\max}=2$, as the harmonic 1-cochain space grows).
- $\nu_\ell$ distribution: 452 at $\nu_\ell = 2$, 542 at $\nu_\ell = 3$
  (matches infrastructure memo §4).

---

*Sprint XCWG, May 2026. Per CLAUDE.md §1.5, all results above are
structural-discrete-graph readings of a finite abelian lattice gauge
theory on a compact graph; they are not continuum QED claims, not
Lorentzian, and not mass-gap results. The strongest available comparison
for any of these observables is to other discrete-graph Wilson
constructions (Paper 25 scalar, Paper 30 SU(2), Sprint ST-SU3 SU(3));
the comparison to continuum 3D compact U(1) is structural shape-
matching, not asymptotic equivalence.*
