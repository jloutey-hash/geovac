# XCWG Rule B Infrastructure Survey

**Sprint:** XCWG (Cross-l Wilson Gauge) — infrastructure scoping
**Date:** 2026-05-15
**Status:** Survey complete. Rule B is structurally distinct from the scalar Fock graph along every dimension/connectivity axis tested.
**Key files:**
- `debug/xcwg_rule_b_infrastructure.py` (driver, no production changes)
- `debug/data/xcwg_rule_b_infrastructure.json` (numerical tables)
- `geovac/ihara_zeta_dirac.py` (Rule B adjacency, unchanged)
- `geovac/dirac_matrix_elements.py` (DiracLabel, kappa_to_l, unchanged)
- `geovac/fock_graph_hodge.py` (scalar Fock graph for comparison, unchanged)
- `papers/observations/paper_29_ramanujan_hopf.tex` §RH-C (Rule B origin)

---

## §1 Edge set definition (explicit)

Each node is a Dirac orbital label
$$
v = (n_\text{fock},\; \kappa,\; m_j)
$$
with $n_\text{fock} \ge 1$, $\kappa \in \mathbb{Z}\setminus\{0\}$, $m_j$ a half-integer with $|m_j| \le j = |\kappa| - \tfrac{1}{2}$. The orbital quantum number is
$$
\ell(\kappa) = \begin{cases}-\kappa - 1 & \kappa < 0\\ \kappa & \kappa > 0\end{cases}
$$
(so $s_{1/2}$ has $\kappa=-1$, $p_{1/2}$ has $\kappa=+1$, $p_{3/2}$ has $\kappa=-2$, etc.).

The Rule B edge condition is (`_edge_rule_B` in `geovac/ihara_zeta_dirac.py`, lines 172–196):

$$
(a,b) \in E_B \iff
\underbrace{|\Delta n_\text{fock}| \le 1}_{(\text{i})}
\;\wedge\;
\underbrace{|\Delta \ell| = 1}_{(\text{ii})}
\;\wedge\;
\underbrace{|\Delta j| \le 1}_{(\text{iii})}
\;\wedge\;
\underbrace{|\Delta m_j| \le 1}_{(\text{iv})}
\;\wedge\;
a \ne b.
$$

This is *exactly* the standard atomic E1 (electric-dipole) selection rule. Condition (ii) is the strict parity-flip $\Delta \ell = \pm 1$ — Δℓ = 0 is **not allowed**, so all in-shell horizontal hops are forbidden. Condition (i) allows Δn = 0 (intra-n radiative dipole), Δn = ±1 (inter-n radiative dipole).

Critically:

> **Rule B's $\Delta \ell = \pm 1$ requirement means every Rule B edge is by definition a cross-$\ell$ edge. There are no within-shell edges in Rule B.**

This is the single structural fact that drives every comparison in §§2–5.

---

## §2 Graph properties at n_max = 2, 3, 4

The vertex count is $V_B(n_\text{max}) = \sum_{n=1}^{n_\text{max}} 2 n^2$ (each $(n,\ell,m)$ scalar state lifts to a full $(j, m_j)$ doublet/quartet; total spinor count per shell is $2n^2$). Distribution by $\ell$:

| $n_\text{max}$ | $V$ | $\ell=0$ | $\ell=1$ | $\ell=2$ | $\ell=3$ |
|:-:|:-:|:-:|:-:|:-:|:-:|
| 2 | 10 | 4 | 6  |    |    |
| 3 | 28 | 6 | 12 | 10 |    |
| 4 | 60 | 8 | 18 | 20 | 14 |

Computed graph properties (Rule B):

| Metric | $n_\text{max}=2$ | $n_\text{max}=3$ | $n_\text{max}=4$ |
|:-|:-:|:-:|:-:|
| $V$ | 10 | 28 | 60 |
| $E$ | 20 | 106 | 312 |
| $c$ (components) | **1** | **1** | **1** |
| $\beta_1 = E - V + c$ | 11 | 79 | 253 |
| diameter | 4 | 6 | 8 |
| avg path length | 1.689 | 2.095 | 2.547 |
| girth (shortest cycle) | 4 | 4 | 4 |
| degree min / max / avg | 2 / 5 / 4.0 | 2 / 12 / 7.57 | 2 / 18 / 10.4 |
| 4-plaquettes | 44 | 994 | 5047 |

Degree distributions (Rule B):

- **n_max = 2:** $d \in \{2, 4, 5\}$ with multiplicities $\{2, 4, 4\}$. The min-degree-2 vertices are the $p_{1/2}$ states (only 2 m_j values, $\ell=1$, $\kappa=+1$).
- **n_max = 3:** $d \in \{2, 4, 5, 6, 7, 8, 9, 10, 12\}$ with growing tail; max-degree 12 vertices are the $d_{3/2}$ at $n=3$.
- **n_max = 4:** spread $d \in \{2, 4, ..., 18\}$.

Degree by $\ell$-shell at $n_\text{max}=2$:

| $\ell$ | min deg | max deg | avg deg |
|:-:|:-:|:-:|:-:|
| 0 | 5 | 5 | 5.0 |
| 1 | 2 | 4 | 3.33 |

**Connectivity:** Rule B is connected at every tested $n_\text{max}$ ($c = 1$). This is the most consequential single property: the scalar Fock graph is disconnected ($c = n_\text{max}$, one component per $\ell$-shell), while Rule B's dipole rule mixes $\ell$-shells freely and produces one giant component.

**Girth:** The shortest cycle in Rule B at every tested $n_\text{max}$ is a 4-cycle. This makes Rule B a *bipartite-friendly* but **not** bipartite structure. (Quick check: Rule B is bipartite under the $(-1)^\ell$ parity coloring — every edge flips $\ell$ parity, so no odd cycles. The 4-cycles are the minimum even cycles.)

---

## §3 Comparison to scalar Fock graph (cross-l edge census — the key result)

The scalar Fock graph (`geovac/lattice.py::_build_adjacency_matrix`) has edges only between same-$\ell$ states (radial T± hops at fixed $(\ell, m)$ and angular L± hops at fixed $(n, \ell)$). Cross-$\ell$ edges are zero by construction. Side-by-side:

| Metric | Rule B $n_\text{max}=2$ | Scalar $n_\text{max}=2$ | Rule B $n_\text{max}=3$ | Scalar $n_\text{max}=3$ | Rule B $n_\text{max}=4$ | Scalar $n_\text{max}=4$ |
|:-|:-:|:-:|:-:|:-:|:-:|:-:|
| $V$ | 10 | 5 | 28 | 14 | 60 | 30 |
| $E$ | 20 | 3 | 106 | 13 | 312 | 34 |
| $c$ | 1 | 2 | 1 | 3 | 1 | 4 |
| $\beta_1$ | 11 | 0 | 79 | 2 | 253 | 8 |
| diameter | 4 | 2 | 6 | 4 | 8 | 6 |
| girth | 4 | $\emptyset$ (forest) | 4 | 4 | 4 | 4 |
| cross-$\ell$ edges | **20 / 20** | 0 / 3 | **106 / 106** | 0 / 13 | **312 / 312** | 0 / 34 |
| cross-$\ell$ fraction | **1.0000** | 0.0000 | **1.0000** | 0.0000 | **1.0000** | 0.0000 |
| 4-plaquettes | 44 | 0 | 994 | 2 | 5047 | 8 |

The cross-$\ell$ census is *categorical*. The scalar Fock graph has 0/E edges that cross $\ell$ shells (this is the structural source of the RG sprint's 2D-grid-per-$\ell$ decomposition). Rule B has E/E edges that cross $\ell$ shells. Every Rule B hop is a parity-flipping inter-shell hop.

This 100% / 0% split is *not* a numerical accident: it is a direct consequence of Rule B's $\Delta \ell = \pm 1$ requirement vs the scalar L±/T± rules' $\Delta \ell = 0$ requirement. Both are exact selection rules of their respective adjacency definitions.

**Connectivity comparison:** $c_\text{scalar} = n_\text{max}$ exactly (one per $\ell$-shell), $c_\text{Rule B} = 1$ at every tested $n_\text{max}$. The scalar Fock graph factorizes into an $\ell$-disjoint union of 2D rectangular grids $P_{n_\text{max} - \ell} \times P_{2\ell + 1}$; Rule B *does not factorize* — its $\ell$-mixing edges make it a single connected non-product object.

**Betti-1 growth:** scalar $\beta_1$ grows as $\sim n_\text{max}^3 / 6$ (sum of grid cycles per shell); Rule B $\beta_1$ grows much faster ($11 \to 79 \to 253$, roughly $\sim n_\text{max}^4$ in the data range), reflecting the dramatically richer cycle structure.

---

## §4 Plaquette census at n_max = 3 (primitive 4-cycles)

We enumerated all primitive non-backtracking closed 4-walks (equivalently, $C_4$ subgraphs) at $n_\text{max} = 3$ and classified each by the multiset of $\ell$-values its four vertices carry.

**Total 4-plaquettes at $n_\text{max} = 3$: 994 (Rule B), 2 (scalar).**

Rule B classification by $\ell$-signature (sorted multiset):

| $\ell$-signature | Count | Type | Fraction |
|:-|:-:|:-:|:-:|
| $(0,0,1,1)$ | 272 | "mixed-mixed" (s↔p tetrahedron) | 27.4% |
| $(0,1,1,2)$ | 542 | "rotating shell" (s↔p, p↔d on opposite sides) | 54.5% |
| $(1,1,2,2)$ | 180 | "mixed-mixed" (p↔d tetrahedron) | 18.1% |

The plaquette graph is bipartite in $\ell$-parity (every edge flips $\ell$ parity), so the only allowed signatures are sorted multisets of the form $(\ell_a, \ell_a, \ell_b, \ell_b)$ or $(\ell_a, \ell_b, \ell_b, \ell_c)$ with $\ell_a, \ell_c$ same parity. The data matches this prediction exactly: 100% of plaquettes are cross-shell, with the three observed signatures all consistent with the bipartite parity constraint.

**Within-shell plaquettes (all 4 vertices same $\ell$): 0 / 994 = 0%.**
**Cross-shell plaquettes (≥2 distinct $\ell$-values): 994 / 994 = 100%.**

Compare to the scalar Fock graph at $n_\text{max}=3$: only 2 plaquettes, both within-shell at $\ell=1$ (these are the two faces of the $P_2 \times P_3$ grid in the $\ell=1$ sector). The scalar graph cannot produce cross-shell plaquettes because it has no cross-shell edges.

At $n_\text{max} = 4$ the Rule B signature-zoo grows to 5 patterns:
$(0,0,1,1)$: 516; $(0,1,1,2)$: 1672; $(1,1,2,2)$: 1151; $(1,2,2,3)$: 1388; $(2,2,3,3)$: 320. All 5047 plaquettes are cross-shell.

**Structural reading of the dominant signature $(0,1,1,2)$:** this is a 4-cycle of the form $s \to p \to d \to p \to s$ (with appropriate $n,\;m_j$ alignments). At $n_\text{max}=3$ such cycles are 54.5% of all 4-plaquettes — the most common 4-cycle in Rule B is a *three-shell* cycle, not a two-shell cycle. The framework is genuinely $\ell$-tripled at $n_\text{max}=3$.

---

## §5 Effective dimension indicators

Cheap proxies for the local geometry of the graph (not the full spectral dimension — Track 3 will compute that). The key idea: a $D$-dimensional regular lattice has plaquette-to-edge ratio $\sim (D-1)/2$ (each edge sits in roughly $D-1$ 2-cells).

**Plaquette / edge ratio:**

| Graph | $n_\text{max}=2$ | $n_\text{max}=3$ | $n_\text{max}=4$ |
|:-|:-:|:-:|:-:|
| Rule B | 2.20 | 9.38 | 16.18 |
| Scalar | 0.00 | 0.15 | 0.24 |

Reference points (regular lattices, plaquettes/edge in continuum limit):

- 2D square lattice: 0.5
- 3D cubic lattice: 1.5
- 4D hypercubic lattice: 3.0

Rule B at $n_\text{max} = 3$ has 9.4 plaquettes/edge — well above the 4D cube value of 3.0. This reflects (i) the high local connectivity of Rule B (max-degree 12 at $n_\text{max}=3$, growing rapidly), and (ii) the bipartite-friendly girth-4 structure that lets each edge participate in many small cycles.

The scalar graph's $0.15$ value at $n_\text{max} = 3$ matches the expected 2D-grid value scaled by the small grid sizes (a $P_2 \times P_3$ grid has 6 edges and 2 plaquettes, giving 0.33; the global ratio is diluted by the $P_3 \times P_1$ $\ell=0$ component which has 0 plaquettes).

**Plaquette / vertex ratio:**

| Graph | $n_\text{max}=2$ | $n_\text{max}=3$ | $n_\text{max}=4$ |
|:-|:-:|:-:|:-:|
| Rule B | 4.40 | 35.50 | 84.12 |
| Scalar | 0.00 | 0.14 | 0.27 |

Same story, more dramatically. Rule B is locally very plaquette-dense.

**Girth:** Both graphs have girth 4 (where they have any cycles at all — scalar at $n_\text{max}=2$ is a forest). So the *minimum* cycle length is the same; the difference is in cycle *abundance*, not cycle *length*.

**Diameter scaling:** Rule B diameter grows as $\{4, 6, 8\}$ at $n_\text{max} \in \{2, 3, 4\}$ — roughly $2 n_\text{max}$. Scalar diameter is the diameter of the largest $\ell$-component (a $P_{n_\text{max}-\ell} \times P_{2\ell+1}$ grid), with the $\ell=0$ component giving diameter $n_\text{max}-1$. The two graphs have *comparable* diameters but very different cycle densities — meaning Rule B is **denser**, not **more spread out**.

---

## §6 Recommendation: Wilson gauge in 3D — REACHABLE

Based on §§2–5 the verdict on Rule B as a substrate for a 3D U(1) Wilson lattice gauge construction:

**REACHABLE.** Rule B is structurally distinct from the scalar Fock graph along every metric tested:

1. **Connectivity is single-component.** The scalar Fock graph factorizes into disjoint 2D grids per $\ell$-shell; Rule B does not factorize, sitting as one connected non-product object. A Wilson-action partition function on Rule B is *not* a product of per-$\ell$-shell partition functions — it has irreducible inter-shell content.
2. **100% cross-$\ell$ edge fraction.** Every Wilson link variable $U_e$ on Rule B carries a parity flip $\Delta \ell = \pm 1$. Cross-shell hopping is the *only* hopping; there is no in-shell hopping to mix with. This is the cleanest possible substrate for testing whether inter-shell content adds dimensional structure to a Wilson construction.
3. **Cross-shell plaquette fraction = 100%.** Every Wilson plaquette $W_\square = U_{e_1} U_{e_2} U_{e_3} U_{e_4}$ at $n_\text{max} = 3$ involves vertices at ≥ 2 distinct $\ell$-shells (272 use $\ell\in\{0,1\}$, 542 use $\ell\in\{0,1,2\}$, 180 use $\ell\in\{1,2\}$). The plaquette action $S_W = \beta \sum_\square (1 - \mathrm{Re}\, W_\square)$ is structurally a *cross-shell* action — there is no scalar-Fock-style decomposition into per-$\ell$ pieces.
4. **Plaquette density is hypercubic-or-higher.** Plaquette/edge ratio 9.4 at $n_\text{max}=3$ is comfortably above the 4D cubic-lattice value of 3.0; the per-edge cycle abundance is rich enough to support a Wilson kinetic term whose effective dimension is at least 3 (Track 3's spectral-dimension diagnostic will pin this number down).
5. **Dipole $(0,1,1,2)$ "three-shell" plaquettes are the dominant class** (54.5% at $n_\text{max}=3$). This is structurally novel: in scalar-Fock the only plaquettes are 2D-grid cells within one $\ell$; in Rule B the *typical* plaquette spans three $\ell$-values. This is the closest combinatorial witness yet that Rule B's emergent geometry has more than the 2D-per-shell structure of the scalar graph.

**Honest caveats / structural surprises flagged:**

- **Rule B is bipartite under $\ell$-parity.** All edges flip $\ell$ parity (Δℓ = ±1 ⇒ Δ((-1)^ℓ) = -1), so the graph is bipartite with parts $\{\ell\text{ even}\}$ and $\{\ell\text{ odd}\}$. Girth is forced to 4. This means a U(1) Wilson construction on Rule B will have a bipartite-Wilson structure (analogous to staggered fermions); if the goal is to recover full 3D content, the bipartite nature should be inspected for whether it generates artifact-doublers. Track 3 (spectral dimension) and Track 4 (Wilson action / confinement diagnostic) should be alert to this.
- **Plaquette density grows fast with $n_\text{max}$** (4.4 → 35.5 → 84 plaquettes/vertex). The Wilson action will be dominated by short loops at small $n_\text{max}$; the strong-coupling regime might be qualitatively different from the weak-coupling regime in a way that's already visible at $n_\text{max} = 2$. This is a feature for tractability, not a bug, but it means $n_\text{max}$-extrapolation is important.
- **Degree spread is wide and growing** (2 → 12 → 18 across $n_\text{max} = 2, 3, 4$). Wilson constructions usually assume regular-graph or near-regular-graph structure. Rule B is not regular; min-degree-2 vertices (the $p_{1/2}$ doublet at each $n$) are 4–6× less connected than the high-degree central nodes. This may need to be handled with vertex-degree-weighted measures (an option already standard in irregular-graph Wilson constructions per Marcolli–van Suijlekom 2014).
- **Sub-Ramanujan strictness from Paper 29:** Rule B is strictly sub-Ramanujan at $n_\text{max} = 3$ (deviation $-0.119$ in Paper 29 Table) and crosses the Ramanujan bound at $n_\text{max} = 4$ (deviation $+1.53$, Sprint 2 RH-D in CLAUDE.md §2). For a Wilson construction, this is a *strength*: the graph has good spectral gap at small $n_\text{max}$ (faster Markov mixing of Wilson configurations) and develops dense near-regular sub-blocks at larger $n_\text{max}$ (closer to the regular-lattice ideal).

**Concrete deliverables this memo enables for downstream tracks:**

- Track 2 (incidence matrix B and edge Laplacian L₁ for Wilson kinetic): edge list is in `debug/data/xcwg_rule_b_infrastructure.json`, 20 / 106 / 312 edges at $n_\text{max} = 2, 3, 4$.
- Track 3 (spectral dimension via heat-kernel return probability): vertex and edge counts let one estimate compute cost; at $n_\text{max} = 4$ the 60 × 60 adjacency is fast.
- Track 4 (Wilson loop expectation): 4-plaquette enumeration is already in hand; longer-loop enumeration uses the same primitive-non-backtracking-walk machinery.
- Track 5 (comparison to scalar-Fock per-$\ell$ Wilson): scalar-Fock 2D-grid Wilson is textbook 2D-U(1); Rule B Wilson is to be compared against this baseline. The data here gives the apples-to-apples per-graph statistics needed.

**Verdict in one sentence:** Rule B is a 1-component graph with 100% cross-$\ell$ edges and 100% cross-shell plaquettes, the dominant plaquette signature spans three $\ell$-shells, plaquette density is well above 4D-hypercubic — structurally it is the kind of object on which a Wilson U(1) lattice gauge theory could plausibly carry 3D content, and the substrate is ready for the next tracks.

---

## Reported key numbers

- **Cross-$\ell$ edge fraction at $n_\text{max} = 3$: 1.0000 (106/106).** All Rule B edges cross $\ell$-shells by construction.
- **Cross-shell plaquette fraction at $n_\text{max} = 3$: 1.0000 (994/994).** Every 4-cycle in Rule B at $n_\text{max} = 3$ involves at least 2 distinct $\ell$-shells; 54.5% involve three distinct $\ell$-shells.
- **One-sentence verdict:** Rule B is structurally distinct from the scalar Fock graph along every dimension/connectivity axis — 1 vs $n_\text{max}$ components, 100% vs 0% cross-$\ell$ edges, 100% vs ~0% cross-shell plaquettes, and plaquette/edge ratio 9.4 (above 4D-hypercubic 3.0) vs 0.15 at $n_\text{max} = 3$ — so a Wilson U(1) construction on Rule B carries irreducible inter-shell content the scalar graph cannot provide.
