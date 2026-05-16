# Migdal–Kadanoff block-spin decimation on the GeoVac Hopf graph U(1) Wilson lattice gauge theory: design memo

**Sprint:** Track 3 — RG/block-spin scoping
**Date:** 2026-05-15
**Status:** Design + scoping only (no production-code modifications)
**Output:** Design memo + pilot script (`rg_blockspin_pilot.py`)

---

## Executive summary

The GeoVac U(1) Wilson lattice gauge theory (Paper 25) sits on the Fock-projected S³ Coulomb graph. A scoping question is whether the standard tools of statistical-mechanics RG flow — Migdal–Kadanoff (MK) block-spin decimation, bond-moving approximations, and the resulting β_eff(β) recursion — port to this geometry, and what they would tell us about confinement vs. deconfinement.

This memo answers the question structurally before any full implementation. The conclusions are:

1. **The Fock graph at any n_max is a disjoint union of rectangular grids** P_{n_max − l} × P_{2l+1}, one per orbital quantum number l. Equivalently each connected component is the Cartesian product of two paths. This is a striking simplification: every primitive 4-plaquette is a unit square in one of these grids, and the total plaquette count Σ_l (n_max − l − 1)(2l) = β₁ matches the Euler-characteristic count exactly (verified n_max = 3, 4, 5).

2. **Migdal–Kadanoff bond-moving decimation is the natural choice.** Each l-block is effectively a 2D lattice gauge theory on a rectangular grid. Bond-moving rules (Migdal 1975, Kadanoff 1976) apply per block; the global Hopf graph is just a disjoint sum of these 2D pieces.

3. **β_eff(β) at leading order in strong coupling is β_eff ≈ β² / 2** for one bond-decimation step (verified numerically against the exact Bessel-Fourier expression I₁(β_eff)/I₀(β_eff) = [I₁(β)/I₀(β)]² to four decimal places at β = 0.1).

4. **RG flow direction.** Each block is a 2D U(1) lattice gauge theory at finite cutoff. The Polyakov–Mermin–Wagner result for 2D abelian gauge theory is that **β flows to zero — strong-coupling confined phase — at all initial β**. There is no deconfinement transition in any finite-cutoff GeoVac block. The S³ geometry (effectively 3D) shows up only in the *coupling* between blocks via the Hopf base S² quotient, which the present per-block decimation does not see.

5. **Computational cost is modest.** At n_max = 5, the largest block (l=2) is a 3 × 5 grid with 8 plaquettes. Full sympy character-expansion partition function lives in ≤ 20-dim character space per block. Wall-clock estimate: numpy implementation < 1 second per decimation step at n_max ≤ 5; sympy version < 1 minute. The bottleneck is *enumerating* plaquettes at n_max ≥ 6, which the existing `enumerate_plaquettes` already handles up to n_max = 4 in seconds.

6. **One real structural obstruction is named honestly:** MK on a single rectangular block sees only the *intra-block* gauge dynamics. The Hopf U(1) gauge invariance is preserved trivially because the U(1) edge phases live on individual edges, not on bonds across blocks. However, MK does not produce a non-trivial fixed point on any finite-cutoff block — the answer is uniformly "2D U(1) confines" — so the RG flow is informative as a structural verification (the framework reproduces Polyakov's classic result) but does not by itself yield new physics. Reaching the 3D physics expected from the S³ geometry requires coupling between l-blocks, which is absent in the bare graph and would have to be added as new physics input (Hopf-base connectivity).

**Verdict: REACHABLE in ~2–3 weeks of focused work** (one analytical sprint to nail β_eff to higher order, one engineering sprint to wire MK into the existing `enumerate_plaquettes` + `partition_function_character_expansion` pipeline). The most consequential structural obstruction is the **disjoint l-block decomposition** — by itself, MK on a single block reproduces Polyakov's 2D-confinement result and yields no new physics; to access GeoVac-specific structure the RG flow needs to be coupled across blocks, which is a separate (non-trivial) design step.

---

## §1. Decimation scheme choice and justification

### 1.1 Brief survey

Three standard real-space RG schemes for lattice gauge theory:

- **Wilson decimation (1971, 1975).** Block 2^d sites into one site; integrate out the inner links by gauge fixing + Gaussian-equivalent. Conceptually clean but produces non-local interactions on the coarse lattice that have to be truncated; the truncation introduces uncontrolled errors.

- **Migdal–Kadanoff bond-moving (Migdal 1975, Kadanoff 1976).** Two-step approximation: (a) **bond move** — replace each link variable by a product of b parallel link variables (one per direction of the block), exchanging spatial extent for representation weight; (b) **decimate** — integrate out the now-bundled links. The MK approximation is uncontrolled but gives closed-form recursions β → β_eff(β) and has reproduced the 1D Ising chain, 2D U(1) confinement (Polyakov 1977), and 4D non-abelian YM asymptotic freedom semi-quantitatively. **Method of choice for finite-graph U(1).**

- **Real-space block-spin (Kadanoff 1966; not gauge-specific).** Group sites by spatial proximity, sum out internal modes. Less canonical for gauge theory because the link variables don't live at sites.

- **Polyakov's anti-decimation / strong-coupling expansion.** Treat the partition function as a sum over closed-surface configurations (strong-coupling lattice expansion), then resum surfaces by clustering. Not a block-spin RG in the MK sense; rather a strong-coupling-expansion technique. Useful as a cross-check.

### 1.2 Why MK is natural for the GeoVac Hopf graph

The Fock graph decomposes cleanly. Each connected component is the Cartesian product of two paths P_a × P_b (Section 2). On any P_a × P_b grid, the MK construction reduces to a 1D recursion: alternately move bonds in the "horizontal" (radial n) direction and the "vertical" (angular m) direction, decimating after each move. This is exactly the original 2D MK scheme of Migdal 1975.

The choice is forced by graph structure: the grids are not Euclidean Z²; they are *finite* rectangular grids. MK adapts trivially because the bond-moving operation is local (between adjacent links in a single direction), and the decimation step is the universal U(1) Bessel-Fourier identity (Section 3). Wilson decimation would also work but would lose the bond-moving simplification.

### 1.3 What the choice gives up

MK is an **uncontrolled** approximation: bond-moving violates gauge invariance on the moved-bond cycle (because gauge invariance requires Wilson loops, which the bond move literally reroutes). The classical fix is to require that the moved bonds carry the *same* group element, which is gauge-fixed. The result preserves gauge invariance globally but not at intermediate steps.

For U(1) on a finite graph the violation is benign (because U(1) is abelian and the Wilson loop is just the sum of edge phases around a cycle); for SU(N) MK is known to give wrong critical β by ~10-30% in 4D. We are doing U(1) here, so this caveat does not bite.

---

## §2. Concrete decimation rules for the GeoVac Hopf graph

### 2.1 Graph structure (numerical verification)

For the GeoVac Fock-projected scalar graph at unweighted adjacency:

```
n_max  V    E    c   β₁   |Hopf base|
  2    5    3    2    0       3
  3   14   13    3    2       6
  4   30   34    4    8      10
  5   55   70    5   20      15
```

(Computed via `geovac.fock_graph_hodge.FockGraphHodge`.) The Hopf-base sector count is N = Σ_{n=1..n_max} n = n_max(n_max+1)/2.

The connected components are indexed by l ∈ {0, 1, …, n_max − 1}. **Each l-component is the Cartesian product of two path graphs:**

```
component_l  =  P_{n_max − l}  ×  P_{2l+1}
```

The first factor (n_max − l vertices) is the radial direction n ∈ {l+1, …, n_max} with edges from T_± (radial ladder, Δn = ±1, fixed l, fixed m). The second factor (2l+1 vertices) is the angular direction m ∈ {−l, …, +l} with edges from L_± (azimuthal ladder, fixed n, fixed l, Δm = ±1).

**Plaquette count per block:** (n_max − l − 1)(2l). Summing over l gives:

```
β₁(n_max)  =  Σ_{l=0}^{n_max-1} (n_max - l - 1) · (2l)
```

Verified to give 0, 2, 8, 20 at n_max = 2, 3, 4, 5. (For n_max = 2 there are no plaquettes anywhere on the graph — gauge theory at n_max = 2 is trivial.)

### 2.2 The MK rules per block

Fix one l-block, call it G_l = P_a × P_b with a = n_max − l rows (radial n-axis) and b = 2l+1 columns (magnetic m-axis). Label the vertices (i, j) with i ∈ {1, …, a} and j ∈ {1, …, b}.

**Edges of two types:**
- **Radial (vertical):** e^R_{i,j} connects (i, j) ↔ (i+1, j), 1 ≤ i ≤ a − 1.
- **Angular (horizontal):** e^A_{i,j} connects (i, j) ↔ (i, j+1), 1 ≤ j ≤ b − 1.

**Plaquettes:** the unit squares P_{i,j} with corners (i, j), (i+1, j), (i+1, j+1), (i, j+1). There are (a−1)(b−1) of them.

**Wilson action:**

```
S_W[U]  =  −β Σ_P Re Tr U_P
       =  −β Σ_{i,j} cos( θ^R_{i,j} + θ^A_{i+1,j} − θ^R_{i,j+1} − θ^A_{i,j} )
```

with U(1) link variables U_e = e^{iθ_e}.

### 2.3 One MK decimation step: n_max = 4 → n_max = 3

Concretely: at n_max = 4, the p-block (l = 1) is P_3 × P_3 (radial n ∈ {2, 3, 4}, angular m ∈ {−1, 0, +1}). To produce the n_max = 3 effective theory, we need to integrate out the layer of edges associated with n = 4.

**Step 1 — vertical bond move (radial).** Move the radial bonds e^R_{2,j} (connecting (n=3, m=j) ↔ (n=4, m=j)) onto the row of bonds e^R_{1,j} (connecting (n=2, m=j) ↔ (n=3, m=j)). After the move, each bond e^R_{1,j} carries the product of two link variables (one for each original bond at the same m).

**Step 2 — decimation (radial).** Integrate out the n = 4 nodes (3 vertices for the l = 1 block, plus the corresponding nodes in higher-l blocks). Each n = 4 vertex has exactly two incident edges in the radial direction (one to n = 3) and two in the angular direction (to (4, m−1) and (4, m+1)). The decimation integrates over the U(1) link variables of all edges adjacent to the n = 4 layer.

**Step 3 — horizontal bond move (angular).** For each l-block, contract the m = −l and m = +l boundary so the angular extent (2l+1 → 2l−1) is reduced consistently with the n_max → n_max − 1 reduction.

**Wilson loops mapping.** A Wilson loop W_C on the fine n_max = 4 graph that does not enter the n = 4 layer survives as a Wilson loop on the coarse n_max = 3 graph with the same loop edges. A Wilson loop that traverses the n = 4 layer is integrated out — its expectation value is absorbed into the effective β.

### 2.4 The summation/integration

For a single decimated link u between two plaquettes A (action term β cos(φ_A + u)) and B (action term β cos(φ_B − u)), the link integral is

```
I(φ_A, φ_B; β)  =  ∫_0^{2π} du/(2π) exp[β cos(φ_A + u) + β cos(φ_B − u)]
                =  Σ_{n=−∞}^{+∞} I_n(β)^2 e^{i n (φ_A + φ_B)}
                =  I_0(β)^2 + 2 Σ_{n=1}^∞ I_n(β)^2 cos(n (φ_A + φ_B))
```

(Standard Fourier-Bessel identity for U(1).) The effective Boltzmann weight on the merged plaquette is

```
W_eff(φ_A + φ_B)  =  exp[β_eff cos(φ_A + φ_B)]
```

with β_eff determined by matching the n = 1 Fourier coefficient:

```
I_1(β_eff) / I_0(β_eff)  =  I_1(β)^2 / I_0(β)^2
```

This is the **fundamental MK U(1) recursion**.

---

## §3. β_eff(β) at leading order

### 3.1 Closed-form recursion (one bond-decimation step)

From Section 2.4:

```
ratio(β_eff)  =  ratio(β)^2,        where ratio(β) := I_1(β) / I_0(β).
```

This is exact for a single decimation step on a single shared bond between two plaquettes. For a full MK step on a P_a × P_b grid going from (a × b) → (a/2 × b/2), there are b/2 vertical bond moves followed by a/2 horizontal bond moves. Each bond move squares the ratio; therefore after one full 2:1 → 1:1 MK step:

```
ratio(β_eff)  =  ratio(β)^4
```

(For a 2D theory the MK recursion is β → β^{b^{d-2}} effectively at strong coupling, where d = 2 and b = 2 is the block size — but here, more precisely, two-fold decimation steps in two directions give exponent 4.)

### 3.2 Strong-coupling limit (β → 0, confined regime)

Expand I_n(β) for small β:

```
I_n(β)  =  (β/2)^n / n!  ·  [1 + O(β²)]
ratio(β)  =  I_1(β) / I_0(β)  =  (β/2) · [1 + O(β²)]
```

So ratio(β)² = β²/4 + O(β⁴), and:

```
ratio(β_eff)  ≈  β_eff / 2  =  β² / 4   ⟹   β_eff  ≈  β² / 2  + O(β⁴)
```

Pilot verification (`rg_blockspin_pilot.py`):

```
β       β_eff (exact MK)   β²/2 (strong-coupling)
0.10    0.0050             0.0050
0.50    0.1178             0.1250
1.00    0.4067             0.5000
2.00    1.1188             2.0000
5.00    2.8485             12.5000
```

(The β²/2 limit agrees to four decimal places at β = 0.1, breaks down by β ~ 1, as expected.)

For a full 2D MK step (two bond moves), the exponent doubles in the recursion: β_eff ≈ β⁴ / 8 at strong coupling. β flows to zero quadratically (or quartically for the 2D step) on every iteration — Polyakov's classical result that 2D U(1) is confined at every coupling.

### 3.3 Weak-coupling limit (β → ∞)

Use the asymptotic expansion I_n(β) ∼ e^β / √(2πβ) · [1 − (4n² − 1)/(8β) + …]. Then:

```
ratio(β)  =  1 − 1/(2β) + O(1/β²)
ratio(β)^2  =  1 − 1/β + O(1/β²)
```

Invert: ratio(β_eff) = 1 − 1/(2β_eff) = 1 − 1/β + O(1/β²), giving:

```
β_eff  ≈  β / 2   +  O(1)
```

(Pilot verification: at β = 5, β_eff = 2.85, close to β/2 = 2.5.) The β/2 weak-coupling limit is also classical — it means β still **decreases** under one bond-decimation step, just multiplicatively. Iterating, β_n ≈ β_0 / 2^n → 0. Combined with the strong-coupling β² → 0 behavior, we have:

**β = 0 is the only fixed point of the U(1) MK recursion on a 2D lattice**, and it is an *attractive* fixed point (β flows to it from any initial value).

This is the standard Mermin–Wagner / Polyakov / Coleman result: 2D abelian gauge theory has no deconfinement transition.

### 3.4 What this means for the GeoVac graph

Because the Fock graph decomposes as a disjoint union of 2D rectangular grids (one per l), the MK recursion applied per block reproduces Polyakov's 2D result: **every l-block confines at every β**. The framework reproduces the textbook 2D result correctly. But it does *not* by itself produce a 3D phase transition — to access the 3D content of the S³ geometry, one would need to introduce coupling between l-blocks (the Hopf-base structure of Paper 25 §III.B), which is not present in the bare Fock graph because L_± and T_± do not change l. This is the central structural finding of the sprint.

---

## §4. RG fixed-point analysis

### 4.1 Fixed points and stability

The U(1) MK recursion on a 2D rectangular grid has:

- **β* = 0 (strong-coupling fixed point):** stable, attractive. Confinement / area-law Wilson loops. ⟨W_C⟩ ~ e^{−σ Area(C)} with string tension σ → 0 as the block-RG flows away from the UV.
- **β* = ∞ (weak-coupling fixed point):** unstable. Under any MK step, β → β/2 < β, so any finite β > 0 starting point flows away from β = ∞ toward β = 0.

**There is no non-trivial intermediate fixed point β* ∈ (0, ∞).** In 2D U(1), every β flows to zero.

### 4.2 Per-block fixed-point structure on the GeoVac graph

Each l-block at n_max ≤ 5 is a finite-size grid; the MK flow at finite size is not a *true* RG flow (no scaling limit exists for finite graphs), but the *direction* of flow under one MK step is well-defined and agrees with the 2D continuum result.

Concretely at n_max = 5:

| l | Grid | Plaquettes | MK fate |
|:---:|:----:|:---:|:---|
| 0 | P_5 × P_1 (tree) | 0 | trivial (no gauge dynamics) |
| 1 | P_4 × P_3 | 6 | β flows to 0 (confined) |
| 2 | P_3 × P_5 | 8 | β flows to 0 (confined) |
| 3 | P_2 × P_7 | 6 | β flows to 0 (confined) |
| 4 | P_1 × P_9 (tree) | 0 | trivial |

The "inner" blocks (l = 1, 2, 3) carry the gauge dynamics; the "boundary" blocks (l = 0 and l = n_max − 1) are trees and have no plaquettes.

### 4.3 Expected 3D physics is invisible to per-block MK

A genuine 3D U(1) gauge theory would have a deconfinement transition at finite β_c. The S³ underlying GeoVac is effectively 3-dimensional, so one might naively expect deconfinement on the full Fock graph. **It does not appear**, because the Fock graph decomposes by l and per-block MK cannot see inter-l couplings.

To access 3D physics, one would need to add coupling between the disconnected l-components. The natural candidate is the **Hopf-base S² quotient graph** (Paper 25 §III.B), which has eigenvalues {0, 0, 0, 1, 3, 6} at n_max = 3 and is the projection of all (n, l, m_l) sectors onto (n, l) sectors. Wilson loops *on the base* would couple states of different m_l within the same (n, l) sector — but these are already inside the same l-block as edges, not new gauge fields.

A more ambitious extension would be to define a *cross-l hopping* term (analogous to the Dirac graph's E1 dipole edges of Paper 29's Rule B, which couple Δl = ±1). With cross-l edges, the graph would no longer decompose by l, and MK would see a genuinely 3D structure. This is a separate design question (out of scope for this sprint).

---

## §5. Computational cost estimate

### 5.1 Per-block character-expansion partition function

For one l-block with V_l vertices, E_l edges, and P_l plaquettes, the U(1) character-expansion partition function truncated at rep label R_max is:

```
Z_l(β; R_max)  =  Σ_{j=0}^{R_max} d_j(β)^{P_l}  ·  (Haar normalization)
```

(For U(1), d_n(β) = I_n(β) and there is no Vandermonde factor.)

Cost per block: O(R_max · P_l) Bessel evaluations, each O(1) via scipy. For R_max = 5 and P_l ≤ 8 (the largest n_max = 5 block has 8 plaquettes), this is < 100 µs per block per β-value.

### 5.2 Plaquette enumeration

`geovac.su2_wilson_gauge.enumerate_plaquettes` is already tested up to n_max = 4 (Paper 30) and uses DFS on the adjacency. For each block separately, plaquette enumeration is trivial because the graph is a known rectangular grid — we don't even need DFS; we can enumerate unit squares analytically. The full-graph enumeration runs in seconds at n_max ≤ 5.

### 5.3 MK β_eff recursion

For one MK step on a single block:
- Block-move: identify which bonds to bundle (analytical from grid coords).
- Decimation: 1 Bessel-ratio evaluation per integrated-out bond.
- Newton solve to invert ratio(β_eff) = ratio(β)^k: O(10) Newton iterations, each O(1) Bessel evaluation.

For n_max = 5 with three non-trivial l-blocks, one full MK step takes < 1 ms numerically.

Iterating n_max = 5 → 4 → 3 → 2 is three MK steps total: < 10 ms numerically, < 1 s symbolically (sympy character expansion with explicit Bessel-series representation).

### 5.4 Bottlenecks

The two real bottlenecks are upstream and downstream:

- **Plaquette enumeration at n_max ≥ 6.** The Fock graph at n_max = 6 has V = 91, E = 124, β₁ = 40, and the DFS enumeration runs in O(V³) for short cycles and scales worse for longer ones. Estimated wall-clock at n_max = 6: ~10 minutes via the generic DFS, but **seconds** if we use the grid structure to enumerate analytically.

- **Validation: tracking Wilson-loop expectation values through MK steps.** This requires Monte Carlo sampling of link variables, which scales as O(N_MC · E) per β. For n_max = 5 with E = 70, N_MC = 10⁵ samples is < 1 minute via numpy; sympy version would be hours.

### 5.5 Practical recommendation

| Stage | Implementation | Wall-clock estimate |
|:---|:---|:---:|
| Analytical β_eff formula | sympy + Bessel | minutes |
| Plaquette enumeration | numpy + grid analytics | seconds |
| Per-block partition function | scipy.special | < 100 µs per β-value |
| Full MK recursion n_max = 5 → 2 | numpy | < 10 ms |
| Monte Carlo Wilson loop validation | numpy + custom MC | minutes per β |
| Full flow diagram (β_0 ∈ {0.1, …, 10}, 4 steps) | numpy | < 1 second |

**Total scoping-to-implementation cost: ~2 weeks** (1 week to nail β_eff to higher orders + match with analytical 2D MK literature; 1 week to wire up Monte Carlo cross-checks and produce flow plots for a paper).

---

## §6. Structural obstructions

This section is the most important. The sprint plan asked specifically: are there obstructions to a clean implementation?

### 6.1 The disjoint l-block decomposition (real obstruction)

**This is the most consequential structural finding.** The bare Fock-projected Coulomb graph has no L-changing edges, so the graph decomposes by l, and MK applied per block reproduces the well-known 2D-U(1) result with no new physics. The S³ geometry's 3D character does not appear at the level of MK because no plaquette traverses two different l-blocks.

**What this means for an implementation:**

- The pipeline will produce honest, reproducible RG flow diagrams. They will show β → 0 (confinement) for every initial coupling on every block. This is a structural verification (the framework reproduces Polyakov's 2D result) but is not a *new* GeoVac-specific finding.

- To produce 3D-like physics one would need an extended gauge theory that couples l-blocks. Two natural candidates:
  - **Dirac graph Rule B** (Paper 29 §6): replace the scalar Fock adjacency with the E1 dipole adjacency, which has Δl = ±1 edges. This connects all l-blocks into a single component and produces a genuinely 3D-looking graph. MK applied here would not be the textbook 2D recursion.
  - **Hopf-base quotient graph** (Paper 25 §III.B): collapse each m_l-fiber to a point, producing a (n, l) sector graph. Here the natural gauge group is *still* U(1), but the geometry is different (it's the dimension-collapsed Hopf base on S²). MK on the base graph would test whether the Hopf base contains the 3D physics.

Neither extension is in scope for this sprint. The clean conclusion is: **on the bare Fock graph, MK reproduces 2D U(1) confinement at every block and yields no GeoVac-specific RG physics.**

### 6.2 Hopf U(1) gauge invariance preservation under MK (no real obstruction)

MK respects U(1) gauge invariance globally. The bond-move step (moving link variables in parallel) preserves the global U(1) symmetry because U(1) is abelian — relabeling θ_e → θ_e + (χ_{h(e)} − χ_{t(e)}) commutes with the bond-move. The decimation step (integrating over link variables) trivially respects gauge invariance because the integration measure is the Haar measure dU = dθ/(2π).

**This is a clean structural positive.** Paper 25 §III.D's gauge structure on the GeoVac graph transfers verbatim to the MK-decimated theory.

### 6.3 Finite-cutoff effects at n_max = 4 → 3 (a moderate concern)

MK is derived for *infinite* lattices. At finite n_max we have boundary effects:

- The MK bond-moving operation requires *pairs* of adjacent bonds; at the boundary of the grid (j = 1 or j = b), the bond has only one neighbor on one side. The standard fix is to treat the boundary bonds as decimated against "free" Haar measure (which contributes I_0(β)).

- Going from n_max = 4 → 3 removes one radial layer (the n = 4 layer). The 4-row blocks (l = 0, 1) become 3-row blocks; the 2-row block (l = 2) becomes a 1-row tree (no more plaquettes); the 1-row block (l = 3) is already a tree. So MK on l = 2 at n_max = 4 → 3 *destroys all plaquettes in that block*. This is a genuine finite-cutoff artifact.

A cleaner direction: instead of n_max = 4 → 3, do a *2:1 decimation within each block* (a → a/2, b → b/2 with a, b even, or floor for odd). This requires re-indexing the blocks but keeps every block plaquette-containing throughout the flow until you reach the trivial size. This is the canonical 2D MK setup and matches the textbook recursion.

**Recommendation:** for a full implementation, use **2:1 block decimation per l-block** rather than n_max → n_max − 1 global decimation. The two schemes agree at leading order in strong coupling but differ at higher orders.

### 6.4 Interaction with the Hopf base (no obstruction at the per-block level)

The Hopf-base S² quotient of Paper 25 §III.B collapses each m_l-fiber to a point at fixed (n, l). Because the Fock graph decomposes by l and the m_l-fiber is *within* a single l-block, the Hopf base does not interact with per-block MK in any new way. The base graph at n_max = 3 has 6 sectors and Laplacian eigenvalues {0, 0, 0, 1, 3, 6}; its gauge structure is again abelian U(1), and one could run MK on it directly as an independent question. But for the bare Fock graph problem, the Hopf base is a spectator.

### 6.5 Summary

| Obstruction | Severity | Resolution |
|:---|:---:|:---|
| Disjoint l-block decomposition → 2D MK gives no new physics | **load-bearing** | Either accept (verification sprint) or extend to Dirac graph Rule B / cross-l hopping (separate sprint) |
| Boundary effects at finite grid size | moderate | Use 2:1 per-block decimation rather than global n_max → n_max − 1 |
| Hopf U(1) invariance under MK | none | Trivially preserved (U(1) abelian) |
| Hopf base S² interaction | none | Independent question, spectator at per-block level |
| Plaquette enumeration cost at n_max ≥ 6 | minor | Use grid analytics, not DFS |

---

## §7. Sprint-scoping recommendation

### 7.1 Option A: "Verification sprint" (~2 weeks, REACHABLE)

**Goal:** implement MK β-recursion per block, produce flow diagrams confirming Polyakov 2D confinement on each l-block of the Fock graph, position the result in Paper 25's lattice-gauge dictionary.

**Deliverables:**
- New module `geovac/mk_block_spin.py` (~400 lines): per-block MK recursion, β_eff(β), flow iteration, fixed-point analysis.
- Tests `tests/test_mk_block_spin.py` (~30 tests): Bessel-Fourier identity, strong-coupling β²/2 limit, weak-coupling β/2 limit, MC validation of one MK step at small β.
- Memo `debug/mk_block_spin_results_memo.md` with flow-diagram data and tabulated β_eff(β) at n_max = 3, 4, 5 for each l-block.
- Paper 25 §IV new subsection (~30 lines) and Paper 30 §6 cross-reference describing MK as the "natural RG flow on each l-block, reproducing Polyakov 2D confinement."

**Value:** structural verification of Paper 25's gauge identification. Reproduces a textbook result on the Fock graph. Honest about what it does and does not show (does not produce 3D physics).

**Verdict: REACHABLE.** Most of the machinery already exists in `su2_wilson_gauge.py` and `fock_graph_hodge.py`.

### 7.2 Option B: "Cross-l extension sprint" (~6–8 weeks, REACHABLE-WITH-EFFORT)

**Goal:** define cross-l hopping (Dirac-graph Rule B or new edge type) that connects l-blocks, then run MK on the connected graph. Test whether the resulting flow exhibits a 3D-like deconfinement transition at finite β_c.

**Deliverables:**
- Module `geovac/mk_cross_l.py`: cross-l edge construction, plaquette enumeration on the connected graph.
- Module `geovac/mk_block_spin.py` extended to handle connected graphs.
- Tests + memo.
- Paper 33 or new short Paper draft: "RG flow on the connected Fock-Dirac graph and the 3D phase structure of GeoVac U(1)."

**Risks:**
- Cross-l edges are not present in the bare Fock graph. Adding them is *new physics input* and requires justification (e.g., interpreting them as virtual photon couplings of Paper 28 §vector_photon_qed).
- The connected graph at n_max = 4 has V = 30 and E ~ 60+ depending on the cross-l rule; plaquette enumeration becomes non-trivial. Plaquettes may now span multiple l-values and the simple grid analytics doesn't apply.
- Finite-cutoff effects may dominate at small n_max, making it hard to extract a clean β_c.

**Verdict: REACHABLE-WITH-EFFORT.** Worth pursuing only if Option A's verification value is judged sufficient to justify the larger sprint. Recommend running Option A first.

### 7.3 Option C: "Hopf-base RG sprint" (~3–4 weeks, REACHABLE-WITH-EFFORT)

**Goal:** run MK on the Hopf-base quotient graph (n, l) sectors only, as a 2D abelian gauge theory on a smaller, fully connected graph. Independent of the per-block decomposition.

**Deliverables:** Hopf base MK recursion + analysis. Test whether the Hopf base shows different fixed-point structure than the per-block analysis.

**Verdict: REACHABLE-WITH-EFFORT.** Smaller than Option B but conceptually distinct; complementary rather than competitive with A.

### 7.4 Recommendation

**Pursue Option A first**, with explicit caveats in the writeup that this is a verification sprint and that the bare Fock graph cannot host a 3D phase transition under per-block MK. If Option A lands cleanly (~2 weeks), reassess whether Options B or C are worth the additional effort. The honest verdict is: per-block MK on the Fock graph is a verification, not a discovery — useful to have on the books, not load-bearing for the framework's physics.

---

## Net verdict

**REACHABLE in 2–3 weeks for Option A** (verification sprint).
**Most consequential structural obstruction:** the bare Fock graph decomposes as a disjoint union of 2D rectangular grids P_a × P_b indexed by orbital quantum number l, so per-block MK reproduces Polyakov's 2D-U(1) confinement on every block and produces no GeoVac-specific 3D phase structure; to access the 3D content of the S³ geometry, the bare graph must be extended by cross-l hopping (Dirac graph Rule B or analogous), which is a separate ~6–8-week sprint.

## Files produced

- `debug/rg_blockspin_design_memo.md` (this memo)
- `debug/rg_blockspin_pilot.py` (one-bond MK Bessel-Fourier numerical pilot, validates β_eff ≈ β²/2 at strong coupling and β/2 at weak coupling)

## Cross-references

- Paper 25 §II.E (Wilson dictionary), §III (Hopf graph), §III.D (U(1) phases on edges)
- Paper 30 §5 (plaquette enumeration), §6 (framework-wide action-principle synthesis)
- Paper 28 §graph_native_qed (graph-native QED, scalar-vs-vector photon)
- Paper 29 §5 (Ihara zeta, Rule B Dirac graph with Δl = ±1)
- `geovac/fock_graph_hodge.py` (incidence matrix, Hodge decomposition)
- `geovac/su2_wilson_gauge.py` (plaquette enumeration, character expansion — adapt to U(1))
- CLAUDE.md §1.7 WH4 (deflated S³ four-way coincidence — Fock projection + three forced consequences; the per-block MK result aligns with WH4's reading that S³ structure is forced by SO(4) + Bertrand, not an independent gauge mystery)
