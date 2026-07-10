# Sprint ST-SU3: SU(3) Wilson Lattice Gauge Theory on the Bargmann-Segal S^5 Graph

**Date:** 2026-05-03
**Status:** Complete
**Verdict:** **POSITIVE WITH STRUCTURAL CAVEAT** — A genuine non-abelian SU(3) Wilson lattice gauge theory exists on the Bargmann-Segal S^5 graph at any N_max. It is gauge-invariant, Haar-normalizable, and reduces to a U(1) × U(1) Wilson theory on the maximal Cartan torus. Sprint 5 Track S5's CG obstruction is **bypassed at the gauge level** but **rediscovered at the matter-coupling level** — Wilson SU(3) is well-defined as a *pure-gauge* theory on Bargmann; coupling it to the (N,0)-shell matter encounters the same Clebsch-Gordan structure Sprint 5 identified.

**Mirrors:** Paper 30's SU(2) Wilson on S^3 = SU(2). The two theories are siblings: same construction (link variables in the gauge group, plaquettes from primitive non-backtracking walks, Wilson action), different graph and different gauge group. The structural similarities and asymmetries are documented below.

---

## §1. Goal and Scope

Paper 25 §VII.A asked: does the S^5 Bargmann-Segal graph carry a non-abelian gauge structure? Sprint 5 Track S5 returned a clean negative for the *adjacency-preserving* SU(3) action: transitions (N,0) → (N+1,0) are CG intertwiners between distinct SU(3) irreps, not group elements, so an SU(3)-action on the (N,0)-tower fails Wilson's "fixed group on every link" requirement.

This sprint tests a structurally different question: **Wilson lattice gauge theory** with link variables U_e ∈ SU(3) as *external* gauge degrees of freedom living on edges. Wilson links are 3×3 matrices in the *fundamental* representation of SU(3), independent of the shell index N. They do not require an action on shell labels and are not reps of SU(3) on the shell tower.

**Outcomes (preview):** All three structural theorems verified.
- Theorem 1 (Cartan reduction to U(1) × U(1)): VERIFIED at machine precision (max |Δ| < 1.5×10⁻¹¹ over random Cartan trials at N_max=2, 3).
- Theorem 2 (weak-coupling kinetic term = 1/12 per plaquette per su(3) component): VERIFIED both symbolically and numerically (relative error scales as O(eps²) consistent with O(A⁴) corrections).
- Theorem 3 (CG obstruction status): Bypassed at gauge level, rediscovered at matter-coupling level. Both verdicts coexist by design.

**Scope statement.** What this is: a gauge-invariant, finite, Haar-normalizable non-abelian lattice gauge theory on the finite Bargmann-Segal S^5 graph. What this is NOT: continuum SU(3) Yang-Mills mass gap (Clay), a confinement diagnostic (graph too small for area/perimeter law families), Lorentzian, or a refutation of Sprint 5 Track S5 (which concerned the irrep-tower action of SU(3), categorically different from Wilson links).

---

## §2. Construction

### 2.1 Graph

We use the Bargmann-Segal lattice from `geovac/nuclear/bargmann_graph.py` (Paper 24, Track NK). At cutoff N_max:

- Nodes: (N, l, m_l), N = 0..N_max, l of same parity as N with 0 ≤ l ≤ N, m_l = -l..+l.
  Total: V = (N_max+1)(N_max+2)(N_max+3)/6.
- Edges: dipole transitions ΔN = ±1, Δl = ±1, Δm_l ∈ {-1, 0, +1}, edge weights nonzero per Wigner 3j × radial r²-matrix-element.

For the gauge construction we only need the underlying connectivity, not the weights. The graph is **bipartite** by N parity, so all primitive cycles have **even length**.

| N_max | V | E | c | β₁ = E−V+c | N₄ | N₆ | N₈ |
|:-----:|:--:|:--:|:-:|:----------:|:---:|:---:|:----:|
| 1 | 4  | 3  | 1 | 0 | 0 | 0 | 0 |
| 2 | 10 | 15 | 1 | 6 | **15** | 21 | 146 |
| 3 | 20 | 42 | 1 | 23 | 88 | 526 | 7151 |
| 4 | 35 | 90 | 1 | 56 | 272 | 2934 | 64021 |

Compare to Paper 30's S^3 Coulomb table (V, E, β₁ at N_max = 2, 3, 4): (5, 3, 0), (14, 13, 2), (30, 34, 8). The Bargmann graph is denser and has many more cycles per shell — at N_max=3 we have β₁ = 23 vs S^3's β₁ = 2.

### 2.2 SU(3) link variables

Each undirected edge {u, v} with chosen orientation u → v carries a link variable U_{u→v} ∈ SU(3). The reverse carries U_{v→u} = U_{u→v}†. SU(3) is parametrized by 8 real parameters via the Gell-Mann basis T^a = λ^a/2, a = 1..8, with conventions tr(T^a T^b) = (1/2) δ^{ab}, [T^a, T^b] = i f^{abc} T^c (standard).

The *Cartan torus* T ⊂ SU(3) is U(1) × U(1):
```
T = { diag(e^{i a}, e^{i b}, e^{-i(a+b)}) : a, b ∈ ℝ }
```
The two phases a, b parametrize T independently; the third is fixed by det = 1.

### 2.3 Plaquettes

Plaquettes are primitive closed non-backtracking walks (same convention as Paper 30, mirrored from Paper 29's Ihara cycles). The bipartite structure means all primitive cycles have even length. The smallest plaquettes are length-4 squares connecting (N, l, m_l) → (N+1, l', m_l') → (N, l'', m_l'') → (N+1, l''', m_l''') → (N, l, m_l).

### 2.4 Wilson action

```
S_W[U] = β · Σ_P (1 − (1/3) Re tr U_P)
```
with U_P = ∏_{e in P} U_e the ordered plaquette holonomy. Vanishes iff all U_P = I (trivial vacuum); non-negative on SU(3). β = 6/g² in Wilson's SU(3) convention.

### 2.5 Gauge transformation

Node-local g: V → SU(3) acts by U_e → g_{src(e)} U_e g_{tgt(e)}†. Plaquette holonomies transform by conjugation U_P → g_{v_0} U_P g_{v_0}†, so tr U_P is gauge-invariant. **Verified to machine precision** at N_max = 3 (`debug/st_su3_theorem3_cg_obstruction.py`):

| N_max | |S_before − S_after| |
|:-----:|:--------------------:|
| 2 | 3.55 × 10⁻¹⁵ |
| 3 | 0 (exact) |

### 2.6 Partition function

Z(β) = ∫ ∏_e dU_e exp(−S_W[U]) with normalized bi-invariant Haar measure on SU(3)^|E|. Finite-dimensional; well-defined for all β ∈ ℝ.

---

## §3. Theorem 1 — Cartan torus reduction to U(1) × U(1)

**Statement.** Restricting the SU(3) Wilson action (§2.4) to link variables in the maximal Cartan torus T = U(1) × U(1) yields a U(1) × U(1) Wilson theory:
```
S_W[U_e = diag(e^{i a_e}, e^{i b_e}, e^{-i(a_e + b_e)})]
= β · Σ_P [1 − (1/3)(cos a_P + cos b_P + cos(a_P + b_P))],
```
where a_P = Σ_{e ∈ P} a_e, b_P = Σ_{e ∈ P} b_e (with reverse-edge phases negated by U_{e^{-1}} = U_e†, which negates both a_e and b_e).

**Proof.** A product of diagonal SU(3) matrices is diagonal: U_e diagonal, U_P diagonal, with phases summing the per-edge phases. The character (1/3) Re tr U_P = (1/3)(cos a_P + cos b_P + cos(a_P + b_P)). Substituting into S_W gives the U(1) × U(1) action.

**Numerical verification.** `debug/st_su3_theorem1_cartan_reduction.py`. At N_max = 2, 3 with random Cartan phases and 5 random β values:

| N_max | n plaquettes (L≤8) | n forward edges | max |S_SU3 − S_U1xU1| |
|:-----:|:------------------:|:---------------:|:---------------------:|
| 2 | 182 | 15 | 5.68 × 10⁻¹⁴ |
| 3 | 7,765 | 42 | 1.46 × 10⁻¹¹ |

Both below 10⁻¹⁰; reduction is exact at machine precision (the small numerical drift at N_max=3 is from accumulated floating-point error over 7,765 plaquette evaluations).

**Comparison to Paper 30 SU(2) → U(1).** Paper 30 §4.1 verified that SU(2) restricted to its U(1) maximal torus reduces to Paper 25's U(1) Wilson action exactly. The SU(3) → U(1) × U(1) reduction here is the rank-2 generalization. The *number of independent U(1) factors equals the rank of the gauge group*: rank 1 for SU(2), rank 2 for SU(3).

**Structural reading.** This makes Paper 30's SU(2) Wilson the *Cartan-subalgebra* sector of an SU(2) theory, and the present Wilson SU(3) has TWO independent U(1) Cartan sectors. Either of them, taken alone, is a Paper 25-style U(1) Wilson theory on the same Bargmann graph. The two U(1) sectors do NOT decouple at full Cartan level: the plaquette character mixes a_P and b_P via the third phase −(a_P + b_P), so even on the Cartan torus the action is a non-trivial *coupled* U(1) × U(1) theory.

---

## §4. Theorem 2 — Edge Laplacian L_1 as weak-coupling kinetic term

**Statement.** The second-order expansion of the SU(3) Wilson action around the trivial vacuum gives
```
S_W[U] = (β/12) · Σ_{a=1..8} ⟨A^a, L_1^{plaq} A^a⟩ + O(A⁴) + (non-abelian),
```
where A^a (a = 1..8) are the eight su(3) components on each link and L_1^{plaq} is the plaquette-incidence bilinear form on 1-cochains. On a triangle-free simplicial 1-complex (such as our bipartite graph), L_1^{plaq} coincides with the edge Laplacian L_1 = B^T B (Paper 25) up to a positive constant determined by the plaquette-length distribution. The eight su(3) components are *kinetically decoupled* at quadratic order; non-abelian self-interactions enter at O(A^4).

**Coefficient derivation.** Expand U_e = exp(i Σ_a A_e^a T^a) ≈ I + i A_e^a T^a + O(A^2). The plaquette holonomy at quadratic order is
```
U_P = I + i (Σ_e A_e^a) T^a − (1/2) (Σ_e A_e^a T^a)(Σ_{e'} A_{e'}^b T^b) + O(A^3) + (commutators).
```
Trace using tr(T^a) = 0, tr(T^a T^b) = (1/2) δ^{ab}:
```
tr U_P = 3 − (1/4) (A_P^a)² + O(A^3),     A_P^a = Σ_{e in P} A_e^a.
```
Hence:
```
1 − (1/3) Re tr U_P = (1/12) (A_P^a)² + O(A^4).      [cubic vanishes by reality]
```
Coefficient: **1/12**.

**General SU(N) formula.** With T^a = λ^a/2 normalized so that tr(T^a T^b) = (1/2) δ^{ab}, the per-plaquette per-component coefficient is
```
1 / (4 · N_c).
```
This gives 1/8 for SU(2) (Paper 30 §5.2), 1/12 for SU(3) (this paper), 1/16 for SU(4), etc. The N_c-dependence is structural: more colors = smaller per-component coefficient because the (1/N_c) Re tr normalization in the Wilson action divides by N_c while the trace normalization holds tr(T^a T^a) fixed at 1/2.

**Symbolic verification.** `debug/st_su3_theorem2_kinetic.py`:
```
sympy_kinetic_coefficient = Rational(1, 12)
tr T^a (a=1..8) = [0, 0, 0, 0, 0, 0, 0, 0]   (all traceless)
tr T^a T^b = (1/2) δ^{ab}                     (verified)
tr T^a T^b off-diagonal max = 0 exactly
```

**Numerical verification on Bargmann at N_max = 2.** With small random A^a coefficients (eps × N(0, 1) per edge per component, 15 edges × 8 components = 120 dof), build U_e = exp(i Σ_a A_e^a T^a), compute the actual S_W and compare to the quadratic prediction (1/12) Σ_P Σ_a (A_P^a)²:

| eps | S_predicted (quadratic) | S_actual (full) | rel error |
|:----|:------------------------:|:----------------:|:---------:|
| 0.05 | 6.95 × 10⁻² | 6.94 × 10⁻² | 1.0 × 10⁻³ |
| 0.01 | 2.78 × 10⁻³ | 2.78 × 10⁻³ | 7.6 × 10⁻⁵ |
| 0.001 | 2.78 × 10⁻⁵ | 2.78 × 10⁻⁵ | 1.4 × 10⁻⁵ (noise floor) |

Relative error scales as O(eps²) at moderate eps, consistent with the O(A⁴) leading correction. At eps = 0.001 the error hits the floating-point noise floor.

**Structural reading.** Paper 25's L_1 = B^T B is the *kinetic* sector of an SU(3) Wilson theory at weak coupling, with the same scalar form for each of the 8 su(3) components. Non-abelian self-interactions (cubic and quartic A-vertices using the SU(3) structure constants f^{abc}) enter only at O(A^3) and O(A^4) and use the [T^a, T^b] = i f^{abc} T^c commutators that vanished from the quadratic expansion. The 8 components are kinetic-decoupled at Gaussian order; SU(3) "color mixing" is purely an interaction effect, not a kinetic one.

---

## §5. Theorem 3 — CG obstruction status

**Question.** Sprint 5 Track S5 (`debug/s5_gauge_structure_memo.md`) showed that the *adjacency-preserving* SU(3) action on the (N, 0) representation tower fails Wilson's "fixed group on every link" requirement: transitions (N, 0) → (N+1, 0) are Clebsch-Gordan intertwiners between distinct SU(3) irreps, not group elements. Does the present Wilson construction bypass this obstruction, or rediscover it?

**Verdict: BOTH.** Two layers, two answers:

### 5.1 Gauge level: BYPASSED.

Wilson link variables U_e are external 3×3 SU(3) matrices in the *fundamental* representation (dimension 3, fixed). They are not representations of SU(3) acting on shell labels. Gauge transformations g_v ∈ SU(3) are also in the fundamental rep, so g_{src(e)} U_e g_{tgt(e)}† is well-defined. The construction works:

| N_max | Links in SU(3) | Plaquette holonomies in SU(3) | Gauge-invariant |
|:-----:|:--------------:|:------------------------------:|:----------------:|
| 2 | True | True (50/50 tested) | True (|S_diff| = 3.6e-15) |
| 3 | True | True (50/50 tested) | True (|S_diff| = 0) |

Sprint 5's negative concerned the *adjacency-preserving* SU(3) action on the (N, 0) tower — which Wilson does not require. The Wilson construction sidesteps the obstruction by putting SU(3) on link variables, not on representations of shell labels.

### 5.2 Matter coupling: REDISCOVERED.

To couple Wilson SU(3) to *matter* on the Bargmann graph nodes — i.e., to define a transformation rule ψ_v → g_v ψ_v — one must specify the SU(3) representation in which ψ_v lives at each node. The natural choice is to use the (N, 0) symmetric SU(3) irrep at each shell N, i.e., the Bargmann construction's intrinsic matter representation:

| N | shell dim (= dim (N,0)) | shell content |
|:-:|:------------------------:|:--------------:|
| 0 | 1 | (l=0, m=0) |
| 1 | 3 | (l=1, m=−1, 0, +1) |
| 2 | 6 | (l=0, m=0); (l=2, m=−2..+2) |
| 3 | 10 | (l=1, m=...); (l=3, m=...) |

The matter representations have **different dimensions** on different shells. Wilson links carry the fundamental 3-dim rep. To compute U_e ψ_v for a link e from shell N to shell N+1, one needs a map (3-dim) ⊗ (dim (N,0)) → (dim (N+1,0)). This is a *Clebsch-Gordan intertwiner* — exactly the structure Sprint 5 identified as the obstruction to a uniform Wilson SU(3) action on the (N, 0)-tower.

So matter coupling on the Bargmann graph requires the same CG-intertwiner structure that Sprint 5 found cannot be a Wilson lattice gauge theory in the strict sense. The Wilson SU(3) of this paper is naturally a *pure-gauge* theory: gauge dynamics without matter coupling to the (N, 0) shells.

### 5.3 Comparison with SU(2) on S^3 (Paper 30).

| Aspect | SU(2) on S^3 (Paper 30) | SU(3) on S^5 (this sprint) |
|:-------|:--------------------------|:----------------------------|
| Manifold = group? | Yes: S^3 = SU(2) | No: S^5 ≠ SU(3) |
| Manifold dim | 3 | 5 |
| Gauge group dim | 3 | 8 |
| Rank | 1 | 2 |
| Cartan torus | U(1) | U(1) × U(1) |
| Matter rep at node | spin-1/2 (dim 2) | (N, 0) symmetric (dim (N+1)(N+2)/2) |
| Link rep | spin-1/2 fundamental (dim 2) | fundamental (dim 3) |
| Matter rep matches link rep? | **YES** (both dim 2 spin-1/2) | NO |
| Natural matter coupling? | YES (Paper 14 Tier 2 Dirac fermions) | NO (CG intertwiner needed) |

The *S^3 = SU(2)* coincidence makes the SU(2) sibling structurally special: the underlying manifold *is* the gauge group, and matter (spin-1/2 fermions) lives in the same rep as the link variables. The S^5 case has neither property: S^5 is the unit sphere in ℂ^3, not SU(3); the matter representations on the (N, 0) shells have dimensions that grow with N, so a uniform link-rep matter coupling is impossible.

This is itself a **new dimension** of Paper 24's Coulomb/HO asymmetry. Paper 24 §IV established that S^3 carries a calibration-π content and a spectrum-computing graph Laplacian, while S^5 carries neither. Paper 30 added the gauge-content dimension: SU(2) Wilson on S^3 is natural and matter-coupled, while SU(2) Wilson would not transfer naturally to S^5. The present sprint shows that **even SU(3) Wilson on S^5 is "pure gauge" only**: the natural matter sector (the (N, 0) Bargmann shells) cannot be coupled in the standard Wilson sense without invoking CG intertwiners — Sprint 5's structural negative re-emerges from the matter side.

### 5.4 Three-tier asymmetry refinement.

Combining Sprint 5 Track S5, Paper 30, and this sprint, the Coulomb/HO asymmetry has now sharpened from one layer to **four layers**:

| Layer | S^3 Coulomb (Paper 7) | S^5 Bargmann (Paper 24) |
|:-----:|:----------------------|:-------------------------|
| (i) Spectrum-computing role of L_0 | YES (κ · (D−A) = Rydberg) | NO (HO spectrum in diagonal) |
| (ii) Calibration π | YES (Fock projection) | NO (linear-affine projection) |
| (iii) Wilson gauge with matter coupling | YES (SU(2) Wilson + Dirac) | NO (CG intertwiner needed) |
| (iv) Universal Wilson-Hodge vocabulary | YES (B, L_0, L_1, β₁) | YES (also well-defined) |

Layer (iv) is the universal content: any finite simplicial 1-complex with chosen orientation has B, L_0, L_1, β₁. Layers (i)-(iii) are Coulomb-specific. Sprint ST-SU3 contributes layer (iii) on the *gauge-with-matter* side: Wilson gauge alone is universal, but coupling Wilson gauge to the natural shell matter is Coulomb-specific.

---

## §6. Plaquette tabulation and Wilson loops

### 6.1 Plaquette counts (verified §2.1 table)

`debug/data/st_su3_plaquettes.json`. The Bargmann graph at N_max = 2, 3, 4 has many more plaquettes than the S^3 Coulomb graph at the same n_max:

| N_max | β₁ (Bargmann) | β₁ (S^3 Paper 30) | N_4 (Bargmann) | N_4 (S^3) |
|:-----:|:-------------:|:------------------:|:---------------:|:---------:|
| 2 | **6** | 0 | 15 | 0 |
| 3 | **23** | 2 | 88 | 2 |
| 4 | **56** | 8 | 272 | 8 |

The Bargmann graph has roughly an order of magnitude more cycles per shell than S^3 (consistent with §2.1 of Paper 24's edge count: 165 edges at N_max=5 vs S^3's 13 at n_max=3, reflecting the dipole multi-channel structure).

### 6.2 Wilson loops at N_max = 2

`debug/data/st_su3_wilson_loops.json`. Metropolis Monte Carlo with local SU(3) updates (auto-tuned scale ~ 1/√β), 600 samples after 200 thermalization, seed = 12, Wilson loop = first primitive 4-cycle:

| β | ⟨W⟩ (MC) | ± stderr |
|:--:|:---------:|:---------:|
| 0.5 | +0.047 | 0.010 |
| 1.0 | +0.100 | 0.010 |
| 2.0 | +0.260 | 0.009 |
| 5.0 | +0.654 | 0.007 |
| 10.0 | +0.885 | 0.002 |

Compare to Paper 30 SU(2) at n_max=3 with the same loop convention:
| β | ⟨W⟩_SU2 | ⟨W⟩_SU3 (this) |
|:--:|:--------:|:----------------:|
| 0.5 | 0.126 | 0.047 |
| 1.0 | 0.158 | 0.100 |
| 2.0 | 0.332 | 0.260 |
| 5.0 | 0.681 | 0.654 |

SU(3) Wilson loops are systematically lower than SU(2) at small-to-moderate β. This is consistent with general SU(N) lattice-gauge-theory expectations: more colors give more degrees of freedom to randomize the trace, suppressing the leading-order character expansion (cf. Drouffe-Zuber 1983 §3.3 for the standard SU(N) leading-order analysis). At large β (continuum limit), both approach 1.

### 6.3 Wilson loops at N_max = 3

| β | ⟨W⟩ (MC) | ± stderr |
|:--:|:---------:|:---------:|
| 0.5 | +0.268 | 0.022 |
| 1.0 | +0.581 | 0.013 |
| 2.0 | +0.851 | 0.005 |
| 5.0 | +0.926 | 0.002 |

⟨W⟩ at N_max=3 is *higher* than at N_max=2 at the same β. This is a **plaquette-count effect**: at N_max=3 the test loop sits inside a graph with 7,765 plaquettes (length ≤ 8) vs 182 at N_max=2, so the action gradient driving links to identity is much stronger per link. The plaquette of the test loop is "frozen near identity" by the surrounding sea of constraint plaquettes.

This is not confinement (no area-perimeter law family available; the Bargmann graph at this size has only one plaquette length class per scale). It is the standard *lattice volume effect*: dense graphs → more constraints → links pinned closer to identity at the same β.

### 6.4 Confinement diagnostic

As in Paper 30, area-perimeter law diagnostics require a continuous family of Wilson loops indexed by area. The Bargmann graph at N_max = 2, 3, 4 has many plaquettes but only a few distinct *cycle lengths* (4, 6, 8 dominate); a meaningful area-perimeter family is not naturally available at these sizes. The MC data at hand are **consistent with deconfined behavior** (⟨W⟩ smoothly interpolates 0 → 1 with β) but *cannot demonstrate* confinement either way at this graph size. This is a structural feature of small graphs, not a statement about the theory.

To probe confinement on the Bargmann graph, one would need either (a) much larger N_max (computationally feasible: N_max=5 has V=56, E=165, β₁=110 from Sprint 5; the MC update scales as |E| per sweep, so N_max=5 with ~165 links is ~10× more expensive than N_max=2's 15 links), or (b) a coarse-graining of the graph that produces a larger family of macroscopic loops.

---

## §7. Computational validation

### 7.1 Test suite

`tests/test_su3_wilson_s5.py`: **39 tests pass, 1 slow Monte Carlo test skipped by default** (collected on first run, ~3.8 seconds total). Test classes mirror Paper 30:

| Class | Tests | What it verifies |
|:------|:------:|:------------------|
| TestGellMann | 7 | λ_a hermitian, traceless, normalized; λ_3 and λ_8 diagonal |
| TestSU3Membership | 5 | Identity, Haar-random, exp(i x_a T^a), Cartan-diag all in SU(3) |
| TestPlaquetteHolonomy | 3 | Trivial, random plaquettes give SU(3); reverse-edge convention |
| TestGaugeInvariance | 2 | S_W invariant under random gauge; identity gauge no-op |
| TestCartanReduction | 4 | Diagonal-link character formula; reduction at machine precision |
| TestBargmannAdapter | 5 | Node counts, symmetry, no self-loops, bipartite |
| TestPlaquetteEnumeration | 5 | N_max=1 has none; bipartiteness; closure; no backtracking |
| TestWilsonActionStructure | 2 | Trivial vacuum gives S=0; non-negativity |
| TestWeakCouplingKinetic | 4 | Coefficient = 1/12 = 1/(4·N_c); SU(2)/SU(3) ratio = 2/3; quadratic expansion |
| smoke (N_max=2, 3) | 2 | End-to-end Bargmann + gauge invariance |
| TestMonteCarlo (slow) | 1 | Wilson loop monotonic in β |

No regression in SU(2) tests (`tests/test_su2_wilson_gauge.py`: 26/26 pass).

### 7.2 Symbolic verifications

`debug/st_su3_theorem1_cartan_reduction.py`, `debug/st_su3_theorem2_kinetic.py`, `debug/st_su3_theorem3_cg_obstruction.py` — three drivers, JSON outputs in `debug/data/st_su3_*.json`.

### 7.3 Production timing

| Operation | N_max | Time |
|:----------|:-----:|:----:|
| Build Bargmann adjacency | 4 | <1s |
| Enumerate plaquettes (L≤8) | 2 | 0.05s |
| Enumerate plaquettes (L≤8) | 3 | 2.2s |
| Enumerate plaquettes (L≤8) | 4 | 22.6s |
| MC Wilson loop (600 samples) | 2 | 0.2s |
| MC Wilson loop (400 samples) | 3 | 0.8s |
| Gauge invariance test (full) | 3 | <1s |

---

## §8. Recommended Paper edits

### 8.1 Paper 30 (priority: HIGH)

Paper 30 is the SU(2) sibling of this work; updating §1 (Introduction) and §7 (Open Questions) to cite the SU(3) sister sprint is the right cadence.

**Recommended edits to Paper 30:**

1. **Abstract or §1.** Add one sentence: "An SU(3) sibling on the Bargmann-Segal S^5 graph (Paper 24) has been constructed in Sprint ST-SU3 [memo: debug/st_su3_wilson_memo.md], confirming Theorem 1's analog (Cartan reduction to U(1) × U(1)) and Theorem 2's analog (kinetic coefficient = 1/12 = 1/(4·N_c)) at machine precision; matter coupling on S^5 encounters the Sprint 5 Track S5 CG-intertwiner obstruction, sharpening Paper 24's Coulomb/HO asymmetry to four layers."

2. **§7 Open questions.** Add a subsection §7.7 "SU(3) Wilson sibling on S^5": brief summary that the Wilson construction works at the gauge level on Bargmann (Theorem 1, 2 verified), matter coupling encounters CG obstruction (Theorem 3), and the four-layer Coulomb/HO asymmetry is the new structural reading (cite §5.4 of the ST-SU3 memo).

3. **§5.1 (current Cartan reduction).** Add a remark: "On a rank-r gauge group with maximal torus T = U(1)^r, the analog of this proposition gives r independent (but coupled, via det-=-1) U(1) sectors. The SU(3) case (r=2) is computed in Sprint ST-SU3."

4. **§5.2 (current kinetic term coefficient 1/8).** Add: "The general SU(N) coefficient is 1/(4 N_c) per plaquette per generator component, derived in Sprint ST-SU3 §4 for N_c = 3."

### 8.2 Paper 25 (priority: MEDIUM)

Paper 25 §VII.A explicitly flagged the SU(3) question as open. This sprint provides a concrete answer.

**Recommended edits to Paper 25:**

1. **§VII.1 (the open question).** Replace or extend with: "The question 'does Paper 25's U(1) Wilson-Hodge structure on the Hopf graph admit a non-abelian extension on the Bargmann S^5 graph?' has the following answer (Sprint ST-SU3, memo debug/st_su3_wilson_memo.md): YES at the gauge level (SU(3) Wilson with link variables in the fundamental rep is well-defined, gauge-invariant, has a Cartan reduction to U(1) × U(1), and a weak-coupling kinetic term L_1 with coefficient 1/12 = 1/(4 N_c)). NO at the matter-coupling level: Wilson SU(3) coupled to the natural (N, 0)-shell matter requires CG intertwiners (Sprint 5 Track S5's negative re-emerges from the matter side). Net verdict: a self-consistent pure-gauge SU(3) Wilson theory on the Bargmann S^5 graph exists; coupling it to matter sends one back to Sprint 5's CG-intertwiner structure."

2. **§VII.1 conclusion.** Note the four-layer Coulomb/HO asymmetry refinement: the universal Wilson-Hodge vocabulary transfers to S^5 (well-known); SU(3) Wilson gauge transfers to S^5 (new from this sprint); matter coupling does not (Sprint 5).

### 8.3 Paper 24 (priority: LOW)

Paper 24 §V (Coulomb/HO asymmetry) is the natural home for the four-layer refinement, but the existing Paper 25/30 cross-references already cover layers (i)–(iii). A brief footnote or pointer is sufficient.

**Recommended edit to Paper 24:**

1. **§V (Coulomb/HO asymmetry).** Footnote at the section opening: "The Coulomb/HO asymmetry has been refined to four layers across Sprint 5 Track S5 (gauge-content), Paper 30 (SU(2) Wilson on S^3), and Sprint ST-SU3 (SU(3) Wilson on S^5): (i) spectrum-computing role of L_0; (ii) calibration π; (iii) Wilson gauge with natural matter coupling; (iv) universal Wilson-Hodge vocabulary. Layers (i)–(iii) are Coulomb-specific; (iv) is universal. See debug/st_su3_wilson_memo.md §5.4."

### 8.4 CLAUDE.md (priority: MEDIUM, PI prerogative)

CLAUDE.md §2 currently has a long bullet on the S^5 sprint with the MIXED verdict from Sprint 5. The natural addition is a one-paragraph update appending the ST-SU3 result.

**Recommended addition (PM should NOT auto-apply this; PI to confirm):**

After the existing "S^5 gauge-structure extension (Sprint 5 Track S5)" bullet in §2:

"**SU(3) Wilson on S^5 Bargmann (Sprint ST-SU3, May 2026):** Tests the question Paper 25 §VII.A flagged as open and Sprint 5 Track S5 partially answered. Wilson lattice gauge theory with link variables in the fundamental SU(3) representation (NOT a representation of SU(3) on the (N,0) tower — categorically different from Sprint 5's adjacency-preserving action). (a) **Construction works at gauge level**: link variables in SU(3), plaquette holonomies stay in SU(3), action gauge-invariant under node-local SU(3), Haar-normalizable partition function. Verified on Bargmann at N_max = 2, 3 to machine precision. (b) **Cartan torus reduction to U(1) × U(1) is exact**: restricting links to T = U(1) × U(1) reduces S_W to the corresponding abelian Wilson action with character (1/3)(cos a_P + cos b_P + cos(a_P + b_P)). Verified at machine precision (max |Δ| < 1.5 × 10⁻¹¹ at N_max = 3, 7,765 plaquettes). (c) **Weak-coupling kinetic coefficient = 1/12 per plaquette per su(3) component**, generalizing SU(2)'s 1/8 via the formula 1/(4 N_c). Eight su(3) components decouple at quadratic order; non-abelian self-interactions enter at O(A^4). (d) **Sprint 5 CG obstruction status**: BYPASSED at gauge level (Wilson links don't act on shell labels), REDISCOVERED at matter-coupling level (matter on Bargmann shells lives in (N, 0) reps of growing dimension; coupling Wilson SU(3) to matter requires CG intertwiners — exactly Sprint 5's structure). Wilson SU(3) is a pure-gauge theory on Bargmann; matter sector is the same CG structure Sprint 5 found cannot be a Wilson theory in the strict sense. (e) **Plaquette table**: V, E, β₁ at N_max = 2,3,4 = (10,15,6), (20,42,23), (35,90,56); N_4 = 15, 88, 272. (f) **Wilson loops** at N_max = 2 with local SU(3) updates: ⟨W⟩(β) = 0.05, 0.10, 0.26, 0.65, 0.89 at β = 0.5, 1, 2, 5, 10; systematically lower than Paper 30's SU(2) at small β (more colors), consistent with deconfined behavior, area/perimeter diagnostic unavailable at this graph size. (g) **Sharpens the Coulomb/HO asymmetry to four layers**: spectrum-computing L_0 / calibration π / Wilson gauge with matter / universal Wilson-Hodge vocabulary; first three Coulomb-specific, fourth universal. **Status**: POSITIVE WITH STRUCTURAL CAVEAT. New module `geovac/su3_wilson_s5.py`, tests `tests/test_su3_wilson_s5.py` (39/39 pass, 1 slow skipped), data `debug/data/st_su3_{plaquettes,wilson_loops,theorem1_cartan,theorem2_kinetic,theorem3_cg_obstruction}.json`, memo `debug/st_su3_wilson_memo.md`."

---

## §9. Summary

| Item | Result |
|:-----|:--------|
| Module | `geovac/su3_wilson_s5.py` (≈ 480 lines) |
| Tests | `tests/test_su3_wilson_s5.py` — 39/39 pass, 1 slow MC skipped |
| Theorem 1 (Cartan → U(1) × U(1)) | VERIFIED at machine precision (|Δ| < 1.5e-11) |
| Theorem 2 (kinetic = 1/12 = 1/(4·N_c)) | VERIFIED symbolically + numerically (rel. err ~ O(eps²)) |
| Theorem 3 (CG obstruction) | BYPASSED at gauge level, REDISCOVERED at matter level |
| Plaquette tabulation | N_max = 1, 2, 3, 4 with β₁ and (N_4, N_6, N_8) |
| Wilson loops at N_max = 2 | β = 0.5, 1, 2, 5, 10 → ⟨W⟩ = 0.05 → 0.89 |
| Wilson loops at N_max = 3 | β = 0.5, 1, 2, 5 → ⟨W⟩ = 0.27 → 0.93 |
| Paper 30 edits recommended | §1 abstract, §7.7 (new), §5.1, §5.2 |
| Paper 25 edits recommended | §VII.1 (rewrite open question with answer) |
| Paper 24 edits recommended | §V footnote (four-layer asymmetry pointer) |
| CLAUDE.md §2 addition | Detailed Sprint ST-SU3 bullet (PI to confirm before applying) |

**Net structural reading.** The SU(3) Wilson construction on the Bargmann-Segal S^5 graph is a clean structural sibling of Paper 30's SU(2) Wilson on S^3 = SU(2). At the gauge level, both work the same way and admit the same three structural theorems (Cartan reduction, L_1 as kinetic term, plaquette enumeration). The asymmetries appear at the *matter-coupling* level: S^3 = SU(2) means matter and gauge live in the same Hilbert space (Paper 14 Tier 2), while on S^5 they are categorically separate, with the Bargmann shell tower carrying CG-intertwiner structure between irreps that the Wilson formalism does not naturally accommodate.

The four-layer Coulomb/HO asymmetry is the consolidated picture: layers (i) spectrum-computing L_0, (ii) calibration π, and (iii) Wilson gauge with matter coupling are all Coulomb-specific; layer (iv) the universal Wilson-Hodge vocabulary applies to both. Sprint ST-SU3 contributes layer (iii) on the gauge-with-matter side.

---

## §10. Files

**Module:** `geovac/su3_wilson_s5.py`

**Tests:** `tests/test_su3_wilson_s5.py`

**Drivers:**
- `debug/st_su3_theorem1_cartan_reduction.py`
- `debug/st_su3_theorem2_kinetic.py`
- `debug/st_su3_theorem3_cg_obstruction.py`
- `debug/st_su3_plaquettes_and_wilson.py`

**Data:**
- `debug/data/st_su3_theorem1_cartan.json`
- `debug/data/st_su3_theorem2_kinetic.json`
- `debug/data/st_su3_theorem3_cg_obstruction.json`
- `debug/data/st_su3_plaquettes.json`
- `debug/data/st_su3_wilson_loops.json`

**Memo:** `debug/st_su3_wilson_memo.md` (this file)
