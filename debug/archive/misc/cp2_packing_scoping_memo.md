# CP² Packing-Axiom Scoping Memo

**Date:** 2026-05-08
**Sprint:** post-Sprint-MH multi-track research launch (Track 3)
**Question:** Does the Paper 0 packing axiom admit a faithful analog on natural compact homogeneous spaces beyond S³ — specifically CP² and the SU(3) (p,q) irrep lattice? Does the construction *close* (well-defined discrete graph + Casimir spectrum + projection)?
**One-line verdict:** **PARTIAL CLOSURE.** The SU(3) (p,q) irrep lattice closes as a graph + Casimir + adjacency structure, and CP² hosts the (p,p)-symmetric sub-lattice with clean Laplace–Beltrami spectrum. **Neither closes from a packing-axiom analog with the same epistemic standing as Paper 0** — Paper 0's "binary distinguishability + isotropy" axioms are rank-1 (SO(3))-specific and do not lift to rank-2 (SU(3)) without becoming representation-theoretic. The construction closes mathematically; the *axiomatic motivation* does not transfer.

---

## §1. Geometry walk

Paper 0's four ingredients (extracted from the .tex, §II):

1. **Reference measure**: σ₀ = π·d₀²/2 (fundamental area per state)
2. **Scaling rule**: r_k = k·d₀ (linear concentric shells)
3. **Induction rule**: N_k = annular_area / σ₀ = 2(2k−1) per shell
4. **Incidence rule**: radial T± edges (Δn = ±1 at fixed l, m) and angular L± edges (Δm = ±1 at fixed n, l)

**Output**: discrete (n, l, m) labeling on a 2D Euclidean plane, with cumulative count 2n² through n shells (factor of 2 from orientability of the compactified S² = CP¹), conformally lifted to S³ via the Fock projection p₀ = √(−2E).

The crucial observation: Paper 0's packing is in **2D** (a plane), produces **rank-1** angular structure (the SO(3) shell labeling l with multiplicity 2l+1), and only acquires **rank-2** structure (the (n, l) tower) via the *cumulative grouping* of shells — which §IV.A of Paper 0 explicitly flags as a structural choice "natural but not uniquely forced by the axioms." S³ enters as the dynamics carrier through Fock, not as the packing manifold itself.

For an SU(3) analog, the relevant rank is **2** (SU(3) has rank 2; its weight lattice is 2D). The dimensional escalation ladder is:

| Group | Rank | Original packing manifold | Dynamics manifold (Fock-class) |
|:------|:-----|:--------------------------|:-------------------------------|
| SO(3) ↔ SU(2) | 1 | R² (= CP¹ compactified) | S³ |
| SU(3) | 2 | R⁴? CP²? | S⁵ Hardy / CP² / unknown |
| SU(4) | 3 | R⁶? CP³? | unknown |

The Bertrand × Hopf-tower truncation theorem (CLAUDE.md §2 "Sprint of 2026-05-07"): U(1)×SU(2)×SU(3) saturate the SM gauge content as a structural upper bound, forced by Bertrand's classical-mechanics closed-orbit theorem (only −Z/r and r² admit all bound orbits closed) crossed with the complex-Hopf S^(2n−1) → SU(n) tower truncated at n ≤ 3. Beyond SU(3) the Hopf bundle structure (S⁷) lacks Lie-group structure (no octonionic group of all bound orbits). This is *gauge-content* truncation; whether it constrains *packing* analogs is a separate question this memo addresses.

### Geometry table

| Candidate | Group | Dim | Rank | Casimir spectrum | Spherical reps | Ingredient closure | Overall verdict |
|:----------|:------|:----|:-----|:-----------------|:---------------|:-------------------|:----------------|
| **CP² = SU(3)/U(2)** | SU(3) | 4 | 2 | λ_p = p(p+2) | (p,p) symmetric | (1)✓ (2)✓ (3)~ (4)✓ | **closes geometrically; axiom motivation absent** |
| **SU(3) (p,q) irrep lattice** | SU(3) | 2D ℤ²₊ | 2 | C₂(p,q) = (p²+q²+pq)/3+(p+q) | all (p,q) | (1)~ (2)✓ (3)✓ (4)✓ | **closes as graph; rep-theoretic, not packing** |
| S⁶ = G₂/SU(3) | G₂ | 6 | 2 | spherical-rep tower | restricted | (1)? (2)✓ (3)~ (4)✓ | partial — exceptional, no SM correspondence |
| HP¹ = S⁴ = Sp(2)/Sp(1)² | Sp(2) | 4 | 2 | λ_k from rank-1 symmetric pair | symmetric (k,0) | (1)~ (2)✓ (3)~ (4)✓ | closes; quaternionic structure, no SM home |
| CP³ = SU(4)/U(3) | SU(4) | 6 | 3 | λ_p = p(p+3) | (p,0,p) symmetric | (1)? (2)✓ (3)~ (4)✓ | closes; Penrose-twistor, SU(4) not in SM |
| Full flag SU(3)/T | SU(3) | 6 | 2 | all (p,q) reps appear | all reps | (1)~ (2)✓ (3)✓ (4)~ | closes; Schubert calc., no clean dynamics |
| OP² = F₄/Spin(9) | F₄ | 16 | 4 | rank-1 symmetric pair | (k,0,0,0) tower | (1)? (2)✓ (3)? (4)✓ | closes; octonionic, magic-square SUGRA |
| Grassmannians SU(N)/(SU(k)×SU(N-k)×U(1)) | SU(N) | 2k(N-k) | min(k,N-k) | Howe duality | hook-shaped reps | varies | varies; CP^{N-1} subcase clean |

**Closure column legend:** ✓ = ingredient extends naturally; ~ = ingredient extends with non-axiomatic input; ? = ingredient requires substantial new work; ✗ = ingredient demonstrably fails.

The pattern: rank-2 (and higher) compact symmetric spaces all admit Casimir spectra and Laplace–Beltrami operators by rep theory (Helgason 1984, Camporesi–Higuchi 1996). The construction *always* closes mathematically. The question is whether any of them admit a **packing-style motivation** — an axiomatic input parallel to "binary distinguishability + isotropy" that produces the construction without invoking representation theory at the start.

The honest answer this memo defends: **none of the candidates do.** Every higher-rank candidate requires representation-theoretic input either at the "scaling rule" step (which reps are allowed) or at the "induction rule" step (how reps grow). Paper 0's axiomatic content is uniquely 2D-isotropy specific.

What this *doesn't* mean: the construction is useless. The (p,q) lattice and CP² spherical-rep spectrum are well-defined, computable, and have content as discrete graph constructions. They are simply not "packings" in Paper 0's strict sense. They are *representation-theoretic discrete frameworks on natural geometries*, which is a weaker but still useful class.

---

## §2. CP² detailed analysis

### §2.1 Manifold and group action

CP² = SU(3)/U(2) (more precisely SU(3)/S(U(2)×U(1))). Real dimension 4. Kähler with Fubini–Study metric — the unique (up to scale) SU(3)-invariant metric. Curvature constant non-zero positive. The pair (SU(3), U(2)) is a compact symmetric pair of rank 1 (the "spherical functions" form a single index family).

### §2.2 Spherical representations

For a symmetric pair (G, K), an irrep V_λ is **spherical** if it admits a non-zero K-invariant vector. For (SU(3), U(2)), the spherical reps are exactly the symmetric (p,p)-reps with both Dynkin labels equal:

V_(p,p) for p = 0, 1, 2, ...

The associated dimensions and Casimirs (numerator-over-3 convention, so C₂ = N/3):

| p | rep label | dim(p,p) = (p+1)³ | C₂ = p(p+2) |
|:-:|:---------:|:-----------------:|:-----------:|
| 0 | (0,0) trivial | 1 | 0 |
| 1 | (1,1) adjoint | 8 | 3 |
| 2 | (2,2) | 27 | 8 |
| 3 | (3,3) | 64 | 15 |
| 4 | (4,4) | 125 | 24 |
| 5 | (5,5) | 216 | 35 |

The eigenvalues p(p+2) are **structurally identical** to the Laplace–Beltrami spectrum of S³ at level n = p (where S³ eigenvalues are n(n+2) with multiplicity (n+1)²). Multiplicities differ: CP² has (p+1)³, S³ has (p+1)². The extra factor of (p+1) is the dimension of the additional U(1) direction in the Hopf-fibration-like CP² structure.

### §2.3 Discrete graph candidate construction

Paper 0's packing on R² produces (l, m) labels on concentric circles. The CP² analog would discretize CP² into "cells" labeled by spherical reps (p,p), with state count (p+1)³ at level p.

**Reference measure (ingredient 1):** SU(3)-invariant Fubini–Study volume on CP², normalized so that the smallest non-trivial cell — the adjoint (1,1) rep, dim 8 — sets the scale. This gives σ₀ = Vol(CP²)/8. Compare to Paper 0 where σ₀ = π·d₀²/2 = Vol(initial annulus)/2. The "8" arises from the dimension of the smallest non-trivial spherical rep, which is structurally analogous to the "2" arising from N_init in Paper 0. **But 8 is not derived from a "minimum to define a metric" argument** — it comes from rep theory. *Ingredient 1 closes with rep-theoretic input but not from a packing axiom.*

**Scaling rule (ingredient 2):** r_p = p (level index, by analogy to k in Paper 0). Cell volume at level p: scales as p³ (since dim(p,p) = (p+1)³ ~ p³ for large p). Compare to Paper 0 where shell-k angular capacity scales as 2k−1 ~ 2k. The CP² version has a **cubic** growth law versus the **linear** growth law of Paper 0 — this is the "rank-2 dimensional upgrade." *Ingredient 2 closes; growth law follows from spherical-rep dimensions.*

**Induction rule (ingredient 3):** Tensor of (p,p) with the adjoint (1,1) decomposes as
(1,1) ⊗ (p,p) = (p+1, p+1) ⊕ (p−1, p−1) ⊕ (p, p) ⊕ (p+1, p−1) ⊕ (p−1, p+1)
(symbolic — exact decomposition by Littlewood–Richardson; the off-diagonal (p+1, p−1) and (p−1, p+1) terms are non-spherical and would not descend to CP²). Restricting to spherical reps, the induction rule (1,1)-tensor sends level p → p+1 and p → p−1 with two surviving channels each. *Ingredient 3 closes if we accept "tensor with adjoint" as the induction operator, but this is not the packing-axiom analog of "annular area / σ₀" — it is rep-theoretic.*

**Incidence rule (ingredient 4):** Two spherical reps (p,p) and (p',p') are adjacent if (p',p') appears in (1,1) ⊗ (p,p). This gives a path graph on the spherical-rep tower (Δp = ±1 only when restricted to spherical reps). *Ingredient 4 closes; the resulting graph is a 1D path on the (p, p) tower, structurally identical to the n-radial chain in S³ packing.*

### §2.4 Critique: CP² closes geometrically but loses packing motivation

Three observations:

(a) **The CP² spherical-rep tower is structurally a 1D chain** (path graph indexed by p). All the "extra" SU(3)-rank-2 structure has been projected away by the spherical-restriction. The (p,p)-tower reproduces the SO(3) (l)-tower of Paper 0's per-shell structure, with eigenvalues p(p+2) instead of l(l+1) and multiplicities (p+1)³ instead of 2l+1. The construction is *more like a deformation of Paper 0's S³ embedding* than a genuinely new packing.

(b) **The eigenvalues p(p+2) match S³ shape but coefficients are tied to rep dimensions, not packing density.** In Paper 0, the eigenvalues come from the graph Laplacian on the (n, l, m) graph at K_vac = −1/16 (one universal kinetic scale). In a CP² analog, the eigenvalues are Casimir values from rep theory; there is no analog of K_vac that emerges from a "fundamental scale" axiom. *The κ = −1/16 structure (CLAUDE.md §2 v2.26.1 derivation) does not lift.*

(c) **The Fock-projection target is unclear.** S³ is the dynamics carrier for the (n, l, m) graph because the Fock projection at p₀ = √(−2E) maps Coulomb hydrogen onto the Laplace–Beltrami operator on unit S³ (Paper 7). For the CP² spherical-rep tower, what is the dynamics carrier? Naively CP² itself, but then we need a **Coulomb-class potential on CP²** whose bound-state spectrum aligns with the (p,p)-tower. The natural candidate is the CP² scalar Laplacian itself (eigenvalue p(p+2) by construction), which describes a free particle on CP² — no analog of "Coulomb attraction." A candidate analog would be the CP² geodesic harmonic oscillator, whose spectrum is presumably p(p+2) shifted, but this is HO-class (Bertrand's other closed-orbit potential), not Coulomb-class.

### §2.5 Physics correspondence candidates

(a) **Adjoint at p=1, dim 8**: This is the **gluon octet** of QCD at the fixed-Z perturbative level (or equivalently the Gell-Mann meson / baryon octet of the eightfold way).
(b) **(2,2) at dim 27**: Symmetric two-adjoint product, appears in gluon-gluon scattering.
(c) Higher (p,p) reps: less obvious physical correspondence.

These correspondences are **not predictions** — SU(3) flavor symmetry of the eightfold way *is by definition* the SU(3) rep lattice, restricted to specific reps that match the observed hadron multiplets. The CP² spherical-rep tower hosts these multiplets only because it inherits SU(3) rep theory, not because the packing construction predicted them. This is a tautological hosting, useful for organizing existing structure but not for predicting new physics.

**One non-tautological observation worth flagging:** the (p+1)³ multiplicity at level p exceeds the (p+1)² multiplicity on S³ at the same eigenvalue. The "extra factor of (p+1)" is the dimension of the U(1) Hopf direction. If this extra factor has dynamical content — e.g., describing a U(1) gauge degree of freedom co-located with the SU(3) structure — then CP² could host the **U(1) hypercharge × SU(3) color** structure of the SM strong sector at the kinematic level. This is speculative and not verified by this scoping; it is flagged as a follow-up question.

---

## §3. SU(3) (p,q) irrep lattice

### §3.1 Lattice structure

Lattice points: ℤ²_{≥0} = {(p, q) : p, q ≥ 0 integers}. Each point carries an SU(3) irrep V_(p,q) of dimension dim(p,q) = (p+1)(q+1)(p+q+2)/2. Casimir C₂(p,q) = (p² + q² + pq)/3 + (p + q) (physicist convention; multiply by 3 for integer numerator).

Explicit table through level k = p+q ≤ 4:

| (p,q) | level k | dim | 3·C₂ | physical name |
|:-----:|:-------:|:---:|:----:|:--------------|
| (0,0) | 0 | 1 | 0 | trivial / vacuum |
| (1,0) | 1 | 3 | 4 | quark triplet (fund) |
| (0,1) | 1 | 3 | 4 | antiquark triplet (anti-fund) |
| (2,0) | 2 | 6 | 10 | symmetric two-quark |
| (1,1) | 2 | 8 | 9 | **adjoint / octet** (Gell-Mann eightfold way) |
| (0,2) | 2 | 6 | 10 | antisymmetric two-antiquark |
| (3,0) | 3 | 10 | 18 | **decuplet** (Δ, Σ*, Ξ*, Ω⁻) |
| (2,1) | 3 | 15 | 16 | mixed |
| (1,2) | 3 | 15 | 16 | mixed |
| (0,3) | 3 | 10 | 18 | anti-decuplet |
| (4,0) | 4 | 15 | 28 | symmetric 4-quark |
| (3,1) | 4 | 24 | 25 | mixed |
| (2,2) | 4 | 27 | 24 | "27" rep |
| (1,3) | 4 | 24 | 25 | mixed |
| (0,4) | 4 | 15 | 28 | symmetric 4-antiquark |

Cumulative state count through level K:

| K | level total Σ dim(p,q) for p+q=K | cumulative through K | formula (K+1)(K+2)²(K+3)/12 |
|:-:|:---------------------------------:|:--------------------:|:----------------------------:|
| 0 | 1 | 1 | 1 |
| 1 | 6 | 7 | 6 |
| 2 | 20 | 27 | 20 |
| 3 | 50 | 77 | 50 |
| 4 | 105 | 182 | 105 |

The closed form Σ_{p+q=K} dim(p,q) = (K+1)(K+2)²(K+3)/12 is computed in `debug/data/cp2_packing_scoping_results.json` (verified by direct sum). Compare to S³ where level-n total is (n+1)² (same dim notation: each (n+1)²-dim level summed gives Σ_{n=0}^{K} (n+1)² = (K+1)(K+2)(2K+3)/6 ~ K³/3).

The SU(3) lattice has level-total ~ K⁴/12, growing one degree faster than S³'s K³/3. This is the rank-2 dimensional upgrade.

### §3.2 Adjacency from fundamental tensor

(1,0) ⊗ (p,q) = (p+1, q) ⊕ (p−1, q+1) ⊕ (p, q−1) (drop reps with negative labels)

(0,1) ⊗ (p,q) = (p, q+1) ⊕ (p+1, q−1) ⊕ (p−1, q) (drop reps with negative labels)

Combining both gives the natural 6-neighbor "hexagonal" connectivity at interior points (p, q ≥ 1). Specifically: (p±1, q), (p, q±1), (p+1, q−1), (p−1, q+1).

This is precisely the **triangular lattice** structure the PI's QCD intuition pointed at. Each point has up to 6 neighbors; boundary points (on p=0 or q=0 axis) have fewer. The lattice IS the SU(3) weight diagram extended to infinity, with the natural "tensor-with-fundamental" adjacency. This is well-known from the eightfold-way construction (see e.g. Georgi 1999 *Lie Algebras in Particle Physics*).

### §3.3 Comparing the four packing ingredients to Paper 0

**Ingredient 1 (reference measure):**
- Paper 0: σ₀ = π·d₀²/2 from Axiom 1 (binary distinguishability fixes N_init = 2 → fundamental area).
- SU(3) lattice: each lattice point (p,q) carries multiplicity dim(p,q). The "fundamental" cell could be the trivial (0,0) (dim 1) or the fundamental (1,0) (dim 3 — analog of "binary 2 → ternary 3" because SU(3) fundamental rep has dim 3, vs SU(2) fundamental dim 2). The latter is the natural "minimum non-trivial" but the **"minimum to define a metric" argument does not transfer** — what's the analog of "two points to define a distance"? Three points to define a 2D simplex? Not quite the same epistemic status.
- *Verdict: closes with rep-theoretic input, but the axiomatic motivation does not lift faithfully.*

**Ingredient 2 (scaling rule):**
- Paper 0: r_k = k·d₀ (linear), giving annular area ∝ 2k−1.
- SU(3) lattice: lattice points at integer (p,q) with no metric scaling on the lattice itself — the "scale" is the irrep dimension, dim(p,q) = (p+1)(q+1)(p+q+2)/2. *Linearity in p, q gives quadratic-in-level state counts; cubic at higher product.*
- *Verdict: closes; the polynomial degree of state counts upgrades from "linear in shell index" (SO(3)) to "quadratic in level index" (SU(3)).*

**Ingredient 3 (induction rule):**
- Paper 0: N_k = annular area / σ₀ = 2(2k−1).
- SU(3) lattice: at level k, the (k+1) distinct (p,q) reps with p+q=k contribute (k+1)(k+2)²(k+3)/12 total states. The "induction" can be rephrased: at level k+1, the new lattice points are obtained by tensor-with-fundamental from level k.
- *Verdict: closes via tensor-with-fundamental; this is structurally a Pieri/Littlewood–Richardson rule rather than an "annular area / unit area" geometric argument.*

**Ingredient 4 (incidence rule):**
- Paper 0: T± and L± edges (radial Δn=±1, angular Δm=±1).
- SU(3) lattice: (1,0) and (0,1) tensor adjacency, giving 6-neighbor hexagonal connectivity.
- *Verdict: closes cleanly; this is the most natural ingredient transfer.*

### §3.4 Comparison to S³ packing structure

| Feature | S³ packing (Paper 0) | SU(3) (p,q) lattice |
|:--------|:---------------------|:---------------------|
| Rank of internal symmetry | 1 (SU(2)) | 2 (SU(3)) |
| Lattice dimensionality | 1D shells (k indexed) + cumulative n | 2D points (p, q) |
| Per-level structure | flat: 2k−1 angular slots, all eigenvalue n²−1 | stratified: (k+1) distinct reps, multiple Casimirs |
| Casimir vs level | flat in level n: λ_n = n(n+2) | not flat in level k: C₂(p,q) varies within level |
| Adjacency type | T± (radial) + L± (angular) | 6-neighbor hexagonal (fund + antifund tensor) |
| Cumulative state count | 2n² (with spin) | (k+1)(k+2)²(k+3)/12 |
| Fock-projection target | S³ (Coulomb dynamics) | unclear (CP², S⁵ Hardy, or unknown) |
| Universal kinetic scale | K_vac = −1/16 (Paper 7 derivation) | no analog; would require new derivation |
| Physical hosts | hydrogen spectrum, periodic table | eightfold way / decuplet / 27 / ... (tautological) |

The structural difference centers on **stratification within levels**. In S³, every state at level n has the same Casimir n²−1 (the angular sub-structure (l, m) has 2l+1 multiplicity but doesn't affect the Casimir). In SU(3), level k = p+q hosts (k+1) distinct (p,q) reps with multiple Casimirs. This means the SU(3) lattice does NOT have the "shell" structure of Paper 0 — its levels are mixtures of representations rather than single irreducible "blocks."

This is consistent with rank: rank 1 has 1D weight diagrams (shells are flat); rank 2 has 2D weight diagrams (shells are stratified by 2nd Casimir).

### §3.5 Closure assessment

The SU(3) (p,q) lattice closes as a graph with:
- ✓ A well-defined infinite vertex set (all (p,q) ∈ ℤ²_{≥0})
- ✓ A well-defined Casimir spectrum (rational values)
- ✓ A well-defined adjacency rule (fundamental + antifundamental tensor)
- ✓ A clean physical labeling that hosts SU(3) rep theory

It does NOT close as a *packing* in Paper 0's sense:
- ✗ The "binary distinguishability" axiom doesn't directly transfer (SU(3) has no "two states define a metric" geometric argument)
- ✗ The "isotropy → uniform density" axiom doesn't directly transfer (SU(3) doesn't have a 2D-isotropy analog; it has SU(3)-equivariance, which is rep-theoretic)
- ✗ The scaling rule (irrep dim formula) is not derived from a geometric area-counting argument; it comes from Weyl's character formula
- ~ The induction rule via tensor-with-fundamental is rep-theoretic, not geometric

**Bottom line: the SU(3) (p,q) lattice is a real, useful, geometric-and-combinatorial object. It is NOT the SU(3) analog of Paper 0's packing in the strong axiomatic sense.** It is the analog in a weaker sense: "natural discrete graph on the rank-2 weight lattice with Casimir spectrum and tensor adjacency."

---

## §4. Cross-comparison and verdict

### §4.1 Summary of closure status

| Construction | Mathematically well-defined | Casimir spectrum | Discrete graph | Fock-projection analog | Packing-axiom motivation |
|:-------------|:---------------------------:|:----------------:|:--------------:|:----------------------:|:------------------------:|
| S³ (Paper 0 baseline) | ✓ | n(n+2) flat | (n,l,m) | S³ via Fock | ✓ (Paper 0 §II) |
| CP² spherical reps | ✓ | p(p+2) flat | (p,p)-tower path | unclear | ✗ (rep-theoretic) |
| SU(3) (p,q) lattice | ✓ | C₂(p,q) stratified | hexagonal triangular | unclear | ✗ (rep-theoretic) |

### §4.2 Where the obstruction sits

Paper 0's axioms are **rank-1 specific**:
- "Binary distinguishability" → 2 states define 1 distance → SU(2) double cover Z₂ structure
- "Geometric indifference" → 2D isotropy → SO(2) rotational symmetry → SO(3)/SU(2) Lie group emergence
- "Concentric shell with area-proportional counting" → 1D radial direction + circular shell structure → integer angular momentum from Bohr–Sommerfeld-like quantization

Each step is *uniquely* aligned with rank 1. The dimensional reductions (2D plane → S² compactification → S³ Fock dynamics) all preserve rank 1 throughout.

The Coulomb/HO asymmetry (Paper 24 §V, four-layer asymmetry of Sprint ST-SU3 May 2026) provides another reading of the same obstruction. The S³ and S⁵ Hardy spaces are the natural "Bertrand-closed-orbit" packings (one for −Z/r, one for r²) at rank 1 and rank 2 respectively. **But neither of those is a packing analog in the SO(3) sense — they are dynamics carriers for closed-orbit central potentials.** A "packing" analog for SU(3) would need a packing-style geometric construction that produces SU(3) labels, and Bertrand-Hopf-tower truncation explicitly says the only natural packings in this family are S¹/S³/S⁵ for U(1)/SU(2)/SU(3) gauge. These are gauge-host manifolds, not packing manifolds.

This is the key structural finding: **gauge-host manifolds (S¹, S³, S⁵) are not packing manifolds** in Paper 0's axiomatic sense. The Hopf-tower truncation confirms the gauge content is saturated; it does NOT produce a packing-axiom analog for SU(3).

### §4.3 Coulomb/HO category mismatch as obstruction

Building on the four-layer Coulomb/HO asymmetry of Paper 24 §V:

| Layer | S³ (SO(3) packing) | S⁵ Hardy (HO Bargmann–Segal) | SU(3) (p,q) lattice |
|:-----:|:------------------:|:---------------------------:|:-------------------:|
| (i) Spectrum-computing role of L_0 | ✓ via Coulomb Fock | ✗ | ✗ (Casimir comes from group theory directly, no L_0 ↔ Casimir conformal equivalence) |
| (ii) Calibration π | ✓ (κ = −1/16, K = π(B+F−Δ)) | ✗ (π-free certificate, Paper 24 §III) | ✗ (rationals from rep theory) |
| (iii) Wilson-Hodge with natural matter coupling | ✓ (Paper 30 SU(2) on S³, matter = electrons) | partial (Sprint ST-SU3: gauge YES, matter NO via CG obstruction) | unclear |
| (iv) Combinatorial vocabulary (universal) | ✓ | ✓ | ✓ |

The (p,q) lattice sits at a fifth category position: it has the universal combinatorial vocabulary (layer iv) but lacks layers i, ii, iii. The construction closes mathematically (layer iv is enough for that) but lacks the Coulomb-projection structure (layers i, ii) that gives Paper 0's S³ packing its axiomatic content and the Wilson-matter coupling (layer iii) that gives the gauge-content saturation its physical content.

**This is a genuine Coulomb/HO-style category obstruction**: the (p,q) lattice and CP² spherical-rep tower are at category positions where no Bertrand-class central potential lives, so they cannot inherit the Fock-projection / κ / calibration-π structure that S³ does.

### §4.4 Verdict

**Primary candidates:**
- **CP²**: closes geometrically (Laplace–Beltrami spectrum p(p+2), multiplicity (p+1)³); does NOT close with packing-axiom motivation. Useful as a discrete graph; not a packing analog.
- **SU(3) (p,q) lattice**: closes as a hexagonal-connectivity graph with Casimir spectrum; does NOT close with packing-axiom motivation. Useful for organizing SU(3) rep theory; not a packing analog.

**Secondary candidates** (all closures rep-theoretic; no genuine packing-axiom transfer):
- S⁶, HP¹, CP³, full flag, OP², Grassmannians: all close mathematically as discrete graphs with Casimir spectra; none provide a packing-axiom analog.

**Negative signal confirmed:** The four-layer Coulomb/HO asymmetry resurfaces at every higher-rank candidate. Calibration π content and Fock-projection structure are S³-specific (Paper 23 §IV Fock rigidity theorem; Paper 24 §V four-layer asymmetry). The packing-axiom MOTIVATION does not transfer beyond rank 1; only the combinatorial structure (layer iv) does.

**Bonus signal partial:** SU(3) (p,q) lattice tautologically hosts the eightfold way (octet at (1,1)), decuplet ((3,0)), and 27-rep ((2,2)). This is structural matching, not prediction — the eightfold way *is* SU(3) rep theory by construction (Gell-Mann 1962). The construction matches because the construction inherits the symmetry whose physics the eightfold way encodes.

**One non-tautological observation flagged for follow-up:** The CP² spherical-rep tower has multiplicity (p+1)³ versus S³'s (p+1)². The extra (p+1) factor is the dimension of the Hopf U(1) direction. If this U(1) direction has dynamical content beyond representation labeling, CP² could potentially host **U(1) hypercharge × SU(3) color** structure at the kinematic level — a candidate "QCD packing" target. This is speculative, would require a separate sprint to verify, and does not affect the primary verdict.

### §4.5 Re-ranking of multi-track sprint priorities

Per the directive: "If Track 3 lands well, the natural next sprint is testing whether the resulting graph hosts known physics ... If Track 3 *doesn't* land well, split focus on the other two going forward."

**Track 3 verdict: PARTIAL CLOSURE — does NOT land well in the strong sense.**

The construction closes mathematically but the packing-axiom motivation does not transfer. Per the directive's exit criteria ("Lands well" requires "the four packing ingredients can be written down for at least one candidate, the construction CLOSES, and small examples can be computed explicitly"), the small-examples and graph-closure conditions are met, but the four ingredients only "extend" via representation-theoretic input rather than packing-axiom input.

**Recommendation:** Do NOT prioritize CP² packing as the next major track. Continue with Tracks 1 (multi-focal precision catalogue) and 2 (chemistry solver re-test) as the active fronts.

The (p,q) lattice and CP² constructions remain available as **secondary tools** for organizing SU(3)-related observations within the existing framework — they are useful catalogue-style discrete graphs even though they don't come from a packing axiom. They could be invoked, e.g., as the natural discrete substrate for any SU(3)-flavor-symmetric extension of GeoVac, in much the way that the second-row frozen-core treatment uses Clementi-Raimondi screening as a calibrated input rather than a derived axiom.

---

## §5. Follow-up sprint plan (conditional)

Per directive: "If a candidate lands: scoped follow-up sprint plan for testing physics correspondence."

The verdict is partial-closure-not-landing, so a follow-up sprint is **not recommended as a top-priority track**. However, two focused sub-sprints are flagged for opportunistic execution if the precision-catalogue (Track 1) or chemistry-solver (Track 2) sprints uncover SU(3)-flavor-symmetry-related questions:

### §5.1 Sub-sprint A: CP² Hopf-direction dynamical content (1-2 weeks, opportunistic)

**Question:** Does the extra (p+1) multiplicity factor on CP² (relative to S³) carry dynamical content, or is it a kinematic redundancy?

**Method:**
- Compute the CP² scalar Laplacian explicitly at small cutoff
- Identify the Hopf-fibration U(1) direction explicitly in the (p,p) eigenspace decomposition
- Check whether projecting onto the U(1)-trivial subspace recovers the S³ structure exactly (in which case the extra factor is redundant) or produces a non-trivial enrichment (in which case it carries dynamical content)
- Cross-reference with Paper 25 (Hopf bundle structure) and the Paper 30 SU(2) Wilson construction

**Exit criterion:** If the Hopf direction is dynamically inert, the construction reduces to S³ packing × U(1) and adds nothing new. If it has content, scope a second sub-sprint to test whether it hosts hypercharge structure.

### §5.2 Sub-sprint B: SU(3) (p,q) lattice as catalogue substrate for QCD-flavor extensions (2-4 weeks, opportunistic)

**Question:** If GeoVac extends to QCD flavor at some future point, does the (p,q) lattice provide the natural discrete substrate for the flavor sector, in the way the S³ (n,l,m) graph provides the substrate for the electromagnetic sector?

**Method:**
- Build a small (k_max=4) (p,q) lattice graph with Casimir-weighted node weights
- Compute the analog of the K_vac kinetic scale (if any)
- Compare to known QCD chiral perturbation theory at the meson sector
- Identify what would have to be calibrated externally vs derivable from the lattice

**Exit criterion:** If the lattice can be used to compute *anything* relevant to QCD flavor (e.g., Gell-Mann–Okubo mass relations, meson mixing patterns) without ad-hoc calibration beyond the SU(3) symmetry itself, this is a positive signal. If everything must be calibrated externally (as expected from this scoping memo), the lattice is a labeling tool, not a predictive framework.

### §5.3 What this scoping rules out

- **Don't pursue** a "second packing axiom" that would generate SU(3) from a higher-rank version of Paper 0's two axioms. Rank-1 specificity of the axioms is structural; this is not a missing piece, it is the identity of Paper 0.
- **Don't pursue** reframing CP² or the (p,q) lattice as the "next" GeoVac geometry. They are not natural geometries in the Paper 7 / Paper 13 / Paper 15 hierarchy sense (no separation of variables of a Coulomb-class Schrödinger equation lands on them).
- **Don't pursue** a Fock-projection analog for CP² unless a Coulomb-class central potential on CP² is identified first. The Bertrand × Hopf-tower truncation suggests no such potential exists in the natural class.

### §5.4 What this scoping leaves open as legitimate frontier

- The "second packing axiom" question (CLAUDE.md §1.7 W3, §13 inner-factor calibration) is unchanged by this scoping. Whatever generates inner-factor calibration data does NOT need to be a packing-axiom analog on a higher-rank geometry. It could be something else entirely (e.g., a refinement of the original Paper 0 axioms, a new inductive structure, or a non-geometric data source).
- The Coulomb/HO asymmetry sharpens: at each layer (spectrum-computing L_0, calibration π, Wilson-with-matter, combinatorial vocabulary), the higher-rank candidates fail at different layers. CP²/(p,q) sit at layer-iv-only — combinatorial vocabulary without dynamics structure. This is a structurally distinct category from S³ (layers i-iv) and S⁵ Hardy (layers iii-partial, iv).

---

## Appendix: Numerical data

Full computed tables in `debug/data/cp2_packing_scoping_results.json`. Highlights:

- SU(3) (p,q) lattice up to level k=4: 15 lattice points, 182 cumulative states
- CP² spherical reps up to p=5: dimensions (1, 8, 27, 64, 125, 216), eigenvalues (0, 3, 8, 15, 24, 35)
- Adjacency tabulated for fundamental and antifundamental tensor neighbors at each (p,q)
- Eightfold way correspondence: octet at (1,1), decuplet at (3,0), 27 at (2,2) — all tautological

---

## TL;DR

**Verdict: PARTIAL CLOSURE.**

(a) The SU(3) (p,q) irrep lattice closes as a graph with hexagonal connectivity, Casimir spectrum, and rep-theoretic adjacency. CP² hosts the (p,p)-symmetric sub-lattice with Laplacian eigenvalues p(p+2) and multiplicities (p+1)³. Both are mathematically well-defined.

(b) Neither closes from a packing-axiom analog with the same epistemic standing as Paper 0. Paper 0's axioms are rank-1 specific — "binary distinguishability + 2D isotropy" → SO(3) → S³. They do not lift to rank 2 without becoming representation-theoretic.

(c) The four-layer Coulomb/HO asymmetry (Paper 24 §V) resurfaces at every higher-rank candidate. Calibration π, Fock-projection structure, and κ = −1/16 are S³-specific; they do not transfer to CP² or the (p,q) lattice.

(d) The eightfold-way correspondence (octet at (1,1), decuplet at (3,0)) is tautological — the SU(3) lattice contains the eightfold way *because* it is SU(3) rep theory. Not a prediction.

(e) **Recommendation:** Do not prioritize this as the next major track. Continue with Tracks 1 (multi-focal precision catalogue) and 2 (chemistry solver re-test). Flag the CP² Hopf-direction question and the (p,q) catalogue substrate as opportunistic sub-sprints if SU(3) flavor questions arise downstream.

The "second packing axiom" question (W3 in CLAUDE.md §1.7 multi-focal-wall taxonomy) is **structurally distinct** from "find a packing on a different geometry" and is not addressed by this scoping. That question remains open.
