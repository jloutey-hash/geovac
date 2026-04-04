# Track AH: Graph-Native Inter-Group Antisymmetry at Level 5

## Investigation Summary

**Question:** Can the Pauli exclusion principle, which is exact and free at Levels 1-4, be extended to inter-group antisymmetry at Level 5 without the Phillips-Kleinman pseudopotential?

**Answer: No.** All three avenues produce structural negative results. Inter-group antisymmetry is intrinsically non-topological in composed geometry because composition breaks the coordinate unification that makes antisymmetry topological at lower levels.

---

## Background: How Antisymmetry Works at Levels 1-4

At each level, antisymmetry is enforced by a different mechanism, but they share a common feature: **all electrons live in a single coordinate system**.

**Level 1 (S3, one electron):** No antisymmetry needed. Single electron.

**Level 3 (Hyperspherical, two electrons):** The two electrons share coordinates (R_e, alpha, theta_12). Exchange symmetry is the map alpha -> pi/2 - alpha, which is a COORDINATE SYMMETRY of the angular Hamiltonian (Paper 13, Sec VII). The singlet constraint selects even grand angular momentum nu = 0, 2, 4, ... This is topological: the symmetry is a property of the coordinate system itself, not an added potential. The angular basis functions sin(2n*alpha) automatically satisfy the constraint when restricted to odd n.

**Level 4 (Mol-frame hyperspherical, two electrons):** Same as Level 3 but in molecule-frame coordinates. The gerade constraint l_1 + l_2 = even (Paper 15, Eq. 14) is the analog of the Level 3 exchange symmetry. It restricts which channels are allowed. This is again topological: the constraint is on the quantum number labels, not an external potential.

**The common feature:** At Levels 3 and 4, both electrons occupy the SAME Hilbert space with a SHARED coordinate system. The antisymmetry is a SYMMETRY OF THE HAMILTONIAN (alpha -> pi/2 - alpha commutes with H_ang). This allows it to be enforced by restricting to a symmetry sector --- no potential needed.

**Level 5 (Composed geometry):** The core electrons live in Level 3 hyperspherical coordinates (R_e^core, alpha_core). The valence electrons live in Level 4 molecule-frame hyperspherical coordinates (R_e^val, alpha, theta_1, theta_2, Phi; rho). These are DIFFERENT coordinate systems with DIFFERENT dimensionalities. There is no shared coordinate in which to define the exchange map "swap a core electron with a valence electron."

---

## Avenue 1: Node Exclusion on the Composed Graph

### Idea
The composed graph G_total = G_core x G_val is a product of two graphs. PK says "valence electrons cannot occupy the same spatial region as core electrons." If we could identify which valence graph nodes correspond to the core region, we could DELETE them, making exclusion topological (a property of the graph, not an added potential).

### Analysis

The core lives on a Level 3 graph with nodes (R_i^core, nu_core). The valence lives on a Level 4 graph with nodes (R_j^val, alpha_k, channel). The core is parameterized by a hyperradius R_e^core = sqrt(r_1^2 + r_2^2) and hyperangle alpha_core = arctan(r_2/r_1), where r_1, r_2 are distances from the Li nucleus.

The valence electrons are parameterized by R_e^val = sqrt(r_3^2 + r_4^2) where r_3, r_4 are distances from the molecular center (or from the nuclei, depending on origin convention).

**Obstruction 1: No coordinate correspondence.** The core coordinate r_1 (distance of core electron 1 from Li nucleus) and the valence coordinate r_3 (distance of valence electron 3 from molecular center) are not the same. Even if r_3 happens to be small (valence electron near Li), r_3 is measured from a different origin and parameterized differently. There is no 1-to-1 map between "core node i" and "valence node j" that would identify which valence nodes to delete.

**Obstruction 2: The core is not a spatial region in valence coordinates.** The PK potential V_PK(r) = A*exp(-Br^2)/r^2 acts on valence electron coordinates (the r in V_PK is the valence electron's distance from the Li nucleus within the Level 4 framework). The "core region" in valence coordinates is the set of points where r is small --- but this is a CONTINUOUS region in the Level 4 angular/radial coordinate space, not a discrete set of nodes.

At Level 4, the angular coordinate alpha parameterizes the ratio r_3/r_4 between the two valence electrons. The radial coordinate R_e parameterizes their joint distance scale. The "core region" r_3 ~ 0 corresponds to alpha ~ 0 AND R_e ~ r_4, which is a continuous submanifold, not a set of graph nodes.

**Obstruction 3: Node deletion changes the spectrum globally.** Even if we could identify which nodes to delete, removing nodes from a graph changes its Laplacian spectrum in a global, non-local way. The PK potential is LOCALIZED (Gaussian decay), but node deletion would change eigenvalues of ALL channels, including high-l channels that should be unaffected. This contradicts the l-dependent PK requirement (Paper 17, Sec IV): the core is pure l=0, so only l=0 channels should feel the exclusion. Node deletion cannot be made l-selective.

**Obstruction 4: Track AD confirmed PK is essential for equilibrium.** Without PK, the PES is monotonically attractive (Paper 17). This means the exclusion must provide a REPULSIVE contribution that GROWS at short R. Node deletion removes degrees of freedom but does not add repulsive energy. The equilibrium requires an actual repulsive potential energy contribution from the orthogonalization kinetic energy.

### Verdict: NEGATIVE

Node exclusion fails because (a) there is no coordinate correspondence between the two graphs, (b) the core region is not a discrete set of nodes in valence coordinates, (c) node deletion is not l-selective, and (d) node deletion cannot provide the repulsive energy needed for equilibrium.

---

## Avenue 2: Fiber Bundle Connection as Antisymmetry

### Idea
The composed geometry is a fiber bundle: G_total = G_nuc ⋉ G_core(R) ⋉ G_val(R, core_state). The Berry connection P_nu_mu encodes how channels change with the hyperradius. Could inter-group antisymmetry be encoded as a gauge condition on this connection?

### Analysis

In single-group problems (Levels 3-4), the Berry connection P_nu_mu(R) = <Phi_nu | d/dR | Phi_mu> is a property of the INTRA-group angular eigenstates. It couples channels within the same electron group. It has nothing to do with antisymmetry --- it encodes the non-adiabatic dynamics when the slow coordinate R changes.

At Level 5, there are TWO fiber bundles stacked:
1. Core bundle: base = R, fiber = Level 3 angular channels
2. Valence bundle: base = R, fiber = Level 4 angular channels (parameterized by core state)

The inter-group coupling occurs through:
- Z_eff(r): core density screens nuclear charge for valence
- PK: core-valence orthogonality enforcement
- Nuclear repulsion: trivial

**Obstruction 1: Gauge transformations preserve physics, not enforce it.** A gauge condition on the connection is a CHOICE OF BASIS for the angular channels, not a physical constraint. Changing the gauge (e.g., from adiabatic to diabatic channels) changes the connection P but not the physical observables. Antisymmetry is a physical constraint (the wavefunction changes sign under electron exchange), not a basis choice.

**Obstruction 2: The connection is within groups, not between them.** The Berry connection P_nu_mu connects channel nu to channel mu within the SAME electron group. There is no "inter-group Berry connection" because the core and valence electrons have different fast coordinates. The composition is a SEMI-DIRECT product, not a direct product: the valence problem depends on the core state, but not vice versa. This one-way dependence means there is no symmetric exchange structure to encode.

**Obstruction 3: Antisymmetry requires a swap operation.** Antisymmetry under exchange of electron i (core) with electron j (valence) requires an operation that maps:

  (core state with electron i) x (valence state with electron j)
  -> (core state with electron j) x (valence state with electron i)

This operation MIXES the fiber spaces: it takes a vector in Fiber_core x Fiber_val and produces a vector in the SAME product space, but with the electron labels exchanged. This is not a connection (which maps fibers at nearby base points); it is a GLOBAL operation on the total Hilbert space.

In the composed geometry, this swap is undefined because "core state with electron j" requires solving the core problem with a different electron, which changes the core coordinates entirely. The Level 3 hyperspherical coordinates (R_e^core, alpha_core) describe TWO core electrons. Replacing one with a valence electron creates a Level 3 system with different Z_eff, different boundary conditions, and a fundamentally different Hamiltonian.

**Obstruction 4: The connection cannot generate repulsion.** Even if a gauge condition could encode some constraint, the Berry connection governs non-adiabatic TRANSITIONS between channels. It does not contribute a REPULSIVE potential energy. The PK pseudopotential's role is to raise the energy when valence electrons penetrate the core --- this is a positive-definite energy contribution. No connection or gauge condition can generate a positive-definite energy contribution from geometric data alone.

### Verdict: NEGATIVE

The fiber bundle connection is an intra-group, basis-dependent quantity that cannot encode inter-group antisymmetry. The swap operation between groups is undefined in composed coordinates, and no connection can generate the repulsive energy required for equilibrium.

---

## Avenue 3: Spectral Exclusion via Graph Laplacian

### Idea
Solve the valence problem in the orthogonal complement of the core eigenstates. If the core occupies eigenstates {phi_1, ..., phi_k} of some operator, constrain the valence solver to the subspace orthogonal to these states. This is spectral exclusion --- projecting out the occupied states.

### Prior art and failure modes
Track Q (v2.0.12) tested spectral PK projectors (rank 1-3) in the Gegenbauer channel basis. All produced monotonically attractive PES --- no equilibrium. The documented conclusion: "Angular-space projectors cannot replace coordinate-space exclusion."

### Analysis

The key question for Avenue 3 is whether GRAPH-LEVEL spectral exclusion (on the discrete nodes of the composed graph) can succeed where ANGULAR-SPACE spectral exclusion (Track Q) failed.

**Obstruction 1: The core and valence graphs have different spectra.** The core graph (Level 3 hyperspherical) has eigenvalues determined by SO(6) Casimir + nuclear/e-e coupling on S^5. The valence graph (Level 4 molecule-frame hyperspherical) has eigenvalues determined by the mol-frame angular Hamiltonian on a different manifold. The core eigenstates live in a DIFFERENT Hilbert space from the valence eigenstates. You cannot project out core states from a space that does not contain them.

To make spectral exclusion work, you would need a SHARED Hilbert space containing both core and valence states. But this is precisely what composed geometry avoids --- the whole point of composition is that each group gets its own natural coordinate system.

**Obstruction 2: Same failure as Track Q, for the same reason.** Track Q's failure was diagnosed as: angular-space projectors cannot replace COORDINATE-SPACE exclusion. The reason is that the core's effect on valence electrons is a SPATIAL phenomenon: the core density occupies a region of physical space near the nucleus, and the valence electrons must be kept out of that region. Spectral exclusion (projecting out eigenstates) operates in FUNCTION SPACE, not coordinate space. A valence eigenstate can be orthogonal to all core eigenstates in function space while still having significant amplitude in the core region of coordinate space.

This is the fundamental distinction:
- PK (coordinate-space exclusion): V_PK(r) is large where the core density is large, regardless of which eigenstate the valence electron is in.
- Spectral exclusion (function-space): projects out specific eigenstates but does not prevent amplitude in any spatial region.

For a concrete example: the 2s orbital of hydrogen is orthogonal to the 1s orbital (spectral exclusion is satisfied), but the 2s orbital has NONZERO amplitude at r=0 (coordinate-space penetration is NOT prevented). The orthogonality comes from the node structure (2s has a node), not from spatial exclusion.

In the PK framework, even a 2s-like valence state pays an energy penalty for having amplitude near the nucleus, because V_PK(r) ~ A*exp(-Br^2)/r^2 is large there. This is the CORRECT physics: the orthogonalization kinetic energy cost of maintaining a node in the core region.

**Obstruction 3: Spectral exclusion in a shared space = full 4-electron problem.** If we wanted to define spectral exclusion properly, we would need to embed both core and valence electrons in a SINGLE coordinate system (e.g., S^11 for 4 electrons at one center) and then solve the 4-electron problem with antisymmetry imposed. But this is exactly what Paper 15 identified as computationally infeasible (~10^6 angular basis functions). The composed geometry was invented precisely to AVOID this.

**Obstruction 4: Graph-level vs. angular-level is a false distinction.** The "graph" in GeoVac IS the discretization of the angular manifold. The graph nodes correspond to angular grid points or spectral basis functions. "Graph-level spectral exclusion" and "angular-space spectral exclusion" are the same thing, just described differently. Track Q's negative result applies directly.

### Verdict: NEGATIVE

Spectral exclusion fails because (a) core and valence states live in different Hilbert spaces with no natural embedding, (b) function-space orthogonality does not prevent coordinate-space penetration (the Track Q obstruction), (c) proper spectral exclusion requires the full N-electron Hilbert space that composition was designed to avoid, and (d) graph-level and angular-level exclusion are the same thing.

---

## Master Obstruction: Why Composition Breaks Antisymmetry

The three avenues all fail for a single underlying reason:

**Antisymmetry is a property of the TOTAL wavefunction, not a property of individual subsystems.**

At Levels 1-4, all electrons share a single coordinate system, so the total wavefunction Psi(r_1, ..., r_N) is defined on a single product space. The exchange operator P_ij that swaps electrons i and j is a well-defined symmetry of the Hamiltonian. Antisymmetry can be enforced by restricting to a symmetry sector (even/odd nu, gerade/ungerade l_1+l_2).

At Level 5, the total wavefunction is FACTORED:

  Psi_total ~ Psi_core(r_1, r_2) * Psi_val(r_3, r_4; core_state)

The exchange operator P_13 (swap core electron 1 with valence electron 3) maps:

  Psi_core(r_1, r_2) * Psi_val(r_3, r_4) -> Psi_core(r_3, r_2) * Psi_val(r_1, r_4)

This is NOT a symmetry of the composed Hamiltonian. The core Hamiltonian H_core acts on (r_1, r_2) in Level 3 coordinates; after the swap, it would need to act on (r_3, r_2), but r_3 is a Level 4 coordinate. The swap mixes the coordinate systems, making the factored form invalid.

The Phillips-Kleinman pseudopotential is the CORRECT resolution: it approximates the effect of inter-group antisymmetry as an effective potential acting within the valence coordinate system. The PK theorem (Phillips & Kleinman, 1959) proves that for any orthogonality constraint |psi_val> _|_ |psi_core>, there exists an effective potential V_PK such that the constrained variational problem is equivalent to an unconstrained problem with V_PK added. The Gaussian/r^2 form is an approximation to this exact (nonlocal) PK potential.

### Classification of the obstruction

In the Paper 18 exchange constant taxonomy:

- **Intrinsic exchange constants** (kappa = -1/16): arise from the graph itself. Present at all levels.
- **Calibration exchange constants** (e^a E_1(a)): arise from projecting discrete spectra to continuous coordinates. Present at Levels 2+.
- **Embedding exchange constants** (1/r_12 on S^5): arise from embedding interacting systems on compact manifolds. Present at Level 4.
- **Flow exchange constants** (mu(R)): arise from adiabatic separation. Present at Levels 3+.

The PK pseudopotential represents a NEW type of exchange constant, not yet classified in Paper 18: the **composition exchange constant**. It is the irreducible cost of approximating inter-group antisymmetry when the groups live in different coordinate systems. It cannot be absorbed into graph topology, bundle connections, or spectral projections. It is intrinsic to the act of composition itself.

---

## Relationship to the 8 PK-Related Failed Approaches (Section 3)

Every prior PK-related failure in Section 3 attempted to MODIFY the form of PK or REPLACE it with a different potential/projector. All failed because they misidentified the problem:

1. **Channel-blind PK at higher l_max**: Failed because PK derived from l=0 core shouldn't affect l>0 channels. l-dependent PK partially fixes this but doesn't address the adiabatic approximation.

2. **Algebraic PK projector (rank-1)**: Rank-1 too weak; valence avoids core in function space while penetrating in coordinate space. This is exactly the function-space vs. coordinate-space distinction from Avenue 3.

3. **Spectral PK projector (rank-k)**: Same obstruction as rank-1, extended to structural conclusion. Angular-space projectors categorically cannot replace coordinate-space exclusion.

4. **Enhanced Z_eff (PK-free)**: 1/r^2 vs 1/r shape mismatch. Core exclusion cannot be encoded in the nuclear potential.

5. **Self-consistent PK**: Reduces drift ~50% but doesn't eliminate it. Adiabatic approximation is the root cause.

6. **Projected PK**: Destroys variational character.

7. **R-dependent PK**: Correct form but R_ref is empirical (Track AD: no ab initio derivation possible).

8. **PK-free (Track AD)**: No equilibrium exists without PK. PK is essential.

Track AH's investigation confirms: the problem is not with PK's FORM but with the STRUCTURE of composed geometry. PK (or some equivalent effective potential) is a necessary feature of any composed approach. The 2D variational solver (which reduces the adiabatic approximation error) is the correct path for improving results, not eliminating PK.

---

## Conclusions

1. **All three avenues produce structural negative results.** Inter-group antisymmetry cannot be made graph-native at Level 5.

2. **The master obstruction is coordinate incompatibility.** Antisymmetry requires a swap operation between electrons, which is well-defined only when all electrons share a coordinate system. Composition deliberately breaks this shared structure for computational tractability.

3. **PK is the correct architecture.** The Phillips-Kleinman pseudopotential is not a hack or approximation that can be replaced by a cleaner topological construction. It is the natural encoding of inter-group antisymmetry in a factored coordinate system.

4. **PK is a composition exchange constant.** In the Paper 18 taxonomy, PK represents a new class: the irreducible cost of working in composed coordinates rather than the full N-electron Hilbert space.

5. **The improvement path is clear.** The l_max divergence is from the adiabatic approximation (confirmed by Tracks A, AC, AD). The 2D variational solver reduces drift 4x. PK form improvements (l-dependent, R-dependent) provide additional gains. But PK itself cannot be eliminated from composed geometry.

6. **A graph-native alternative would require the full S^(3N-1) manifold.** The only way to make all inter-electron antisymmetry topological is to put all electrons in a single coordinate system --- which is the S^(3N-1) approach of Paper 7 Sec VI. For 4 electrons (LiH), this is S^11, estimated at ~10^6 angular basis functions. This is the cost of making antisymmetry free; PK is the price of avoiding it.
