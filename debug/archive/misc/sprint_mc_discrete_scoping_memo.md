# Quantum Monte Carlo × GeoVac's discrete structure: a scoping pass (2026-09-25)

**Status:** scoping and diagnostics only. No production code, no paper edits. Probe: `debug/mc_discrete_probe.py` → `debug/data/mc_discrete_probe.log` (single-threaded, 34 s, exit 0).
**Seed:** the PI's daydream. Monte Carlo samples along a probability axis give a distribution that "sounds like a wavefunction"; can we "find the closed form of that integration over time"?

---

## 0. The answer in one sentence

In GeoVac's discrete machinery the sign problem enters at **the two-body V_ee matrix elements written in the orbital basis**. It does **not** enter through the one-body graph, which is bipartite and therefore sign-free under the gauge (−1)^(n+m). It does not enter through antisymmetry first either: the He singlet is nodeless in real space, yet its configuration space already carries frustrated cycles. And it does not enter through the qubit encoding, because in a fixed-N sector Jordan–Wigner (JW) is a diagonal gauge. The Gaunt/M_L selection rules do not cure the sign problem. They *contain* it: they split the off-diagonal graph into small symmetry sectors with 26× fewer independent cycles and a 28× smaller stoquastic gap than the same physics in a dense basis. GeoVac's actual signs are no more coherent than random signs placed on the same graph.

---

## 1. What the PM's quick read got right and wrong

| # | PM claim | Verdict |
|---|---|---|
| 1 | VMC samples \|Ψ\|²; the sample histogram is the Born density; for an exact Ψ the local energy collapses to a spike | **Correct.** The zero-variance property holds for any eigenstate, not only the ground state. This is textbook material and not new. |
| 2 | "Integration over time" is imaginary-time projection e^{−Hτ} (DMC / Feynman–Kac) | **Correct, but the closed form is circular.** A closed-form e^{−Hτ} needs the spectrum of H, and for the interacting problem the spectrum is exactly the unknown. The closed form exists only where the spectrum is already known, i.e. for the one-body graph. |
| 3 | The Fock graph has an exactly known integer spectrum, and its heat trace is Paper 18's Mellin object | **Wrong in two ways.** (a) The spectrum −(n²−1) belongs to the *continuum* S³. The lattice Laplacian L = D − A is positive-semidefinite, with spectrum {2−2cos(jπ/(n_max−l)) + 2−2cos(kπ/(2l+1))}. Those eigenvalues are algebraic, not integer (`geovac/lattice_spectrum.py`; `memory/kappa_observation_not_derived.md`). (b) The M1/M2/M3 Mellin engine is built on continuum Dirac/Laplace–Beltrami spectra, not on this lattice heat kernel. The two agree only in the continuum limit. The lattice heat kernel *is* closed form (P2), but as a product of path-graph kernels. |
| 4 | A bare graph Laplacian is stoquastic, so the sign enters somewhere specific | **Half right.** With κ = −1/16, GeoVac's physical one-body operator has **100% positive** off-diagonals: h1_hybrid = diag + κ(−A) = diag + A/16. It is *not* stoquastic as written. It becomes stoquastic under a diagonal gauge only because the graph is bipartite (P1a). The claim that the sign enters at antisymmetry is **wrong for He**: the He singlet has no fermion nodes at all, yet its orbital configuration space has Δ > 0. The sign enters at the basis representation of V_ee. |

---

## 2. Framework: the sign problem as Z₂ holonomy on cycles

Take a real symmetric H. Its off-diagonal graph G has an edge (i,j) wherever H_ij ≠ 0, and each edge carries the sign σ_ij = sgn(H_ij). A diagonal ±1 gauge S makes every off-diagonal ≤ 0 **iff** every cycle C satisfies ∏_{(ij)∈C}(−σ_ij) = +1. In other words the signed graph must be *balanced* in Harary's sense: the Z₂ "flux" through every cycle is pure gauge.

- Balance can be decided in linear time. BFS a spanning tree, fix the gauge along it, then check the β = E − V + C non-tree edges.
- The number of edges that must be removed to restore balance is the **frustration index**. Computing it is NP-hard, and that is the discrete face of Troyer–Wiese.
- The gauge-invariant severity measure used here is the **stoquastic gap** Δ = E₀(H) − E₀(H_st), where H_st = diag(H) − |offdiag(H)|. In a projector or FCIQMC walk without annihilation, the average sign decays as e^{−Δτ}. Spencer–Blunt–Foulkes (2012) tie FCIQMC's annihilation plateau to the gap between H and this sign-stripped Hamiltonian.

**None of this framework is new** (see §5). Its value for GeoVac is that "frustrated cycles on a graph" is the corpus's native vocabulary. Paper 29 counts cycles on the Fock graph with Ihara zeta functions, and Paper 30 computes gauge holonomies (Wilson loops) on the same cycles. In this language the sign problem is a **Z₂ Wilson loop on the determinant graph**, and a sign-curable Hamiltonian is one whose Z₂ connection is flat with trivial holonomy.

---

## 3. Probe results (P1: frustration map; P2: heat kernel)

### P1a. The one-electron Fock graph

| n_max | nodes | edges | β | raw positive off-diag | bipartite under (n+m) parity | frustrated edges after BFS gauge | Δ |
|---|---|---|---|---|---|---|---|
| 3 | 14 | 13 | 2 | **1.000** | yes | 0 | 1e-16 |
| 5 | 55 | 70 | 20 | **1.000** | yes | 0 | −2e-16 |

Every edge flips the parity of n+m: radial edges change n by 1, angular edges change m by 1. So every cycle has even length and the all-positive signs are pure gauge.

- Without the gauge, e^{−τ h1} has 11% negative entries.
- With the gauge, its minimum entry is 0.

(The GS-sign-gauge column in the log reports frustration here. That is an artefact: the graph is disconnected across l, and the ground state is zero on the l > 0 blocks.)

### P1b. He 2e singlet (hybrid h1 + exact Slater V_ee, `LatticeIndex slater_full`)

All numbers refer to the connected sector containing the ground state. E₀ is identical across bases, which confirms each rotation is a pure basis change.

| n_max | basis | sector dim | edges | β | Δ (Ha) | Δ with random signs on the same graph |
|---|---|---|---|---|---|---|
| 2 | native | 9 | 16 | 8 | 3.5e-3 | 3.8e-3 |
| 2 | exact-h1 (V_ee is the only off-diagonal) | 5 | 10 | 6 | 9.1e-4 | 1.5e-3 |
| 2 | radial-rotated (Gaunt kept) | 9 | 16 | 8 | 1.4e-1 | 5.7e-2 |
| 2 | fully rotated (Gaunt lost) | 15 | 105 | 91 | 2.2e-1 | 1.9e-1 |
| 2 | determinant space (FCIQMC) | 13 | 35 | 23 | 1.7e-2 | 1.1e-2 |
| 3 | native | 57 | 337 | 281 | 4.2e-2 | 2.0e-2 |
| 3 | exact-h1 | 19 | 153 | 135 | 2.3e-2 | 7.5e-3 |
| 3 | radial-rotated | 57 | 337 | 281 | 4.8e-2 | 3.6e-2 |
| 3 | fully rotated | 105 | 5460 | 5356 | 9.7e-1 | 1.5 |
| 3 | determinant space | 100 | 1018 | 919 | 5.7e-2 | 4.6e-2 |
| 4 | native | 241 | 4298 | 4058 | 5.8e-2 | 4.7e-2 |
| 4 | exact-h1 | 57 | 1418 | 1362 | 3.4e-2 | 1.9e-2 |
| 4 | radial-rotated | 241 | 4298 | 4058 | 4.7e-2 | 7.0e-2 |
| 4 | fully rotated | 465 | 107,880 | 107,416 | **1.61** | 2.46 |
| 4 | determinant space | 452 | 14,944 | 14,493 | 7.3e-2 | 5.2e-2 |

What the table shows:

1. **The first frustrated cycle appears at n_max = 2 in the V_ee-only (exact-h1) singlet.** The GS sector there has 5 configurations and 10 edges, frustration index ≤ 2, and Δ = 0.9 mHa. The one-body part is balanced, so the frustration must come from V_ee's two-body couplings.
2. **The He singlet is nodeless in real space.** Two electrons in a spatially symmetric state have a positive ground state (Perron–Frobenius on the continuum). So every mHa of Δ above is *basis-induced*: it is the orbital representation of V_ee, not fermionic antisymmetry.
3. **Determinant space costs more than configuration space.** At n_max = 4, Δ is 0.073 in determinant space against 0.058 in the spin-adapted singlet space. The extra is spin-coupling signs, and the known cure is spin adaptation (GUGA-FCIQMC; Dobrautz–Smart–Alavi 2019).
4. **Selection rules contain the sign problem; they do not cure it.** At equal dimension, removing the Gaunt/M_L structure raises the cycle rank 26× (107,416 vs 4,058) and Δ 28× (1.61 vs 0.058 Ha). A radial rotation that keeps the Gaunt structure changes Δ by only 0.8–1.2× at n_max = 3–4. It gives 40× at n_max = 2, where mixing 1s with 2s is simply a bad basis. So the protection is **angular** (sector fragmentation plus Gaunt zeros). The hydrogenic radial basis is not special.
5. **GeoVac's signs are not coherent.** The native Δ sits at or above the random-sign control on the same graph (n_max = 3: 0.042 vs 0.020; n_max = 4: 0.058 vs 0.047; 5 seeds). The benefit is entirely topological (fewer cycles), with no favourable sign correlation.
   - Curve-fit audit: this is 3 sizes × 5 seeds. It is not a law, and no scaling claim is made.
6. **Side fact.** The hybrid h1's m±1 graph edges do not conserve M_L (40 such entries at n_max = 4). They are, however, energetically inert for the ¹S ground state: E₀ is identical to 1e-6 with all M_L, restricted to M_L = 0 (as `build_graph_native_fci` enumerates), or with the m-hops deleted. They enlarge the connected sector (57 → 241 at n_max = 4) and add cycles without changing the energy. This is not a defect in the headline, just a place where the graph adds frustration it does not need.

### P1c. Composed LiH (Q = 30, `build_composed_hamiltonian(lih_spec())`, JW operator restricted to N = 4)

- The N = 4 sector has 27,405 determinants and 118,047 nonzeros, splitting into **9,190 connected components**.
- The physical sector (core 2e, bond 2e, Sz = 0) has 3,400 determinants in 924 components.
- **The ground state lives in a 14-determinant component:** 49 edges, β = 36, frustration index ≤ 10, **Δ = 4.3 mHa** (random-sign control 11–13 mHa).
- The spin-ordering convention was checked: under the interleaved convention, block occupations are constant on every component.
- Energies are reported as the builder's sector eigenvalues (E₀ = −14.144, which includes PK). They are not interpreted.

The composed architecture's block-diagonal ERIs conserve particle number and M_L block by block, and that shrinks the ground state's determinant graph to 14 nodes. **The sign problem is present but irrelevant:** a 14×14 block is diagonalised exactly, so QMC has nothing to do here.

### P1d. The qubit encoding

In a fixed-N sector, the JW computational basis *is* the determinant basis, and changing the orbital ordering multiplies basis states by ±1, which is a diagonal gauge. Frustration and Δ are gauge-invariant, so **JW adds no sign problem** (P1c was run on the JW operator itself). Two caveats:

- Bravyi et al.'s *term-wise* stoquasticity of the Pauli sum is a stronger, basis-local notion.
- Z₂ tapering is a Clifford conjugation, i.e. a genuine basis change, and could move Δ either way. It was not tested.

### P2. The Fock-graph heat kernel is closed form (n_max = 5, 55 nodes)

| t | max\|expm(−tL) − closed form\| | deviation of row sums from 1 | minimum entry |
|---|---|---|---|
| 0.1 | 6.9e-16 | 2e-16 | 0 |
| 1 | 3.9e-16 | 4e-16 | 0 |
| 10 | 4.4e-16 | 4e-15 | 0 |

The closed form is K_t = ⊕_l [K^path_{n_max−l}(t) ⊗ K^path_{2l+1}(t)], with path kernels V diag(e^{−tλ_j}) Vᵀ and λ_j = 2−2cos(jπ/m). This is an exact continuous-time Markov kernel.

**What it buys a projector method:** an exact, Trotter-free propagator for the *graph-kinetic* piece, usable in an interaction-picture or continuous-time expansion (the CT-INT / worm-type idea). **What it does not buy:**

- It does not commute with the diagonal potential. Both diag(−Z²/2n²) and −Z/n² vary along the n-edges, so the potential still needs the interaction picture.
- GeoVac's accurate one-body operator is the *exact* hydrogenic h1, which is **diagonal**. Its propagator is trivial, e^{−τε_n}, and needs no random walk at all.

So the closed form is real, and FORCED, but its computational payoff is near zero.

**Transcendental tag.** The eigenvalues λ are algebraic, lying in ℚ(cos π/m). For rational τ ≠ 0 and algebraic λ ≠ 0, e^{−τλ} is transcendental by Lindemann–Weierstrass. This is Euclidean-time evolution: Paper 18's *flow* class, Paper 34's temporal (imaginary-time) projection. π enters only on temporal compactification to a finite β (Matsubara; Paper 35, WH7). The lattice Mellin transform Σλ^{−s} is a finite sum of algebraic numbers and contains no π, consistent with the π-free graph principle.

---

## 4. Verdicts: FORCED / FREE / WALL

| Item | Verdict | Mechanism |
|---|---|---|
| One-body Fock graph is sign-free | **FORCED** | Bipartite by (n+m) parity, so the gauge (−1)^(n+m) makes it stoquastic. Δ = 0 bit-level. |
| Sign of κ = −1/16 (all off-diagonals positive) | **FREE → gauged away** | A convention choice. Only bipartiteness makes it harmless. |
| Closed-form lattice heat kernel / random walk | **FORCED** | Path-product kernel, exact to 7e-16. Algebraic spectrum, transcendental only through e^{−τλ} (flow). |
| Identifying the lattice heat trace with the Paper 18 Mellin engine | **Not supported** | The continuum-limit relation only. |
| V_ee in the orbital basis produces frustrated cycles | **WALL for diagonal gauge; FREE w.r.t. basis** | Δ > 0 is gauge-invariant, so no ±1 gauge cures it. It appears at n_max = 2 even for the nodeless He singlet, so a basis change (real space) removes it for 2e. |
| Fermionic antisymmetry (N ≥ 3, same-spin pairs) | **WALL (known)** | NP-hard in general (Troyer–Wiese). GeoVac structure offers no cure. Not isolated separately in this probe. |
| Spin coupling (determinants vs CSFs) | **FREE** | Removable by spin adaptation (known). |
| JW encoding in a fixed-N sector | **FORCED neutral** | A diagonal gauge, so no new frustration. Tapering is untested. |
| Gaunt/M_L selection rules | **FORCED containment, not cure** | 26× lower cycle rank and 28× lower Δ vs a dense basis. Signs no better than random. |
| Composed block factorization | **FORCED containment** | The LiH ground state lives in a 14-determinant component, which makes QMC moot there. |

---

## 5. Prior-art table (verified against primary records 2026-09-25 unless marked)

| Result | Source | Status |
|---|---|---|
| Imaginary-time path integral / diffusion ↔ Schrödinger | Kac, Trans. AMS 65, 1 (1949) (Feynman–Kac) | standard; not web-checked this pass |
| Random-walk (DMC) solution of the Schrödinger equation | Anderson, J. Chem. Phys. 63, 1499 (1975) | verified |
| Fixed-node/release-node QMC for fermions | Ceperley & Alder, PRL 45, 566 (1980) | verified |
| Walkers with signs in determinant space; annihilation | Booth, Thom & Alavi, JCP 131, 054106 (2009) | verified |
| Sign problem ↔ population dynamics / plateau in FCIQMC | Spencer, Blunt & Foulkes, JCP 136, 054110 (2012) | verified |
| Lattice-regularized DMC | Casula, Filippi & Sorella, PRL 95, 100201 (2005) | verified |
| Stoquastic definition; complexity class of stoquastic LH | Bravyi, DiVincenzo, Oliveira & Terhal, QIC 8, 361 (2008) | verified |
| Sign problem is NP-hard | Troyer & Wiese, PRL 94, 170201 (2005) | verified |
| Curing by local basis change: decision problem | Klassen & Terhal, Quantum 3, 139 (2019); Marvian, Lidar & Hen, Nat. Commun. 10 (2019) (NP-complete for single-qubit Clifford/orthogonal cures) | verified |
| Easing the sign problem by optimized basis changes; non-stoquasticity measures | Hangleiter, Roth, Nagaj & Eisert, Sci. Adv. 6, eabb8341 (2020) | verified |
| Balanced signed graphs; cycle criterion; frustration index | Harary (1953) and the signed-graph literature (e.g. Zaslavsky, arXiv:1701.07963) | criterion verified in secondary; Harary primary not fetched |
| Spin-adapted (GUGA) FCIQMC | Dobrautz, Smart & Alavi, JCP 151, 094104 (2019) | verified earlier in corpus (`debug/lit_search_algebraic_methods_compatibility_2026_06_05.md`) |

**What is not new:** everything in §2, the stoquastic gap as a severity measure, symmetry blocking as a sign-problem reducer, the nodelessness of the 2e singlet, and the interaction-picture use of an exact reference propagator.

**What is ours:** only the *measurements* on GeoVac's Hamiltonians. These are the bipartite (n+m) gauge of the Fock graph, the level-by-level location table, the 26×/28× containment ratio, the "signs no better than random" control, and the 14-determinant composed-LiH ground-state component. These are corpus facts, not contributions to QMC.

---

## 6. Ranked leads

**L1: fixed-node DMC on the existing LiH VMC-over-FCI trial function.** Highest value; not GeoVac-native.
- The corpus frontier is the LiH energy: −8.047 variational, with −8.062 needing 3-body correlation.
- DMC on Ψ_CI × Jastrow (from `debug/lih_vmc.py`) projects out *all* non-nodal error, which covers the missing 3-body and cusp correlation, leaving only the fixed-node error of Ψ_CI's nodes. This is the "integration over time" of the seed, applied where it pays.
- Standard algorithm: Umrigar–Nightingale–Runge drift-diffusion-branching (this citation was not verified in this pass).
- **Falsifier:** FN-DMC fails to land more than 5 mHa below the VMC −8.047, which would mean the Ψ_CI nodes are poor.
- **Cheapest next step:** add a DMC branch to a *copy* of the VMC driver, with small τ and 2–3 time steps for τ → 0 extrapolation. This must run after the current VMC jobs finish, because it must not touch the running files.
- **Novelty:** none. It is a standard method applied to a GeoVac trial function. Its value is accuracy.

**L2: minimise the stoquastic gap Δ within the Gaunt-preserving basis class.** GeoVac-flavoured; idea known.
- Optimise per-l radial orthogonal rotations, which keep every selection rule (the §4 prime directive), to minimise Δ.
- A random radial rotation already lowered Δ at n_max = 4 (0.058 → 0.047).
- **Falsifier:** the optimised Δ fails to fall below 50% of native at n_max = 3–4.
- **Cheapest next step:** scipy minimisation of Δ over exp(antisymmetric) per-l generators for He n_max = 3 (57-dim sector; minutes).
- **Novelty:** the method is Hangleiter et al. (2020), and natural orbitals reduce FCIQMC plateaus, which is known. The only new element is the selection-rule-preserving constraint. It matters only if GeoVac Hamiltonians ever reach FCIQMC-scale spaces, which composed ones do not (L3).

**L3: library-wide ground-state component size.** A positioning check.
- **Claim:** across the 37-system library the ground state's connected determinant component stays small enough for exact diagonalisation, making QMC unnecessary for composed Hamiltonians.
- **Falsifier:** some library molecule has a ground-state component above ~10⁵ determinants at max_n = 2.
- **Cheapest next step:** loop the P1c sector analysis over `hamiltonian()` names.
- **Novelty:** a direct corollary of block-diagonal ERIs (Paper 14/22). Worth one line in a paper at most.

**Killed in Phase B:**
- "The lattice heat kernel is the Paper 18 Mellin object." It is continuum-only.
- "Gaunt sparsity gives sign-coherent Hamiltonians." The random-sign control refutes it.
- "The graph heat kernel as a useful reference propagator." Production h1 is diagonal, so the payoff is ~0.
- "The sign problem enters at antisymmetry." False for He; the orbital basis enters first.

---

## 7. Honest overall verdict

This is **known physics seen from a new angle**. There is no GeoVac-native contribution to QMC here.

The angle is still worth having. Phrased as Z₂ holonomy on the off-diagonal graph, the sign problem is a cycle question, and it plugs straight into the corpus's cycle and holonomy machinery (Papers 29/30). The measured map is clean: FORCED sign-free for the one-body graph, a WALL (gauge-invariant Δ) at V_ee in the orbital basis, FORCED containment by the selection rules, and neutrality of the encoding.

The one lead with real payoff is L1, standard fixed-node DMC on the LiH trial function the corpus already has. That is where "integrate over (imaginary) time" meets the actual frontier.
