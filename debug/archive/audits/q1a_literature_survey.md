# Q1-A: Literature Survey — QED Photon Field on S³

**Track:** Q1-A | **Date:** April 15, 2026 | **Status:** Complete

---

## 1. What's Known

### Lattice gauge theory on S³

Standard lattice gauge theory (Wilson 1974) places matter fields on lattice sites and gauge connections on links (edges). The formulation is inherently flat-space: the lattice is a regular hypercubic grid in R⁴ (Euclidean). Luscher and collaborators developed lattice gauge theory on *finite-volume* manifolds (including spheres) primarily as an infrared regulator. Luscher (1982, "Topology of Lattice Gauge Fields," Comm. Math. Phys. 85) showed that the topological charge is well-defined on a discretized S⁴, but this was a compactification tool, not a claim that S³ is the photon's natural home. Berg and Luscher (1981) studied lattice gauge theory on S⁴ for instanton physics. **No paper was found formulating dynamical lattice QED on S³ as a natural geometry in the Fock sense.** Lattice QCD simulations sometimes use S³ × R as a spatial geometry to study finite-volume effects and confinement (e.g., Aharony et al. 2004, "Hagedorn/Deconfinement Phase Transition in Weakly Coupled Large N Gauge Theories"), but this is a computational choice, not a natural geometry argument.

### Chern-Simons theory on S³

Witten (1989, "Quantum Field Theory and the Jones Polynomial," Comm. Math. Phys. 121) formulated Chern-Simons gauge theory on S³ and showed it computes the Jones polynomial of knots. This is a **topological** field theory: the action S = (k/4π) ∫ Tr(A ∧ dA + (2/3)A ∧ A ∧ A) depends only on the topology of S³, not on any metric. The partition function Z(S³) is exactly computable and gives knot invariants. Key distinction: Chern-Simons is a 3D theory with no local propagating degrees of freedom (no photons). QED is a 4D theory with local dynamics. The mathematical frameworks share the language of connections on S³ (gauge fields as 1-forms), but the physics is categorically different. Chern-Simons on S³ does NOT give dynamical photon propagation.

### Conformal invariance of Maxwell in 4D

The source-free Maxwell action ∫ F_μν F^μν √g d⁴x is conformally invariant in exactly 4 spacetime dimensions (and only in 4D). This is a classical result going back to Cunningham (1909) and Bateman (1910). Under conformal compactification, Minkowski R^{3,1} maps to (a portion of) S³ × R (where R is conformal time). This is the **conformal boundary of AdS₅** — the same S³ that appears in AdS/CFT.

Specifically: the conformal group of Minkowski space is SO(4,2), which is also the isometry group of AdS₅. The conformal boundary of AdS₅ is S³ × R. Under this compactification, the spatial slices of the conformal boundary are S³. Maxwell's equations on S³ × R are the conformally mapped Maxwell equations from flat Minkowski space. This is well-established in the AdS/CFT literature (Maldacena 1997, Witten 1998).

Eastwood and Singer (1985, "A Conformally Invariant Maxwell Gauge," Phys. Lett. A 107) showed that the Maxwell field on conformally compactified spacetime has a well-defined formulation on S⁴ (Euclidean) or S³ × S¹ (thermal). The photon propagator on S³ × R can be expanded in vector spherical harmonics on S³. This IS a natural "photon on S³" formulation — but it's the BOUNDARY of the story, not a Fock-type momentum-space projection.

### Fock projection: why it doesn't work for photons

Fock's 1935 projection works because:
1. The hydrogen Schrodinger equation in momentum space has the form (p² + p₀²)ψ = V*ψ
2. The constraint p₀² = -2E (energy shell) maps {p₀, p₁, p₂, p₃} to a 3-sphere
3. The Coulomb potential 1/r in position space becomes a convolution kernel on S³
4. The resulting equation is the eigenvalue problem for the Laplace-Beltrami operator on S³

For photons (massless spin-1):
- The energy-momentum relation is E = |p|c (a light cone, not a sphere)
- There is no p₀ parameter to define a 3-sphere in momentum space
- The photon has no bound states (no discrete spectrum to index nodes)
- The photon is a gauge field (1-form), not a scalar field (0-form)

**There is no direct "Fock projection for photons."** The obstruction is structural, not technical: the Fock projection uses the massive particle's energy-shell (a sphere) as the compactification target. Massless particles have a cone (null surface), which does not compactify to a sphere.

### Hodge Laplacian on S³

The spectrum of the Hodge Laplacian Δ_p on p-forms on S³ is exactly known (Ikeda and Taniguchi 1978, "Spectra and Eigenforms of the Laplacian on S^n and P^n(C)"):

- 0-forms (scalars): eigenvalues n(n+2), degeneracy (n+1)², n = 0, 1, 2, ...
- 1-forms (vectors/gauge fields): eigenvalues n(n+2), degeneracy 2n(n+2), n = 1, 2, 3, ...
- 2-forms: eigenvalues n(n+2), degeneracy 2n(n+2), n = 1, 2, 3, ...
- 3-forms (pseudo-scalars): eigenvalues n(n+2), degeneracy (n+1)², n = 0, 1, 2, ...

Notable: The eigenvalues n(n+2) are the SAME for all form degrees — this is special to S³ (a consequence of the round metric and maximal symmetry). The degeneracies differ. For 1-forms, the degeneracy 2n(n+2) decomposes into exact (d-exact) and coexact (δ-exact) pieces by the Hodge decomposition. On S³ (which has β₁ = 0), there are no harmonic 1-forms, so every eigenspace splits cleanly into exact and coexact halves.

The coexact 1-forms on S³ are the transverse vector harmonics — these are precisely the modes that a photon field on S³ would occupy (the "Coulomb gauge" modes). Their degeneracy is n(n+2), matching the scalar degeneracy at level n+1 (shifted by one unit).

### Graph Hodge theory

Lim (2020, "Hodge Laplacians on Graphs," SIAM Review 62) and Schaub et al. (2020, "Random walks on simplicial complexes and the normalized Hodge 1-Laplacian," SIAM Review 62) developed the combinatorial Hodge Laplacian on simplicial complexes. The edge Laplacian L₁ = B₁ᵀ B₁ + B₂ B₂ᵀ (where B₁ is the node-edge incidence matrix and B₂ is the edge-triangle incidence matrix) is the discrete analog of Δ₁. For a graph without triangles (or ignoring higher simplices), L₁ = B₁ᵀ B₁. This has been applied to signal processing, network analysis, and data science. **No paper was found applying graph Hodge theory to lattice gauge theory or to quantum field theory on graphs.** This appears to be a genuine gap in the literature.

---

## 2. What's NOT Known / The Gap

1. **No one has formulated a "Fock projection for photons."** The massless dispersion relation (light cone vs energy shell) is the structural obstruction. The conformal invariance of Maxwell provides an S³ home for the photon via conformal compactification of spacetime, but this is a different mechanism than Fock's momentum-space projection.

2. **No one has connected graph Hodge theory to lattice gauge theory.** Lattice gauge theorists use Wilson loops and link variables; graph theorists use incidence matrices and edge Laplacians. These are the same mathematical objects (gauge field on edges = 1-cochain on graph) but the communities have not cross-pollinated on this specific point.

3. **The relationship between the S³ Hodge-1 spectrum and the GeoVac graph is unexplored.** The continuous Hodge Laplacian on 1-forms on S³ has known spectrum. The GeoVac graph is a discretization of S³. The edge Laplacian of the GeoVac graph is computable. Whether these converge (as the graph refines) is a well-posed mathematical question that no one has asked.

4. **The Hopf fibration's role in QED on S³ is suggestive but undeveloped.** The U(1) fiber of the Hopf bundle S³ → S² is the gauge group of QED. In lattice gauge theory, gauge transformations live on vertices and gauge fields live on edges. In the Hopf bundle, the U(1) acts on each S¹ fiber. The connection between these two pictures has not been made explicit in the literature for the specific case of the GeoVac S³ graph.

---

## 3. Key Question: Is There a "Fock Projection for Photons"?

**No, not in the direct (momentum-space stereographic) sense.** The obstruction is the massless dispersion relation: E = |p|c gives a cone, not a sphere. There is no energy-shell parameter p₀ to define a stereographic map p-space → S³.

**However, there IS a natural "photon on S³" via a different route:** conformal compactification of spacetime. Maxwell's equations are conformally invariant in 4D. Under conformal compactification, the spatial slices of the conformal boundary are S³. The photon modes on S³ are the coexact vector spherical harmonics (transverse 1-forms), with spectrum n(n+2) and degeneracy n(n+2).

The distinction matters:
- **Fock projection (matter):** momentum space → S³. The S³ is the space of quantum states. Nodes = electron states.
- **Conformal compactification (gauge):** spacetime → S³ × R. The S³ is the spatial manifold. 1-forms on S³ = photon modes.

These are **different S³s playing different roles.** The electron's S³ is in momentum space; the photon's S³ is in configuration space (after conformal compactification). Whether these can be unified into a single mathematical object is an open question — but it is NOT the naive expectation that both fields live on the same S³ in the same way.

---

## 4. Assessment: Is QED-on-S³ a Well-Posed Question?

**Yes, but it requires careful formulation.** The question is NOT "does the photon have a Fock projection to S³?" (answer: no, structural obstruction from masslessness). The well-posed questions are:

1. **(Graph Hodge theory question)** Does the edge Laplacian L₁ of the GeoVac S³ graph encode photon-like degrees of freedom? This is purely combinatorial and answerable by computation (Track Q1-B). The edge Laplacian is the discrete Hodge-1 Laplacian. If its spectrum converges to the continuous Hodge-1 spectrum on S³ as n_max → ∞, then the GeoVac graph naturally discretizes both matter (nodes = 0-forms) and gauge fields (edges = 1-forms).

2. **(Conformal field theory question)** In the AdS/CFT picture, the hydrogen atom's SO(4,2) symmetry maps to the conformal boundary S³ × R. The photon on S³ × R is a well-defined CFT₄ operator. Does the GeoVac graph encode this boundary CFT structure? This connects to Paper 3's holographic framework.

3. **(Hopf gauge structure question)** The Hopf fibration S³ → S² has U(1) fiber = QED gauge group. In the GeoVac graph, the Hopf fiber structure is already computed (Paper 2, Phase 4B). Do the U(1) holonomies along the Hopf fibers reproduce lattice Wilson loops? This is a concrete, computable question.

**The category error to avoid:** treating the photon as if it should have a Fock projection (it shouldn't — it's massless). The correct framing is that matter and gauge fields discretize differently on the same graph: matter on nodes (0-forms, Laplace-Beltrami), gauge fields on edges (1-forms, Hodge Laplacian). This is exactly the lattice gauge theory paradigm (Wilson 1974), applied to the GeoVac S³ graph rather than a flat hypercubic lattice.

**Bottom line:** The photon does NOT have a "Fock projection" analogous to the electron's. But it may have a natural S³ discretization through the edge Laplacian / Hodge-1 structure of the GeoVac graph. Track Q1-B (edge Laplacian computation) is the decisive test.

---

## References

- Fock, V. (1935). Z. Phys. 98, 145. [S³ stereographic projection for hydrogen]
- Wilson, K. (1974). Phys. Rev. D 10, 2445. [Lattice gauge theory]
- Witten, E. (1989). Comm. Math. Phys. 121, 351. [Chern-Simons on S³]
- Luscher, M. (1982). Comm. Math. Phys. 85, 39. [Topology of lattice gauge fields]
- Cunningham, E. (1909). Proc. London Math. Soc. 8, 77. [Conformal invariance of Maxwell]
- Bateman, H. (1910). Proc. London Math. Soc. 8, 223. [Conformal invariance of Maxwell]
- Eastwood, M. & Singer, M. (1985). Phys. Lett. A 107, 73. [Conformally invariant Maxwell gauge]
- Ikeda, A. & Taniguchi, Y. (1978). Osaka J. Math. 15, 515. [Spectra of Laplacians on spheres]
- Lim, L.-H. (2020). SIAM Review 62, 685. [Hodge Laplacians on graphs]
- Schaub, M. et al. (2020). SIAM Review 62, 353. [Random walks on simplicial complexes]
- Maldacena, J. (1998). Adv. Theor. Math. Phys. 2, 231. [AdS/CFT]
- Aharony, O. et al. (2004). Adv. Theor. Math. Phys. 8, 603. [Gauge theory on S³]
- Berg, B. & Luscher, M. (1981). Nucl. Phys. B 190, 412. [Lattice gauge fields on S⁴]
