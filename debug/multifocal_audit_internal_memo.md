# Multi-Focal Composition Internal Audit (Phase A, Track 1)

**Date:** 2026-05-07
**Scope:** Catalogue every multi-focal composition in the GeoVac project record (works and failures), classify failures into walls W1 / W2, and propose sub-walls.
**Method:** Reading-only audit. Sources: CLAUDE.md §1.6–§3, CHANGELOG v2.7+, Papers 17 / 19 / 23 / 32 §VIII.B–§VIII.C, Sprint H1 / LS-8a / HF-1..5 memos, multi-focal-wall + structural-skeleton memory files.
**Author:** PM (audit fork; no production files modified).

The "focal length" $p_0 = \sqrt{-2E}$ is the Fock projection scale that maps a discrete sector onto $S^3$. Distinct binding energies $\Rightarrow$ distinct focal lengths. A "multi-focal composition" is any construction that combines two or more sectors at different focal lengths into a single physical observable. Discrete-label tensor products (e.g. proton spin $\otimes$ electron spin) are *not* multi-focal in this sense — they live on combinatorial structure that exists prior to any Fock projection, on a single common graph. The audit treats those as the working baseline and looks at what happens when *spatial / energy-scale* compositions are attempted.

---

## Section A: Compositions that WORK

### A1. Composed natural geometry (PK-mediated): Level 5 core / valence (Paper 17)

- **Source:** Paper 17, $G_{\text{total}} = G_{\text{nuc}} \ltimes G_{\text{core}}(R) \ltimes G_{\text{val}}(R, \text{core state})$. Verified for LiH (R$_\text{eq}$ 5.3 % at $l_\text{max}=2$, l-dependent PK), BeH$_2$ (R$_\text{eq}$ 11.7 %), H$_2$O (R$_\text{eq}$ 26 %).
- **Mechanism:** Fiber bundle. Each electron group is solved at its own natural geometry/focal length (Level 3 hyperspherical for the heavy-atom core at charge $Z_A$; Level 4 mol-frame hyperspherical for the valence pair at $R$), and the two are coupled *by classical scalar information passed between layers*: the core solve produces $Z_\text{eff}(r)$ and a Phillips-Kleinman pseudopotential $V_\text{PK}(r) = A e^{-Br^2}/r^2$ whose parameters $(A, B)$ are derived from core wavefunction expectation values and atomic ionization energy. The valence solver then runs on a modified single-electron potential. Inter-group antisymmetry is enforced by the PK barrier acting *as a one-body operator on the valence Hamiltonian*, not by Slater determinants spanning core + valence.
- **Bridge type:** *Other* — classical-scalar parametric coupling. The bridge is a number ($A$, $B$, $Z_\text{eff}(r)$ as a function), not an operator on a joint Hilbert space.
- **Scope:** Atoms with closed-shell cores (Li $[\text{He}]2s$, Be $[\text{He}]2s^2$, etc.), core $\subset$ Be$^{2+}$ scale, valence $\sim$ bonding scale. Works for hydrides through second / third row via frozen-core extension.
- **Limitations:** PK overcorrects systematically. R$_\text{eq}$ accuracy degrades from 5.3 % $\to$ 11.7 % $\to$ 26 % as block count grows ($1 \to 2 \to 4$ valence-style blocks). $l_\text{max}$-divergence: $R_\text{eq}$ drifts $+0.15$–$0.22$ bohr per added angular channel and the single-point energy violates the variational bound at $l_\text{max} \ge 2$ — confirmed structural to the composed/PK architecture (Track BQ, v2.0.32). Six PK modifications attempted, all failed (CLAUDE.md §3). Polyatomic lone-pair coupling fails at $Z_\text{eff} > 4$ (Paper 17, Sec VIII).
- **Computational signature:** Within-molecule $N_\text{Pauli} \sim O(Q^{2.2})$ across LiH/BeH$_2$/H$_2$O/HF/NH$_3$/CH$_4$. Uniform 11.11 $\times Q$ across 28 molecules in the composed library. ERI tensor block-diagonal (cross-block ERIs identically zero). 51 $\times$ – 1,712 $\times$ Pauli advantage over Gaussian baselines.

### A2. Balanced-coupled architecture (cross-center $V_{ne}$ via multipole expansion, Paper 19)

- **Source:** Paper 19, Sprint CD (v2.0.39–42). Balanced coupled LiH: 0.20 % energy error at $n_\text{max}=3$, only bound 4-electron FCI configuration among four tested coupling schemes. Verified at LiH, BeH$_2$, H$_2$O.
- **Mechanism:** PK is removed entirely; replaced by *explicit two-center one-body integrals* $\langle \psi_{nlm}^A | -Z_B / |\mathbf{r} - \mathbf{R}_B| | \psi_{n'l'm'}^A \rangle$ via Legendre multipole expansion of $1/|\mathbf{r}-\mathbf{R}_B|$ in the orbital frame on center $A$. The expansion *terminates exactly* at $L_\text{max} = 2 l_\text{max}$ by Gaunt selection rules (no truncation error). For non-collinear geometry (H$_2$O), the one-center matrix is computed in a $z$-axis-aligned frame and rotated into the lab frame via block-diagonal real-spherical-harmonic Wigner $D$-matrices.
- **Bridge type:** *Separable-multipole* (genuinely cross-focal in the sense that orbitals on $A$ now feel the potential of nucleus $B$, but mediated by an algebraic multipole expansion rather than by a joint two-body operator).
- **Scope:** First-row diatomics and small polyatomics with single-center valence orbital basis on each block. Works because (i) the Coulomb potential admits a separable multipole expansion, (ii) Gaunt selection rules truncate the expansion exactly at finite $L$, (iii) the orbitals on each center are still single-center single-focal-length objects — only the *potential* spans both centers.
- **Limitations:** Inflates Pauli count ($2.63 \times$ for LiH, $7.45 \times$ for H$_2$O — $O(B^2)$ in block count). Frozen-core failure mode for second-row hydrides: Sprint 7 NaH/MgH$_2$ balanced FCI shows monotonic over-attraction with no equilibrium (CLAUDE.md §3, Track CB-extension; root cause: frozen [Ne] core means valence sees full $Z=11/12$ from cross-center $V_{ne}$ without adequate core screening). $R_\text{eq}$ drifts structurally ($+0.053$ bohr per $n_\text{max}$, three times smaller than PK's drift but still nonzero). The architecture cannot reach beyond first row without a deeper fix to the screening of cross-center $V_{ne}$ by the frozen core.
- **Computational signature:** LiH 878 Pauli (Q=30), BeH$_2$ 2,652 (Q=50), H$_2$O 5,798 (Q=70). $\lambda$ (1-norm): LiH 74.1 Ha, BeH$_2$ 304.7 Ha, H$_2$O 1,509.3 Ha. Energy errors at $n_\text{max}=3$: LiH 0.20 %, BeH$_2$ 10.7 %.

### A3. Two-species nuclear-electronic tensor product (Paper 23 Track NI)

- **Source:** Paper 23, §V; Sprint NI 2026 (deuterium PoC at 26 qubits, 614 Pauli terms). Hyperfine validation: singlet–triplet gap 1.62 $\times 10^{-7}$ Ha = 21 cm line.
- **Mechanism:** Two qubit registers — proton on harmonic-oscillator $(n_r, l, m_l, m_s)$, electron on hydrogenic $(n, l, m, m_s)$ — tensored together. The Hamiltonian is $H = H_\text{nuc} \otimes I_e + I_\text{nuc} \otimes H_e + V_\text{fs} + H_\text{hf}$. Cross-register couplings are only spin–spin: $H_\text{hf} = A_\text{hf} \mathbf{I} \cdot \mathbf{S}$, four-qubit operator on (proton 0s spin up/down, electron 1s spin up/down). The finite-size correction $V_\text{fs}$ is a one-body *electronic* operator with $R_\text{nuc}$ as a classical scalar parameter.
- **Bridge type:** *Tensor-product-bilinear* on discrete spin labels only. The two registers' spatial / focal-length structure is spectator-additive (zero cross-coupling between proton spatial and electron spatial); only the spin labels couple.
- **Scope:** The 21 cm hyperfine constant (singlet-triplet splitting). Architectural validation. Demonstrates that a quantum register can host two species at vastly different binding scales (~13 orders of magnitude coefficient ratio between nuclear and hyperfine terms) within a single Pauli decomposition.
- **Limitations:** The $\sim 10^{13}$ coefficient ratio prevents single-pass variational diagonalization on near-term hardware; must solve in block-partitioned form (Born-Oppenheimer analog). Cross-register architecture handles spin–spin coupling but no spatial-spatial coupling exists (HF-3 verdict). Cannot derive recoil, finite-size charge corrections to $A_\text{hf}$, or Zemach.
- **Computational signature:** 16 qubits nuclear + 10 qubits electronic = 26 qubits. 592 nuclear + 10 electronic + 12 cross-register hyperfine = 614 Pauli terms. Coefficient ratio $\sim 2 \times 10^{13}$.

### A4. Balanced-coupled cross-center $V_{ne}$ in non-collinear geometry via Wigner D-matrix rotation (H$_2$O, Track CD v2.0.42)

Catalogued separately because Wigner-D rotation is a genuinely distinct piece of machinery from the multipole expansion.

- **Source:** Track CD (v2.0.42).
- **Mechanism:** $V_{ne}$ matrices on a non-collinear bond axis are obtained by computing the matrix in a $z$-aligned frame (where the multipole expansion is diagonal in $m$) and applying a block-diagonal Wigner D-matrix rotation in the real-SH basis. l=2 Wigner D extension applied (Track CM, v2.1.1) unblocks $n_\text{max}=3$ for all molecules.
- **Bridge type:** *Geometric-rotation*. A coordinate change *within* one focal length, not across focal lengths. The rotation acts on the angular labels of orbitals on a fixed center, which already exist on a single graph.
- **Scope:** Any composed/balanced architecture with non-collinear bond axes (H$_2$O, NH$_3$, CH$_4$, ...). Reaches general-l via the recursive complex Wigner d-matrix converted to real SH basis.
- **Limitations:** Strictly a coordinate change inside one focal length; provides nothing genuinely new at the multi-focal level. Its working-ness is a special case of A2 once you grant that A2 works.
- **Computational signature:** Identity for $l=0$, $3 \times 3$ permuted Cartesian for $l=1$, $5 \times 5$ for $l=2$.

### A5. Composed Pauli encoding for multi-center molecules (Tracks CU–CX, v2.3.0)

- **Source:** v2.3.0; LiF, CO, N$_2$, F$_2$, NaCl, CH$_2$O, C$_2$H$_2$, C$_2$H$_6$.
- **Mechanism:** Each block is mapped to a specific nucleus via `OrbitalBlock.center_nucleus_idx`; multi-center geometry is stored in `MolecularSpec.nuclei`. Composed Hamiltonian uses sub-block positions per nucleus. The bridge between blocks is again the PK pseudopotential, applied per-block, so the coupling is the same parametric type as A1.
- **Bridge type:** *Other* (classical-scalar parametric, same family as A1).
- **Scope:** All 28 molecules in the current library across both single-center and multi-center geometries with closed-shell cores.
- **Limitations:** Inherits A1's PK overcorrection. $R_\text{eq}$ accuracy is no better than the PK-mediated single-center pattern.
- **Computational signature:** Universal $N_\text{Pauli} = 11.11 \times Q$ across all 28 molecules, depends only on Q (not on block count, structure, or species).

### A6. Two-species JW encoding within nuclear shell model (Paper 23, deuteron and He-4)

- **Source:** Paper 23 §III–§IV.
- **Mechanism:** Proton register and neutron register tensored on a single qubit set with intra-species antisymmetrization plus a Minnesota NN potential acting on the joint register. Both species use the *same* HO basis at the same $\hbar\omega$, so this is *single-focal-length two-species* rather than multi-focal.
- **Bridge type:** *Tensor-product-bilinear*, but on a single shared HO focal length. Genuine many-body Coulomb (proton-proton) and Minnesota nucleon-nucleon are operators on the joint register.
- **Scope:** Closed-shell light nuclei in the HO basis. He-4 has $12.25 \times$ larger Hilbert space than deuteron but only $1.20 \times$ more Pauli terms.
- **Limitations:** Both species share a focal length by construction. Does NOT generalize to systems where proton and neutron have different binding scales without the same composition wall A3 hits.
- **Computational signature:** Deuteron 16 qubits / 592 Pauli, He-4 16 qubits / 712 Pauli, both $\sim 227$ MeV 1-norm class.

### A7. Multi-discrete-label gauge tensor products (Paper 32 §VIII.B; Papers 25 / 30 / ST-SU3)

- **Source:** Paper 25 (U(1) on $S^3$ Hopf graph), Paper 30 (SU(2) on $S^3$ = SU(2)), ST-SU3 (SU(3) on Bargmann $S^5$). Each is an independent Wilson lattice gauge construction on a *single* manifold.
- **Mechanism:** Gauge group element $U_e \in G$ on each edge of a fixed graph; plaquette / Wilson loop machinery; weak-coupling kinetic coefficient $1/(4 N_c)$ universal across $G \in \{U(1), SU(2), SU(3)\}$.
- **Bridge type:** *Algebraic-discrete-label*. The gauge element is a discrete label, not a spatial focal length. The "composition" is between matter labels and gauge labels on the same graph, all at one focal length.
- **Scope:** Works on each of the three gauge manifolds independently. *Does not work* across manifolds (G4b — see Section B).
- **Limitations:** Three separate constructions on three sub-manifolds. No cross-manifold unification; that is G4b's wall.
- **Computational signature:** Maximal-torus reduction connects U(1) and SU(2) on $S^3$; weak-coupling coefficient identical.

### A8. Sprint H1 almost-commutative extension (electroweak fiber on truthful CH triple)

- **Source:** Sprint H1 (May 2026), Paper 32 §VIII.C addendum.
- **Mechanism:** $\mathcal{T}_\text{AC} = \mathcal{T}_\text{GV} \otimes \mathcal{T}_F$ with $\mathcal{A}_F = \mathbb{C} \oplus \mathbb{H}$ on doubled $\mathcal{H}_F = \mathbb{C}^4_\text{matter} \oplus \mathbb{C}^4_\text{anti}$. KO-dim $3 + 6 = 9 \equiv 1 \pmod 8$, $J^2 = -I$, $JD = +DJ$ exact. Inner fluctuations $D \mapsto D + \omega + \epsilon' J\omega J^{-1}$ produce gauge sector (recovers Papers 25 / 30 at operator level) and Higgs sector (non-trivial iff $D_F$ has non-zero Yukawa).
- **Bridge type:** *Tensor-product-bilinear* in the Connes-Marcolli sense. Both factors are spectral triples; the tensor product is a textbook NCG construction at finite $n_\text{max}$.
- **Scope:** The construction is constructively well-defined at the architectural level. Higgs sector exists.
- **Limitations:** GeoVac does not autonomously select the Yukawa $Y$. With $Y = 0$: $\Phi = 0$ exactly. With $Y > 0$: $\| \Phi \|_\text{max} \in [0.05, 0.27]$. The Yukawa is a free input — Marcolli-vS-without-Higgs at the structural level.
- **Computational signature:** $\| \omega \| \sim 1.4$ for pure Higgs at $n_\text{max}=2$, $\| \omega_\text{gauge} \| \sim 1.3$ for pure gauge.

---

## Section B: Compositions that FAIL

### B1. Two-body coordinate operator across electron and nuclear registers (HF-3 recoil)

- **Source:** Sprint HF-3 (May 7, 2026).
- **Failure mechanism:** Track NI's V_fs / V_ne machinery is one-body on the electronic register with $R_\text{nuc}$ a *classical scalar parameter*. There is no operator $V_{ne}(\mathbf{r}_e, \mathbf{R}_n) = -Z/|\mathbf{r}_e - \mathbf{R}_n|$ that quantum-couples electron position to a *dynamical* nuclear position. Spatial-spatial cross-register Pauli count: 0 of 71. Therefore variational diagonalization of $H_e + H_p + V_{ne}$ on the joint register cannot exist as a single object — only the spectator-additive $H_e + H_p$ exists, plus the spin-spin and parametric-finite-size cross-pieces. Reduced-mass / recoil corrections like $(1+m_e/m_p)^{-3}$ on $|\psi(0)|^2$ are external, applied by hand using textbook Bethe–Salpeter.
- **Wall classification:** **W1** — cross-register two-body spatial coordinate operator missing.
- **What's missing structurally:** A two-body spatial coordinate operator on the joint (proton spatial $\otimes$ electron spatial) tensor product, of the form $V_{ne}(\hat{\mathbf{r}}_e, \hat{\mathbf{R}}_p)$ with both arguments operators on their respective registers. Equivalently, a center-of-mass / relative-coordinate split applied to the joint Hamiltonian that promotes $\mu_\text{red}$ to the conformal scale of the Fock projection on the relative-coordinate sector.
- **Status of failure:** *Empirically demonstrated* by direct register inspection (HF-3). No no-go theorem; the architecture *permits* the right operator to be added but the operator has not been built. Recasting the Fock projection to use $\mu_\text{red}$ instead of $m_e$ would reproduce $(1+m_e/m_p)^{-3}$, but that is no different in content from the textbook reduced-mass replacement. Demonstrably absent from the current code; structurally addable as a multi-week sprint.

### B2. Magnetization-distribution operator (HF-4 Zemach radius)

- **Source:** Sprint HF-4 (May 7, 2026).
- **Failure mechanism:** `hyperfine_coupling_pauli` is structurally point-like in the proton spatial coordinate. Repository-wide grep: zero hits on `Zemach`, `magnetization`, `r_Z`, `R_M`, `magnetic moment distribution`. The single proton-structure parameter `R_PROTON_BOHR` is the proton's *charge* radius, wired to the Foldy binding-energy correction (one-body on the electron register), not to the hyperfine coupling. There is no quantum operator for the proton magnetization distribution on the nuclear register and no separate $r_Z$ input.
- **Wall classification:** **W1** — cross-register two-body operator missing, but a *different sub-flavor* than B1: the missing operator here is not a spatial coupling, it is a *magnetization-density-weighted* coupling (weighted convolution of $\rho \star m$). Charge radius $R_p$ and Zemach radius $r_Z$ are *two different Layer-2 focal lengths*, both QCD-internal, both absent from the framework's current operator inventory. They are structurally distinct: $\rho(r)$ is the SU(3) quark color-electric distribution, $m(r)$ is the quark color-magnetic-moment distribution.
- **What's missing structurally:** A magnetization-distribution operator on the proton spatial register, plus a coupling of that distribution to the electron spin and electron contact density. Or equivalently, a smearing of the contact-density delta-function in the hyperfine coupling weighted by a separate $r_Z$ focal length.
- **Status of failure:** *Empirically demonstrated* by repository inspection. No no-go theorem; operator is in principle addable. The sub-wall is W1-magnetization (vs W1-coordinate of B1).

### B3. Multi-loop QED renormalization counterterms ($Z_2 / \delta m / Z_3$): LS-8a, HF-5

- **Source:** Sprint LS-8a (May 7) for hydrogen 2S two-loop SE; Sprint HF-5 (May 7) for the two-loop $a_e$ (Petermann coefficient).
- **Failure mechanism:** The bare iterated CC spectral action on Dirac-$S^3$ reproduces every structural feature of multi-loop QED — selection rules at all four vertices, $(\alpha/\pi)^2$ prefactor, three-topology decomposition (rainbow / SE / crossed for $a_e$, rainbow / crossed for SE), correct UV degree, correct sign at every $n_\text{max}$. Bare spectral sums are power-law divergent: $\sim N^{3.43}$ for two-loop SE, $\sim N^{3.79}$ for two-loop vertex. Two natural regularizations available within the GeoVac toolkit fail: (i) subtracting $[\Sigma_{1L}]^2$ removes only $<0.1$ % of the divergence (it lives in the connected diagram); (ii) Drake-Swainson asymptotic subtraction does not apply to power-law divergences (only logarithmic). Finite extraction of $C_{2S} = +3.63$ or $a_2 \approx -0.328$ requires field-theoretic $Z_2$ (wave-function), $Z_3$ (photon), and $\delta m$ (mass) renormalization counterterms not produced by the bare CC spectral action.
- **Wall classification:** **W2** — vertical UV / IR composition missing. The "focal length" scale here is energy-cutoff scale, not a bound-state focal length: composition runs from UV to IR through the renormalization group, and the framework supplies the bare integrand but not the counterterms.
- **What's missing structurally:** A renormalization machinery built from GeoVac internals — for example, an on-shell $Z_2^{(1)}$ extracted from the framework's own one-loop electron self-energy as a function of the cutoff $n_\text{max}$, used to subtract the disconnected piece at the spectral level — and a corresponding $\delta m$ from the on-shell mass renormalization. Equivalent in flat-space QED but not yet built on the discrete S^3 spectrum.
- **Status of failure:** *Empirically demonstrated* by power-law fits of the bare spectral sum across $n_\text{max} \in [2, 6]$. The agent flags this as the "LS-8a-renorm extension," a multi-sprint scope item explicitly deferred. No no-go theorem — the renormalization machinery is in principle implementable on the discrete spectrum, but doing so would import flat-space conventions that the structural-skeleton scope memo flags as importing data the framework does not natively produce.

### B4. Genuine cross-manifold spectral-triple unification: G4b ($\mathcal{T}_{S^3} \otimes \mathcal{T}_{S^5}$)

- **Source:** Paper 32 §VIII (G4 four-gap analysis), G4 scoping memo (May 6, 2026).
- **Failure mechanism:** Connes' tensor product of spectral triples requires both factors to be Riemannian spectral triples in the same sense. $\mathcal{T}_{S^3}$ carries the Camporesi–Higuchi Riemannian spinor Dirac (KO-dim 3). $\mathcal{T}_{S^5}$ as built in Paper 24 carries a *first-order complex Euler operator* on the holomorphic Hardy sector $H^2(S^5)$ — not a Riemannian Dirac. Two candidate $\mathcal{T}_{S^5}$ exist: Option A (round-$S^5$ Riemannian Dirac, KO-dim 5) gives a textbook tensor product but discards Paper 24's $(N,0)$-tower content and reintroduces calibration $\pi$ via standard Weyl asymptotics. Option B (Bargmann-Euler-based) keeps the $(N,0)$ content and the $\pi$-free certificate but is *not* a Riemannian spectral triple in the published sense; KO-dim is not classically defined; and there is no published prescription in the NCG literature for tensoring a Riemannian spectral triple with a Hardy-sector first-order complex object.
- **Wall classification:** **W2** — but a structurally different sub-flavor than B3. Here the wall is *categorical-mismatch between two types of geometric structure* (Riemannian-spinor vs Hardy-sector-complex-analytic), not a renormalization gap. This sits at the spectral-triple level rather than at the perturbation-theory level. The four-layer Coulomb / HO asymmetry of Paper 24 §V resurfaces here as the fourth and deepest layer.
- **What's missing structurally:** Either an extension of NCG to handle Riemannian-Hardy mixed tensor products, or a unifying ambient manifold whose GeoVac sub-structure contains both $S^3$ Coulomb and $S^5$ Bargmann as natural sub-objects.
- **Status of failure:** *Documented structural obstruction*, recorded as an open question rather than a sprint target. The closest published analogs (Hawkins 2000 on the unit disk, Hekkelman 2022 / 2024, Hekkelman-McDonald 2024) stay within either Berezin-Toeplitz on a single manifold or Riemannian on a single manifold — never crossing.

### B5. Single-center molecular encoding (Sturmian shared-$p_0$ basis): GUARDRAIL Papers 8–9

- **Source:** Papers 8–9 Sturmian Structural Theorem; Track DF Sprints 4–5 (re-derived this at significant cost).
- **Failure mechanism:** When all electrons of a heteronuclear molecule are placed in a single Sturmian basis with a shared exponent $k$ (single $p_0$), the matrix elements satisfy $H_{ij} \propto S_{ij}$ in the shared-$p_0$ basis, so the eigenvalues are $R$-independent: a generalized eigenvalue problem with $H = \lambda S$ trivially gives the same spectrum at every $R$. Equivalently, no PES exists: the Hamiltonian's spectrum cannot depend on internuclear distance because the basis already carries the only scale.
- **Wall classification:** Pre-W1 — this is the single-focal-length wall. A single $p_0$ cannot represent *two* binding scales (heavy core, light bond). It is the natural-geometry hierarchy's *foundation*: composition is necessary because of B5, not because of B1–B4.
- **What's missing structurally:** *Multiple* focal lengths. The composed framework (A1) is the response.
- **Status of failure:** **PROVEN no-go theorem** (Sturmian Structural Theorem, Papers 8–9). Track DF spent six sprints re-deriving this in three molecular variants (single-center 33.7 % R$_\text{eq}$ error, charge-center 48.2 % energy error, heterogeneous Löwdin destroys Gaunt sparsity at $1{,}711$ vs $120$ Pauli). All negative. CLAUDE.md §3.5 made this a guardrail paper.

### B6. Graph-concatenation molecular Hamiltonian: GUARDRAIL FCI-M

- **Source:** Paper FCI-M (29-version diagnostic arc, v0.9.8–0.9.37).
- **Failure mechanism:** When two atomic graphs (each at its own scale) are concatenated by adding cross-edges to form a "molecular graph," the graph Laplacian's intra-atom kinetic energy is *R-independent* — the atomic lattice structure does not change as atoms approach. The PES is monotonically attractive with no equilibrium.
- **Wall classification:** A different flavor than W1/W2: this is a *naive composition* failure (treat two single-focal-length objects as if they could be glued without re-projecting). The natural geometry hierarchy's response (Levels 2–4 with prolate-spheroidal / molecular-frame hyperspherical coordinates) sidesteps this by re-doing the projection on a coordinate system that handles both centers.
- **What's missing structurally:** A coordinate system spanning both centers' scales. Levels 2 and 4 supply this for one- and two-electron systems. The composed framework (A1) supplies it for core+valence.
- **Status of failure:** **PROVEN no-go theorem** (R-independent kinetic operator).

### B7. Frozen-core balanced PES for second-row hydrides (NaH, MgH$_2$ at $n_\text{max}=2$)

- **Source:** Sprint 7 balanced second-row (v2.19.4, April 2026).
- **Failure mechanism:** With a frozen [Ne] core absorbing 10 electrons into a screening function, the cross-center $V_{ne}$ from the multipole expansion couples valence orbitals to the *full* nuclear charge ($Z=11$ Na, $Z=12$ Mg) on the other center, not to a screened effective charge. The valence is therefore over-attracted by the cross-center $V_{ne}$ and the PES has no equilibrium. LiH succeeds because it has an *explicit* core block providing dynamical screening.
- **Wall classification:** Sub-wall of W1 — the missing piece is *two-body screening of a cross-center potential by a frozen core*. The frozen-core architecture provides *one-body* screening (via its $Z_\text{eff}(r)$ on the same-center potential) but not *cross-center* screening of the bare $V_{ne}$ from the other nucleus. Once you cross centers, the frozen-core Z_eff(r) machinery doesn't reach.
- **What's missing structurally:** Cross-center screening: the frozen-core electronic density on center $A$ attenuating the bare $V_{ne}$ from center $B$ that the *valence* orbital on $A$ sees. Sprint 7b partially addressed this for SO splitting (screened $\xi_\text{SO}$ via FrozenCore $Z_\text{eff}(r)$) but not for the PES.
- **Status of failure:** *Empirically demonstrated.* Documented as a structural limitation of balanced + frozen-core for PES.

### B8. Lone-pair coupling at $Z_\text{eff} > 4$ (Paper 17 §VIII)

- **Source:** H$_2$O composed solver, Track CD' diagnostic.
- **Failure mechanism:** Bond-bond Slater-integral coupling validated at $\sim 0.5$ Ha for BeH$_2$. Lone pair coupling with $Z_\text{eff}=6$ (oxygen) produces $-28$ Ha bond-lone and $-15$ Ha lone-lone matrix elements — exceeding the total electronic energy. The Slater integrals on a single high-$Z_\text{eff}$ orbital block are physically unphysical at this scale.
- **Wall classification:** Sub-wall of A1's classical-scalar parametric coupling: the *parameter* (Slater integral on a single block) is not the physically relevant quantity at high $Z_\text{eff}$. The correct treatment requires orbital-level exchange, which the composed framework currently doesn't have.
- **What's missing structurally:** Inter-fiber exchange that distinguishes the orbital structure across blocks (rather than treating a block as a single $Z_\text{eff}$-scaled scalar).
- **Status of failure:** *Empirically demonstrated.* Disabled in production for $Z_\text{eff} > 4$.

### B9. Coupled composition without cross-center $V_{ne}$ (Track CB negative result)

- **Source:** Track CB (v2.0.37).
- **Failure mechanism:** Replace PK with explicit cross-block ERIs (cross-block electron-electron repulsion) without adding cross-center electron-nucleus attraction. Result: 29 % FCI energy error, *worse* than PK's 15 %. Reason: cross-block ERIs add inter-group $V_{ee}$ repulsion without the balancing $V_{ne}$ attraction. Hamiltonian is energetically unbalanced.
- **Wall classification:** Sub-W1 — the missing operator is again cross-center $V_{ne}$, exactly the class of B1, but here within the *balanced-coupled architecture* before the multipole-expansion fix landed. Track CD then closed it (positive, A2).
- **What's missing structurally:** Same as A2's mechanism — multipole expansion of cross-center $V_{ne}$. Adding it converted B9 (negative) into A2 (positive).
- **Status of failure:** *Empirically demonstrated*, *resolved by A2*.

### B10. Single-constant graph-to-continuum QED projection (C × F_2 → α/(2π))

- **Source:** Sprint GN-QED (April 2026), Spectral-Zeta Projection Sprint.
- **Failure mechanism:** Graph-native QED (scalar photon) on the Fock graph computes $F_2$ values that are π-free in $\mathbb{Q}[\sqrt{2}, \sqrt{3}, \sqrt{6}]$ and decrease as $F_2 \sim n_\text{max}^{-0.57} \to 0$. The continuum target is $\alpha / (2\pi)$. Multiplying by a single graph-to-continuum projection constant $C$ does not work: $C \times F_2$ *grows* with $n_\text{max}$ (0.053 at $n=3$, 0.075 at $n=4$), diverging from $\alpha/(2\pi)$ by 46–65×. Three independent projection constants $C_\text{VP}$, $C_\text{SE}$, $C_{F_2}$ have *different* power-law exponents (1.48 vs 1.48 vs 0.57); the projection is topology-dependent.
- **Wall classification:** Sub-W2 with a discrete-graph flavor: the wall is between graph-level and continuum-level QED. The graph cleanly produces the algebraic content but the projection back to the continuum requires per-diagram-topology-specific factors, not a single multiplicative rescaling.
- **What's missing structurally:** Either a vector photon promotion at $1/(4\pi)$ per loop (S² Weyl exchange constant of the Hopf base, Paper 33's 1+6+1 partition) or recognition that no single scalar $C$ can reconcile graph-native QED with continuum $a_e$.
- **Status of failure:** *Resolved as scope statement* — Paper 33 (1+6+1 partition) explains why the projection is per-diagram. Demoted from "negative" to "structurally clarified": vector-photon promotion (Paper 33's six recovered selection rules) closes 7 of 8 with $1/(4\pi)$ calibration per loop, and the eighth comes from Dirac spinor phase constraint at zero cost.

### B11. Co-exact mode q-labeling for SO(4) channel count, Ward identity, charge conjugation (transverse photon analysis)

- **Source:** Mode-resolved transverse self-energy analysis (May 2026).
- **Failure mechanism:** Co-exact eigenvectors of the plaquette Laplacian $d_1 d_1^T$ have correct *eigenvalues* (approaching $n(n+2)$) but wrong *support* — confined to adjacent shell-pairs ($|\Delta n| \le 1$) by graph connectivity. Triangle inequality $|n_1 - n_2| \le q \le n_1 + n_2 - 2$ fails at 63–86 % for high-q modes. Ward identity diagnostic $\| [\Sigma_T, H_0] \| / \| \Sigma_T \| \approx 0.58$.
- **Wall classification:** Sub-W2 in the photon sector. The 5/8 graph-intrinsic ratio of selection-rule recoveries is maximal without promoting 1-cochains to vector bundle sections indexed by $(q, m_q)$.
- **What's missing structurally:** Vector bundle structure on the photon sector with explicit $(q, m_q)$ angular momentum labels.
- **Status of failure:** *Documented structural negative*, resolved by Paper 33's vector-photon-promotion analysis.

### B12. Yukawa selection from GeoVac internals (Sprint H1 falsifier; G2/G3 collapse)

- **Source:** Sprint H1 (May 6); G2/G3 collapse in Paper 32 §VIII.B/C.
- **Failure mechanism:** The AC inner-fluctuation construction is well-defined and produces a Higgs sector for any non-zero Yukawa $Y$, but no GeoVac data couples to the off-diagonal $\mathbb{C} \leftrightarrow \mathbb{H}$ block of $D_F$ that would select $Y$. Sprint G3 sharpens: $\gamma_\text{GV}$ (the chirality grading on $\mathcal{H}_\text{GV}$) and $\gamma_F$ (KO-dim 6 grading on $\mathcal{H}_F$) are independent commuting $\mathbb{Z}_2$'s with operator-norm residual $\| \gamma_\text{GV} \otimes I_F - I_\text{GV} \otimes \gamma_F \|_\text{op} = 2$ exactly at every $n_\text{max} \in \{1, 2, 3\}$. The Yukawa lives on the $\gamma_F$ flip; no GeoVac-side data couples to that flip.
- **Wall classification:** Inner-factor parameter selection — the *inner-factor input data* tier of Paper 18 (sixth tier, added 2026-05-07). This is a *fifth wall* class, distinct from W1 (cross-register) and W2 (UV/IR), at the level of *what data the framework's input contains*.
- **What's missing structurally:** A separate axiom or structure that supplies inner-factor data ($A_F$ choice, Yukawa values, hypercharges, generation count, renormalization counterterms). The "second packing axiom" question.
- **Status of failure:** *PROVEN structural absence* via G3's commuting-$\mathbb{Z}_2$ theorem ($\| \cdot \|_\text{op}$ residual exactly 2). The two Z$_2$'s cannot be identified by tensor-product factorization.

---

## Section C: Are the working cases special cases of each other?

The working cases organize naturally into a small lattice of bridge types. Here is the reduction map.

### C1. A1 (PK-mediated composed) and A2 (balanced-coupled) are *competing* responses to the same problem

Both A1 and A2 address the multi-block problem (single-center valence orbitals on each block, bond between centers). A2 is *not* a special case of A1: A2 *eliminates* PK (the central bridge of A1) and replaces it with explicit cross-center $V_{ne}$. They are two architectures, neither subsumed by the other. Reduction map: A1 = "PK as parametric bridge"; A2 = "multipole-expanded $V_{ne}$ as separable-multipole bridge"; bridges are categorically different. The negative result B9 establishes that A1's PK is *not* algebraically removable from A1 without adding A2's mechanism (otherwise the Hamiltonian is unbalanced).

### C2. A4 (Wigner D rotation in non-collinear $V_{ne}$) is a coordinate change inside A2, not a separate bridge

Wigner D-matrix rotation acts on angular labels of single-center orbitals. The orbitals already exist on a single graph at a single focal length. The rotation produces the cross-center matrix in the lab frame from the diagonal-$m$ frame. A4 reduces cleanly to: "A2 + a coordinate change on the single-center orbital basis." Strictly speaking A4 is *not* genuinely cross-focal; it is necessary for A2 to handle non-collinear bond axes but is not multi-focal in the audit's sense.

### C3. A5 (multi-center composed) is A1 with `OrbitalBlock.center_nucleus_idx` plumbing

A5 generalizes A1 from single-center to multi-center by adding nucleus-indexing to each block. The bridge is still PK (parametric scalar). Same bridge type as A1. Reduction: A5 ⊆ A1 + multi-center plumbing.

### C4. A3 (nuclear-electronic tensor product, Paper 23 NI) and A6 (proton + neutron in same HO basis) are different

A3 uses *different focal lengths* (proton in HO basis at MeV scale, electron in hydrogenic at Ha scale) but couples them only by spin-spin (a discrete-label bridge), so it is single-focal-bridge despite multi-focal substrate. A6 uses a single shared HO focal length for both species and couples them by genuine many-body operators (Coulomb, NN). Different bridge types: A3 is *tensor-product-bilinear-on-spin-only*; A6 is *tensor-product-bilinear-on-full-spatial-and-spin*. A6 is genuinely many-body; A3 is many-body only on the spin sector.

### C5. A7 (multi-discrete-label gauge tensor products) is *single-focal-bridge across multi-focal substrate*

The three Wilson constructions (Papers 25 / 30 / ST-SU3) live on three sub-manifolds. The "tensor product" between them is the orthogonal-direct-sum sense ("we have three independent constructions") rather than the Connes-Marcolli sense. Genuine cross-manifold spectral-triple unification (G4b) does not exist (B4). So A7 in its current form is *three independent constructions*, not one composed object. The closest thing to a unified A7 is the U(1) × SU(2) co-located on $S^3$ (Papers 25 + 30), which is a *single-focal* electroweak target (Paper 32 §VIII.B's "most concrete near-term target") and is not yet completed.

### C6. A8 (Sprint H1 AC extension) is the *one* genuinely-cross-focal working composition in the spectral-triple sense

A8 is a Connes-Marcolli tensor product of *two* spectral triples ($\mathcal{T}_\text{GV}$ on truncated $S^3$, $\mathcal{T}_F$ on a finite-dimensional inner factor). KO-dim arithmetic works ($3+6=9 \equiv 1 \pmod 8$). The $J = J_\text{GV} \otimes J_F$ structure is canonical. This is the *gold standard* of multi-focal composition in the framework, and its limitation (Yukawa undetermined) is recorded as the inner-factor-input-data tier (B12), not as a structural failure of the composition itself. A8 reduces no further: it's the canonical NCG construction.

### C7. Net reduction lattice

The bridge types reduce to four genuinely distinct classes:

1. **Classical-scalar parametric** (A1, A5; PK, $Z_\text{eff}$, $A_\text{hf}$ scalar). Used when the inter-block coupling can be summarized by a scalar function passed between layers.
2. **Separable-multipole** (A2, A4 as a coordinate-rotation extension). Used when the cross-block potential admits an algebraic multipole expansion that terminates by Gaunt selection rules.
3. **Tensor-product-bilinear** on discrete labels (A3 spin-only, A6 spatial-and-spin within shared focal length, A7 gauge labels, A8 Connes-Marcolli inner factor). The "bilinear" is in the operator sense — two registers, one operator coupling them.
4. **Geometric-rotation** (A4). Coordinate change inside one focal length, not a multi-focal bridge per se.

Among these, only (3) — tensor-product-bilinear — is a genuinely-multi-focal bridge in the spectral-triple sense, and only A8 (Connes-Marcolli on two spectral triples) is a *full* multi-focal composition. (1) and (2) are *partial* bridges that hide the multi-focal-ness inside a parameter or expansion that is computed from one focal length and consumed in another. The walls in Section B are about what (3) *cannot* do natively: B1 / B2 want spatial-spatial bilinear couplings that don't exist; B3 wants UV/IR bilinear couplings that don't exist; B4 wants cross-manifold bilinear couplings that don't exist; B12 wants Yukawa-data selection that doesn't exist.

---

## Section D: Sub-walls within W1 and W2

After cataloguing the failures, the W1 / W2 split is correct as a first cut but compresses three or four genuinely distinct sub-walls. The proposed refinement:

### W1 splits into at least three sub-walls

- **W1a — Spatial-coordinate cross-register (B1 recoil).** Missing operator: $V(\hat{\mathbf{r}}_e, \hat{\mathbf{R}}_n) = -Z/|\mathbf{r}_e - \mathbf{R}_n|$ on the joint register. Tensor product of *spatial* operators across registers at different focal lengths. The architecture *permits* this; the operator has not been built.
- **W1b — Magnetization-distribution cross-register (B2 Zemach).** Missing operator: a magnetization-density operator on the proton register convolved with the electron contact density. *Different* from W1a: not simply a coordinate operator but a density-weighted convolution. Charge radius $R_p$ and Zemach radius $r_Z$ are categorically distinct focal lengths (color-electric vs color-magnetic). Charge radius is wired (Foldy correction); Zemach is not.
- **W1c — Cross-center screening of $V_{ne}$ by frozen core (B7 NaH/MgH$_2$).** Missing operator: the frozen-core electronic density on center $A$ attenuating the *bare* $V_{ne}$ from center $B$ that the valence orbital on $A$ feels. The frozen-core $Z_\text{eff}(r)$ machinery handles *same-center* screening; it does not reach across centers without explicit work.

These three are unified by: *the framework lacks two-body operators that couple spatial registers at different focal lengths*, but they break into three sub-cases by what kind of two-body operator (coordinate vs density-weighted vs cross-center-screening).

### W2 splits into at least two sub-walls

- **W2a — UV / IR vertical composition / renormalization (B3 LS-8a, HF-5).** Missing structure: counterterms ($Z_2$, $Z_3$, $\delta m$) that compose a UV cutoff prescription with an IR observable. The bare iterated CC spectral action reproduces the divergent integrand correctly but lacks the subtraction machinery. This is *vertical* in the sense of running between scales rather than between manifolds.
- **W2b — Cross-manifold spectral-triple composition (B4 G4b).** Missing structure: a published NCG framework that handles tensor products between Riemannian-spinor and Hardy-sector-complex-analytic spectral triples. The four-layer Coulomb / HO asymmetry of Paper 24 §V is the deeper root. This is *horizontal* in the sense of running between two different geometric manifolds at the same level of QFT.

W2a and W2b have categorically different mathematical content: W2a is about renormalization-group running; W2b is about category-theoretic compatibility of spectral-triple tensor products. They share only the labeling "vertical-or-horizontal composition that the framework cannot do natively."

### A potential fifth wall: W3 — inner-factor input data

B12 (Yukawa selection from GeoVac internals; G2/G3 collapse) does not fit cleanly into either W1 or W2. The failure mode is *not* that an operator coupling registers is missing; the operator (off-diagonal $D_F$ block) is well-defined and produces a Higgs sector once you specify it. The failure is that *no GeoVac-internal data selects the value of that operator*. This is a different wall: *parameter selection within a structurally permissible operator class*. CLAUDE.md memory `inner_factor_mellin_engine.md` calls this the "inner-factor input data" tier, sitting alongside Yukawa, hypercharge, generations, renormalization counterterms, and the choice of $A_F$ itself.

Whether to call this a third wall (W3) or a sub-wall of W2 is a framing choice. The argument for W3: B12's failure mode is structurally distinct from B3/B4 — B12 is about *selecting a value within a permissible class*, B3 is about *generating counterterms that close a divergent series*, and B4 is about *category-mismatch between two geometric structures*. They share only the negative form ("the framework supplies the structure but not the calibration"). The argument for W3 ⊂ W2: every multi-loop renormalization counterterm *is* an inner-factor input in this sense (Paper 18 fourth tier puts $Z_2 / \delta m$ in the same tier as Yukawa). Under that reading W2a *is* W3 restricted to RG counterterms.

The cleanest taxonomy may be:

- **W1** = cross-register two-body spatial operator (sub-walls a, b, c).
- **W2** = cross-scale composition (sub-walls a UV/IR-renormalization, b cross-manifold).
- **W3** = inner-factor parameter selection (Yukawa, hypercharge, generation count, $A_F$ choice).

with W2a and W3 partially overlapping (RG counterterms are inner-factor data).

### Other candidate sub-walls flagged but not load-bearing

- **W1d — Lone-pair coupling at high $Z_\text{eff}$ (B8).** Resembles W1c but is mediated by classical-scalar coupling rather than cross-center potential. Currently mitigated by disabling the coupling at $Z_\text{eff} > 4$. Could be subsumed into W1c if framed as "cross-block coupling lacks orbital-resolution-aware screening at high effective charges."
- **B10 / B11 (graph-to-continuum QED projection, transverse photon mode-labeling).** These are now resolved as scope statements via Paper 33's 1+6+1 partition (vector-photon promotion at $1/(4\pi)$ per loop closes them). They are *not* live walls; they are clarified scope.

---

## Section E: Honest scope

What this audit covered: the project's record of multi-focal compositions through v2.30.0, including all five-track HF sprint memos, H1, LS-8a, Paper 32 §VIII.B/C, and the four "skeleton" memos. Cataloguing is exhaustive within the scope of the W1/W2 framing question.

What this audit did *not* cover deeply:

1. **The Sprint MR-A/B/C master Mellin engine domain partition** (M1/M2/M3 sub-mechanism partition). This is foundational for understanding *which* transcendentals appear where in compositions, but is taxonomic-of-transcendentals rather than taxonomic-of-compositions, so it sits orthogonally to the W1/W2 question.
2. **Paper 17 §VIII rho-collapse cache** as a working multi-focal mechanism. Listed under A1 (it is a sub-mechanism of A1's PK pipeline) but not separately catalogued.
3. **Paper 18 §IV six-tier exchange-constant taxonomy** as a frame for the failures. This audit references the "inner-factor input data" tier (sixth tier, added 2026-05-07) but does not attempt to map every failure to a Paper 18 tier; that's a separate exercise.
4. **The G4a sprint scoping** (Paper 32 §VIII.D candidate addendum, $A_F = \mathbb{C} \oplus \mathbb{H} \oplus M_3(\mathbb{C})$, predicted positive-thin). Documented but not yet attempted as a sprint, so no pass/fail data exists. Listed neither in Section A nor Section B.
5. **Sprint MR-B's closed form $\epsilon(t)$ for the modular residual** as evidence that the M2 mechanism is predictive within its own domain. This refines the Paper 35 Prediction 1 graduation but is not itself a multi-focal composition.

What I am uncertain about:

- **Whether B7 (frozen-core balanced PES negative) should be a fully separate sub-wall (W1c) or whether it is downstream of W1a's missing two-body operator.** Argument for sub-wall: the architectural fix would be different — adding a *screened* cross-center $V_{ne}$ from a Z_eff(r) generated by the frozen core, rather than a full quantum two-body coordinate operator. Argument for downstream-of-W1a: the deeper fix would be the same as W1a (give cross-center $V_{ne}$ a dynamical second register).
- **Whether C5's "A7 in current form" should count as a working multi-focal composition or as a placeholder.** It has the form of multi-focal (three sub-manifolds), but each sub-construction is single-focal and the bridges between them don't exist. A more honest reading: A7 is "three single-focal compositions tabulated together," not one multi-focal construction. I have classified it as A7 ("works" in the sense that each piece is a well-defined Wilson construction) but flagged the absence of cross-manifold integration via B4 (G4b).
- **Whether B10 and B11 should be removed from Section B entirely** since they are now resolved by Paper 33. I have left them as "documented sub-walls that are now structurally clarified" because the *original* compositions they describe (single-constant projection, co-exact mode labeling) are still inaccessible — Paper 33 introduces a *different* mechanism (vector-photon promotion) rather than fixing the original.

What I would want to follow up on:

- **Whether a "second packing axiom" candidate has emerged anywhere in the project record** (CLAUDE.md memory mentions it as the speculative frontier; the audit didn't find any concrete proposal). If yes, that would suggest a route into W3 (inner-factor input data).
- **Whether the W3 framing in Paper 18's six-tier taxonomy is exhaustive** — i.e., are there input-data items that don't fit into any tier? This is a question for a separate Paper 18 audit.
- **Whether A2's success (separable-multipole) generalizes to other cross-focal interactions besides $V_{ne}$** — for example, whether a similar multipole expansion would close W1a (recoil) by giving an algebraic decomposition of $V(\mathbf{r}_e, \mathbf{R}_p)$ in proton-CoM coordinates. The HF-3 memo flags this as a possible direction but does not attempt it.

No factual errors found in CLAUDE.md while reading. One framing observation: CLAUDE.md §3.5 GUARDRAIL papers are presented as a "two-paper" set (Papers 8–9 + FCI-M); this audit confirms they are independent failure modes, both load-bearing, but distinguishable as B5 (single-focal-length wall, Sturmian) vs B6 (R-independent kinetic, graph-concatenation). Worth keeping separate in the failure taxonomy.

---

**End of Phase A Track 1 audit.**
