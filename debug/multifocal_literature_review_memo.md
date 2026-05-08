# Multi-Focal Composition: Literature Review

**Date:** 2026-05-07
**Author:** Sub-agent (literature review track)
**Frame:** GeoVac is an almost-commutative spectral triple in the Marcolli–van Suijlekom 2014 lineage (Latrémolière propinquity convergence proven on the Camporesi–Higuchi spectral triple, Loutey 2026, Paper 38). Two walls block the multi-focal frontier:
- **W1 — Cross-register spatial coupling:** the framework couples discrete labels but has no native two-body coordinate operator across distinct Fock projections (e.g. quantum-mechanical nucleus × electron, proton magnetization × electron 1s).
- **W2 — Vertical UV/IR composition:** the framework reproduces UV-divergent integrands of multi-loop QED but cannot autonomously generate $Z_2/\delta m$ counterterms.

This memo surveys published machinery that could be tested against either wall.

---

## 1. Multi-λ Coulomb-Sturmian basis sets for multi-particle systems  *(targets W1)*

### State of the art

The single-λ Coulomb-Sturmian school is mature: Avery, J. (1989) *Hyperspherical Harmonics: Applications in Quantum Theory* (Kluwer); Avery, J. and Avery, J. E. (2006) *Generalized Sturmians and Atomic Spectra* (World Scientific). The seminal multi-particle paper is **Avery, J., "Many-particle Sturmians," J. Math. Chem. 21, 285 (1997)**, generalized to atoms in Avery & Avery (2003) *J. Phys. Chem. A* 107, 7218 [doi:10.1021/jp040225m]. Many-center extension via hyperspherical harmonics: Avery J. & Avery J. (2004) *Int. J. Quantum Chem.* 100, 121 [doi:10.1002/qua.10820]. Recent active line: **Mitin, A. V. et al., "Quantum chemistry with Coulomb Sturmians: Construction and convergence at HF level," Phys. Rev. A 99, 012512 (2019)** [arXiv:1811.05777]; **Frapiccini, A. L., et al., "Sturmian expansions for two-electron atomic systems," Phys. Rev. A 82, 042503 (2010)**; Abdouraman et al. "Sturmian bases for two-electron systems in hyperspherical coordinates" (Few-Body Systems, 2016).

**The crucial structural fact:** in every published "many-particle Sturmian" formulation (Avery 1997 onward), all members of the basis set are *isoenergetic* — chosen so that all basis functions correspond to the same energy of the state being represented. The "weighting factors $\beta_\nu$" in the generalized Sturmian construction are weights on a *single* zeroth-order potential $V_0$, not independent exponents per particle. The method optimizes orbital exponents jointly, but the underlying $\lambda$ is a one-parameter family across all particles. For three-body calculations with arbitrary masses (Frapiccini, Granados-Castro, Gasaneo, Mitnik et al.), the *coordinate system* uses different reduced masses per pair, but the Sturmian generating potential is still chosen to give one dominant scale.

**Honest assessment:** I could not locate a published $N$-particle Sturmian formalism in which particle $i$ uses Sturmian functions at exponent $\lambda_i \neq \lambda_j$ in the literal sense W1 needs (electron at $\lambda_e = Z_e/n$, nucleus at $\lambda_N = Z_N/n_N$ with $\lambda_e \neq \lambda_N$). The tensor-product two-particle constructions of Frapiccini–Granados-Castro–Gasaneo (e.g. Granados-Castro & Gasaneo 2014, Few-Body Systems 55, 1029) discretize each coordinate independently with its own scale, but the angular machinery is built for the same scale per coordinate. The closest published thing to a literal "multi-λ" object is the *multi-channel* Sturmian basis used in scattering, where different channels carry different λ — but that is a labelling of asymptotic states, not a basis for the bound-state Hamiltonian. Porting GeoVac's discrete cross-register (Track NI deuterium PoC) to a multi-λ Sturmian language would therefore require: (a) explicitly defining the cross-particle two-body operator, (b) constructing $\langle Y_{n_e l_e m_e}^{(\lambda_e)} Y_{n_N l_N m_N}^{(\lambda_N)} | V_{eN}(r_e - R_N) | \cdot \rangle$ via Shibuya–Wulfman-like multipole expansion **with mismatched exponents**, (c) confirming it terminates at $L_{\max} = 2 \max(l_e, l_N)$ as the matched-exponent case does. This is the *mathematical* counterpart of GeoVac's open W1, not a published bypass.

**Speculative class: Moderate** — the math is published only at the single-λ level; multi-λ is a natural-but-unpublished extension. Sprint scope: 4–8 weeks to draft the integral structure, longer if a clean Gaunt-style termination doesn't survive.

---

## 2. Tensor products of metric spectral triples (Farsi–Latrémolière)  *(targets W1 in NCG language)*

### State of the art

The relevant published thread: **Latrémolière, F., "The Gromov–Hausdorff propinquity for metric Spectral Triples," Adv. Math. 404 (2022), 108393 [arXiv:1811.10843]**; **Farsi, C., Landry, T., Latrémolière, F., Packer, J., "Convergence of inductive sequences of spectral triples for the spectral propinquity," Adv. Math. 437 (2024), 109442 [arXiv:2301.00274]**; **Farsi, C., Latrémolière, F., "Collapse in Noncommutative Geometry and Spectral Continuity," arXiv:2404.00240 (2024)**; **Farsi, C., Latrémolière, F., "Continuity for the spectral propinquity of the Dirac operators associated with an analytic path of Riemannian metrics," arXiv:2504.11715 (April 2025)**; **Latrémolière, F., "Spectral continuity of almost commutative manifolds for the C¹ topology on Riemannian metrics," arXiv:2603.19128 (March 2026)**. Closely related: Aguilar, K. (arXiv:1907.07357) on quantum metrics on $C(X) \otimes A$ for $A$ AF.

### What is and is not proved

Latrémolière's 2022 propinquity paper introduces the metric spectral triple framework but, on a careful reading of the abstract and the public material, focuses on the metric structure of a *single* spectral triple, not products. The 2024 Farsi–Latrémolière "Collapse" paper explicitly mentions "collapse of product of spectral triples with one Abelian factor" as an example, indicating products *appear* but are not the headline. Aguilar's 2019 paper handles $C(X) \otimes A$ where $A$ is AF; matrix compact quantum metric spaces are stable under minimal tensor products (cf. Latrémolière's Bures-distance line). The Latrémolière 2026 paper (arXiv:2603.19128) is the most relevant: it proves spectral continuity for **almost-commutative manifolds**, i.e. products $T_M \otimes T_F$ where $T_M$ is the canonical Riemannian spin spectral triple and $T_F$ is finite-dimensional. This is the Connes-SM paradigm and is exactly the structure GeoVac sits in (per Sprint H1 and §VIII.C of Paper 32).

**What is NOT in the published literature:** a constructive theorem for the Latrémolière propinquity of two infinite-dimensional non-flat metric spectral triples — e.g. the GeoVac S³ Camporesi-Higuchi triple tensored with a hypothetical S⁵ Bargmann triple. The 2026 paper handles "almost-commutative" with a finite second factor. Two infinite second factors is not handled. Latrémolière's 2018 paper says matrix spectral metric spaces are stable under "external product" but this is a finite-tensor statement.

**Honest assessment:** the tensor product of two infinite metric spectral triples is at the edge of what is published. For GeoVac's W1, the most useful theorem would be: given $(\mathcal{A}_1, \mathcal{H}_1, D_1)$ and $(\mathcal{A}_2, \mathcal{H}_2, D_2)$ both Latrémolière-propinquity-converging, does $(\mathcal{A}_1 \otimes \mathcal{A}_2, \mathcal{H}_1 \otimes \mathcal{H}_2, D_1 \otimes I + \gamma_1 \otimes D_2)$ converge? For the *almost-commutative* case (one factor finite) Latrémolière 2026 says yes. For two infinite factors, no published theorem covers it. The framework's G4b cross-manifold blocker (Paper 32 §VIII.C) is structurally the same wall.

**Speculative class: Moderate-to-Weak** — Latrémolière 2026 covers the almost-commutative slice cleanly; this directly serves Connes-style internal $D_F$ Higgs constructions but does NOT serve true cross-manifold composition. The infinite-infinite case appears genuinely open and would itself be a substantial NCG theorem.

---

## 3. Berezin–Toeplitz quantization on non-flat manifolds (Hekkelman, McDonald, Connes–vS)  *(targets W2 in NCG language)*

### State of the art

The closest named published works:

- **Hekkelman, E.-M. (2022 MSc thesis), "Truncated Geometry on the Circle"**, Radboud University.
- **Hekkelman, E.-M., "Truncated Geometry on the Circle"** (published as Hekkelman, J. Geom. Phys. 2023; arXiv preprint).
- **Hekkelman, E.-M., McDonald, E. A., "A noncommutative integral on spectrally truncated spectral triples, and a link with quantum ergodicity," J. Funct. Anal. 289 (2025); [arXiv:2412.00628]** — proposes a simple approximation of the noncommutative integral in the Connes–vS paradigm, with a Szegő limit formula and a quantum ergodicity definition for compact spectral triples.
- **McDonald, E. A., et al. (arXiv:2506.21950, 2025), "Trace Formulas in Noncommutative Geometry"** — extends Hekkelman's Connes–vS-paradigm approximation.
- **van Suijlekom, W. D., "A generalization of K-theory to operator systems," arXiv:2409.02773 (Sep 2024)** — generalizes K-theory to truncated spaces.
- **Leimbach, M., van Suijlekom, W. D., "Gromov–Hausdorff convergence of spectral truncations for low-dimensional tori," Adv. Math. 439 (2024), 109496** — the workhorse torus theorem, the model that GeoVac's Paper 38 transports to S³.
- **Connes, A., van Suijlekom, W. D., "Tolerance relations and operator systems," Indagationes Math. 34 (2023), 606–621**.
- **Connes, A., van Suijlekom, W. D., "Spectral truncations in noncommutative geometry and operator systems," CMP 383 (2021), arXiv:2004.14115** — the originating paper.

### What this gives W2 and what it does not

The Hekkelman–McDonald 2024b NC integral approximation (arXiv:2412.00628) is the cleanest tool for cutoff-dependent observables in the Connes–vS paradigm: it converts the noncommutative integral into a sum that depends on the cutoff $\Lambda$ in a well-controlled way, and connects this to Widom's theory of Szegő-style asymptotics in pseudodifferential operator theory. **However, it is single-cutoff.** It provides a clean asymptotic expansion as one cutoff $\Lambda \to \infty$ but does not provide a structure for *composing two cutoffs* (which is what RG / EFT matching needs for W2).

For the W2 wall (multi-loop QED renormalization on Dirac-S³), what is wanted is a Berezin–Toeplitz analog of multi-cutoff RG: given cutoffs $\Lambda_1 > \Lambda_2$, a "matching map" $T_{\Lambda_1} \to T_{\Lambda_2}$ that absorbs the difference into a counterterm structure. **This does not appear in any published work I located.** The published Connes–vS / Hekkelman–McDonald program is uni-directional UV cutoff; renormalization-style matching at multiple cutoffs is open.

**Honest assessment:** the available machinery handles the *truncation* end of W2 (Sprint TS-E1 / Paper 32 §VIII master Mellin engine sits naturally inside this framework) but does not address counterterm generation. The most concrete near-term step would be to verify whether the Hekkelman–McDonald NC integral on Dirac-S³ at finite $n_{\max}$ reproduces GeoVac's already-computed Sprint TS-E1 results for $\pi$-source observables; if it does, GeoVac is "Hekkelman–McDonald-compatible," which would tighten Paper 32's structural alignment claim. This is a sprint-scale check (~2 weeks).

**Speculative class: Moderate** — strong for cataloguing GeoVac inside the published Connes–vS/Hekkelman–McDonald ecosystem; weak for providing renormalization machinery for W2.

---

## 4. Higher Racah–Wigner recoupling (9j, 12j, 3nj symbols)  *(targets W1 only via discrete-label coupling)*

### State of the art

Standard references: Yutsis, Levinson & Vanagas (1962) *Mathematical Apparatus of the Theory of Angular Momentum*; **Varshalovich, D. A., Moskalev, A. N., Khersonskii, V. K., (1988) *Quantum Theory of Angular Momentum*** (World Scientific). The computational core is **Johansson, H. T. & Forssén, C., "Fast and accurate evaluation of Wigner 3j, 6j, and 9j symbols using prime factorization and multi-word integer arithmetic," SIAM J. Sci. Comput. 38, A376 (2016) [arXiv:1504.08329]**, implemented in `wigxjpf` and `fastwigxj` libraries. SymPy's `sympy.physics.wigner` covers 3j, 6j, 9j, Gaunt, Racah; 12j is sparsely supported. Recent: Bonatsos & Gnezdilov (2025), "Coupling and recoupling coefficients for Wigner's U(4) supermultiplet symmetry," EPJ Plus 140 (2025); semiclassical 12j: Yu (2011, arXiv:1104.3275) and the Anderson–Aquilanti school.

**N-particle "focal-length recoupling" formulation.** Spin-network 3nj symbols are the standard machinery, but their use in atomic physics is restricted to recoupling *quantum-number labels* at a single energy scale (LS↔jj transformation, intermediate-coupling). I was unable to find a published formulation that uses 3nj symbols to recouple states at *different* focal lengths $\lambda_e, \lambda_N$ — i.e. a generalization in which the 3j or 6j coefficients carry an additional scale label that breaks the orthogonality structure. Closest: nuclear shell model uses Talmi–Moshinsky brackets to transform between center-of-mass and relative coordinates with mass-dependent kinematics, which is *exactly* the GeoVac Track NF / Paper 23 mechanism; this transforms two HO bases at the same $\hbar\omega$ but different reduced masses — not multi-λ in the Coulomb sense.

**Honest assessment:** higher 3nj symbols are mature as a *single-scale* recoupling tool. They offer no published bridge across focal lengths. They are useful inside GeoVac for the angular sector (already used in Tracks NE/NF, Paper 22 angular sparsity, Paper 14 §III) but they don't directly attack W1.

**Speculative class: Weak** — well-developed but doesn't address the wall.

---

## 5. Operator-product expansion and scale composition in CFT  *(targets W2)*

### State of the art

The OPE is mature in Lorentzian/Euclidean flat-space CFT (Polchinski, Mack, Rychkov; modern bootstrap revival of Rattazzi–Rychkov–Tonni–Vichi 2008+). Recent reformulations in non-CFT or curved settings:
- **Penedones, J., et al., "Renormalization group flows in AdS and the bootstrap program," (review/work) arXiv:2305.11209 (2023)** — sum rules from two-point functions, RG flow tracked via varying the AdS radius.
- Heat-kernel side: **Vassilevich, D. V., "Heat kernel expansion: user's manual," Phys. Rep. 388 (2003), 279 [hep-th/0306138]**; **Codello, A. & Zanusso, O., "Heat kernel coefficients for massive gravity," J. Math. Phys. 65, 082301 (2024)**.
- **Marcolli, M., "Spectral action gravity and cosmological models," CRP Caltech (2014)**; rationality of spectral action for Robertson–Walker (Marcolli & van Suijlekom 2014, JHEP 12 (2014) 064).

### What is and is not relevant

OPE provides a clean structural template — composition of two operators at nearby scales as a sum of local operators with scale-dependent coefficients. A direct OPE on the Fock-projected S³ graph is not literally available because OPE requires a continuous Lorentzian neighborhood structure that the graph lacks. However, the *Mellin-transform structure* (Paper 18 §III.7 master Mellin engine; Sprint TS-E1) is a discrete cousin of OPE: $\mathrm{Tr}(D^k e^{-tD^2})$ is the small-$t$ asymptotic of a heat kernel that *would* have an OPE if a flat continuum limit were taken. The natural test is to compare GeoVac's master Mellin output against the published heat-kernel coefficients on $S^3$ (Vassilevich 2003, Camporesi–Higuchi 1996) — this is partially done and corresponds to the Paper 36 LS-1..LS-7 sub-percent Lamb shift closure.

For **multi-scale operator-system composition** in the sense W2 needs (an analog of OPE-style $\Lambda_1 \to \Lambda_2$ matching with explicit Wilson coefficients), I found no published theorem. The closest existing thing is the *spectral action* itself: Chamseddine–Connes 1996 [hep-th/9606001] writes a single regularized trace $\mathrm{Tr}\,f(D/\Lambda)$, and the RG flow is added by hand (Marcolli–vS 2014; recent reviews in Connes' *Noncommutative Geometry, Quantum Fields and Motives*). The *rationality of spectral action for Robertson–Walker metrics* (Marcolli–vS, JHEP 12 (2014) 064) is interesting because it shows the spectral action expansion has rational coefficients in a setting structurally similar to S³ — which echoes Paper 18 §IV's "composition" tier.

**Honest assessment:** OPE itself is not directly applicable — GeoVac is not a CFT and does not have an obvious continuum dual (Paper 4 archived 2026-05-03, this is a deliberate scoping decision). The genuinely interesting connection is between GeoVac's master Mellin engine and the Marcolli–vS rationality theorem: both observe that spectral-action-like coefficients on a homogeneous compact space sit in clean rational/transcendental rings. **A targeted comparison sprint between Sprint TS-E1's case-exhaustion theorem and Marcolli–vS rationality is a concrete sprint-scale test** (~3 weeks).

**Speculative class: Moderate** — reframes existing GeoVac results inside the Marcolli–vS lineage; does NOT close W2.

---

## 6. Sub-questions while searching

### 6a. 2024–2026 work on discrete spectral triples and multi-scale composition

The Connes–van Suijlekom 2021 framework remains the dominant published thread. Significant 2024–2026 advances:
- **van Suijlekom 2024 [arXiv:2409.02773]** generalizes K-theory to operator systems — opens the door to topological invariants at finite truncation (potentially relevant to the Sprint TS-E3 discrete $c_1$ on Hopf S³).
- **Hekkelman–McDonald 2024 [arXiv:2412.00628]** NC integral approximation in the Connes–vS paradigm — directly applicable to GeoVac.
- **Leimbach–van Suijlekom 2024** — torus GH convergence; transported to S³ in Loutey 2026 (Paper 38).
- **Latrémolière 2026 [arXiv:2603.19128]** spectral continuity for almost-commutative manifolds in the C¹ topology — directly relevant to Sprint H1 stability under varying Yukawa.

**No published paper I found directly attempts a two-cutoff (multi-resolution) composition theorem in the Connes–vS framework.** This is the W2 frontier, and it appears genuinely open in the published literature.

### 6b. Recent QED / atomic-physics on recoil + hyperfine + finite-size composition

- **Pachucki, K., Patkóš, V., Yerokhin, V. A., "Recoil Corrections to the Energy Levels of Hydrogenic Atoms," PRL 130, 023004 (2023)** — completes pure-recoil corrections of order $(Z\alpha)^6$ to Coulombic bound states without approximation in the particle masses. *This is the closest published thing to a Hamiltonian-level (rather than perturbative-Feynman-graph) recoil treatment.*
- **Eides, M. I., Shelyuto, V. A., Tomalak, O., et al., "New spin structure constraints on hyperfine splitting and proton Zemach radius," Phys. Lett. B (2024) [doi:10.1016/j.physletb.2024.139049]** — improves proton polarizability uncertainty in hyperfine splitting.
- The standard reference book **Eides, M. I., Grotch, H., Shelyuto, V. A., "Theory of Light Hydrogenic Bound States" (Springer, 2007; updates in 2019)** and its "Theory of Light Hydrogenlike Atoms" review [hep-ph/0002158] remain canonical for the layered-EFT decomposition (Coulomb / Breit / NRQED / pNRQED).
- **Effective field theory for muon conversion**, JHEP 05 (2025) 171 — combines HQET, NRQED, pNRQED, SCET, boosted HQET into a single matched chain. This is a published *production-grade* multi-scale composition pipeline; it works *because* every matching step is at a single physical scale and uses standard $\overline{MS}$ counterterms. **This is the structural archetype W2 lacks in GeoVac.**

The Pachucki–Patkóš–Yerokhin 2023 result is significant: their non-perturbative-in-mass-ratio Hamiltonian recoil structure is, in spirit, the W1 wall's correct address. Their construction goes through a Foldy–Wouthuysen reduction to an effective two-particle Hamiltonian where the nucleus is treated quantum-mechanically — exactly the cross-register coupling W1 names. **Porting Pachucki–Patkóš–Yerokhin 2023's Hamiltonian construction onto Dirac-S³ is a concrete, sprint-scale W1 attack** (~6–10 weeks).

### 6c. Bargmann transform as focal-length bridge between Coulomb and HO sectors

The Coulomb–HO duality via the Kustaanheimo–Stiefel (KS) transformation is well-established (KS 1965; Duru–Kleinert 1979 for path integrals; Cornish 1984). The KS map is a 4D-3D dimensional bridge, and the Bargmann transform sits naturally in the 4D oscillator setting. This is precisely Paper 24's Bargmann–Segal lattice mechanism: the 3D HO discretizes on $S^5$ Hardy space via $SU(3)$ $(N,0)$ irreps, and the Coulomb on $S^3$ via Fock projection — but the *two are not joined into a single bigger spectral triple* in any published work. **The Coulomb/HO asymmetry of Paper 24 §V already names this as a 4-layer structural asymmetry.**

I found **no published paper that frames the Bargmann transform as a literal focal-length bridge in the GeoVac sense** (i.e. as a unitary conjugation between two spectral triples whose Dirac operators have different scaling laws). This appears to be a genuine GeoVac-internal observation. The closest published thing is Andrianov & Cannata's "PT-symmetric Coulomb–HO correspondence" (Phys. Lett. A 2000), which is a unitary statement on a single Hilbert space, not a triple-level statement.

### 6d. Spectral triples of multi-particle quantum mechanics beyond the SM paper

The Connes–Chamseddine SM construction (Connes 1996; Chamseddine–Connes 1997; Chamseddine–Connes–Marcolli 2007) is the canonical "almost-commutative" multi-sector spectral triple, but it treats the *internal* (matrix) sector as the multi-particle structure, not real-space multi-particle quantum mechanics. Type III σ-spectral triples (Marcolli, Caltech notes) treat statistical mechanical systems but are infinite-temperature objects, not literal multi-particle wavefunctions. **Spectral noncommutative geometry, standard model and all that** (arXiv:1906.09583, review) catalogues the standard tensor-product almost-commutative ecosystem.

A **genuine spectral triple of $N$-particle non-relativistic QM** (e.g. of a helium atom or a deuterium atom) does **not** appear to have been written down in the literature in the form I expected. The Connes ACG is the closest analog, but it is for internal SM degrees of freedom, not for spatial multi-particle. **GeoVac's Track NI deuterium PoC is, as far as I can tell, the first explicit construction of a spectral-triple-style multi-particle QM Hamiltonian in the Connes lineage** (i.e. Hilbert space tensor product, Dirac operator with cross-register coupling, finite operator system at each truncation).

This is a genuine gap-in-the-published-literature finding.

---

## Top 3 candidates to test first

### Candidate 1 — Pachucki–Patkóš–Yerokhin 2023 Hamiltonian recoil on Dirac-S³  *(W1)*
**Mechanism:** the PRL 130, 023004 (2023) construction gives an effective two-particle Hamiltonian, exact in mass ratio at $(Z\alpha)^6$, where the nucleus is treated quantum-mechanically (not classically as in Track NI). Porting this to Dirac-S³ would mean: (a) define a two-register spectral triple $T_e \otimes T_N$ where $T_N$ has its own focal length $\lambda_N \neq \lambda_e$, (b) implement the Pachucki–Patkóš–Yerokhin reduction at the qubit-Hamiltonian level on the cross-register Pauli strings, (c) verify that the resulting recoil shifts on hydrogen 1S match published data. This directly attacks W1 and is *not* an arbitrary construction — it is calibrated against an actively-validated atomic-physics literature.
**Sprint scope:** 6–10 weeks. Risk: the spectral triple structure may force a different decomposition than Pachucki et al.'s, in which case the discrepancy itself is informative.

### Candidate 2 — Hekkelman–McDonald NC integral on Dirac-S³ at finite $n_{\max}$  *(W2 framing, not W2 closure)*
**Mechanism:** arXiv:2412.00628 gives a clean cutoff-dependent NC integral approximation that fits the Connes–vS paradigm GeoVac already lives in. Sprint test: verify whether the Hekkelman–McDonald approximation, applied to the Camporesi–Higuchi $D$, reproduces Sprint TS-E1's case-exhaustion classification of $\pi$-sources. If yes, GeoVac inherits a published asymptotic theory and the master Mellin engine becomes a special case of Hekkelman–McDonald. If no, the discrepancy is itself a structural finding.
**Sprint scope:** 2–4 weeks. This does not close W2 — it puts GeoVac inside a published asymptotic theory and helps frame what a counterterm structure would look like.

### Candidate 3 — Multi-λ Sturmian three-body integral on hydrogen with quantum nucleus  *(W1 in atomic-physics language)*
**Mechanism:** explicitly construct $\int Y_{n_e l_e m_e}^{(\lambda_e)}(\vec{r}_e)^* Y_{n_N l_N m_N}^{(\lambda_N)}(\vec{R}_N)^* V_{eN}(|\vec{r}_e - \vec{R}_N|) Y_{n'_e l'_e m'_e}^{(\lambda_e)}(\vec{r}_e) Y_{n'_N l'_N m'_N}^{(\lambda_N)}(\vec{R}_N) \, d^3r_e \, d^3R_N$ for $\lambda_e \neq \lambda_N$, then check whether the Shibuya–Wulfman multipole expansion still terminates at $L_{\max} = 2 \max(l_e, l_N)$ as in the matched case. If yes, this is the algebraic backbone of cross-register two-body coordinate operators and W1 closes structurally. If no (most likely outcome — the multipole tail will probably not terminate cleanly across mismatched exponents), the failure mode itself characterizes W1 in algebraic terms.
**Sprint scope:** 4–6 weeks. This is the GeoVac-internal analog of Avery's matched-λ Shibuya–Wulfman work but at multi-λ; it is a derivable mathematical question.

---

## Surprises

**S1. The published "many-particle Sturmian" school does NOT do multi-λ.** I expected to find Avery (1997) defining $N$-particle Sturmians at $N$ different exponents. It does not — the generalized Sturmian basis is by construction isoenergetic at a single energy, and the weights $\beta_\nu$ are weights on $V_0$, not exponents. The Avery school's "different exponents" are *between basis members* of a single particle, not *across particles*. **Implication: the GeoVac literature claim "Sturmian work is well-developed" is correct but covers a narrower domain than I would have guessed — it does not address W1 directly.**

**S2. Tensor-product convergence theorems for two infinite metric spectral triples appear to be open.** The almost-commutative case (one infinite + one finite) is now closed (Latrémolière 2026 [arXiv:2603.19128]), but two infinite factors is not in the published literature. **Implication: the G4b cross-manifold blocker (Paper 32 §VIII.C) is not a GeoVac-specific limitation — it is an open question in NCG with no published answer**. Paper 38's S³ result + a parallel S⁵ result + a tensor-product theorem would be a substantial NCG contribution.

**S3. The closest published thing to a "spectral triple of an atom" is GeoVac itself.** I expected to find Connes-style spectral-triple constructions for non-relativistic helium or hydrogen with quantum nucleus. They do not appear to exist in the form GeoVac has constructed (Track NI deuterium PoC, composed nuclear-electronic at 26 qubits). The Connes SM almost-commutative geometry treats *internal* gauge degrees of freedom, not real-space multi-particle. **Implication: GeoVac's nuclear-electronic composed sector is more original than the "we live in the Marcolli–vS lineage" framing suggests. Paper 23 / Track NI may be publishable as a stand-alone contribution (analog of Paper 38) once the multi-focal compositionwall is partially addressed.**

**S4. Pachucki–Patkóš–Yerokhin 2023 is the published-physics-side answer to W1.** Their recoil Hamiltonian, exact in mass ratio at $(Z\alpha)^6$, is a Foldy–Wouthuysen-style effective Hamiltonian with quantum-mechanical nucleus. **GeoVac has no native version of this** — the Track NI cross-register couples spin labels but evaluates spatial operators at the classical proton position. Implication: the multi-focal wall has a *concrete, calibrated, atomic-physics target* against which a GeoVac extension can be tested. This is the strongest single sprint candidate.

**S5. Hekkelman–McDonald + Marcolli–vS rationality together form a tighter framing for Paper 32.** The Marcolli–vS rationality of spectral action for Robertson–Walker metrics (JHEP 12 (2014) 064) shows clean rational/$\pi$ structure in spectral-action coefficients on homogeneous compact backgrounds. Combined with Hekkelman–McDonald 2024, this nearly is the Sprint TS-E1 master Mellin engine. **If Paper 32 §VIII cited both, the case-exhaustion theorem would land inside an already-published ecosystem rather than as a stand-alone GeoVac claim.** Recommend adding both citations next time §VIII is touched.

**S6. The multi-loop QED W2 wall has a published archetype: matched EFT chains.** JHEP 05 (2025) 171 (muon conversion EFT chain) shows what a *production* multi-scale composition pipeline looks like. It is structurally everything W2 lacks: a chain of single-scale matchings with $\overline{MS}$ counterterms generated by hand at each step. **The honest scope statement is therefore: W2 is not a NCG-specific wall — it is the standard QFT renormalization story translated into the spectral-action setting, and no current NCG framework solves it autonomously.** GeoVac's "framework reproduces UV-divergent integrand but cannot generate counterterms" is structurally identical to the position Connes–Chamseddine 1996 was in, and Marcolli–vS 2014 partially addresses but does not close. **W2 may be a permanent open question of the spectral-action program, not just of GeoVac.**

---

## Honest scope and uncertainties

**What I did not have time to cover:**
- Detailed reading of Pachucki–Patkóš–Yerokhin 2023 for the explicit Hamiltonian form (only the abstract and PRL summary). Candidate 1 sprint scope is therefore conservative.
- A serious survey of Lattice QCD's multi-scale composition machinery (Symanzik improvement, gradient flow, etc.) which is the most mature *discrete* W2 ecosystem in physics. Lattice QCD has decades of published work on counterterm generation in discrete settings; some of it could be ported to GeoVac. *Recommend a follow-up sprint specifically on Symanzik improvement and gradient flow as templates for GeoVac counterterms.*
- Mick Gielen + van Suijlekom (2023) "Operator systems for tolerance relations on finite sets" *Indagationes Math.* 34, 606. Could be relevant to the discrete-truncation side of W2 but I did not read it.
- van Suijlekom 2024 [arXiv:2409.02773] generalization of K-theory to operator systems — quickly catalogued, not read in depth.
- Fermi liquid + Schrieffer–Wolff transformations as multi-scale composition templates (Sci. Rep. 15, 29243 (2025) on closed-form spin-relativistic corrections). Could be a structural archetype for cross-register coupling; not pursued.
- Quantum-graph spectral truncation work (de Jong, Gielen) post-2024.

**Uncertainties:**
- I cannot fully verify whether Latrémolière's 2018 propinquity paper [arXiv:1811.10843] does or does not contain an external-product theorem without the full text. The abstract and 2019 follow-ups do not headline products.
- The Avery school's most recent work (2018+) is largely in books and Few-Body Systems proceedings that may contain multi-λ extensions I missed.
- Hekkelman–McDonald's NC integral approximation may not extend cleanly to the Camporesi–Higuchi spinor sector — the published result is for scalar Laplacian / single Dirac, and the Dirac sector at half-integer Hurwitz (Paper 28 / Sprint MR-A) may behave differently. This is exactly Candidate 2's open question.

**What I am highly confident about:**
- The published "many-particle Sturmian" school does not provide multi-λ basis sets.
- Tensor products of two infinite metric spectral triples are open in the published NCG literature.
- Pachucki–Patkóš–Yerokhin 2023 is the best-calibrated atomic-physics target for W1.
- W2 is structurally the standard renormalization story and is unlikely to close inside the spectral-action paradigm without importing field-theoretic counterterm machinery from outside.
