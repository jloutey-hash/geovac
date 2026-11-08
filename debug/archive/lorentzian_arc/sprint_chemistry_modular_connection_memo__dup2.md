# Sprint Chemistry-Modular Connection — diagnostic memo

**Date:** 2026-05-25
**Sprint:** Diagnostic-only test of whether the modular / spectral-triple / master-Mellin-engine machinery (Papers 32, 35, 38, 42, 43, 45, 46, 47, 48, 49) has anything substantive to say about chemistry correlation-energy walls, specifically W1e (NaH inner-region overattraction wall named 2026-05-23).
**Mode:** Diagnostic-before-engineering. No code, paper, or memory modifications.
**Decision gate:** (i) clean negative / (ii) partial positive / (iii) substantive positive.
**Cross-references:** CLAUDE.md §1.7 WH1; §2 Sprint TS, Sprint MR-A/B/C, Paper 35 Prediction 1, Paper 49 OSLPLS, Sprint W1c F3 closure, Sprint W1e Schmidt + core correlation closeout; memory/{mellin_taxonomy_engine, mellin_engine_domain_partition, multi_focal_wall_pattern, geovac_structural_skeleton_scope_pattern, pythagorean_orthogonality, modular_flow_regularization_negative, paper27_entropy_thesis, sprint_modular_alpha_arc}; Paper 18 §III.7; Paper 34 §III.28; Paper 35 §VI; Paper 42; Paper 27; Paper 19 §sec:w1c_residual + §sec:w1e_refinement_f4_f6.

---

## §0. Executive summary

The chemistry W1e wall (NaH inner-region overattraction at max_n=2, currently ~3.9 Ha PES well depth vs 0.072 Ha experimental, after F4+F5+F6+Schmidt+core-correlation ruled out at sprint-scale) lives in a region of the framework that the modular / spectral-triple machinery developed in Papers 32, 38–47, and 48–49 does NOT structurally reach.

Across the five probes below, the cumulative verdict is **(i) clean negative with one structural sharpening**:

- **Probes 1, 2, 3, 4 are clean negatives.** The modular flow / Tomita-Takesaki machinery is signature-Lorentzian and KMS-thermal; chemistry FCI ground states are pure, non-thermal, time-independent, and have no boost generator. The master Mellin engine's three sub-mechanisms M1/M2/M3 are all spectral-side and transcendental-class; chemistry binding errors are state-side and dimensional. Paper 35 Prediction 1's "π enters iff continuous integration over temporal/spectral parameter" gives no chemistry foothold because the W1e wall is at the eigenvalue level, not at any continuous integration. Paper 42's four-witness theorem requires a wedge KMS state with rapidity-rate generator, none of which chemistry FCI has.
- **Probe 5 is the one structural sharpening.** The "structural-skeleton scope" pattern — already observed (CLAUDE.md §1.7) — extends cleanly across both the QED-side LS-8a wall and the chemistry-side W1e wall via the multi-focal-composition wall pattern. The cross-domain extension is not new content; it was observed at the W1c full arc F3 closure (2026-05-23) and reaffirmed at W1e Schmidt + core-correlation closeout. The substantive content is the **state-side / spectral-side disjointness theorem** (Sprint TD Track 5, PSLQ negative): the framework's state-side observables (correlation entropy, mutual information, Wasserstein, Uhlmann deficit) are categorically distinct from the master Mellin engine ring. This is empirically established and applied here as the structural boundary statement for the modular/chemistry question.

The one substantive new content of this diagnostic, beyond reaffirming the existing boundary statements, is the **inverse direction**: Paper 49's twin-paradox-as-quantum-information identification (Uhlmann monotonicity bounds reverse-triangle deficits along KMS-state stacks) is the **operator-algebraic dual** of a structurally analogous chemistry phenomenon — the FCI correlation energy being lower-bounded by the entropy of the natural-orbital occupation deviation from {0, 1}. This is a structural correspondence (chemistry analog of TICI), not a quantitative reading. It does not close W1e but it sharpens *why* W1e is not closable by single-particle or mean-field approximations: those are zero-entropy density-matrix classes (idempotent 1-RDM), and the correlation energy deficit is bounded below by an information-theoretic quantity that vanishes precisely on those classes.

**Final verdict classification: (i) clean negative on transfer, with (ii)-grade partial-positive structural reading on the entropy-bound side that is consistent with but does not add to the existing structural-skeleton-scope statement.**

The value of this diagnostic is the boundary crystallization — we now know not to chase modular machinery as a chemistry tool, and the entropy-bound reading clarifies *why* the W1e wall (and the broader multi-focal-composition wall) is intrinsic to the FCI correlation-energy structure, not contingent on the framework's specific architecture.

---

## §1. Probe 1 — Uhlmann-relative-entropy reading of correlation energy

**Question.** Paper 49 identifies the operator-algebraic twin paradox: along thermal-time stacks across distinct KMS states, the reverse-triangle deficit on OSLPLS is bounded by Uhlmann's relative-entropy monotonicity. The Triple-Intersection Cocycle Identity (TICI) + Uhlmann's inequality give strict super-additivity. The chemistry side has Paper 27's entropy structure on FCI 1-RDMs (one-body entanglement-inert + V_ee generates all entropy + universal $S_B = A(\tilde{w}_B/\delta_B)^\gamma$ scaling with $\gamma_\infty \approx 1.96$ matching 2nd-order Rayleigh-Schrödinger). Is FCI correlation energy similarly bounded by an Uhlmann deficit?

### 1.1 Literature anchor for correlation-energy-as-entropy

The question has a clean literature anchor that I should record honestly. The information-theoretic lower bound on the multi-electron correlation energy via density-matrix entropy is **a real published direction** with multiple anchors:

- **Levy 1979, Lieb 1983** (universal density functional construction): the exact $F[\rho]$ is convex, and the correlation contribution $E_c[\rho] = F[\rho] - T_s[\rho] - J[\rho]$ has well-defined sign on the density manifold.
- **Mazziotti** (Phys. Rev. Lett. 108, 263002, 2012; Phys. Rev. A 85, 062507, 2012): variational 2-RDM theory with N-representability constraints; the correlation energy is bounded below by the gap between the variational 2-RDM ground state and the closest single-determinant.
- **Esquivel, Karwowski, Liu, et al.**: a now-large series of papers (mid-2010s onward) on the **complexity / entropy of electron correlation** measured via the single-particle 1-RDM (Shannon entropy of natural-orbital occupations) and the two-particle cumulant. The relevant statement: *correlation energy is bounded below by an entropic quantity that vanishes exactly when the 1-RDM is idempotent (i.e., a single Slater determinant)*.
- **Tropp / Kollmar / Heßelmann** (correlation-as-entanglement): the two-particle reduced density matrix's von Neumann entropy in the natural-orbital basis is a known correlation quantifier.

None of these published bounds is *exactly* the Paper 49 twin-paradox-style Uhlmann deficit between two KMS states. The chemistry analog is Uhlmann monotonicity on the manifold of **N-representable density matrices** (or 2-RDMs) rather than across the modular flow of KMS states. The two scenarios differ structurally:

| | Paper 49 (Uhlmann/twin) | Chemistry (correlation/entropy bound) |
|:---|:---|:---|
| State family | KMS states $\omega_\beta$ on $\mathcal{A}$ across wedge boost orbits | N-representable density matrices on $\mathcal{H}^N$ at fixed $N$, ground state of $H$ |
| Flow | Modular flow $\sigma_t^\omega$, Connes 1973 cocycle | None (single ground state; or RG-style flow on $E_c[\rho]$) |
| Triangle inequality | TICI thermal-time stack, reverse-triangle super-additivity | Subadditivity of von Neumann entropy under 1-RDM / 2-RDM trace |
| Bound type | Uhlmann's relative-entropy monotonicity | Lieb-Robinson, N-representability, $f$-divergence |
| Variable | Thermal time $\tau_\mathrm{mod}^\omega$ | Correlation energy $E_c$ vs an entropic functional of 1-RDM |

The two scenarios are **structurally analogous but not the same theorem**. Paper 49's TICI is specifically thermal-time-cocycle additivity across distinct KMS states; the chemistry analog is N-representability + Lieb-style convexity of $F[\rho]$.

### 1.2 Does Paper 49 give a quantitative bound for W1e?

No. Paper 49's framework operates on the manifold of KMS states of a C*-algebra and on the Tomita-Takesaki modular flow, neither of which exists for the FCI ground state in a non-thermal, time-independent chemistry problem. The chemistry side has no boost generator, no temperature, no wedge, no modular cocycle. Attempting to import Paper 49's framework would require either:

- Treating $e^{-\beta H}/Z$ at finite $\beta$ as a KMS state of the FCI Hamiltonian (the **apparatus identity**, Paper 34 §III.28). This is well-defined — but the ground state limit $\beta \to \infty$ destroys the modular flow's nontriviality.
- Identifying a non-thermal modular structure on the FCI ground state (e.g., a Wright-Ouwehand modular operator on a non-tracial state). This is an open mathematical question and does not have a published chemistry application.

**Probe 1 verdict: clean negative on quantitative transfer; structural-analog observation only.**

The structural analog observation (correlation energy lower-bounded by an entropic functional that vanishes on the single-Slater-determinant submanifold) is **already in the chemistry literature** and is **already aligned with Paper 27's structural finding** (single-determinant / one-body operators are entanglement-inert; V_ee generates all correlation entropy). This is consistent with W1e being structurally distinct from the F4/F5 mean-field-class closures, which all operate on density matrices in the zero-entropy class (idempotent 1-RDMs).

### 1.3 Why W1e specifically isn't closable by PK or Hartree (entropic reading)

Paper 27 §sec:cusp finding: in the eigenbasis of the one-body operator $H_1$, $V_{ee}$ is overwhelmingly diagonal (Frobenius ratio 0.892 at $n_\text{max}=4$), but the **relative commutator** $\|[H_1, V_{ee}]\|_F / \|H_1\|_F = 5.3\%$ at $n_\text{max}=4$ does **not vanish with basis size**; it saturates. The off-diagonal residual is concentrated on the (1s,1s) pair-state (cusp hot node).

This is the structural reason W1e is not closable by mean-field. The off-diagonal residual generates the correlation entropy ($\sim 10$–$39\%$ of correlation energy is carried by it per Paper 26). PK (rank-1 single-particle) and Hartree (mean-field $J - K$) project onto the zero-entropy submanifold. The Uhlmann-class lower bound on correlation energy (Esquivel-Karwowski-Liu style entropy bounds, NOT Paper 49's KMS-state Uhlmann deficit) vanishes by construction on that submanifold. So the F4/F5 architectural closures *cannot* close W1e: they live on the wrong submanifold of the N-representable density matrices.

This is a structural reading that the chemistry literature already supports, and that Paper 27's framing already covers. It's not new content — but it explains *why* the F4 (43% saturated) / F5 (25.7% ceiling) / Schmidt (0.2% closure) results all saturated below the wall. Each closure target was a zero-entropy projection of the FCI correlation channel.

**This sharpening is at the (ii) partial-positive level for *structural understanding*** of why prior closures didn't work, but it is **not load-bearing as a quantitative tool**, and it does not predict a new sprint-scale closure mechanism for W1e.

---

## §2. Probe 2 — Master Mellin engine on chemistry observables

**Question.** The master Mellin engine $\mathcal{M}[\mathrm{Tr}(D^k e^{-tD^2})]$ at $k \in \{0, 1, 2\}$ classifies transcendental signatures M1/M2/M3 (Hopf-base measure / Seeley-DeWitt / vertex-parity Hurwitz). Does it have any transcendental signature in chemistry?

### 2.1 Chemistry observables are dimensional, not transcendental

The framework's chemistry observables (binding energy $D_e$ in Hartree, equilibrium bond length $R_\text{eq}$ in bohr, ionization potential, polarizability) are dimensional quantities. They sit in Layer 2 of Paper 34 with variables (mass, length, energy) introduced by specific projections (rest-mass §III.14, multipole expansion §III.21, adiabatic §III.24).

There is no continuous integration over a temporal or spectral parameter in computing the FCI ground-state energy. The FCI matrix is built from exact-rational Slater integrals on the discrete S^3 graph (per CLAUDE.md §12 Algebraic Registry); diagonalized to machine precision; the ground-state eigenvalue is a finite-precision real number in Ha. There is no Mellin transform anywhere.

**Per Paper 35 Prediction 1, chemistry binding energies should be π-free at the eigenvalue level.** This is verified empirically: graph-native CI matrix elements are pure rational (Paper 7 §VI), the FCI matrix is rational, and the eigenvalues are algebraic in the appropriate sense. The W1e wall is at the eigenvalue level — and **the eigenvalues are not in any transcendental ring**, they are real algebraic numbers (irrational because diagonalization in higher-than-quartic-degree polynomial spaces is needed, but not in any structured transcendental class).

### 2.2 But π *does* appear in the cross-block ERI kernel...

The cross-block V_ne kernel ($1/|\mathbf{r} - \mathbf{R}_B|$ Coulomb singularity) is treated via multipole expansion with Gaunt-coefficient termination at $L_\text{max} = 2 l_\text{max}$ (Roothaan, §III.21). The Gaunt integrals contain $\pi$ via the spherical-harmonic normalization, and the spatial integration over $r$ generates exponential integrals $E_1(a)$ etc. — i.e., the *embedding* of the discrete Coulomb graph onto the metric continuum introduces embedding-class transcendentals.

These transcendentals exist (the 1/r₁₂ embedding exchange constant per Paper 18), but they are **integrated out** into rational R^k integrals in the framework's algebraic Slater-integral implementation. By the time the FCI matrix is built, every π has been absorbed into the rational matrix elements. The eigenvalues are π-free.

This is exactly the Paper 35 framing: there's no continuous integration over a temporal/spectral parameter in computing FCI eigenvalues; therefore no π appears at the eigenvalue level; therefore the master Mellin engine has no transcendental signature to attach to W1e.

**Probe 2 verdict: clean negative.** The master Mellin engine $M_1 \cup M_2 \cup M_3$ is structurally tied to the spectral-side trace formula at integer operator orders $k = 0, 1, 2$. Chemistry observables (and W1e specifically) live in the algebraic-eigenvalue regime — no Mellin transform, no transcendental ring, no $k$-domain to occupy.

### 2.3 Could W1e have a projection-mechanism signature?

W1e is the residual after F1-F6 + Schmidt + core correlation. The mechanism is documented as "deep multi-determinant FCI correlation that resists basis-truncation, single-particle PK, and mean-field Hartree-class closure at scales reachable within sprint-scale compute" (Paper 19 §sec:w1e_refinement_f4_f6).

Per Paper 34's projection dictionary, this maps onto **multiple projections operating simultaneously**:
- §III.20 Phillips-Kleinman / core-valence orthogonality (the H-Na orthogonality piece)
- §III.21 Multipole expansion / Gaunt termination (the cross-block V_ne kernel)
- §III.7 Camporesi-Higuchi spinor lift (the natural-orbital basis structure)
- §III.24 Born-Oppenheimer adiabatic separation (the implicit frozen-nuclear approximation)

W1e is therefore in the **multi-focal-composition wall regime** of Paper 34. It is NOT a single-projection mechanism that the engine could classify; it's a residual that depends on multiple projections composing in a way the framework doesn't handle natively (CLAUDE.md §1.7 multi-focal-composition wall, six instances).

This recovers the existing structural-skeleton-scope statement, not new content.

---

## §3. Probe 3 — Paper 35 Prediction 1 chemistry analog

**Question.** Prediction 1: π enters a GeoVac observable iff its evaluation includes a continuous integration over a temporal/spectral parameter promoted from the discrete graph spectrum. Chemistry has no temporal compactification, but it has several "projection-like" steps. Does Prediction 1 give W1e a Layer-1 / Layer-2 reading?

### 3.1 Where does π enter chemistry?

Three places, all consistent with Prediction 1:

1. **Cross-block ERI kernel / multipole expansion (§III.21).** Spherical-harmonic normalization carries $\pi$ at the Gaunt-coefficient step. This is an integration over the angular-momentum spectrum (effectively over the SO(3) Haar measure), i.e., a continuous integration over a spectral parameter. **Consistent with Prediction 1.**

2. **Adiabatic curve $U(R)$ (§III.24 + §III.25).** The Born-Oppenheimer potential energy curve is the result of diagonalizing $H_\text{elec}$ at each R-point. This is point-by-point algebraic. But evaluating $D_e = U(R_\text{eq}) - U(\infty)$ involves taking a limit, which is mathematically a continuous integration over R. **Marginally consistent with Prediction 1.**

3. **Spectral-action regularization at Mellin level (§III.6).** $\zeta_{D^2}(s)$ derivatives at $s = 0$ involve $\pi$. NOT used in chemistry FCI — used in QED (Paper 28). Not relevant for W1e.

### 3.2 Does W1e have a temporal/spectral integration in its evaluation?

No. W1e is measured by:

```
D_e^F3 = U(R_\text{min}) - U(R = 10\,\text{bohr})
```

where $U(R)$ at each R-point is the FCI ground-state energy at fixed nuclear geometry. Each FCI eigenvalue is a discrete algebraic quantity. The R-grid is discrete (8 points in F3). No continuous integration over a temporal/spectral parameter.

**Per Prediction 1, W1e is therefore π-free at the eigenvalue level.** Empirically confirmed: $D_e^\text{F3} = +4.374$ Ha is a real algebraic number with no clean rational or transcendental signature.

The L1/L2 boundary statement: **W1e lives entirely in Layer 1** (bare combinatorial graph + algebraic FCI on Slater integrals). The Layer-2 projections that *do* operate (cross-block V_ne via §III.21, BO via §III.24) introduce embedding-class transcendentals that get absorbed into the rational matrix elements before the eigenvalue is computed.

**Probe 3 verdict: clean negative.** Prediction 1 is consistent with W1e being π-free, and that's all it tells us. It does not give a structural mechanism for closing W1e — only a confirmation that W1e is not in the framework's "where π enters" sector.

The boundary statement: **the W1e wall is at the algebraic-eigenvalue level, NOT at a transcendental-injection level**. This is the same reason the master Mellin engine has no chemistry signature (Probe 2). Chemistry binding errors live below the framework's transcendental-classification machinery.

---

## §4. Probe 4 — Modular flow / KMS structure on chemistry

**Question.** Paper 42's four-witness theorem is structurally tied to Bisognano-Wichmann / Unruh / signature-Lorentzian content (wedge KMS state, rapidity-rate boost generator, σ_{2π}(O) = O closure). Is there a non-relativistic chemistry analog? Or a Pythagorean HS-orthogonality analog (Paper 43 §10.2 cor:pythagorean_orthogonality)?

### 4.1 Direct application: clean negative

Chemistry has:
- No causal structure / no light cone / no Bisognano-Wichmann wedge
- No boost generator $K_\alpha$ / no rapidity rate
- No KMS state at the FCI ground state (which is pure, not thermal)
- No modular flow $\sigma_t^\omega$ on the ground state (modular flow is trivial on pure states)

The four-witness theorem closure $\sigma_{2\pi}(O) = O$ requires a wedge KMS state with rapidity-rate generator, none of which the chemistry FCI ground state has. **Direct transfer is clean negative.**

This is reinforced by the **modular-flow-regularization-negative** memo (memory/modular_flow_regularization_negative.md, 2026-05-18): the analogous question on the QED side (could the Krein wedge modular flow regulate LS-8a counterterms?) was also a clean negative. The Bostelmann-Cadamuro-Minz 2025 obstruction (modular flow is non-local for interacting theories) is cited as the structural reason why modular flow doesn't function as a regulator. The chemistry analog is even stronger: there's no flow to begin with.

### 4.2 Could a thermal reading of the FCI ground state via $\beta \to \infty$ Gibbs help?

The apparatus identity (Paper 34 §III.28) is well-defined on the FCI Hamiltonian: $\rho_\beta = e^{-\beta H}/Z$ is a Gibbs state of the FCI Hamiltonian at any $\beta > 0$. The Sprint TD Track 2 verification (2026-05-08) confirms the apparatus identity at machine precision (max residual $4.8 \times 10^{-14}$) across H/He/Li⁺ at 30 log-spaced temperatures.

But:
- The $T \to 0$ limit (which is the FCI ground state) makes $\rho_\beta$ a pure projector onto the ground state. Modular flow on this pure state is trivial.
- At any finite $\beta$, the modular flow $\sigma_t^{\omega_\beta}(a) = e^{i\beta^{-1} H t} a e^{-i\beta^{-1} H t}$ exists, but this is just the Heisenberg flow rescaled by $\beta^{-1}$. The σ_{2π} closure becomes $a \mapsto e^{i 2\pi/\beta \cdot H} a e^{-i 2\pi/\beta \cdot H}$ which is *not* identity for generic $a$ (it's identity only if $H$ has spectrum in $\beta \mathbb{Z}/2\pi$, which it doesn't).

So the natural chemistry analog of the four-witness theorem **doesn't close**: there's no integer-spectrum boost generator. The chemistry Hamiltonian has a generic real-valued spectrum, not the integer two_m_j spectrum that makes the Paper 42 BW-α closure bit-exact.

### 4.3 Pythagorean HS-orthogonality analog?

Paper 43 §10.2: on the Krein hemispheric wedge, $\langle H_\text{local}, D_W^L \rangle_\text{HS} = 0$ bit-exact, with closed form $r^2(n; \kappa_g) = \kappa_g^2 S(n)/(4\pi^2) + D(n)$. The $1/\pi^2$ prefactor is the M1 Hopf-base measure signature.

Is there a chemistry analog of the H_local / D_W^L decomposition? Specifically, can FCI Hamiltonian decompose as $H = H_0 + V_\text{corr}$ with $\langle H_0, V_\text{corr} \rangle_\text{HS} = 0$ in some natural inner product?

The Hartree-Fock decomposition $H = H_\text{HF} + V_\text{corr}$ where $V_\text{corr}$ is the fluctuation potential. Under the Frobenius / Hilbert-Schmidt inner product on operators, $\langle H_\text{HF}, V_\text{corr} \rangle$ is **generically nonzero** — there's no chemistry analog of the BW-wedge $\mathbb{Z}_2$ chirality parity argument (Paper 43 thm:pythagorean_formal, established by Sprint Pythagorean Orthogonality 2026-05-23) that would force orthogonality.

The Z₂ wedge-chirality parity grading on the Lorentzian Krein space (Π_W) is a specific structural feature of the four-witness construction that has no chemistry analog: chemistry has no chirality grading, no signature switch, no involution that makes $H_0$ even and $V_\text{corr}$ odd under a natural parity operation.

**Probe 4 verdict: clean negative on direct transfer; clean negative on Pythagorean analog.**

The structural reason: the modular and Pythagorean structures in Papers 42, 43 rest on **signature-Lorentzian content** (wedge, boost, chirality grading) plus **bit-exact integer spectrum of the boost generator**. Neither exists for the FCI Hamiltonian. The chemistry side is structurally non-relativistic, time-independent, and has continuous-spectrum content that the modular machinery doesn't engage.

---

## §5. Probe 5 — Structural-skeleton-scope theorem for chemistry

**Question.** The multi-focal-composition wall pattern (5 physics + 1 chemistry = 6 instances) is cross-domain. Is there a unified theorem-statement of "GeoVac is a structural-skeleton framework" that covers both LS-8a-style physics walls and W1e-style chemistry walls? What would it look like (analog of η-trivialization for inner factor)? Worth crystallizing as Paper 50?

### 5.1 The pattern as documented

From memory/multi_focal_wall_pattern.md (2026-05-08): five instances + the chemistry extension (Sprint W1c arc 2026-05-23) gives six independent observables with the same structural reason. The unifying statement:

> The graph couples discrete labels. The Fock projection couples space. The framework has no native composition theorem for multiple Fock-style projections at once.

This is already stated as a paper-worthy candidate in the memory file. The question is whether to crystallize it formally.

### 5.2 What would a formal theorem look like?

The structural pattern is:

| Sector | What the framework determines | What the framework does NOT determine |
|:-------|:------------------------------|:--------------------------------------|
| QED (LS-8a) | UV-divergent integrand structure (sign, prefactor) | Renormalization counterterms ($Z_2$, $\delta m$) |
| QED (HF-3/4/5) | Discrete-label hyperfine couplings | Multi-projection composition (nuclear motion + spatial × magnetization) |
| Higgs (H1) | AC factor architecture, $J_F$, $\gamma_F$ | Yukawa values |
| Chemistry (W1e) | Selection rules, sparsity, eigenvalue algebraicity | Multi-determinant correlation channels resisting single-particle / mean-field projection |

The unifying observation: **the framework's natural output is the algebraic skeleton + selection rules + scaling laws + transcendental classification. The empirical input is parameter values + multi-projection composition + multi-determinant correlation + renormalization counterterms.**

A formal theorem-statement would look something like:

> **Structural-skeleton-scope theorem (informal).** Let $T_\text{GeoVac}$ be the GeoVac spectral triple at finite cutoff. The framework's natural map $T_\text{GeoVac} \to \text{Observables}$ is functorially complete for single-projection Layer-2 observables (selection rules, scaling laws, transcendental classes) and for two-discrete-label compositions (Wigner-symbol couplings, isostructural invariance). For multi-projection compositions that involve at least one *spatial* projection beyond the first (e.g., recoil = mass-ratio composition; cross-register V_eN = two spatial registers; multi-determinant correlation = density-matrix submanifold beyond zero-entropy class), the framework requires an external composition rule that is NOT generated by the spectral-triple structure alone.

This is the structural-skeleton-scope pattern (memory/geovac_structural_skeleton_scope_pattern.md, 2026-05-07) sharpened.

### 5.3 Is this worth crystallizing as a Paper 50?

Arguments for:
- Six independent observables converging on the same reason is past the point of coincidence.
- Cross-domain (QED + chemistry) consolidation is a natural framing move.
- The η-trivialization theorem for inner factor (Sprint H1 / Angle 2, memory/inner_factor_mellin_engine.md) is the precedent: a structural theorem unifying multiple empirical observations.
- The boundary statement is sharper than the current Paper 32 §VIII.C / Paper 34 §VI Prediction 1' / CLAUDE.md §1.7 formulations.

Arguments against:
- The framework's spectral-triple machinery is structurally tied to the QED-side / single-particle-physics analysis. Extending the same theorem rigorously to chemistry FCI would require new mathematical content (N-representability, density-matrix submanifold geometry) that the math.OA arc (Papers 38–47, 48–49) does NOT cover.
- The CLAUDE.md §1.7 WH-register formulation is sufficient for the framework's operational discipline (when dispatching sprints, the PM uses it to triage; the §1.5 rhetoric rule keeps the formal claims appropriately scoped).
- The W1e wall is at a precision-physics level (recovering 75 mHa binding from a 4000 mHa wall) that no NCG theorem will help with; closure requires *chemistry methods* (Schmidt + multi-zeta + correlation), not NCG machinery.
- Paper 27 already contains the structural-skeleton-scope-chemistry side at the entropy level (one-body entanglement-inert + V_ee generates all correlation entropy). Paper 32 §VIII.C contains the QED side. The cross-paper cross-reference is sufficient.

**Probe 5 verdict: (ii) partial positive at the framing level; (i) clean negative at the formal-theorem level.**

The cross-domain unification IS real and IS paper-worthy, but the appropriate venue is a **synthesis / observation paper** that **cross-references existing structural results** rather than producing new mathematical content. The current CLAUDE.md §1.7 + Paper 34 §VI + Paper 32 §VIII.C network already does this informally; a Paper 50 would just consolidate.

**Recommendation: not worth crystallizing as a separate Paper 50 right now.** The discipline is to use the existing §1.7 framing as the operational tool, cite Papers 27 / 32 / 34 / 35 / 38 for the structural pieces, and reach for the formal theorem only if a new cross-domain phenomenon emerges that the current framing doesn't capture.

If the PI decides otherwise, the right shape is: ~15-page synthesis observation paper, audience = NCG community + quantum-chemistry community, content = (1) the six-wall pattern; (2) the η-trivialization theorem for inner factor extended to "Mellin-factorization + entropy-bound theorem for multi-determinant correlation"; (3) the boundary statement at framework level; (4) open questions (W1e correlation closure, LS-8a renormalization, multi-projection composition). This would be the same role for the multi-focal-composition arc that Paper 49 plays for the OSLPLS strong-form arc.

### 5.4 The substantive new content from the diagnostic

The probes above reveal one substantive new structural reading:

**Chemistry FCI correlation energy and QED renormalization counterterms are in different "non-skeleton" classes.**

QED-side LS-8a / H1 / HF-3/4/5 walls are all in the **calibration-data** class — the framework determines the structure (counterterm structure, integrand form, divergence sign) and external physics supplies the values (Yukawas, $Z_2$, $\delta m$).

Chemistry-side W1e is in the **multi-determinant-correlation** class — the framework determines the FCI matrix structure (Slater integrals, selection rules, eigenvalues) and external chemistry-methods (Schmidt, multi-zeta, CCSD core correlation) supply the correlation channels not reachable from zero-entropy density matrices.

These are categorically different "non-skeleton" classes. The framework's structural-skeleton scope is a single boundary, but what sits **outside** that scope partitions into:

1. **Calibration parameter values** (Yukawas, $Z_2$, $\delta m$ — supplied by external physics / UV completion)
2. **Multi-focal composition** (HF-3/4/5 spatial × spatial — supplied by external composition rules)
3. **Multi-determinant correlation** (W1e — supplied by chemistry-methods on the N-representable density-matrix manifold)

The cross-domain unification at the boundary is real (multi-focal-composition wall pattern). The cross-domain unification *outside* the boundary is illusory — the three "outside" classes are structurally different.

**This is the one durable structural reading of the diagnostic.** It refines the structural-skeleton-scope statement from "GeoVac is a structural-skeleton framework" to "GeoVac is a structural-skeleton framework whose external-input space partitions into at least three categorically distinct classes". Worth recording as a §5.4 finding in this memo; not worth a paper edit until the partition is empirically confirmed at additional walls (e.g., the LS-8a-renorm sprint if it eventually runs).

---

## §6. Honest verdict

**Outcome (i) clean negative on transfer**, with one structural sharpening at (ii)-partial-positive level.

The modular / spectral-triple / master-Mellin-engine machinery (Papers 32, 35, 38, 42, 43, 45, 46, 47, 48, 49) does NOT transfer to chemistry correlation-energy walls, specifically W1e. Four out of five probes returned clean negatives:

1. **Probe 1 (Uhlmann-relative-entropy reading of correlation energy):** clean negative on quantitative transfer; structural analog observation (correlation energy bounded below by entropic functional vanishing on zero-entropy submanifold) is already in chemistry literature and Paper 27, not new content. Sharpens why F4/F5/Schmidt saturate (all project onto zero-entropy submanifold) but doesn't close W1e.

2. **Probe 2 (master Mellin engine on chemistry observables):** clean negative. Chemistry observables are at the algebraic-eigenvalue level, not the transcendental-ring level. No Mellin transform, no $k$-domain.

3. **Probe 3 (Paper 35 Prediction 1 chemistry analog):** clean negative. W1e is π-free at the eigenvalue level per Prediction 1, but Prediction 1 doesn't give a closure mechanism — only a sectoral classification (W1e is in Layer 1).

4. **Probe 4 (modular flow / KMS structure on chemistry):** clean negative on direct transfer (chemistry has no boost, no wedge, no KMS state on ground state); clean negative on Pythagorean analog (no $\mathbb{Z}_2$ chirality grading on chemistry Hamiltonian); reinforced by modular-flow-regularization-negative on the QED side.

5. **Probe 5 (structural-skeleton-scope theorem for chemistry):** (ii) partial positive at the framing level. The cross-domain multi-focal-composition wall pattern is real and paper-worthy at a synthesis level, but a formal Paper 50 is not justified — the CLAUDE.md §1.7 + Paper 27 + Paper 32 + Paper 34 network already covers the structural content. The one substantive new reading is the **partition of the framework's external-input space into three categorically distinct classes** (calibration parameters / multi-focal composition / multi-determinant correlation), which is worth recording but not promoting to a paper.

### Why this is the right outcome

The modular / spectral-triple machinery is structurally tied to:
- Signature-Lorentzian content (wedge, boost, KMS, modular flow)
- Spectral-side trace formula (Mellin transform, heat kernel, master engine)
- Transcendental rings (π-bearing observables via temporal/spectral integrations)
- KMS-state geometry (Connes 1973 cocycle, TICI, Uhlmann monotonicity)

Chemistry FCI is structurally:
- Non-relativistic, time-independent, no boost
- Pure ground state (not Gibbs)
- Algebraic-eigenvalue level (no transcendental ring at eigenvalue level)
- N-representable density-matrix manifold (no KMS state structure)

The two regimes are categorically separated by the same boundary that separates the framework's transcendental-classification machinery from its eigenvalue-computation machinery: **W1e lives below the structural-skeleton scope** that the modular machinery operates on. Crossing that boundary requires chemistry methods (Schmidt orthogonalization, multi-zeta basis, CCSD or DMRG for correlation), not NCG methods.

### Value of this diagnostic

The boundary crystallization. We now know:
- The modular machinery is the right tool for QED-side / signature-physics questions.
- It is NOT the right tool for chemistry binding / correlation-energy questions.
- The cross-domain multi-focal-composition wall is real but doesn't motivate a unified formal theorem; the existing §1.7 + paper network suffices.
- The framework's external-input space partitions into at least three categorically distinct classes (record-only finding, not paper-worthy on current evidence).

This averts:
- A wild-goose chase trying to apply Paper 49's twin-paradox machinery to FCI correlation energy.
- An ill-motivated Paper 50 synthesis paper that would not add new mathematical content.
- Future re-derivation of the same boundary statement under different framing.

### Cross-references and follow-on

The relevant boundary statements are recorded in:
- CLAUDE.md §1.7 (structural-skeleton-scope, multi-focal-composition wall, six instances)
- Paper 27 (one-body entanglement-inert; V_ee generates all correlation entropy; cusp as entropy concentration)
- Paper 32 §VIII.C (Sprint H1 verdict; AC-extension admits Higgs structurally but Yukawa is external)
- Paper 34 §III.28 + §VI (apparatus identity; state-side / spectral-side disjointness; Prediction 1' Layer-2-presence bound)
- Paper 34 §V.B (off-precision matches with error-source classification)
- Paper 19 §sec:w1c_residual + §sec:w1e_refinement_f4_f6 (chemistry-side wall documentation)
- memory/modular_flow_regularization_negative.md (QED-side parallel negative)
- memory/multi_focal_wall_pattern.md (cross-domain pattern)
- memory/geovac_structural_skeleton_scope_pattern.md (boundary statement)

**No follow-on sprint is recommended.** The chemistry arc is paused at W1e (per Sprint W1e Schmidt + core-correlation closeout), and the modular arc is paused at Paper 49 (per math.OA arc consolidation). This diagnostic crystallizes the boundary between them. The natural next direction is whichever PI-selected track is the highest-value forward direction at the current state — but it is NOT a modular-chemistry hybrid sprint.

---

**End of memo.**
