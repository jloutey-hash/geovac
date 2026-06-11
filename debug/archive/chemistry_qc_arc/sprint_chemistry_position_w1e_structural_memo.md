# Sprint Chemistry Position — W1e structural reading vs the multi-focal-composition wall pattern

**Date:** 2026-05-25
**Sprint:** Diagnostic-only structural reading
**Mode:** No code modifications, no paper edits, no memory writes. Memo only.
**Predecessor inputs:** `memory/multi_focal_wall_pattern.md`, `memory/geovac_structural_skeleton_scope_pattern.md`, `memory/sprint_w1c_arc_f4_f6_extension.md`, `debug/sprint_w1c_full_arc_synthesis_memo.md`, `debug/sprint_beta_update_f4f6_memo.md`, `debug/sprint_w1e_schmidt_diagnostic_memo.md`, `debug/sprint_w1e_core_correlation_estimate_memo.md`, Paper 19 §sec:w1e_schmidt_core_correlation, Paper 34 §sec:conv_w1e_cross_domain_wall, Paper 36 §sec:ls8a_result, CLAUDE.md §1.7 WH register.

---

## §0. Executive summary

The W1e wall (inner-region overattraction in NaH binding after W1c+multi-zeta+cross-block-h1 closes W1d at the FCI level) has been documented in the catalogue (Paper 34 §sec:conv_w1e_cross_domain_wall) as the **sixth instance and first cross-domain instance** of the multi-focal-composition wall pattern. The five physics-side instances (H1 / LS-8a / HF-3 / HF-4 / HF-5) and W1e share a unifying structural reason — **the framework determines the structural skeleton but does not autonomously generate calibration data** — but they fall into structurally distinct sub-classes once the unifying statement is unpacked. The honest classification of W1e is **closure-prospect category (a-with-caveat): sprint-scale exhausted but multi-month tractable for the *correlation* sub-mechanism and partial closure of *Pauli-pressure* sub-mechanism, with the deepest open question being whether multi-month efforts can stack to ~98% wall closure or whether the wall is structurally bounded at the ~10–35% combined ceiling of the sprint-scale mechanisms.** W1e is structurally **closer to H1 than to LS-8a**: it admits closable architectural structure (the framework CAN, in principle, do self-consistent HF or strict Schmidt at basis-construction), but the framework does not autonomously do these things from the existing axiom set. The chemistry side's distinguishing feature is that the relevant external input (atomic structure data) is **internally computable in principle** rather than UV-physics-tied as in LS-8a/HF-5.

---

## §1. Wall comparison matrix

| Wall | What's structurally missing | Paper 34 projection slot | Framework-internal vs calibration-external |
|:-----|:---------------------------|:------------------------|:------------------------------------------|
| **H1** (Yukawa selection) | Selection rule for $Y_F$ on the AC inner factor. Inner fluctuations admit Higgs structurally (off-diagonal $D_F$ block between $\mathbb{C}$ and $\mathbb{H}$ summands of $\mathcal{A}_F$), but no GeoVac axiom picks a specific $Y_F > 0$. | Paper 18 §IV sixth tier ("inner-factor input data" / Yukawa Dirichlet ring; added 2026-05-07 from η-trivialization theorem). NOT a Paper 34 §III projection — sits one tier higher in the exchange-constant taxonomy. | **Calibration-external.** No autonomous GeoVac mechanism selects $Y_F$; this is the structural-skeleton-scope statement at the AC inner factor. |
| **LS-8a** ($Z_2/\delta m$ counterterms) | Renormalization counterterms required to extract finite $C_{2S} = +3.63$ from the bare iterated CC spectral action's UV-divergent two-loop SE integrand (divergence $\sim N^{3.43}$, sign and prefactor correct, structure faithful). | §III.6 spectral action — but at the renormalization layer the framework's bare action does NOT autonomously generate. Catalogued as the LS-8a wall in Paper 34 §V.B + Paper 35 §VII.3 (Refined Prediction 1: π in finite parts, divergences from underlying field theory). | **Calibration-external (UV-physics-tied).** Counterterms parameterize an OPEN-ENDED UV problem (Landau pole is a real feature); structurally bounded. |
| **HF-3** (recoil cross-register $V_{ne}(\hat r_e, \hat R_n)$) | Two-body coordinate operator coupling electron position to quantum-mechanical nucleus position. Track NI's cross-register Hamiltonian couples 12 spin Pauli strings across (proton spin) ⊗ (electron spin) but exactly 0 spatial qubit pairs. | Sits in the gap between §III.14 (rest-mass projection, scalar mass-ratio) and §III.17/18 (nuclear charge-density / magnetization-density operators at the electronic-coordinate position). NOT a current §III slot. Catalogued as a multi-focal sub-wall in the chain. | **Framework-internal-missing.** The framework has the discrete-label algebra to couple spin labels cleanly; the missing piece is the multi-focal SPATIAL-coordinate composition operator. Mitigability: W1a class (recoil-mixing extension to §III.18 closed at leading order 2026-05-09; full two-body operator is sprint-scale W1a refinement, NOT W2a). |
| **HF-4** (Zemach magnetization-distribution) | Magnetization-density operator at the proton spatial coordinate. The framework's `hyperfine_coupling_pauli` is structurally point-like in proton position; $r_Z$ must enter as Layer-2 scalar input. | §III.18 (nuclear magnetization-density / Zemach, the 18th projection added 2026-05-09). | **Framework-internal-extendable.** Operator-level Zemach landed in W1b operator-level extension (May 2026, `geovac/magnetization_density.py`). NLO recoil-mixing reproduces Eides Tab 7.3 at 0.012% of LO shift. The "wall" here is largely closed at leading and NLO operator level; remaining residual is calibration-external ($r_Z$ value from QCD-internal nucleon structure). |
| **HF-5** (multi-loop $a_e$ on Dirac-S³) | Same renormalization wall as LS-8a, vertex sector (two-loop topologies VP-on-photon / SE-on-electron / crossed). Bare topologies UV-divergent $\sim N^{3.79}$; framework reproduces (α/π)² prefactor and SO(4) selection at every vertex but no autonomous $Z_2/Z_3/\delta m$ extraction. | §III.6 spectral action (same as LS-8a). Catalogued in Paper 34 §V.B as the cleanest known A-class error code. | **Calibration-external.** Same wall as LS-8a; the structural-skeleton-scope statement at the renormalization layer. |
| **W1e** (inner-region overattraction at NaH after W1d closure) | Pauli repulsion at the F3-constructed bonding orbital level + deep multi-determinant FCI correlation in the bonding/antibonding subspace + sub-shell correlation at the heavy-atom valence/core boundary. Three sub-sub-mechanisms (basis-truncation 10.2% / Hartree-pressure 25.7% / deep-correlation-residual 65–90%) at sprint scale. | §III.20 Phillips-Kleinman / core-valence orthogonality (the 20th projection added 2026-05-15 Sprint 1) is the named slot for the Pauli-pressure sub-sub-mechanism. §III.25 coupled-channel / adiabatic-curve (the 25th projection) is the slot for the deep-correlation sub-sub-mechanism in spirit, though the W1e case sits in single-center FCI rather than nuclear adiabatic. | **Mixed.** Framework-internal-extendable for basis-truncation (W1e-basis: bigger basis is mechanical) and Hartree-pressure (W1e-Hartree: explicit-core FCI is multi-week but well-defined). Multi-determinant correlation in the bonding subspace at the heavy-atom valence/core boundary is the deepest sub-component; whether it's framework-internal-extendable (via CCSD-on-12-electrons-class infrastructure) or calibration-external (analog of inner-factor input data — atomic structure information from outside the bare GeoVac axiom set) is an open question this memo addresses below. |

### §1.1 The two-bucket sub-classification within the wall pattern

The six walls split into two structural sub-classes that the unifying statement does not initially separate:

- **Calibration-external (UV-tied / inner-factor-tied):** H1, LS-8a, HF-5. These tie to physics outside the bare GeoVac axiom set in a structurally permanent way — UV physics for LS-8a/HF-5, AC inner-factor data for H1. No autonomous extension of the framework's existing machinery can close these without importing external physics. Mitigability is structurally bounded.

- **Framework-internal-extendable (architectural-extension-tractable):** HF-3, HF-4, and the basis-truncation/Hartree-pressure components of W1e. These tie to data that the framework's own internal physics CAN in principle compute — two-body spatial operators (HF-3), magnetization-density operators (HF-4 partially landed), atomic-structure data for second-row binding (W1e). Each is a multi-week to multi-month architectural extension, but the closure path is well-defined within the framework's existing computational primitives.

- **The W1e-deep-correlation sub-sub-mechanism (65–90% of wall by sprint-scale measurement)** is the ambiguous member of this taxonomy. The chemistry-arc literature (Iron-Martin 2003, Sylvetsky-Martin 2018) places binding-relevant CV correlation at sub-percent of $D_e$, far below the structural wall. The sub-mechanism's true magnitude must therefore be different from what the F4/F5/F6/Schmidt/core sprint scales measured — it's plausibly a **single-determinant FCI variational physics deficit at the active-space level**, not a core-correlation effect per se. This places it closer to the basis-truncation sub-component than the Day-1 core-correlation memo's literature framing suggested; the sub-component is plausibly closable by extending the FCI to include all 12 NaH electrons (lifting [Ne] out of the screening potential).

### §1.2 The unifying statement re-stated more precisely

CLAUDE.md §1.7 / `memory/multi_focal_wall_pattern.md` give the unifying statement: *"The graph couples discrete labels. The Fock projection couples space. The framework has no native composition theorem for multiple Fock-style projections at once."* This is correct but **too coarse to separate calibration-external from framework-internal-extendable walls.** A sharper statement after W1e:

> The framework's compact-graph axiom (Paper 0 packing) outputs the algebraic skeleton (n,l,m,s and operators on it) at a single focal length. Calibration data spans two structural classes:
>
> **(α) Calibration data tied to physics outside the GeoVac scope** (UV-completion counterterms, AC inner-factor selection) is structurally bounded — no extension of the framework's primitive operations can generate it autonomously. The "second packing axiom" question (could a different structural axiom generate this data?) is open and unsolved.
>
> **(β) Calibration data computable from the framework's existing computational primitives, but not generated by the bare axiom set** (multi-Fock-projection spatial operators; atomic-structure data for heavy-atom binding) is architecturally tractable — the framework's own physics can compute it given a multi-week to multi-month engineering commitment, but the bare axiom set does not autonomously execute these computations.

W1e is **predominantly class (β)** at the sprint-scale-empirical level, with a residual fraction of class (α)-flavored deep correlation that remains to be tested at the multi-week+ scale. The five physics-side walls split: H1 / LS-8a / HF-5 are (α); HF-3 / HF-4 are (β) with HF-4 substantially closed already.

---

## §2. W1e's distinguishing features

### §2.1 Cross-domain transition: physics → chemistry

W1e is the first chemistry-domain instance of the wall pattern. The structural-skeleton-scope statement now generalizes cross-domain: the framework determines the algebraic skeleton across both physics and chemistry observables, with calibration data being external input in both domains. This is the substantive Paper 34 §sec:conv_w1e_cross_domain_wall content added 2026-05-23.

The cross-domain transition reveals what was implicit in the physics-side framing: **the wall pattern is not specifically about Fock-S³ projections or QED.** It's about the relationship between the bare axiom set and calibration data more generally. Chemistry adds a domain where the calibration data is **atomic structure** (frozen-core potentials, correlation energies, basis-shape parameters), which is structurally distinct from UV counterterms or Yukawa values.

### §2.2 Structural reading: W1e is closer to H1 than LS-8a

Three independent comparisons place W1e in the H1 sub-class rather than the LS-8a sub-class:

**(a) The framework CAN build the right structure.** H1's AC extension admits a Higgs structurally (inner-fluctuation off-diagonal $D_F$ block); W1e's cross-block h1 extension constructs the bonding orbital structurally (FCI naturals $[1.9991, 0.0007]$, four-orders-of-magnitude advance over W1c+mz alone). Both walls have the same shape: structural construction works, autonomous numerical preference does not. LS-8a/HF-5 by contrast cannot autonomously build the right finite parts — the UV divergence is genuinely present, not a numerical-preference issue.

**(b) The closure mechanism is "external preference data," not "external counterterms."** H1 needs an external $Y_F$ value supplied to the Yukawa Dirichlet ring. W1e needs an external prescription for *which* configuration is energetically preferred at $R_\text{eq}$ — equivalently, an external mechanism that supplies the Pauli pressure (Schmidt orthogonalization, explicit-core FCI, or self-consistent HF) that the bare two-electron FCI architecture lacks. LS-8a needs renormalization counterterms parameterizing UV physics; these have a different structural character (open-ended Landau-pole class vs. closed-form atomic-structure class).

**(c) The "framework-internally computable in principle" property.** H1's $Y_F$ is computable from any UV completion that fixes the Higgs Yukawa coupling — the framework's structural-skeleton role is to provide the AC architecture, the Yukawa value comes from a UV completion outside GeoVac. W1e's atomic-structure data (CR67 exponents, self-consistent HF energies, sub-shell correlation) is computable from atomic-physics machinery the framework already partially has (e.g., `geovac/atomic_classifier.py`, `geovac/neon_core.py` for screened radial Schrödinger). The chemistry-side gain over H1 is that the relevant atomic-physics machinery is *closer at hand* — multi-week to multi-month architectural commitments rather than a fundamental UV completion.

### §2.3 The "active-space deficit" structural reading

The cleanest single statement of W1e is **the framework's two-electron active-space treatment lacks the multi-determinant correlation channels that the full 12-electron problem provides at NaH binding geometries.** The frozen-[Ne] screening potential captures core-electron mean-field at the cross-V_ne kernel level (W1c), but the 2-electron FCI cannot represent the multi-configurational correlation in the bonding/antibonding subspace that would supply the Pauli pressure preventing inner-region collapse. The Phillips-Kleinman family of methods is the standard atomic-physics response to this deficit; F4's rank-1 PK saturated at 43% closure ceiling because rank-1 is not enough — full PK or Schmidt at the basis-construction level should do better in principle but Schmidt was ruled out empirically at H 1s scale.

**The honest reading is therefore:** W1e is not a wholly new wall class but is the chemistry-analog of "the structural skeleton is right; calibration data must be supplied externally." The form of the calibration data is **multi-determinant correlation information at the heavy-atom valence/core boundary** — which is in principle computable from atomic-physics infrastructure the framework's neighbors (modern quantum chemistry codes) routinely deliver.

---

## §3. Paper 18 / Paper 34 readings

### §3.1 Transcendental taxonomy reading (Paper 18)

W1e is **dimensionful Ha-scale energy wall**, not a transcendental coefficient with structural significance. The W1e magnitude (4.37 Ha at NaH max_n=2; ~3.9 Ha at max_n=4) is an FCI residual against experimental $D_e = 0.072$ Ha — pure numerics, no transcendental signature.

That said, the cross-block-h1 kernel that *constructs* the bonding orbital (and that surfaces W1e by enabling FCI bonding without supplying repulsion) is in the algebraic-implicit tier of the Paper 18 hierarchy: closed-form analytic two-center one-body integrals via 2D axial Gauss-Legendre + radial-grid quadrature in `geovac/cross_block_h1.py`. The wall's *appearance* depends on a Paper 18 mechanism (the framework's per-block algebraic h1 architecture extended to inter-block coupling), but the wall *itself* is a numerical-energetic deficit not a transcendental classification.

**Verdict:** W1e has no transcendental-taxonomy reading. This is consistent with the broader observation that chemistry-side walls in the multi-focal pattern tend to manifest as energetic deficits rather than transcendental-classification mismatches; the QED-side walls (LS-8a / HF-5) manifest as both (divergent integrand AND missing finite parts requiring renormalization).

### §3.2 Paper 34 projection slot

W1e is currently catalogued in Paper 34 as §sec:conv_w1e_cross_domain_wall (in the §V.D living catalogue of convention exposures and architectural limitations, class iv "framework-architecture limitation"). It does NOT have a dedicated §III projection slot like §III.20 PK or §III.25 coupled-channel.

**Does W1e correspond to a missing §III projection?** Candidates:

- **§III.20 Phillips-Kleinman / core-valence orthogonality** is the closest existing slot. It captures the *operator form* of the Pauli-repulsion mechanism but is implicitly tied to the rank-1 / single-particle architecture that F4 ruled out at 43% saturation. Either (a) §III.20 should be re-read to include the full-rank Schmidt analog and the operator-level explicit-core FCI extension (with W1e as an "off-axis" application where the single-particle reading fails), or (b) a new §III.29 projection slot should be added: "multi-determinant Pauli pressure / explicit-core correlation" — an operator-level extension of §III.20 to the case where rank-1 PK cannot supply enough repulsion. The latter is the cleaner taxonomic option.

- **§III.25 coupled-channel / adiabatic-curve** is the analogous slot in spirit (multi-channel correlation), but it's tied to nuclear-frame adiabatic separation rather than electron-correlation channels. Not a direct match.

- **Paper 18 §IV sixth tier (inner-factor input data)** is the structural-skeleton-scope tier for H1; the analogous chemistry-side tier (if added) would be "atomic-structure input data" — frozen-core exponents, correlation energies at the heavy-atom valence/core boundary, etc. This is the cleanest taxonomic location for the calibration-external residual of W1e, leaving the calibration-internal-extendable sub-components in §III slots.

**Recommendation (not a paper edit, just a flag):** The cleanest taxonomy refinement would be a new §III.29 projection (multi-determinant Pauli pressure / explicit-core correlation), absorbing F4/F5/Schmidt into a single operator-family entry, plus a Paper 18 §IV seventh tier (atomic-structure input data) holding the multi-month-architectural-extension-tractable but currently-external-input residual.

### §3.3 Two structural-skeleton tiers

After W1e the structural-skeleton-scope statement reads two-tiered:

- **Tier 1 (current Paper 18 § IV.6 inner-factor input data):** Calibration-external, UV-physics-tied or AC-inner-factor-tied. H1 / LS-8a / HF-5 sit here. Structurally bounded; the "second packing axiom" is the speculative path forward.

- **Tier 2 (proposed Paper 18 §IV.7 atomic-structure / multi-focal input data):** Calibration-internal-extendable but external to the bare axiom set. HF-3 / HF-4 (operator-level largely closed) / W1e (partially closed) sit here. Mitigability is via multi-week to multi-month architectural extensions; closure path is well-defined within the framework's existing computational primitives.

This refinement is the **substantive taxonomic content** the W1e arc has produced beyond the H1 / LS-8a / HF-3/4/5 set.

---

## §4. Closure-prospect classification

**Honest verdict: (a) Sprint-scale exhausted but multi-month tractable, with a calibration-external residual of unclear magnitude.** 

Specifically: W1e decomposes into three sub-sub-mechanisms (per `memory/sprint_w1c_arc_f4_f6_extension.md` §1):

1. **W1e-basis-truncation (~10.2% of wall, F6 measured at max_n=4):** Mechanically closable by basis enlargement — but F6 already paid the algorithmic exponent and got 10.2%, suggesting diminishing returns at max_n=5+. Closure path is multi-week implementation at max_n=5 + extrapolation analysis. Expected outcome: 15–25% ceiling.

2. **W1e-Hartree-pressure (~25.7% of wall, F5 predicted by Hartree J-K analysis):** Mechanically closable by **explicit-core FCI** — lift the [Ne] core out of the screening potential and into the FCI sector. This is a multi-month engineering commitment (12-electron NaH FCI at max_n=2 has FCI dim that exceeds sprint-scale compute by orders of magnitude) but well-defined: existing `geovac/neon_core.py` infrastructure can be inverted from screening-potential mode to explicit-orbital mode. Expected outcome: ~25–35% closure on this sub-mechanism alone.

3. **W1e-deep-correlation-residual (65–90% of wall, untested at sprint scale):** This is the dominant component and the deepest open question. Two readings:
   - **Reading A (framework-internal-extendable, multi-month):** It's CCSD/CCSD(T)-class multi-determinant correlation that the FCI active space cannot supply because the active space is too small. Closable by post-FCI extensions (CCSD on 2-electron active space + the 10 [Ne] electrons; or perturbation-theoretic corrections post-FCI). Multi-month sprint, well-defined. Expected outcome: substantial closure (modern CCSD(T) on NaH recovers $D_e$ to sub-kcal/mol against experimental).
   - **Reading B (calibration-external residual):** The "active-space deficit" cannot be cleanly closed without lifting the framework's fundamental two-electron active-space architecture, which would amount to abandoning the composed-architecture decomposition that the framework's qubit-encoding advantage rests on. The deep-correlation residual is structurally tied to the framework's architectural choice (composed/qubit-friendly two-electron active space vs explicit-correlation 12-electron treatment). In this reading, the residual is calibration-external in a chemistry-specific sense: the framework has chosen a structural scope (two-electron active space + frozen-core screening) that doesn't reach atomic-structure-grade precision at the heavy-atom valence/core boundary.

**Most likely actual closure path (Reading A weighted):** Multi-month combined sprint of explicit-core FCI (W1e-Hartree-pressure, ~30% closure) + max_n=5 basis enlargement (W1e-basis, ~15% closure) + post-FCI CCSD-class correlation (W1e-deep, ~40–60% closure) = combined ~85–95% closure, leaving a small calibration-external residual that the modern quantum-chemistry community handles via empirical core-polarization potentials (CPP) or basis-set extrapolations beyond CBS.

**Most likely structural answer (Reading B weighted):** The framework reaches ~30–50% closure via the combined multi-month sprint, hits a true architectural-scope boundary, and the remaining ~50% is calibration-external in the sense that recovering it requires abandoning the composed-architecture decomposition.

The two readings cannot be distinguished without running the multi-month explicit-core FCI sprint. The named diagnostic that *could* distinguish them is **a Hartree-Fock benchmark on the explicit-core 12-electron NaH problem**: if HF gives binding (D_e in [0.05, 0.10] Ha range, internal minimum present), then the W1e wall is genuinely an active-space deficit (Reading A correct) and CCSD-class post-HF corrections recover sub-percent residual. If HF still gives monotone descending PES, then the wall is structurally tied to the composed-architecture architectural choice (Reading B correct).

This is the **single most important next diagnostic** for placing W1e structurally. It's a ~2-4 week sprint (implement explicit-core 12-electron HF in the framework's basis machinery) and would decisively distinguish (a) multi-month tractable from (b) structurally bounded.

---

## §5. Specific named follow-ons

### §5.1 Top priority: explicit-core HF diagnostic (~2–4 weeks, distinguishes Reading A from Reading B)

**Mechanism:** Implement a Hartree-Fock-on-12-electron-NaH solver in the framework's existing basis machinery. Use single-zeta CR67 + multi-zeta for valence (machinery exists in `geovac/multi_zeta_orbitals.py`), 12-electron HF at the production-feasible matrix size (~120 orbitals at max_n=2 → 120×120 Fock matrix; routine for HF). Run 8-point PES sweep across [2, 10] bohr.

**Decision gate:** If 12-electron HF gives internal minimum at $R_\text{eq} \in [3.0, 4.5]$ bohr with $D_e \in [0.04, 0.12]$ Ha → **Reading A confirmed** (W1e is active-space deficit, multi-month CCSD-class post-HF closure path is well-defined, framework can reach sub-percent NaH binding). If HF still gives monotone-descending PES → **Reading B confirmed** (W1e is structurally tied to composed-architecture choice, framework's scope limit is at chemistry-grade precision for second-row hydride binding).

**Expected outcome (PM judgment):** Reading A is more likely (modern HF-on-NaH gives ~80% of experimental $D_e$ before correlation corrections; the framework's machinery on 12 explicit electrons should not be structurally worse). But Reading B is the deeper structural finding if it surfaces — it would be the chemistry-side analog of H1's "no autonomous Yukawa selection" verdict and worth its own paper.

### §5.2 Conditional follow-on (if Reading A): multi-month explicit-core CCSD/CCSD(T) (~3–6 months)

Conditional on Reading A confirmation. Implement post-HF correlation corrections at CCSD or CCSD(T) level on the 12-electron NaH active space. Sub-percent NaH binding is the deliverable. Expected outcome: ~85–95% closure of the W1e wall, with the residual being Layer-2 polish (relativistic corrections, core polarization potentials).

This is a substantial engineering commitment and may compete with other priorities (math.OA arc, Paper 48-49 polish, precision-catalogue extension). The PI judgment is whether second-row binding accuracy is load-bearing for the framework's positioning — if the framework's primary computational advantage is quantum-resource estimation (Paper 14/20), the answer is plausibly "no, residual structural-skeleton-scope at chemistry precision is acceptable."

### §5.3 Conditional follow-on (if Reading B): document the architectural-scope boundary (~1 week, paper edit only)

Conditional on Reading B confirmation. The W1e wall is then the cleanest cross-domain instance of the multi-focal-composition wall pattern in its **calibration-external (architecturally-scoped) sub-class**, structurally analogous to H1 / LS-8a in the physics domain. Paper 19 / Paper 32 / Paper 34 updates document the architectural-scope boundary explicitly:

- Paper 18 §IV gains a seventh tier (atomic-structure / multi-focal input data) holding the chemistry-side architecturally-scoped residual.
- Paper 32 §VIII.D extends with a "chemistry-domain architectural-scope" subsection paralleling the §VIII.B/VIII.C SM-gauge / inner-factor verdicts.
- Paper 34 §sec:conv_w1e_cross_domain_wall is sharpened from "first cross-domain wall instance" to "first chemistry-domain architecturally-scoped wall instance."

The structural-bound theorem (analog of η-trivialization for inner-factor calibration data) would read: **on the framework's composed-architecture active space, multi-determinant correlation channels at the heavy-atom valence/core boundary are bounded above by the active-space size; the binding-energetic residual at second-row hydrides scales with the active-space deficit and is structurally tied to the framework's architectural choice of two-electron active space + frozen-core screening.**

This is a positive structural-theorem candidate IF the explicit-core HF diagnostic returns Reading B — and it would be a Paper-48-class deliverable in its own right.

### §5.4 Lower priority: self-consistent HF iteration in `geovac/neon_core.py`

The PK cross-center synthesis memo (chemistry_arc_paused_w1c_residual) named self-consistent HF iteration as "the principled long-term alternative that closes the cliff for ALL Z" — at the heavy-atom screening problem in Cs HFS (Z>20 cliff, Paper 34 §V.C.6). This is a different sprint than the W1e closure (Cs HFS is electronic structure on a heavy atom; W1e is bond formation on a light-heavy diatomic), but the underlying machinery (self-consistent $Z_\text{eff}$ depending on cross-center density) is the same. A combined sprint could in principle attack both walls with one implementation.

Multi-month, well-defined, lower priority than the explicit-core HF diagnostic.

### §5.5 Lowest priority: pause chemistry, return when atomic-physics infrastructure is needed elsewhere

The natural alternative is to pause the W1e closure attempt and pivot back to math.OA / precision-catalogue / Paper 48-49 polish, leaving the W1e architectural-scope verdict at the empirically-grounded "sprint-scale exhausted, multi-month tractable for components" status. This is consistent with the F6 §5 recommendation that closed the F4-F6 arc (PI default option: documentation batch + pivot back to math.OA).

The PI choice between §5.1 (do the diagnostic) vs §5.5 (pause and pivot) is a strategic-direction question rather than a technical one. The diagnostic is cheap enough (~2-4 weeks) that running it before any deeper pivot is the recommended default.

---

## §6. Honest scope and confidence

| Claim | Confidence | Reason |
|:---|:---:|:---|
| W1e is the sixth instance and first cross-domain instance of the multi-focal-composition wall | HIGH | Paper 34 §sec:conv_w1e_cross_domain_wall and Paper 19 §sec:w1e_schmidt_core_correlation both document this explicitly; structural unifying statement holds |
| W1e is structurally closer to H1 than LS-8a | MEDIUM-HIGH | Three independent comparison arguments (§2.2 (a)-(c)); not yet structurally proved as a theorem |
| W1e's three sub-sub-mechanisms (basis / Hartree / deep-correlation) are independent | MEDIUM | Sprint-scale measurements (10.2% / 25.7% / ~65-90%) give partial-additivity ceiling ~35-45%; whether they stack to the full wall under multi-month implementation is an open question |
| W1e's deep-correlation sub-sub-mechanism is structurally analogous to H1's Yukawa Dirichlet ring | MEDIUM-LOW | Speculative; rests on the "active-space deficit" reading which is testable via explicit-core HF; if reading B (architectural-scope boundary), then closer to H1 than current literature CV-correlation framing |
| The explicit-core HF diagnostic distinguishes Reading A vs Reading B at decision-gate-grade rigor | MEDIUM-HIGH | Standard chemistry-physics: HF-on-12-electron-NaH recovers ~80% of experimental D_e in standard calculations; if framework HF doesn't, that's a strong Reading B signal |
| The natural Paper 18 §IV.7 / Paper 34 §III.29 taxonomic addition (atomic-structure input data; multi-determinant Pauli pressure) is the right framing | MEDIUM | Recommended as a flag but not a paper edit; depends on PI judgment about whether the cross-domain pattern crystallization is paper-worthy or stays as catalogue cross-reference |

The high-confidence claims are the structural-classification verdict (W1e is the sixth-and-first-cross-domain wall instance, with the structural-skeleton-scope reading). The medium-confidence claims are the closure-prospect details, which depend on the explicit-core HF diagnostic outcome.

---

## §7. Net position

**Closure-prospect classification: (a-with-caveat) Sprint-scale exhausted but multi-month tractable for components, with a calibration-external residual whose magnitude (Reading A small ~5-15%, Reading B substantial ~50%) cannot be determined without running the explicit-core HF diagnostic.**

**The structural-skeleton-scope framing extends from "five physics observables" (CLAUDE.md §1.7) to "five physics + one chemistry observable, with the chemistry instance distinguishing calibration-internal-extendable from calibration-external sub-residuals."** The substantive new content beyond the existing cataloguing is the **two-tier refinement** of the calibration-data classification (Paper 18 §IV.6 inner-factor input data; proposed §IV.7 atomic-structure / multi-focal input data) and the **named diagnostic** (explicit-core HF on 12-electron NaH) that would distinguish Reading A from Reading B at decision-gate-grade rigor.

**Most likely actionable next step:** If the PI continues the chemistry arc, run §5.1 (explicit-core HF diagnostic, ~2-4 weeks) before any further structural commitment. If the PI pivots back to math.OA / precision catalogue, the current W1e characterization (sprint-scale exhausted, multi-month tractable for components, cross-domain wall instance documented) is a clean pause point and can be left at this maturity without further sprint commitment.

---

**End of memo.**
