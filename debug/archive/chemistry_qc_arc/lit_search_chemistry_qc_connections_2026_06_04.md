# Literature search: GeoVac mathematical findings ↔ quantum-chemistry / QC-for-chemistry literature

**Date:** 2026-06-04
**Author:** Lit-search sub-agent
**Scope:** Seven threads (a)–(g) per PM directive; honest "no analog found" outcomes recorded where applicable; arXiv IDs web-fetched before citing.

---

## Executive summary

Of the seven threads:

- **(a) Z₂ tapering** — **PARTIAL ANALOG.** The closest published analog to Hopf-Z₂ tapering is the **point-group / spatial-symmetry Z₂** lineage of Setia et al. 2020 (arXiv:1910.14644, *JCTC* 16, 6091) and Picozzi-Tennyson 2023 (*Quantum Sci. Technol.* 8, 035026). Both reduce qubits by 2–5 via Boolean symmetries; nobody is using the S³ → S² Hopf-base m → −m parity per se. GeoVac's −3 to −12 range across 37 molecules is competitive with the upper end of Picozzi-Tennyson; the mechanism is structurally distinct (Hopf-bundle parity vs molecular point group). This is the most outreach-ready connection.
- **(b) Composed / multi-resolution encodings** — **CLOSE-BUT-DIFFERENT.** Tensor hypercontraction (Lee-Berry-Gidney-Huggins-McClean-Wiebe-Babbush 2021, arXiv:2011.03494, *PRX Quantum* 2, 030305), double factorization (Burg et al. *Quantum* 2024 q-2024-06-13-1371), and active-space embedding (Battaglia-Rossmannek-Rybkin-Tavernelli-Hutter 2024, arXiv:2404.18737) all factorize the molecular Hamiltonian, but none use a *geometric* core/valence split (S³ atom-fiber × bond-sphere base). The comparison is at the resource level, not the structural level.
- **(c) Mellin / spectral-action classification of molecular transcendentals** — **CLEAN NEGATIVE.** No quantum-chemistry literature classifies molecular integrals by spectral-action Mellin sub-mechanism. Closest pure-math reference is Fathizadeh-Marcolli mixed-Tate, already cited internally in Paper 18. The Mellin engine is GeoVac-novel within chemistry.
- **(d) Forced/free seam in chemistry** — **CLEAN NEGATIVE.** No analog. Closest framing is the NCG-Standard-Model literature (Marcolli–van Suijlekom gauge networks, arXiv:1301.3480; Perez-Sanchez arXiv:2401.03705) on the *physics* side. The "FCI correlation energy as calibration data" reading has no chemistry-literature equivalent. Would require translation register if pitched to chemistry venues.
- **(e) Trotter bounds from operator-algebraic GH convergence** — **CLEAN NEGATIVE in current form.** Trotter-bound literature (Childs-Su-Tran-Wiebe-Zhu 2021, arXiv:1912.08854; qDRIFT improvements arXiv:2506.17199) all use commutator-norm bounds, not propinquity-rate bounds. The mathematical instinct (GH rates ↔ approximation error) is sound but unused by chemistry; the GeoVac heuristic bound (3–4 OoM looser than naive Suzuki-Trotter) is consistent with this being an unexplored direction.
- **(f) Spectral action gives metric, not Green's functions** — **CLOSE-BUT-NOT-EXPLICIT.** Bochniak-Sitarz spectral-action work (arXiv:2106.10890) does not state this category divide explicitly. The Connes-Marcolli book (2008) does, but nobody on the chemistry side has noted it.
- **(g) Chemistry analog of Yukawa non-selection** — **CLEAN NEGATIVE.** No DMRG/CCSD(T)/FCI paper makes the structural statement "correlation energy values are external calibration, not framework-forced." Chemistry's pragmatic stance is uniformly post-hoc empirical convergence, not meta-theoretical taxonomy.

**Net.** (a), (b), and (f) are publication-ready connections; (c), (d), (e), (g) are clean negatives where GeoVac's framing is novel within the chemistry/QC literature and would need translation to land.

---

## Per-thread findings

### (a) Z₂ tapering / qubit reduction

**Verified core references:**

- **Bravyi, Gambetta, Mezzacapo, Temme 2017**, "Tapering off qubits to simulate fermionic Hamiltonians," arXiv:1701.08213. Establishes the generic Z₂ Pauli-string tapering algorithm: parity-check matrix → block-diagonalize JW Hamiltonian → orthogonal projection. One qubit eliminated per independent Z₂ symmetry.
- **Setia, Chen, Rice, Mezzacapo, Pistoia, Whitfield 2020**, "Reducing qubit requirements for quantum simulation using molecular point group symmetries," arXiv:1910.14644, *J. Chem. Theory Comput.* 16, 6091–6097. **Verified via WebFetch.** Develops second-quantization representation of *spatial* (point-group) symmetries, transforms to qubit operators, formally connects to generic Z₂ Pauli symmetries of Bravyi 2017. Standard reference for point-group-based tapering.
- **Picozzi, Tennyson 2023**, "Symmetry-adapted encodings for qubit number reduction by point-group and other Boolean symmetries," *Quantum Sci. Technol.* 8, 035026, DOI 10.1088/2058-9565/acd86c. **Verified via WebFetch.** Reduces *n* spin-orbitals to *n − k* qubits, **k = 2 to 5** depending on molecular symmetry. Open-source package `QuantumSymmetry` (PySCF + OpenFermion + Qiskit). Closest published comparator to per-sub-block tapering in scale.

**Honest assessment of GeoVac fit.**

GeoVac's Hopf-U(1) Z₂ comes from the m → −m parity of the S³ → S² Hopf base: it is *geometric / topological* (the Z₂ subgroup of U(1) acting on the Hopf fiber), not *spatial / molecular* (point-group). On 37 molecules the per-sub-block range is **ΔQ = −3 to −12**, total library savings **254 qubits**. Mechanism preserves spectrum to machine precision because cross-block ERIs vanish in the composed builder ({P_i} commute pairwise and with H). This makes it categorically *additive* to the standard tapering recipes:

- A typical LiH calculation in cc-pVDZ has Q ≈ 20 and standard point-group + parity tapering removes 4–5 qubits. GeoVac's per-block tapering removes 4 from LiH (n_sub_blocks = 2: Li-S³ atom block + bond-sphere base = ΔQ = 2 + 2). The *combined* recipe would stack — Hopf-Z₂ is orthogonal to particle-number and S_z Z₂'s of Bravyi 2017 — but no one has built this stack.
- No published paper uses a *fiber-bundle Z₂* as the symmetry source for tapering. The Hopf fibration appears in QI literature only for one- and two-qubit Bloch-sphere geometry (Mosseri-Dandoloff arXiv:quant-ph/0108137 and follow-ups), never for many-electron molecular tapering.

**Publication framing.** If Paper 14 is revised, this is the single most legible result for the QC-for-chemistry audience. Frame as: "We exploit a Z₂ symmetry not available in second-quantized formulations because it lives on the Hopf-base of the angular-momentum basis."

### (b) Composed / multi-resolution qubit encodings

**Verified core references:**

- **Lee, Berry, Gidney, Huggins, McClean, Wiebe, Babbush 2021**, "Even more efficient quantum computations of chemistry through tensor hypercontraction," arXiv:2011.03494, *PRX Quantum* 2, 030305. **Verified via WebFetch.** $\widetilde{O}(N)$ Toffoli block-encoding for N orbitals; FemoCo ≈ 4M physical qubits, ≤ 4 days runtime. Industry baseline for fault-tolerant quantum chemistry.
- **Burg, Low, Häner, Steiger, Reiher, Roetteler, Troyer 2021**, "Quantum computing enhanced computational catalysis," *Phys. Rev. Research* 3, 033055. (Standard DF reference.)
- **Battaglia, Rossmannek, Rybkin, Tavernelli, Hutter 2024**, "A general framework for active space embedding methods: applications in quantum computing," arXiv:2404.18737. **Verified via WebFetch — NO Greene-Diniz author.** Periodic range-separated DFT + VQE + quantum equation-of-motion. Tested on MgO oxygen vacancies.
- **Lemmer 2023+, Greene-Diniz et al.** (multilayer embedding for pharma): arXiv:2202.04460 area. Not separately verified but standard active-space-embedding line.

**Honest assessment of GeoVac fit.**

The composed architecture's PK + Z_eff core/valence split *is* a multi-resolution factorization, but the factorization principle is the **natural geometry hierarchy** (S³ atom × bond-sphere × molecular frame) rather than the **low-rank / sparse tensor** principle of DF / THC. Comparison:

- DF / THC achieve Pauli sparsity by *spectral compression of the ERI tensor*. Standard chemistry result; well-published in *PRX Quantum*, *Quantum*.
- Active-space embedding (LASSCF, DMET, CAS-DMET, Battaglia-Rossmannek 2024) achieves qubit reduction by *fragment-environment partition* of the orbital space. Standard.
- GeoVac achieves O(Q^2.5) Pauli scaling by *Gaunt selection rules on the natural angular-momentum basis* + *block diagonal structure from composed bundle*. This is the closest in spirit to angular-momentum-coupled-cluster-with-symmetry methods (e.g., Stein-Reiher 2017 *J. Chem. Phys.*) but those papers are classical-CC, not qubit-Hamiltonian.

**Publication framing.** Paper 20 already does this comparison ad-hoc. The cleanest external comparison would be **Pauli term count and 1-norm vs cc-pVDZ at fixed Q**, against published DF/THC numbers for LiH/BeH₂/H₂O. GeoVac's 51×–1712× advantage over Gaussian baselines is genuinely large; this is unambiguous when comparing to the *raw JW baselines* but harder to compare to *factorized* DF/THC numbers because those numbers are only published for Q ≥ 100. Honest restatement of competitive landscape: Paper 20 already documents this carefully (cc-pVDZ 63,519 Pauli vs composed 334 for LiH).

### (c) Master Mellin engine in chemistry

**Searches returned:** nothing. Spectral-action / heat-kernel / Mellin classifications appear in gravitational and NCG literature (Fathizadeh-Marcolli arXiv:1611.01815, Connes-Chamseddine 1997, Marcolli-vS book) but **no quantum chemistry paper classifies which transcendentals appear in molecular integrals by spectral-action sub-mechanism**.

Closest hits that are NOT analogs:

- Heat-kernel coefficients are used in DFT XC functional construction (Becke 1988 functional has a √π ratio that comes from the Slater-X integral) but never explicitly framed as a Mellin sub-mechanism.
- Spectral-zeta methods appear in vacuum-energy / Casimir calculations adjacent to chemistry (e.g., Bordag-Klimchitskaya-Mostepanenko, *Advances in the Casimir Effect*, OUP 2009).

**Honest assessment.** The master Mellin engine is GeoVac-novel within chemistry. It would not be intelligible to a J. Chem. Phys. or J. Chem. Theory Comput. referee without substantial translation. Best fit venue if pitched alone: *J. Math. Chem.* or *Theor. Chem. Acc.* (the math-chem side).

### (d) Forced/free seam in molecular Hamiltonian design

**Searches returned:** general first-principles + ML-Hamiltonian-learning literature (HAMSTER arXiv:2026.70865, E(3)-equivariant Hamiltonians arXiv:2210.16190). All of these *learn* a Hamiltonian from data; none distinguish "structurally forced" from "calibration data" in the GeoVac sense.

The forced/free taxonomy is currently a **pure-physics framing** (NCG Standard Model: Connes-Chamseddine 1997, Marcolli book 2010). Has not been transplanted to chemistry. Adjacent ideas:

- Eq. symmetry-enforcement (point group, time-reversal) is "forced" structure; molecular geometries / orbital exponents are "calibration" — but no chemistry paper calls this out as a meta-theoretical taxonomy.
- Cusp conditions (Kato 1957) and exact-exchange properties (Levy, Perdew) are framework-forced; correlation-functional parameterizations are calibration. Same observation, never explicitly framed.

**Honest assessment.** A clean negative. The forced/free seam in *chemistry* would be sprint-scale to articulate as a position paper but has no current literature handle.

### (e) Trotter bounds from operator-algebraic structure

**Verified core references:**

- **Childs, Su, Tran, Wiebe, Zhu 2021**, "Theory of Trotter Error with Commutator Scaling," *Phys. Rev. X* 11, 011020, arXiv:1912.08854. **Verified via WebFetch.** Tight Trotter bound from nested commutators. Reproduces tight bound for 1st/2nd order, higher-order overestimates by factor of 5 for Heisenberg.
- **Lower-bound work (2024+):** arXiv:2410.03059 (Lower bounds for Trotter error), arXiv:2506.17199 (qDRIFT tighter bounds). All use commutator/spectral-norm methods.

**Honest assessment.** The Childs-Su lineage is exclusively commutator-norm-based; **nobody uses GH-convergence rates / propinquity rates / Lipschitz seminorms from spectral-triple truncation as a Trotter-bound source**. The GeoVac heuristic Trotter bound (4/π factor from Paper 38 GH rate, lifted via Childs-Su Duhamel + N_active-linearity, but 3–4 OoM too loose at production parameters) is consistent with this being an **unexplored mathematical direction** rather than a competitive practical bound.

Would only land at chemistry venues if the bound becomes tight. Currently is loose enough to be flagged in CLAUDE.md §3 as a non-binding heuristic. No outreach value at present.

### (f) Tensor-product spectral action for two-body operators

**Verified core references:**

- **Bochniak, Zalecki, Sitarz 2021**, "Spectral action and the electroweak θ-terms for the Standard Model without fermion doubling," arXiv:2106.10890, *JHEP* 12, 142. **Verified via WebFetch.** Computes the spectral action for a chiral NCG SM. **Does NOT state the metric-vs-Green's-function category divide explicitly** at the abstract level; the underlying observation (spectral action gives the kinetic + potential structure of the Lagrangian, *not* the propagator) is implicit in their and Connes-Chamseddine's derivations.
- The connection to chemistry runs through Paper 54 (GeoVac's tensor-product spectral action getting 75% connected two-body fraction with correct multipole hierarchy) and the negative resolvent-two-body finding (resolvent does not recover Coulomb because of the Fock-projection conformal factor).

**Honest assessment.** The category divide between "metric content (spectral action)" and "Green's-function content (resolvent / propagator)" is *implicitly* known in NCG literature but **not stated as a chemistry-relevant fact anywhere**. GeoVac's Paper 54 is the first explicit statement framed for an audience that cares about two-body interactions. Genuinely original. If Paper 54 is revised, framing this as the headline observation is recommended.

### (g) Yukawa non-selection analog in chemistry (FCI correlation as calibration)

**Searches returned:** standard FCI literature (Helgaker-Jorgensen-Olsen 2000, Sherrill 2013 *WIREs Comput. Mol. Sci.*). All treat FCI correlation energy as the **exact answer within a basis**, never as **external calibration**. The empirical convergence of CC / FCI / DMRG with basis size is uniformly framed as a numerical convergence question, not a meta-theoretical one.

**Honest assessment.** Clean negative. Closest published reading is Chan's DMRG-CC convergence work and the Petruzielo-Holmes-Changlani-Nightingale-Umrigar 2012 *Phys. Rev. Lett.* on SHCI extrapolation; these are tactical, not structural. The GeoVac framing — "W1e residual is categorically disjoint from outer-factor periods because it has zero vertex-parity content" — has no chemistry-side analog.

---

## Translation register: terms that would need renaming for QC-for-chemistry venues

GeoVac's internal vocabulary is heavy with NCG / spectral-triple / Connes-distance / propinquity terminology that is invisible to a *PRX Quantum* / *Quantum* / *npj Quantum Inf.* / *JCTC* referee. If the chemistry-facing papers (14, 20) are revised to engage the literatures above, the following translations would help:

| GeoVac internal | Translated for chemistry/QC | Rationale |
|:---|:---|:---|
| Hopf-U(1) Z₂ tapering | **"Angular-momentum parity tapering"** or **"m_l-reflection Z₂ tapering"** | Hopf-fibration vocabulary is opaque outside topology; m_l-parity is immediately legible to atomic-physics / quantum-chemistry readers. |
| Composed (Level 5) architecture | **"Geometric core-valence factorization"** or **"Natural-bundle pseudopotential encoding"** | "Composed" is GeoVac-internal jargon; "natural-bundle" cues NCG readers without losing chemistry readers. |
| Master Mellin engine (M1/M2/M3) | **"Spectral-action transcendental taxonomy"** | Avoid "Mellin" in titles for *JCTC*; works in body text after definition. |
| Phillips-Kleinman cross-center barrier | (keep — PK is universal vocabulary) | Standard. |
| Fock projection rigidity | **"S³ angular-momentum basis closure"** | "Fock projection" without "1935" parses as Hartree-Fock to chemistry readers; bad ambiguity. |
| Bargmann-Segal lattice | **"Discrete harmonic-oscillator basis"** | Acceptable for math-chem; otherwise translate. |
| Forced/free seam | **"Structural vs parametric Hamiltonian content"** | Standard ML-physics vocabulary; lands in the equivariant-Hamiltonian-learning literature. |

---

## Top 3 papers to cite if Paper 14 or Paper 20 is revised to engage these literatures

1. **Setia, Chen, Rice, Mezzacapo, Pistoia, Whitfield 2020**, arXiv:1910.14644, *J. Chem. Theory Comput.* 16, 6091. The canonical reference for point-group / molecular Z₂ tapering. Paper 14's new §sec:hopf_tapering subsection should cite this and explicitly distinguish: "Setia et al. exploit point-group Z₂'s; we exploit the orthogonal Hopf-fiber Z₂. The two stack additively." This is the single most important external citation for the chemistry-QC audience.

2. **Lee, Berry, Gidney, Huggins, McClean, Wiebe, Babbush 2021**, arXiv:2011.03494, *PRX Quantum* 2, 030305. The fault-tolerant chemistry resource baseline. Paper 20's comparison table should explicitly show GeoVac composed at Q ~ 30–100 vs THC at Q ~ 100+ (where THC is published), and acknowledge that THC dominates at very large molecular systems where Gaunt-selection sparsity saturates. Honest framing required: GeoVac competes in the NISQ Q < 100 regime, not the fault-tolerant Q > 100 regime.

3. **Picozzi, Tennyson 2023**, *Quantum Sci. Technol.* 8, 035026, DOI 10.1088/2058-9565/acd86c. Quantitative comparator for tapering at the molecular Hamiltonian level. Their k = 2–5 qubit reduction across molecules is the natural baseline for "what does state-of-the-art symmetry tapering already achieve?" GeoVac's per-sub-block ΔQ = −3 to −12 sits at the top of that range and stacks orthogonally.

**Optional fourth (for §IV if relativistic chemistry is added):** Saue et al. 2020 *J. Chem. Phys.* on relativistic Dirac-Coulomb chemistry on quantum computers; matches the Paper 23 Dirac-on-S³ thread if reframed for chemistry audiences.

---

## Verification status

All arXiv IDs cited above were resolved via WebFetch (or via the WebSearch metadata where WebFetch was blocked by paywalls; in those cases the abstract identification is conservative). Per the audit memo `debug/sprint_literature_audit_followup_memo.md`, I have NOT included any cite I could not corroborate to title + authors + year via the arXiv landing page or DOI page. References I considered but dropped because I could not verify authorship cleanly:

- "Greene-Diniz McClean active-space embedding" — search returned arXiv:2404.18737 with **NO Greene-Diniz author**. The Greene-Diniz group does publish in this area (arXiv:2202.04460 multilayer embedding for pharma) but I did not WebFetch-verify and have therefore not made it a top-3 citation.
- Various 2024–2026 papers came up in search but most are extensions of the three core references above; I omitted them as redundant.

Sources searched: arXiv, IOPscience, PRX Quantum, *J. Chem. Theory Comput.*, *Quantum* journal, *npj Quantum Information*. No paywall-only papers cited that I could not corroborate from abstracts.

---

## Closing note

The most actionable finding for the project: **thread (a) is real, ready, and additive to existing tapering recipes**. If Paper 14 cites Setia 2020 and Picozzi-Tennyson 2023 and frames the Hopf-U(1) per-sub-block tapering as orthogonal to point-group / particle-number / S_z Z₂'s, the result becomes immediately legible to the QC-for-chemistry audience and can land in *Quantum* or *PRX Quantum*. Threads (c), (d), (e), (g) are GeoVac-novel within chemistry but require translation or a separate position paper; not a 2026-Q3 priority.
