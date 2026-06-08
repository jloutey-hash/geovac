# Sprint Hybrid-Pipeline Scoping — canonical memo

**Date:** 2026-06-07.
**Sprint position:** Multi-month-roadmap scoping for a hybrid pipeline that pairs GeoVac (skeleton: sparse integrals + selection rules + qubit encoding) with an external correlation engine (CCSD(T) / DMRG / AFQMC / VQE) that consumes them to handle multi-determinant correlation — i.e. the W1e calibration tier (Class 1, Paper 18 §IV.6 chemistry-side analog) that GeoVac structurally cannot generate.
**Verdict line:** **GO** at Phase 1 scale. Top-2 first integrations are **(1) DMRG via Block2 / pyscf-DMRG** and **(2) VQE via Qiskit-Nature / OpenFermion-PySCF** in that order. The first-system pick is **LiH at the established `R = 3.015` benchmark** (production-stable, sub-percent reference, eight-orbital active space is the DMRG-sweet-spot), with NaH at `R_eq = 3.566` as the second-system stress test of the W1e tier specifically. Phase 1 is sprint-scale (~2 months) and is bounded by interface work, not by physics research.
**Cross-references:** `geovac/ecosystem_export.py` (current external interface), `papers/group4_quantum_computing/paper_14_qubit_encoding.tex` (qubit interface, O(Q^2.5) Pauli scaling claims), `papers/group4_quantum_computing/paper_20_resource_benchmarks.tex` (chemistry-consumer audience framing), `papers/group3_foundations/paper_31_universal_coulomb_partition.tex` (A/D split universal vs. Coulomb-specific), `debug/sprint_w1e_period_class_memo.md` (W1e structural classification), `memory/external_input_three_class_partition.md` (three-class boundary), CLAUDE.md §1.7 WH5 (α as projection constant, cosmic-Galois categorization), CLAUDE.md §1.5 (positioning), CLAUDE.md §3 (failed approaches: PK, screened-Schrödinger valence, multi-zeta, kernel-shape, bonding-PK, basis enlargement, explicit-core Hartree — all six W1e closure attempts).

---

## §0. Executive summary

The cosmic-Galois categorization (Paper 18 §III.7 master Mellin engine + §IV.6 inner-factor input-data tier) has classified W1e as Class 1 calibration data — categorically external to GeoVac's outer-factor mechanism (Sprint W1e period-class, 2026-06-04: 0/11 outer-factor identifications across all eleven NaH correction terms). Six W1e closure attempts (CLAUDE.md §3 PK cross-center, screened-Schrödinger valence, multi-zeta, kernel-shape, bonding-PK rank-1, basis enlargement, explicit-core Hartree) all failed at sprint scale for the same structural reason: they attempt outer-factor mechanisms on Class-1 calibration content.

Under this classification, the principled architecture is to STOP attempting W1e closure inside GeoVac and START a hybrid pipeline where GeoVac supplies the skeleton (sparse one-/two-electron integrals + Gaunt-enforced selection rules + qubit encoding with O(Q^2.5) Pauli scaling and per-block Hopf-U(1) tapering Δ_Q = 2 + n_sub_blocks) and an external correlation engine handles the multi-determinant N-electron problem on top.

This was always a possible move; under the cosmic-Galois reading it is now the *correct* move because it respects the structural-skeleton-scope boundary.

**Top-2 first integrations** ranked by interface cleanliness × chemistry-community-size × headline-publishability:
1. **DMRG (Block2 / pyscf interface).** Native consumer of (h1, eri) tuples; FCIDUMP file format is a 50-year-old standard; tested at the active-space sizes GeoVac produces; sub-percent reference accuracy achievable; venue-ready in JCP or JCTC.
2. **VQE (Qiskit-Nature / OpenFermion-PySCF).** Already half-integrated via `ecosystem_export.to_qiskit()` / `.to_pennylane()`; NISQ audience large; headline benchmark = pair-natural-orbital VQE on the GeoVac-tapered Hamiltonian; sub-block tapering paper (v3.52.0, Paper 14 §sec:hopf_tapering) is a natural lead-in.

**Deprioritized but kept on register:** CCSD(T) (interface clean but headline less novel — sub-percent achievable but only for closed-shell singletons that DMRG also handles), AFQMC (interface needs an FCIDUMP exporter; community size smaller than DMRG; phase-problem-handling adds risk; archived for Phase 2).

**Roadmap:** ~2 months Phase 1 (one consumer, one molecule, end-to-end pipeline + sub-percent reference); ~3-4 months Phase 2 (production sweep across 28 molecules with both consumers); ~2 months Phase 3 (chemistry paper, target JCTC / JCP, with arxiv mirror).

The roadmap is bounded by mechanical interface work and convention-pinning (orbital ordering, integral normalization, spin convention, CAS/active-space construction, Pauli/qubit-tapering integration), not by physics-research uncertainty. Phase 1 launches as soon as the PI signs off on first-consumer choice.

---

## §1. Interface inventory

### §1.1 What `ecosystem_export.py` currently exposes

`ecosystem_export.hamiltonian(system, R, max_n, core_method, tapered)` returns a `GeoVacHamiltonian` wrapper with:

- `.to_openfermion()` → OpenFermion `QubitOperator` (native; canonical)
- `.to_qiskit()` → Qiskit `SparsePauliOp` (with reversed qubit convention)
- `.to_pennylane()` → PennyLane `Hamiltonian`
- `.n_qubits`, `.n_terms`, `.one_norm`, `.one_norm_full` (Pauli sum metrics)
- `.h1_pk` — the PK one-electron integral matrix (M×M spatial) for composed systems
- `.pk_classical_energy(one_rdm)` — algebraic classical evaluation of the PK contribution from a measured 1-RDM
- `.metadata` (system, R, M, Q, N_pauli, tapering parameters if applied)

The system registry covers 28 molecules across three rows (He, H₂, LiH..HF, NaH..HCl, KH..HBr, ScH..ZnH, SrH, BaH, LiF, CO, N₂, F₂, NaCl) via three internal builders (`_build_hydride`, `_build_multi_center`, `_build_tm_hydride`) plus heavy-atom alkaline-earth monohydrides and `_build_he` / `_build_h2`.

The Hopf-U(1) `m → -m` Z₂ tapering (v3.52.0, Paper 14 §sec:hopf_tapering) is integrated as the `tapered={None, 'global', 'per_block'}` keyword. Per-block tapering removes ΔQ = 2 + n_sub_blocks across 37/37 molecules (range −3 for He/H₂, −12 for CO/N₂/F₂), saves 254 qubits library-wide, with bit-exact spectrum preservation on the 6 directly-verifiable cases.

### §1.2 What the four candidate consumers need (gap analysis)

| Consumer | Native input | Currently provided | Gap | Gap type |
|:-----------|:---------------|:---------------------|:----|:---------|
| **CCSD(T)** | (h1, eri) tuple in chemist or physicist notation + n_electrons + n_orbitals + symmetry labels; usually via FCIDUMP or pyscf `mf.mo_coeff/mo_energy/mo_occ` | Pauli `QubitOperator` only (after JW encoding) | **Add a `to_fcidump(filename)` exporter** that emits the un-JW'd integrals. The integrals exist internally inside `build_composed_hamiltonian` / `build_balanced_hamiltonian` before JW transformation; they just are not surfaced. | **Mechanical** (~2 days). The integrals exist; need to plumb them out. |
| **DMRG (Block2 / pyscf-DMRG)** | FCIDUMP or pyscf `SCF` object + bond dimension + active space | Same as CCSD(T) — needs (h1, eri) tuple, not Pauli | Same as CCSD(T): `to_fcidump`. Block2 reads FCIDUMP directly; no further conversion. | **Mechanical** (subsumed in `to_fcidump`). |
| **AFQMC (ipie / pauxy / QMCPACK)** | (h1, eri) tuple, often Cholesky-decomposed; needs trial wavefunction (HF or CASSCF) and walker count | Same as CCSD(T) + Cholesky of eri | `to_fcidump` + `to_cholesky_eri` (Cholesky is 1-2 days; existing `composed_qubit` uses analytical multipole decomposition which IS Cholesky in disguise — Paper 14 §sec:hopf_tapering's predecessor §df_multipole_lift_memo confirms double-factorization rank = multipole channel count). Also needs a trial-state generator (HF MO or CASSCF). | **Mechanical** for FCIDUMP+Cholesky; **moderate** for trial-state choice. |
| **VQE (Qiskit-Nature / OpenFermion-PySCF / PennyLane-Catalyst)** | `SparsePauliOp` (Qiskit), `QubitOperator` (OpenFermion), or `qml.Hamiltonian` (PennyLane) + initial state (HF or CASSCF) + ansatz (UCCSD / hardware-efficient / ADAPT) + classical optimizer | Already exposed: `.to_qiskit()`, `.to_openfermion()`, `.to_pennylane()` | **Trivial gap on Pauli side.** Real gap: VQE needs (a) an initial-state preparation circuit (Slater determinant on the GeoVac orbitals) and (b) an ansatz that respects the GeoVac sparsity (UCCSD on the active space + Gaunt-allowed singles/doubles). | **Mechanical** (initial-state) + **moderate** (ansatz selection). |

**Key finding:** the cleanest gap-to-close is **`to_fcidump` exporter**, which immediately unlocks CCSD(T), DMRG, AFQMC, and any pyscf-compatible solver. The interface for VQE is already mostly there; the gap is on the ansatz side, not the Hamiltonian side. Both first integrations (DMRG and VQE) are within sprint-scale reach.

### §1.3 Convention-pinning checklist (the load-bearing piece)

Before any consumer can verify its energy against GeoVac's internal FCI, the following conventions must be pinned and documented:

1. **Orbital ordering.** GeoVac internal: per-block (Z, n, l, m) lexicographic. Pyscf/FCIDUMP standard: (energy-sorted, but spatial only). Need a `geovac_to_pyscf_orbital_map` dict and a CI sanity check.
2. **Integral normalization.** GeoVac: chemist's notation `(pq|rs) = ∫ φ_p(1) φ_q(1) (1/r12) φ_r(2) φ_s(2)`. FCIDUMP standard: same (chemist's). Confirm both use Hartree, not eV / a.u.²/Ha.
3. **Spin convention.** GeoVac qubits: α spin orbital indices = 2k, β = 2k+1 (OpenFermion standard). pyscf RHF: spatial only (spin is implicit in occupation 0/2). UHF/UCCSD: separate α/β blocks. The map from GeoVac's `n_spin_orbitals = 2M` to pyscf's `n_orbitals = M` + `nelec = (n_α, n_β)` is one-liner but must be tested per molecule.
4. **Active space / CAS.** GeoVac builds a fixed active space (max_n=2 valence: 1s, 2s, 2p for first-row; matching for higher rows via frozen-core analytical Z_eff(r)). The frozen-core energy enters as a constant (the identity term in the Pauli operator). DMRG / AFQMC / CCSD(T) need this constant passed alongside the integrals, or it can be absorbed into the orbital energies + nuclear repulsion. Pinning: pass as `frozen_core_energy` scalar; CCSD(T) result + `frozen_core_energy` = total energy.
5. **Hopf-U(1) tapering and the post-tapering integral interpretation.** Tapered Hamiltonians lose qubits but the spectrum is preserved. CCSD/DMRG/AFQMC are agnostic to tapering (they consume the un-JW'd integrals); VQE consumes the tapered Pauli operator natively. Pinning: for non-VQE consumers, pass un-tapered (h1, eri); for VQE, pass tapered Pauli.
6. **Reference state for correlation energy.** CCSD/DMRG report correlation energies relative to HF; GeoVac internal FCI reports absolute energies. Need a `geovac_hf_reference_energy(spec)` helper that runs a Hartree-Fock on the GeoVac integrals and returns the reference. (This can use pyscf's HF on the FCIDUMP, which is the most consistent route.)

Items 1-6 are mechanical but cumulative: the Phase 1 first sprint must close ALL six on the chosen first system, with documented invariants.

---

## §2. Three-class partition for chemistry — which consumer addresses what

Per `memory/external_input_three_class_partition.md` and the W1e diagnostic memo, the W1e wall splits into three sub-categories: (A) external input data — Clementi-Raimondi exponents for the frozen core; (B) algebraic-implicit content — the FCI eigenvalues themselves; (C) experimental observation — measured D_e for benchmarking. The four candidate consumers address these sub-categories differently:

| Consumer | Class 1 (calibration / W1e data) | Class 2 (multi-focal composition) | Class 3 (multi-determinant correlation) |
|:---------|:--------------------------------|:--------------------------------|:---------------------------------|
| **CCSD(T)** | Does **not** improve the [Ne] core fits (consumes the same Z_eff(r) GeoVac provides). | N/A — single-reference perturbative; not a multi-focal compositor. | **Yes** for systems where the reference is single-determinant (closed-shell ground states); the "(T)" perturbative triples capture dynamical correlation that exceeds FCI on the active space when the active space is large enough. Fails for strongly-correlated systems (multi-reference character). |
| **DMRG** | Same — does not improve calibration data. | N/A. | **Strongest of the four** for the W1e tier: handles arbitrary multi-determinant correlation up to FCI within a chosen bond dimension; converges to FCI in the bond-dimension → ∞ limit; pyscf-DMRG and Block2 both production-stable. Sweet spot is exactly the GeoVac active-space sizes (15-50 spatial orbitals). |
| **AFQMC** | Same. | N/A. | Strong for multi-determinant correlation but adds phase-problem risk and statistical noise; less deterministic than DMRG. Better for very large active spaces (DMRG bond dimension prohibitive). |
| **VQE** | Same — calibration data passed through unchanged. | N/A. | Approximate-correlation handler; quality depends entirely on ansatz. UCCSD on GeoVac is the obvious first test; hardware-efficient ansatz is the NISQ-ready option. **In principle reaches FCI, in practice noise-limited.** |

**The W1e tier (Class 1 calibration data) is NOT addressed by ANY of the four consumers.** This is correct per the cosmic-Galois reading: Class 1 is external to both GeoVac AND the consumers, by design. The frozen-core Z_eff(r) profiles, Clementi-Raimondi exponents, and multi-zeta coefficients are atomic-physics calibration data that ALL methods consume.

**Class 3 (multi-determinant correlation) is what the consumers actually do.** This is the W1c-residual correlation physics that GeoVac's spectral-triple machinery does not natively engage (because spectral triples live on the algebra / Dirac side; correlation lives on the N-representable density-matrix manifold). DMRG and CCSD(T) are the gold standards; AFQMC and VQE are the second-tier options.

**Class 2 (multi-focal composition — recoil cross-register, Zemach magnetization-density) is NOT chemistry's domain.** None of the four consumers handle Class 2 walls. Class 2 closures live in GeoVac itself via architectural extension (W1b operator-level Zemach extension is the precedent).

**Net partition reading:** the hybrid pipeline targets Class 3 correlation cleanly. Class 1 calibration is pass-through (consumers and GeoVac share the same atomic-physics inputs). Class 2 is orthogonal.

---

## §3. First-integration prioritization

### §3.1 Ranking criteria

(a) **Interface cleanliness given current `ecosystem_export.py`** — how much new plumbing is needed?
(b) **Chemistry-community size** — who reads the resulting paper?
(c) **Headline benchmark possible** — what number can we publish?

### §3.2 Per-consumer scores

| Consumer | (a) Interface cleanliness | (b) Community size | (c) Headline benchmark | Composite |
|:---------|:--------------------------|:--------------------|:----------------------|:----------|
| **DMRG (Block2/pyscf)** | Needs `to_fcidump` (~2 days). Block2 then reads it directly with zero further integration. Cleanest non-Pauli interface. | Large (~10⁴ practitioners) and growing. Standard reference method in modern QC. | LiH at R_eq sub-percent vs. FCI/cc-pVTZ reference; NaH W1e tier closure quantification; library sweep across 28 molecules with DMRG-vs-GeoVac-FCI cross-validation. **Strong headline:** "GeoVac integrals + DMRG correlation = production-grade quantum chemistry at 10× fewer Pauli terms than Gaussian DMRG inputs." | **#1** |
| **VQE (Qiskit-Nature)** | Already ~80% integrated. Needs initial-state circuit + UCCSD ansatz. Real work but no new exports. | Very large (NISQ + quantum hardware audience). | "First demonstration of sub-block-tapered GeoVac VQE on LiH with N_qubits-3 = 27 active qubits"; sub-percent on small systems via shot-budgeted noise model; demonstrates compatibility with all major quantum-hardware SDKs. **NISQ-publishable headline.** | **#2** |
| **CCSD(T)** | Needs `to_fcidump` (subsumed in DMRG work). | Very large (CCSD(T) is the gold-standard wall in quantum chemistry). | Sub-percent on LiH/H₂O/etc., but CCSD(T)-on-GeoVac-integrals is not as novel as DMRG-on-GeoVac-integrals (CCSD(T) is widely available, headline is harder to land). Useful as a cross-validation reference for DMRG. | **#3 (consolation track)** |
| **AFQMC** | Needs `to_fcidump` + `to_cholesky_eri` + trial-state generator. Additional 1-2 weeks of work over DMRG. | Smaller community (~10³ practitioners) but growing. | Strong on large-active-space molecules where DMRG bond dimension is prohibitive (TM hydrides, F₂); phase-problem-handling adds risk to first-Phase-1 attempt. **Phase-2 target, not Phase-1.** | **#4** |

### §3.3 Top-2 picks with reasoning

**Pick #1: DMRG.** The interface gap is the smallest among non-Pauli consumers; the community is large; the headline is strong (production-grade quantum chemistry via the structurally-sparse interface). DMRG is also the cleanest way to validate the W1e structural classification: if DMRG-on-GeoVac-integrals reaches sub-percent on NaH, this *proves* the W1e wall was correctly attributed to multi-determinant correlation (Class 3) rather than GeoVac's calibration data (Class 1). This is a Paper 18 §IV.6 falsifiability check at production scale.

**Pick #2: VQE.** Already mostly integrated; the per-block Hopf-U(1) tapering work (v3.52.0) is a natural lead-in (tapered VQE on real hardware is publishable on its own); the NISQ audience is the natural reader for the Paper 14 / Paper 20 narrative. Lower headline-accuracy ambition than DMRG (VQE on current hardware is shot-noise-limited at ~1% on small molecules), but high *positioning* impact (demonstrates that GeoVac is hardware-ready out of the box).

**Sequence:** start with DMRG, layer in VQE in Phase 2. Phase 1 deliverable is a sub-percent DMRG-on-GeoVac-integrals number on LiH with full convention-pinning; Phase 2 deliverable is the same plus VQE-on-tapered-Hamiltonian on the same system.

---

## §4. Multi-month roadmap

### §4.1 Phase 1 (~2 months): prototype against DMRG on LiH

**Goal:** end-to-end pipeline from `ecosystem_export.hamiltonian('LiH', R=3.015)` → Block2 / pyscf-DMRG input → sub-percent energy. Bounded by interface plumbing and convention-pinning, not by physics-research uncertainty.

**Sprint sequence (suggested):**

- **Sprint H1-Interface (2-3 weeks).** Implement `to_fcidump(filename)` on `GeoVacHamiltonian` and a paired `from_geovac_to_pyscf(ham) → pyscf.scf.RHF` helper. Pin convention items 1-5 of §1.3. Sanity check: pyscf RHF on the FCIDUMP recovers the GeoVac HF reference to machine precision.
- **Sprint H1-DMRG (1 week).** Run Block2 / pyscf-DMRG on the LiH integrals at increasing bond dimension D ∈ {50, 100, 200, 500, 1000}. Expect convergence to GeoVac internal FCI at D ~ 200 for LiH-sized active space. Quantify the energy at each D.
- **Sprint H1-Reference (1 week).** Run pyscf cc-pVTZ FCI (or DMRG with cc-pVTZ + 8-electron active space) for the gold-standard external reference. Compare GeoVac-FCI, GeoVac-DMRG, and external cc-pVTZ-FCI on the same LiH benchmark.
- **Sprint H1-NaH stress test (1 week).** Repeat the entire pipeline on NaH at R_eq = 3.566. Quantify the W1e wall closure when DMRG handles the multi-determinant correlation. Expected outcome: significant W1e closure (consistent with the Class-3 attribution). Falsifier: NaH remains overattracted even with DMRG → wall is NOT pure Class 3 → diagnostic implication.
- **Sprint H1-Memo (~3 days).** Standard canonical memo + CHANGELOG entry; Paper 20 §sec:hybrid_pipeline new subsection drafted (do not apply until Phase 2 close).

**Exit criterion:** sub-percent LiH energy via DMRG-on-GeoVac-integrals at production basis (n_max=2, max_n=2), and a structurally-clean NaH wall-closure result (either positive — confirming W1e = Class 3 — or negative — flagging Class-3 as itself partitioned).

### §4.2 Phase 2 (~3-4 months): expand to library sweep + VQE second integration

**Goal:** production pipeline across the 28-molecule library; VQE integration; first hardware-aware benchmarks.

**Sprint sequence (suggested):**

- **Sprint H2-Sweep (~6 weeks).** Run the Phase-1 DMRG pipeline across 28 molecules in the library. Tabulate energy error vs. external cc-pVTZ-FCI reference where available, vs. CCSD(T) where FCI is prohibitive. Identify systematic patterns (does the W1e tier scale linearly with electron count? are there outliers — N₂ multi-reference character?).
- **Sprint H2-CCSD(T) (~2 weeks).** Add CCSD(T) cross-validation track. Easy at this point (`to_fcidump` already exists); produces an independent reference for closed-shell systems where DMRG and CCSD(T) should agree.
- **Sprint H2-VQE (~4-6 weeks).** Phase-2 VQE integration: (a) initial-state preparation circuit (Slater determinant on GeoVac orbitals via OpenFermion + Qiskit-Nature), (b) UCCSD ansatz on the active space with Gaunt-allowed singles/doubles only, (c) noise model + shot budget for the per-block tapered Hamiltonian. Headline: tapered-VQE on LiH at N_qubits=27 with sub-1% accuracy under realistic noise. Pursue PennyLane integration in parallel for the differentiable-programming audience.
- **Sprint H2-W1e-Catalogue (~2 weeks).** With Phase-1 + sweep data, quantify the W1e tier across all 28 molecules: how much of the residual is closed by DMRG/CCSD(T) (the Class-3 correlation), how much remains as calibration mismatch (Class-1)?

**Exit criterion:** complete benchmark library for the chemistry paper; W1e tier characterized across the molecular library; VQE first hardware-aware benchmark on LiH.

### §4.3 Phase 3 (~2 months): chemistry paper

**Goal:** publishable chemistry paper integrating Phase 1+2 results.

**Suggested structure:**

- §I Introduction (sparse skeleton + external correlation engine architecture)
- §II GeoVac integrals (brief — pointer to Paper 14, Paper 17, Paper 22)
- §III Interface design (FCIDUMP exporter, convention pinning, tapered VQE)
- §IV Benchmarks (28-molecule sweep, LiH/NaH stress tests, VQE hardware-aware)
- §V W1e characterization (the structural-skeleton-scope finding made empirical: DMRG closes the multi-determinant correlation, GeoVac supplies the calibrated skeleton, ALL methods inherit the Class-1 calibration data)
- §VI Comparison to alternative architectures (Gaussian DMRG, Gaussian VQE, GeoVac-FCI)
- §VII Open questions

**Estimated venue:** JCTC (J. Chem. Theory and Comput.) for the main chemistry-community pitch. JCP (J. Chem. Phys.) is a fallback. Quantum or PRX Quantum for a VQE-focused spinoff (Phase-2's H2-VQE result is potentially its own paper if it lands clean). Per CLAUDE.md §1 ("no journal-track language unless asked"), positioning is for the project's institutional record + arxiv mirror; the venue is the PI's call.

### §4.4 Total timeline and budget

- Phase 1: ~2 months wall (DMRG plumbing + LiH/NaH benchmarks)
- Phase 2: ~3-4 months wall (sweep + CCSD(T) cross-validation + VQE)
- Phase 3: ~2 months wall (chemistry paper)
- **Total: ~7-8 months end-to-end** to publishable chemistry paper.

This is sprint-scale at the Phase-1 level (any Phase-1 sub-sprint is ~1-3 weeks); Phase 2 is multi-sprint accumulation; Phase 3 is a single sustained writing sprint.

---

## §5. Open questions for PI

Listed in order of decision-urgency. PI input needed before Phase 1 launches:

1. **First-consumer choice: DMRG or VQE?** The memo recommends DMRG (#1 in §3.3) because the interface gap is smallest and the headline is strongest. VQE (#2) is the alternative if NISQ-positioning is more important than sub-percent accuracy. **PM default if no input: DMRG.**

2. **First-molecule choice: LiH or NaH?** LiH (R=3.015) is the established production benchmark; NaH (R_eq=3.566) is the W1e stress test. Memo recommends LiH first, NaH second within Phase 1. **PM default if no input: LiH first, NaH as Sprint H1-NaH.**

3. **Add `block2` as an optional dependency, or shell out to a separate pyscf environment?** Block2 is pip-installable but ~50MB binary; pyscf is ~100MB. Adding both as optional dependencies follows the OpenFermion/Qiskit/PennyLane precedent. **PM default if no input: add as optional dependencies, mirror the `_HAS_OPENFERMION` sentinel pattern.**

4. **Where does the chemistry paper sit in the corpus?** Paper 20 (resource benchmarks, chemistry-consumer audience) is the natural home for the integration writeup; the W1e structural finding could either go in Paper 20 §sec:hybrid_pipeline or as a standalone "GeoVac + DMRG" methods paper. Memo recommends extending Paper 20 (no new paper needed; the integration story is the natural Paper 20 §V or §VI). **PM default if no input: Paper 20 extension.**

5. **Sprint cadence: Phase 1 as one continuous arc, or staged with PI approval gates between H1-Interface, H1-DMRG, H1-Reference, H1-NaH?** The full Phase 1 is ~2 months which is long for a single autonomous sprint; the natural gate is after H1-Interface (FCIDUMP works, conventions pinned) before launching DMRG. **PM default if no input: gate after H1-Interface; the LiH-via-DMRG sub-sprint is launched only after PI signs off on the interface design.**

6. **Falsifier handling for the W1e structural classification.** If Sprint H1-NaH returns DMRG-on-GeoVac-integrals NaH energy that is *not* sub-percent (i.e. W1e wall is NOT closed by external correlation), the cosmic-Galois reading is partially falsified — W1e would not be pure Class 3, but a mix of Class 3 + something else (perhaps a hidden Class 2 cross-shell composition wall that GeoVac is also missing). PI direction: should the NaH stress test be a Phase-1 *exit* criterion (block Phase 2 until resolved), or a Phase-1 *finding* (record the partial falsification, proceed to Phase 2 with explicit caveat)? **PM default if no input: finding-not-exit; record in Paper 18 §IV.6 chemistry-side analog and proceed with caveat. The structural classification is sharper than its first numerical test.**

---

## §6. Honest scope

What this scoping demonstrated:

- The current `ecosystem_export.py` interface is ~80% ready for VQE (`to_qiskit`, `to_pennylane`, `to_openfermion` already exist with metadata) and ~40% ready for DMRG / CCSD(T) / AFQMC (needs `to_fcidump` to surface the un-JW'd integrals, which exist internally inside `build_composed_hamiltonian` and `build_balanced_hamiltonian`).
- The W1e cosmic-Galois classification (Sprint W1e period-class, 2026-06-04) implies that external correlation engines are the structurally-correct closure path for the W1e tier. Six prior in-framework attempts failed at sprint scale; this is the principled architectural move.
- DMRG is the top-ranked first consumer by interface cleanliness × community size × headline-publishability. VQE is the second-ranked first consumer due to existing partial integration and NISQ-audience size. CCSD(T) and AFQMC are Phase-2 / cross-validation tracks.
- Phase 1 (~2 months) is sprint-scale at the sub-sprint level; full Phase 1 is a multi-sprint arc. Phase 2 (~3-4 months) is sweep-scale work. Phase 3 (~2 months) is sustained writing. Total ~7-8 months to publishable chemistry paper.

What this scoping did NOT demonstrate:

- **No code was written.** This is a pure scoping memo per the sprint mandate. `to_fcidump` is sketched but not implemented; convention-pinning items are listed but not pinned; DMRG is not yet tried.
- **No FCIDUMP convention has been verified.** The integrals exist internally but the orbital ordering, normalization, and frozen-core handling have not been tested against external consumers. The first Sprint H1-Interface action item is exactly this.
- **No NaH W1e closure was attempted.** The claim "DMRG handles the multi-determinant correlation that GeoVac structurally cannot generate" is a hypothesis grounded in the W1e cosmic-Galois classification, not an empirical result. Sprint H1-NaH is the load-bearing falsifiability test.
- **No Paper 20 §sec:hybrid_pipeline draft was written.** Drafting is Phase 3 work; the section is named here as the planned home but not drafted in this sprint.
- **No CCSD(T) / AFQMC roadmap details.** These are Phase-2 / consolation tracks; the scoping addresses them only at the level of "interface gap is the same as DMRG (`to_fcidump`)." Detailed Phase-2 roadmap is a follow-on scoping sprint once Phase 1 lands.

Decision-gate outcome (memo's own verdict): **GO — Phase 1 is sprint-scale and bounded by interface plumbing, not physics-research uncertainty. PI input needed on the six items in §5 before launch.**

---

## §7. Follow-on register

| Priority | Item | Cost | Status |
|:--:|:-----|:----|:-------|
| 1 | PI sign-off on §5 items 1-6 (first consumer, first molecule, dependency policy, paper home, sprint gating, falsifier handling) | ~30 min PI time | NAMED gate |
| 1 | Sprint H1-Interface: `to_fcidump(filename)` + convention-pinning items 1-5 + pyscf RHF cross-check | ~2-3 weeks | NAMED Phase 1 launch sprint |
| 2 | Sprint H1-DMRG: Block2 / pyscf-DMRG on LiH at increasing bond dimension; sub-percent vs FCI external reference | ~1 week | NAMED follow-on (gated on H1-Interface) |
| 2 | Sprint H1-NaH: stress test of the W1e structural classification | ~1 week | NAMED falsifier sprint (gated on H1-DMRG) |
| 3 | Sprint H2-Sweep: 28-molecule DMRG library sweep | ~6 weeks | NAMED Phase 2 target |
| 3 | Sprint H2-VQE: tapered-VQE on LiH with UCCSD ansatz and realistic noise model | ~4-6 weeks | NAMED Phase 2 target (parallel to H2-Sweep) |
| 4 | Sprint H3-Paper: Paper 20 §sec:hybrid_pipeline + benchmark tables + W1e characterization | ~2 months | NAMED Phase 3 deliverable |
| 5 | (Phase 2 follow-on) AFQMC integration via `to_cholesky_eri` + trial-state generator | ~3-4 weeks | DEFERRED to post-Phase-2 |
| 5 | (Phase 2 follow-on) CCSD(T) cross-validation track via pyscf | ~2 weeks | DEFERRED but trivial once `to_fcidump` exists |

---

## §8. Files

### Created (this sprint)
- `debug/sprint_hybrid_pipeline_scoping_memo.md` (this memo).

### NOT modified
- Production `geovac/` modules — pure scoping sprint per sprint mandate.
- `geovac/ecosystem_export.py` — no `to_fcidump` implementation yet; this is Sprint H1-Interface work.
- Papers — Paper 20 §sec:hybrid_pipeline subsection is named as the planned home but not drafted.
- CLAUDE.md — no edits; scoping memo only.
- Tests — no code changes; no tests added.

---

**End of Sprint Hybrid-Pipeline scoping memo. Verdict: GO at Phase 1 scale; top-2 first integrations are DMRG (Block2/pyscf) and VQE (Qiskit-Nature/OpenFermion); first system LiH at R=3.015, second system NaH at R_eq=3.566 as W1e structural-classification stress test; total ~7-8 month roadmap to publishable chemistry paper; six PI-input items in §5 gate the launch.**
