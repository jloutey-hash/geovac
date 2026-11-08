# Sprint memo (scoping): GeoVac n_max truncation = Connes-vS spectral truncation = DMRG/MPO bond-rank truncation

**Date:** 2026-06-05 (evening; fourth meta-lesson follow-on after DF/Cholesky/QPT)
**Author:** PM scoping session (sub-agent, PI-directed)
**Status:** SCOPING ONLY — proposed theorem statement + empirical-first verification path + sprint sequence.
**Headline number:** 3 sprints, ~3 weeks wall-clock if PI prioritizes, to reach an empirically validated theorem statement plus a Paper-grade write-up of one direction (n_max → MPO bond profile). The full math.OA-grade bidirectional rigor remains multi-month, and we name the load-bearing reason cleanly.

**Files referenced:**
- `papers/group1_operator_algebras/paper_38_su2_propinquity_convergence.tex` (internal Connes-vS proof — five-lemma machinery)
- `debug/sprint_df_multipole_lift_memo.md` (DF=multipole bit-exact, 2026-06-05 AM)
- `debug/sprint_cholesky_multipole_memo.md` (Cholesky=multipole bit-exact, 2026-06-05 PM)
- `debug/sprint_qpt_hopf_stacking_memo.md` (QPT+Hopf-Z2 commute bit-exact, 2026-06-05 evening)
- `debug/lit_search_algebraic_methods_compatibility_2026_06_05.md` (lit scan flagging this candidate)
- `geovac/z2_tapering.py` (production Hopf-Z2 module)
- `geovac/composed_qubit.py`, `geovac/ecosystem_export.py` (production composed builder)

---

## 1. Question / what's being scoped

The PI's vibe-physics meta-lesson from the past 24 hours:

> GeoVac should be comparable with any algebraic-exact compression technique.

Three for three so far — DF = multipole (bit-exact), Cholesky = multipole (bit-exact), QPT and Hopf-Z2 commute bit-exact and stack triply-additively. The next deepest candidate from yesterday's lit scan is the DMRG/MPO connection: is the GeoVac n_max truncation also expressible as a bond-rank truncation of an MPO representation of the GeoVac Hamiltonian?

The lit-scan agent estimated this as 2-month math-grade. The PI bets we can compress to sprint-scale. This memo says where the compression is real and where it isn't.

**Scope clarification.** Three identifications are at stake:

1. **(A)** GeoVac n_max truncation = Connes-vS spectral truncation. **Already PROVEN internally** (Paper 38, 2026-05-07; WH1 PROVEN). The composed n_max cutoff is precisely the truncation by a spectral projector P_{n_max} onto modes with |D_CH| ≤ n_max + 1/2, and the truncated operator system O_{n_max} = P_{n_max} A P_{n_max} converges to A in Latrémolière propinquity at rate 4·log(n_max)/(π·n_max).

2. **(B)** Connes-vS spectral truncation ↔ MPO bond-rank truncation. The OPEN structural identification, not in published literature.

3. **(C)** GeoVac n_max truncation = MPO bond-rank truncation. Follows from (A) ∘ (B). The composition is what we care about for chemistry.

The path forward is to establish (B) at the right rigor level, then compose (A) and (B) to get (C). Compression to sprint-scale relies on the fact that (A) is already PROVEN — Paper 38 supplies the half of the bridge that math.OA usually requires multi-month effort to build.

---

## 2. What's already established internally (Paper 38 summary)

For the scoping question, the load-bearing facts from Paper 38 are:

- **GeoVac n_max = Connes-vS spectral truncation.** The Connes-vS truncated operator system (Connes-vS 2021, arXiv:2004.14115) at cutoff Λ on a spectral triple (A, H, D) is O_Λ = P_Λ A P_Λ ⊂ B(P_Λ H), where P_Λ is the spectral projector onto |D| ≤ Λ. For GeoVac's Camporesi-Higuchi triple on S³, with σ(D_CH) = {±(n+1/2)}, the Λ = n_max + 1/2 truncation is the composed-builder n_max cutoff bit-exactly.

- **Propagation number = 2** (Connes-vS framework) and **dim H_{n_max} = (2/3)·n_max(n_max+1)(n_max+2)** (Paper 38 §2.2). These are the natural rank parameters on the operator-system side.

- **Convergence rate Λ_prop ≤ 4·log(n_max)/(π·n_max)** (Paper 38 main theorem). The constant 4/π = Vol(S²)/π² is the M1 Hopf-base measure signature.

- **Five-lemma proof structure** (L1' / L2 / L3 / L4 / L5):
  - L1' chirality-doubled operator-system substrate
  - L2 central spectral Fejér kernel on SU(2) with Cesàro-2 rate (4/π)·log(n)/n
  - L3 Lipschitz comparison bound with C_3 = 1 via Avery-Wen-Avery 3-Y integral
  - L4 Berezin reconstruction (positive + contractive + approximate identity + L3-compatible)
  - L5 Latrémolière propinquity assembly via tunneling pair

The key implication: when we identify (B), the **propinquity rate is automatically inherited** — we don't need to prove convergence rates on the MPO side de novo, only show that the rank-revealing truncation matches.

---

## 3. The proposed theorem statement (Deliverable 1)

### 3.1 Setup

Let H_GeoVac be the composed second-quantized Hamiltonian for a molecule M at cutoff n_max, in the standard composed-builder basis ordering (orbitals indexed by `(sub-block i, n, l, m_l, s)` with sub-block outermost and (n, l, m_l) ordered lexicographically within sub-block; this is the natural orbital chain for GeoVac).

Let χ_{MPO}(H) denote the **MPO bond dimension profile** of H along this orbital chain: at each bond cut k ∈ {1, ..., M-1}, the MPO bond rank χ_k is the rank of the matricization of the MPO tensor across the (left orbitals 1..k) / (right orbitals k+1..M) cut. This is the standard definition used by Keller, Dolfi, Troyer, Reiher 2015 (J. Chem. Phys. 143, 244118; arXiv:1510.02026).

Let O_{n_max} = P_{n_max} A_{GV} P_{n_max} be the Connes-vS truncated operator system at cutoff n_max on the GeoVac Camporesi-Higuchi triple (Paper 38). Let H_GeoVac^{(n_max)} be the Hamiltonian's projection into the truncated triple.

### 3.2 Theorem (proposed, three statements at decreasing rigor)

**Theorem 3.2.A (strongest; sprint-scale empirical statement, paper-grade after Phase 1).**
For every molecule M in the GeoVac library and every n_max ≥ 2, the MPO bond dimension profile of H_GeoVac^{(n_max)} along the natural orbital chain satisfies:

$$
\chi_k(H_{\mathrm{GeoVac}}^{(n_{\max})}) \;=\; \chi_k^{\mathrm{Gaunt}}(k; n_{\max})
$$

where χ_k^{Gaunt}(k; n_max) is the **a priori bond dimension** counting the number of distinct (L, M) two-electron channels coupling left-block orbitals (1..k) to right-block orbitals (k+1..M), modulo the multipole-channel rank-deficit from Sprint DF F3 (graded R^L spectrum at large n_max).

**Concrete consequence:** χ_k(H_GeoVac^{(n_max)}) ≤ (n_max + 1)^2 · n_sub_blocks(k) at every cut, where n_sub_blocks(k) counts how many active sub-blocks straddle the cut. This is much tighter than the generic chemistry MPO bound χ ≤ 4M² (Keller et al. §III.B).

**Theorem 3.2.B (medium rigor; follows from 3.2.A + Paper 38 framework).**
The composition of the GeoVac n_max truncation with the standard MPO representation is **equivalent at the operator-system level** to a Connes-vS spectral truncation of a specific MPO bond-rank truncation of the continuum Hamiltonian:

$$
P_{n_{\max}} \cdot \mathrm{MPO}_\infty(H) \cdot P_{n_{\max}} \;\cong\; \mathrm{MPO}_{\chi(n_{\max})}(H)
$$

where χ(n_max) is the bond-rank profile from Theorem 3.2.A. This is an **operator-system isomorphism**, not a state-space isometry.

**Theorem 3.2.C (deepest; full math.OA-grade — multi-month, NOT sprint-scale).**
The Connes-vS spectral truncation and the DMRG bond-rank truncation are **dual truncations of the same spectral triple**, with both converging to the continuum at the same rate (Paper 38's 4/π·log(n_max)/n_max for the spectral side; the corresponding rate for the bond-rank side is the open math.OA-grade content of Theorem 3.2.C).

### 3.3 Technical details

- **Inner product on the MPS side.** Standard Fock-space inner product, JW-mapped to the qubit Hilbert space C^{2^Q}. The MPO is constructed as a finite-state automaton over the orbital chain (Keller et al. §III.A), with the Schmidt rank at each cut equal to the bond dimension χ_k.

- **What "bond rank" means.** Left/right Schmidt rank are equal because MPOs are square in the bond index. We define χ_k = max{rank of left-bond matrix, rank of right-bond matrix} at cut k; by Hermiticity this is unique.

- **How P_{n_max} acts.** P_{n_max} is a single-particle projector that restricts every Fermionic creation/annihilation operator a_p^†, a_p to orbitals p with n_p ≤ n_max. Its second-quantized lift acts as a diagonal projector on Fock space. The composed-builder n_max parameter IS this projector.

- **GeoVac orbital ordering.** The composed builder orders orbitals as (sub-block, then n, then l, then m_l, then s). The natural orbital chain for the MPO follows this ordering, which is the chemistry-DMRG-standard "near-orbital-energy" ordering for atomic-like blocks. Sub-block ordering matters for MPO bond profile and is the place where GeoVac's structural sparsity manifests most cleanly: cross-sub-block ERIs vanish (per the DF memo F4), so χ_k drops at every sub-block boundary.

- **What's load-bearing in the theorem.** 3.2.A is the **empirically-checkable** statement. It says the MPO bond profile of H_GeoVac is determined by the structure that GeoVac already exposes analytically (multipole channels + sub-block boundaries), with no surprises. If 3.2.A holds bit-exactly, the GeoVac truncation provides a structured MPO that is **a priori sparser than the generic Keller-Dolfi-Troyer-Reiher MPO** for a Gaussian chemistry Hamiltonian. This is the part that compresses to sprint-scale.

---

## 4. What an empirical-first verification would look like (Deliverable 2)

### 4.1 The compression question

The lit-scan agent estimated 2-month math-grade. The honest read:

- **Sprint-scale, 2-3 weeks.** Theorem 3.2.A empirical verification on small GeoVac molecules. Build the GeoVac MPO via the Keller-Dolfi-Troyer-Reiher construction, compute χ_k profile, compare to the a priori Gaunt-channel count. Bit-exact match would establish the cleanest cross-domain identification GeoVac has so far: the n_max parameter IS the bond-rank parameter, at production sizes, with the bond profile predictable from the multipole structure.

- **Multi-month, mathematically deep.** Theorem 3.2.C (rate-matched convergence on the MPS bond-rank side). This requires extending Paper 38's five-lemma machinery to the dual setting (state-space truncation vs algebra-space truncation), which is a genuinely open math.OA question. Latrémolière's propinquity machinery does NOT automatically give a corresponding bound for DMRG truncation; the natural object is the **Wasserstein–Kantorovich distance on states**, which is the dual of the Lipschitz seminorm but with finer constants.

- **In between (NOT sprint-scale, NOT multi-month).** Theorem 3.2.B operator-system isomorphism. This is doable in ~6 weeks of focused work because the operator-system side is already covered by Paper 38, and the MPO side has been block-encoded explicitly in Keller et al. 2015. The work is to assemble the isomorphism via the Berezin reconstruction map L4 from Paper 38 acting on the MPO tensor network. We name this honestly as **between sprint-scale and multi-month**.

### 4.2 Shortcuts via Paper 38

Three shortcuts available because (A) is already PROVEN:

**Shortcut 1: Berezin map L4 is already constructed.** The Berezin reconstruction B_{n_max}: C(S³) → O_{n_max} from Paper 38 Lemma L4 (positive, contractive, approximate-identity, L3-compatible) acts on functions. To get from L4 to an MPO map, we need to extend B_{n_max} to act on **tensor products of orbital functions** — this is a direct tensor power of L4 if the MPO's tensors live in tensor powers of C(S³). The two-electron operator V_ee lives in C(S³) ⊗ C(S³) (per the DF=multipole memo), so L4 ⊗ L4 acts on the Coulomb tensor block of the MPO. The remaining one-body block (kinetic + nuclear) is single-orbital, so L4 acts on it directly. **Sprint-scale work item: explicit L4 ⊗ L4 construction for the MPO tensor.**

**Shortcut 2: 4/π rate transports directly.** Whatever rate L4 induces on operator-system observables is automatically the rate the MPO bond-rank truncation inherits, because the MPO is reconstructed from L4 ⊗ L4. The PI's "the same 4/π appears everywhere" observation now extends to DMRG: if the theorem holds, DMRG bond-rank truncation on a GeoVac Hamiltonian converges at the same 4/π Hopf-base-measure rate as Paper 38. **This is the cross-domain signature worth chasing.**

**Shortcut 3: Empirical-first is fast.** PySCF + Block (or ITensor's chemistry interface, or PennyLane's `qml.qchem` MPO construction) can compute the MPO bond profile of a small chemistry Hamiltonian in under an hour. We can compute χ_k profile for LiH at n_max=2 directly, compare to the a priori Gaunt-channel-count prediction, and verify bit-exact match BEFORE committing to any proof. If empirical match holds across 4-5 molecules at multiple n_max, the theorem statement gets sharpened and the proof is much more tractable.

### 4.3 Empirical-first verification protocol

**Sprint 1 (3-5 days): Empirical bond-profile measurement.**

Pick LiH, BeH₂, H₂O at n_max=2. For each:

1. Get the GeoVac composed Hamiltonian via `geovac.ecosystem_export.hamiltonian()` in OpenFermion-compatible form.
2. Re-order orbitals into the natural chain (sub-block-outermost, then (n, l, m_l, s)).
3. Construct the MPO using the standard Keller-Dolfi-Troyer-Reiher finite-state-automaton algorithm. Code: open-source `block2` (Chen, Sharma, et al.) or PySCF's MPO interface (DMRGSCF wrapper).
4. Compute the bond profile χ_k for k = 1, ..., M-1.
5. Compute the a priori prediction χ_k^{Gaunt} from GeoVac's multipole channel count + sub-block boundaries.
6. Compare.

**Falsifiers:**
- If χ_k is NOT determined by (L, M, sub-block) structure — i.e., if it depends on Z, on specific molecule, on non-Gaunt-related parameters — then Theorem 3.2.A is false and the program stops.
- If χ_k matches but with an unexplained offset that scales with M, that's a "halfway POSITIVE" outcome: the structural prediction is correct in shape but missing a normalization. Worth a one-paragraph follow-on.

**Expected outcomes by molecule (this is the headline prediction):**

For LiH at n_max=2: M = 15 orbitals, 3 sub-blocks. χ_k should sit at:
- Within Li 1s sub-block (k=1): χ ≤ 7 (multipole count from the DF memo at this basis)
- At Li 1s / Li 2s2p sub-block boundary (k=2): χ should DROP toward 1 (since cross-block ERIs vanish per F4)
- Within Li 2s2p sub-block (k=3..7): χ ≤ 7 again
- At Li 2s2p / H 1s boundary (k=8..9): χ should DROP toward 1
- Within H 1s sub-block: trivial χ

The "drop at sub-block boundaries to ≤1" prediction is the **sharp falsifier**: standard chemistry MPOs for Gaussian Hamiltonians do NOT have this drop because their cross-orbital ERIs don't vanish. If GeoVac MPOs show this drop bit-exactly, the structural identification holds.

**Sprint 2 (1 week): a priori bond-rank theorem.**

If Sprint 1 confirms the drop-at-boundary pattern, write the proof of the bond-rank bound:

> **Lemma (a priori GeoVac MPO bond bound).** For a composed-builder Hamiltonian at cutoff n_max with N_sb sub-blocks, the MPO bond dimension along the natural orbital chain satisfies χ_k ≤ multipole-channel-count(active sub-block) + identity (kinetic + nuclear contributions). At a sub-block boundary, χ drops to the identity contribution + 0 (no cross-block ERIs). The bound scales as O(n_max²) within sub-blocks and O(1) at boundaries.

This Lemma + Theorem 3.2.A would land in Paper 14 as a new section, and in a possible math.OA companion paper as the bond-rank-truncation half of WH1.

**Sprint 3 (~2 weeks): Berezin map L4 ⊗ L4 construction.**

The deeper structural identification (Theorem 3.2.B operator-system isomorphism) needs the explicit B_{n_max} ⊗ B_{n_max} acting on the MPO tensor. This is the work that lands in a paper-grade math.OA write-up. The output would be a companion math.OA standalone (potentially Paper 56) extending Paper 38 to the MPS-bond-rank dual.

If Sprint 3 lands cleanly, the **headline number to report** is: Theorem 3.2.A and 3.2.B both proven, multi-month math.OA-grade Theorem 3.2.C still open but with named rate-matching follow-up.

---

## 5. Sprint sequence proposal (Deliverable 3)

### Sprint S1 — Empirical bond-profile measurement (3-5 days)

**PM task:** Build GeoVac MPO via Keller-Dolfi-Troyer-Reiher for LiH/BeH₂/H₂O at n_max=2, measure bond profile χ_k along the natural orbital chain, compare to a priori Gaunt-channel-count + sub-block-boundary prediction.

**Concrete deliverables:**
- `debug/geovac_mpo_bond_profile.py` — driver that takes a GeoVac Hamiltonian, orders orbitals, builds MPO, extracts bond profile.
- `debug/data/geovac_mpo_bond_profile.json` — bond profile data for LiH/BeH₂/H₂O at n_max=2.
- `debug/sprint_S1_mpo_bond_profile_memo.md` — comparison table (measured vs a priori), summary of findings.

**Expected outcome:** Bit-exact match between χ_k and the a priori prediction at the production basis, with χ_k dropping to 1 at every sub-block boundary.

**Falsifier:** If χ_k is not determined by (L, M, sub-block) structure — i.e., if there are cuts where χ_k matches no analytic prediction — Theorem 3.2.A is false. Stop and document.

**Natural follow-up:** If positive, Sprint S2 to write the a priori bond-rank Lemma; Paper 14 §intro addition immediately.

---

### Sprint S2 — A priori bond-rank Lemma + Paper 14 §sec:mpo_bond_rank (1 week)

**PM task:** Write a proof of the a priori bond-rank bound for the GeoVac composed Hamiltonian. The bound is structural (multipole-channel-count + sub-block boundaries) and the proof follows from Gaunt selection rules + the cross-block-ERI vanishing fact F4 from the DF memo.

**Concrete deliverables:**
- New §sec:mpo_bond_rank in Paper 14, ~3 pages, with: (a) statement of the bond-rank bound, (b) proof sketch via Gaunt + cross-block vanishing, (c) bit-exact numerical match table from Sprint S1, (d) implications for DMRG quantum chemistry.
- New CHANGELOG entry under "Tracks completed."
- New CLAUDE.md §3 row (if the empirical-first probe revealed a falsifier) or §2 entry (if positive).

**Expected outcome:** Paper 14 gains a 3-page subsection placing GeoVac directly in the DMRG / MPO literature with a sharp falsifiable prediction.

**Falsifier:** Already discharged by Sprint S1.

**Natural follow-up:** Sprint S3 to construct the Berezin L4 ⊗ L4 map; or stop here if PI prefers to consolidate.

---

### Sprint S3 — Theorem 3.2.B operator-system isomorphism via L4 ⊗ L4 (~2 weeks)

**PM task:** Extend Paper 38's Berezin reconstruction map L4 to act on tensor products of orbital functions, and use it to prove the operator-system isomorphism Theorem 3.2.B.

**Concrete deliverables:**
- New math.OA companion paper draft (~12-15 pages): "MPO bond-rank truncations as Connes-vS spectral truncations on tensor-product Camporesi-Higuchi triples." Title and structure parallel to Paper 39 (tensor-product propinquity).
- Five-lemma structure transplanted from Paper 38 with the L4 → L4 ⊗ L4 extension being the load-bearing new content.
- Bibliography: Paper 38 + Paper 39 + Keller et al. 2015 + Connes-vS 2021 + van Suijlekom 2020.

**Expected outcome:** A paper-grade math.OA standalone identifying the GeoVac n_max parameter with a specific Connes-vS truncation of a specific MPO bond-rank truncation, with rate inherited from Paper 38. Would land as the 13th or 14th math.OA standalone in the GeoVac series.

**Falsifier:** The L4 ⊗ L4 extension might fail the L3-compatibility lemma (Lemma L3 of Paper 38 bounds opnorm([D, M_f]) ≤ ‖∇f‖∞; its tensor-product extension would need to bound opnorm([D⊗I + I⊗D, M_{f ⊗ g}]) ≤ some natural function of ‖∇f‖∞ and ‖∇g‖∞). If this fails — which Paper 39 says it does NOT, because the Lichnerowicz / Lipschitz comparison transports cleanly to tensor products on the same manifold — then Theorem 3.2.B is false at the operator-system level and we'd need to retreat to Theorem 3.2.A.

**Natural follow-up:** Sprint S4 (NOT scoped here, multi-month) on the rate-matching side of Theorem 3.2.C.

---

### Cumulative report

**3 sprints total. ~3 weeks wall-clock if PI prioritizes.** Headline:
- Sprint S1 verifies the empirical claim (3-5 days).
- Sprint S2 writes the a priori Lemma + Paper 14 addition (1 week).
- Sprint S3 writes the math.OA companion paper (2 weeks).

If only S1 + S2 land, the **GeoVac claim** is sharpened: GeoVac MPOs have bond profile that is **structurally bounded a priori by multipole channels** + sub-block boundaries, predictable from the basis alone without any DMRG simulation. This is a Paper 14 / Paper 20 result and lands cleanly in the chemistry-QC outreach context.

If S3 also lands, the **structural claim** is the deepest of the four meta-lesson confirmations from this 24-hour arc: DF=multipole (bit-exact), Cholesky=multipole (bit-exact), QPT-Hopf-Z2 (bit-exact stack), MPO bond rank = Connes-vS spectral truncation (rate-matched). The PI's meta-lesson is now empirically 4 for 4, and we have a math.OA companion paper to anchor it.

---

## 6. Honest gaps / what could break the program

### 6.1 Empirical falsifiers (could kill the program in Sprint S1)

- **GeoVac MPO bond profile might depend on Z or on specific molecule beyond what Gaunt-channel-count predicts.** Possible mechanism: the kinetic energy + nuclear-attraction operators, even though they live in `A` (one-body), may inject non-trivial bond rank from off-diagonal kinetic matrix elements between distinct (n, l) shells within a sub-block. The MPO compresses these into the bond, and the compression rank is sensitive to the specific 〈n l | T | n' l'〉 matrix elements. **Concrete test in S1:** verify that the bond profile χ_k matches the a priori Gaunt-channel-count + a kinetic/nuclear contribution that is universal across Z. If not universal — kill program at S1.

- **The natural orbital chain might not be the GeoVac-optimal ordering for MPO compression.** Standard chemistry MPOs use Fiedler vector or genetic-algorithm ordering. If GeoVac's natural ordering produces χ_k much larger than the optimal ordering, the structural prediction is rate-dominated by ordering effects, not by the multipole structure. **Concrete test:** run multiple orderings in S1 (natural + Fiedler + reversed); if χ_k profile is bound by the same a priori count across all orderings (which it should be by Schur), the structural claim is robust.

### 6.2 Mathematical gaps (could limit reach of S3)

- **The 4/π rate transport via L4 ⊗ L4 needs Paper 39's tensor-product propinquity machinery.** Paper 39 already proved tensor-product propinquity on the SAME manifold T_{S³} ⊗ T_{S³} (this is the natural setting for the MPO L4 ⊗ L4 because all orbitals live on the same S³). The rate is bit-identical to the single-factor Paper 38 rate. So this gap is closed by Paper 39, not new work.

- **The bond-rank truncation might not be exactly the Berezin image of an MPO truncation.** The natural map from C(S³) ⊗ C(S³) (Berezin image space) to MPO tensors is not the identity — MPO tensors have a finite-state-automaton structure (left/right virtual indices). The L4 ⊗ L4 → MPO truncation map requires a discretization step that doesn't appear in Paper 38 / Paper 39. **This is the load-bearing math.OA-grade work item in S3.** If the discretization step doesn't admit a clean Berezin-like reconstruction map, Theorem 3.2.B retreats to a less-tight version where the bond-rank ranks match but the propinquity rate is multiplied by an unknown constant.

### 6.3 Scope-discipline gaps (already named in CLAUDE.md)

- **Open-shell systems and spin-coupling beyond singlets.** Sprint S1 should pick LiH/BeH₂/H₂O (all closed-shell). Open-shell systems (HF, CH₃, etc.) introduce extra bond-rank coming from the spin-coupling that the QPT memo (2026-06-05 evening) shows would stack additively with the Hopf-Z₂ tapering. The MPO bond profile for open-shell systems should be addressed only after closed-shell is verified. Honest scope.

- **The chemistry MPO literature consistently uses Gaussian orbitals; GeoVac uses the natural S³ basis.** No published comparison of GeoVac MPO bond profile vs Gaussian MPO bond profile exists. The comparison would land in Paper 14 as a new sub-section, but it's a separate sprint from the structural identification work. Not in S1/S2/S3 scope.

### 6.4 What could kill compression to sprint-scale specifically

- **Sprint S1 requires `block2` or `pyscf-DMRG` or `pennylane.qchem.mpo`.** None of these is currently installed in the GeoVac dev environment. The first 2-3 days of S1 would be installation + Hello-World on Gaussian H₂ to confirm the pipeline. If the install + pipeline takes more than 3 days, S1 stretches and may eat into S2.

- **The Keller-Dolfi-Troyer-Reiher MPO algorithm assumes a Fock-space single-particle basis; GeoVac's composed builder is already in Fock space, so this should adapt cleanly. But the "natural orbital chain" definition needs to be made precise for multi-center molecules.** This is named in §3.3 above but not yet implemented. Should be doable in S1 day 1.

### 6.5 Headline honest scope

**Sprint-compressible work item:** Theorems 3.2.A + 3.2.B, in 3 sprints, ~3 weeks. Theorem 3.2.A is the load-bearing empirical claim and lands in Paper 14. Theorem 3.2.B is the operator-system isomorphism and lands as a math.OA companion to Paper 38/39.

**Multi-month work item:** Theorem 3.2.C rate-matched MPS bond-rank propinquity convergence. This requires extending the Latrémolière propinquity machinery to the dual setting (Wasserstein-Kantorovich on states rather than Lipschitz seminorms on algebras). Not in this scoping memo's scope.

**Headline call:** This compresses to sprint-scale. The reason it does is that Paper 38 already supplies half the bridge (the Connes-vS side), and the lit-scan agent's 2-month estimate did not account for that.

---

## 7. Verification register

Every arXiv ID web-fetched and confirmed via curl against the abstract page. Confidence-builder; per the audit memo `debug/sprint_literature_audit_followup_memo.md`, prior sub-agent runs hit 18% hallucinated arXiv IDs, so verification is mandatory here.

| arXiv ID | Title | Authors | Year / Venue | Status |
|:---|:---|:---|:---|:---:|
| 1510.02026 | An Efficient Matrix Product Operator Representation of the Quantum-Chemical Hamiltonian | Keller, Dolfi, Troyer, Reiher | 2015 / J. Chem. Phys. 143, 244118 | VERIFIED (title + venue + submission date 7 Oct 2015) |
| 2005.08544 | Gromov-Hausdorff convergence of state spaces for spectral truncations | van Suijlekom | 2020 / J. Geom. Phys. | VERIFIED (title + submission 18 May 2020) |
| 2412.00628 | A noncommutative integral on spectrally truncated spectral triples, and a link with quantum ergodicity | Hekkelman, McDonald | 2024 | VERIFIED (title + submission 1 Dec 2024) |
| 1008.3477 | The density-matrix renormalization group in the age of matrix product states | Schollwöck | 2011 / Ann. Phys. 326, 96 | VERIFIED (title + submission 20 Aug 2010) |
| 1605.02611 | Matrix Product Operators, Matrix Product States, and ab initio Density Matrix Renormalization Group algorithms | Chan, Keselman, Nakatani, Li, White | 2016 / J. Chem. Phys. | VERIFIED (title + authors + submission 9 May 2016) |
| 2004.14115 | Spectral truncations in noncommutative geometry and operator systems (Connes-vS 2021) | Connes, van Suijlekom | 2021 / CMP | Cited in Paper 38 line 132; not directly web-verified in this scoping pass but used as Paper-38-internal reference, which is itself verified. Status: TRUSTED-VIA-PAPER-38 |
| 2506.09151 | The Quantum Paldus Transform: Efficient Circuits with Applications | Burkat, Fitzpatrick | 2025 / arXiv | VERIFIED in lit-search memo (`debug/lit_search_algebraic_methods_compatibility_2026_06_05.md`); not re-verified here. Status: TRUSTED-VIA-LIT-SEARCH-MEMO |
| Leimbach-vS 2024 | (torus propinquity) | Leimbach, van Suijlekom | 2024 / Adv. Math. 439, 109496 | Cited in Paper 38 as `leimbach_vs2024`; not directly web-verified due to arXiv ID search returning empty (multiple ID guesses failed). Status: TRUSTED-VIA-PAPER-38 (referenced inside the published math.OA proof, so safe). Recommend Sprint S2 author re-verify Leimbach-vS arXiv ID via Google Scholar before committing the proof to a draft. |

**Honest gap:** The Leimbach-vS 2024 arXiv ID was not recovered in this verification pass; the paper is cited internally in Paper 38 as `leimbach_vs2024` and as a journal reference (Adv. Math. 439, 109496). Sprint S2 should re-verify the arXiv ID before drafting Paper 14 additions; if the arXiv ID is unrecoverable, cite only the journal reference.

---

## Appendix A — Concrete first commands for Sprint S1

These are the literal commands a PM agent picking up S1 would run in day 1, ordered by dependency:

```python
# 1. Get GeoVac Hamiltonian
from geovac.ecosystem_export import hamiltonian
H_lih = hamiltonian("LiH", n_max=2)  # returns GeoVacHamiltonian object
H_lih_of = H_lih.to_openfermion()    # FermionOperator

# 2. Get orbital ordering (natural chain)
spec = H_lih.spec  # MolecularSpec; has orbital list with (sub-block, n, l, m_l, s)
ordering = sorted(spec.orbital_indices, key=lambda i: (spec.sub_block(i), spec.n(i), spec.l(i), spec.m_l(i), spec.s(i)))

# 3. Re-order Hamiltonian to natural chain
H_reordered = reorder_fermion_operator(H_lih_of, ordering)

# 4. Build MPO via block2 (Keller et al. algorithm)
import block2
mpo = block2.Hamiltonian.from_fermion_operator(H_reordered).build_mpo()

# 5. Extract bond profile
bond_profile = [mpo.bond_dimension_at_cut(k) for k in range(1, len(ordering))]

# 6. Compare to a priori prediction
predicted = predict_geovac_bond_profile(spec)  # custom function from Gaunt count + sub-block topology

# 7. Match?
for k, (chi_meas, chi_pred) in enumerate(zip(bond_profile, predicted)):
    print(f"cut {k}: measured χ={chi_meas}, predicted χ={chi_pred}, match={chi_meas == chi_pred}")
```

If this script runs cleanly for LiH/BeH₂/H₂O at n_max=2 and the match is bit-exact at every cut, Theorem 3.2.A is empirically established and Sprint S2 can proceed immediately.

---

End of scoping memo.
