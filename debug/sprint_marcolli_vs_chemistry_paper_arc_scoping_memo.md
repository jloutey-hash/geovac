# Sprint M-vS Chemistry — Paper-arc scoping memo

**Date:** 2026-06-07
**Triggering event:** H₂ Bratteli pilot returned with bit-exact (9.7×10⁻¹⁸ residual) match between Track CD's H₂ one-body Hamiltonian and the Marcolli-vS 2014 gauge-network construction on the H₂ bond quiver. See `debug/bratteli_h2_pilot_memo.md`.

**Decision target:** Should we open a math.OA paper arc claiming GeoVac molecular Hamiltonians are Marcolli-vS gauge networks?

**Verdict:** YES. The H₂ pilot provides empirical seed evidence at bit-exact precision. Scale-up sprints (LiH, NaH) are well-defined. Target venue is J. Geom. Phys. or J. Noncommut. Geom. Effort estimate: 4-6 months end-to-end including LiH/NaH scale-up + Bratteli combinatorial machinery + arXiv submission.

---

## §1 Paper title and one-line pitch

**Working title:** *"Molecular Hamiltonians as gauge networks: bit-exact correspondence between graph chemistry construction and Marcolli-van Suijlekom spectral triples"*

**One-line pitch:** GeoVac's chemistry composition (atomic spectral triple per atom + multipole bimodule per bond) is a Marcolli-vS 2014 gauge network on the molecular bond quiver, verified bit-exactly at the one-body Hamiltonian level on H₂, with structural extension to multi-atom and two-body content.

The paper would be the 15th math.OA standalone in the GeoVac series (siblings: Papers 32, 38-49, 53, 55, 56).

## §2 Main theorem (target)

**Theorem (informal target).** *Let M be a molecule with nuclear set N and bond set E, encoded in GeoVac at angular cutoff n_max ≥ 1. There exists a Marcolli-vS gauge network (Q, {(A_v, H_v, D_v)}, {(φ_e, L_e)}) on the bond quiver Q = (N, E) such that:*

*(i) Vertex prespectral triples carry atomic Camporesi-Higuchi Dirac data: A_v = M_d(ℂ) with d = d(n_max, Z_v) the orbital count, H_v = ℂ^d, D_v = atomic kinetic + same-side V_ne.*

*(ii) Edge bimodules carry multipole-V_ne cross-center data: (φ_e, L_e) = (identity, h1_cross[v, w]).*

*(iii) The assembled global Hamiltonian H = ⊕_v D_v + ⊕_e L_e coincides bit-exactly with Track CD's `build_balanced_hamiltonian` output one-body sector h1.*

*(iv) The spectral action S(H) = Tr f(H/Λ) is well-defined at every n_max and converges to the continuum spectral action in the Paper 38 propinquity sense as n_max → ∞.*

**Corollary (target).** *The two-body sector — Track CD's cross-block ERI tensor — admits a "Bratteli self-loop / plaquette" representation in the Perez-Sanchez 2024 sense, giving the full molecular Hamiltonian (h1 + eri + ecore) as a discrete quiver spectral-action evaluation.*

Status of corollary: speculative. Pilot did one-body only. Corollary needs a separate test sprint.

## §3 Section outline

| § | Section | Content |
|---|---|---|
| 1 | Introduction | Motivation — chemistry composition as graph operation. Cite Marcolli-vS 2014, WH1 PROVEN lineage, Track CD multipole. State main theorem informally. |
| 2 | Marcolli-vS gauge network preliminaries | Define quiver, vertex prespectral triple, edge bimodule, assembled Dirac, spectral action. Follow MvS 2014 §2 closely. |
| 3 | GeoVac chemistry composition recap | Define atomic Camporesi-Higuchi triple, Track CD balanced builder, multipole cross-V_ne, cross-block ERI. Cite Paper 32, 23. |
| 4 | Main theorem | H₂ pilot result formalized. Statement, proof for H₂ (bit-exact by construction), proof outline for general M (induction on bond count, requires LiH and NaH scale-up data). |
| 5 | Spectral action and continuum limit | Cite Paper 38 propinquity. Spectral action of finite-cutoff Hamiltonian converges to round-S³ spectral action in the truncated spectral triple limit. |
| 6 | Numerical verification panel | Bit-exact tables: H₂ (done), LiH (scale-up sprint), NaH (scale-up sprint). Eigenvalue match, spectral action match across Λ panel. |
| 7 | Open questions | Two-body ERI Bratteli reading (corollary). Unitarity gap with strict Perez-Sanchez 2024. Continuum limit at the chemistry-observable level. |
| App A | Bratteli combinatorics | Path-algebra picture for multi-bond molecules. Required for LiH (3-vertex with Li_core + LiH_bond_center + LiH_bond_partner = 3 sub-blocks). |
| App B | Reproducibility | Driver paths, JSON data, package versions. |

Estimated length: 18-22 pages, three-pass clean, 30-40 bibitems.

## §4 What's already done vs what needs new work

### Done (today, 2026-06-07)
- H₂ pilot at n_max=2: bit-exact match Marcolli-vS construction vs Track CD h1.
- Eigenvalue match: max residual 0.0 across all 10 eigenvalues.
- Spectral action match: bit-exact at Λ ∈ {1, 2, 4} for the full-Bratteli (vertex-Dirac-restored) variant.
- Driver: `debug/bratteli_h2_pilot_driver.py` (~440 lines, runs in ~5s)
- Construction notes: `debug/bratteli_h2_pilot_construction_notes.md` (~3000 words)
- Verdict memo: `debug/bratteli_h2_pilot_memo.md` (~3500 words)

### Needed (sprint-scale, in order)

**Sprint M-vS-1 (1-2 weeks): LiH scale-up.**
Goal: verify bit-exact match on LiH at n_max=2, which has 3 sub-blocks (Li_core, LiH_bond_center, LiH_bond_partner) and tests the multi-vertex Bratteli combinatorics. The 3-vertex case is the first non-trivial multi-bond instance — H₂ is a 2-vertex special case.

Decision gate: max residual ≤ 1e-10 between Marcolli-vS assembled H and Track CD h1.

Files: extend `debug/bratteli_h2_pilot_driver.py` to handle 3-vertex; new `debug/bratteli_lih_pilot_memo.md`.

If FAIL: structural surprise. Stop arc, write a "structural negative" note explaining what doesn't match.

**Sprint M-vS-2 (1-2 weeks): NaH scale-up.**
Goal: verify bit-exact match on frozen-core NaH (default frozen path; the explicit-core path is B.1's territory). NaH is structurally interesting because the frozen [Ne] core lives in `nuclear_repulsion_constant` as a CLASSICAL ENERGY OFFSET, not as quantum-encoded orbitals. The Bratteli reading should account for this: maybe the [Ne] core is encoded as a non-trivial Higgs term on the Na vertex (the "vertex self-loop Higgs" picture in Perez-Sanchez 2024 §4)?

Decision gate: max residual ≤ 1e-10 for the quantum-encoded part; honest treatment of the classical frozen-core energy term.

If PASS: the frozen-core treatment fits the Bratteli framework cleanly. If FAIL: the classical energy offset needs a different structural reading — possibly an additive constant outside the Bratteli framework.

**Sprint M-vS-3 (3-4 weeks): two-body ERI Bratteli reading.**
Goal: test the corollary. Does Track CD's cross-block ERI tensor admit a Bratteli self-loop / plaquette reading?

This is the load-bearing test for the FULL paper. If it PASSES, the molecular Hamiltonian (h1 + eri + ecore) is a discrete quiver spectral-action evaluation — a much stronger claim. If it FAILS, the paper is one-body-only and we cite the two-body gap as future work.

Approach: the cross-block ERI tensor at chemistry order (pq|rs) couples four orbitals across (up to) two atoms. In Perez-Sanchez §4 the spectral-action expansion includes plaquette terms ∝ Tr(L_e L_e' L_e'' L_e''') = closed-walk traces of length 4. The empirical test: does the cross-block ERI tensor match such a closed-walk-4 reading?

Decision gate: structural identification with empirical numerical match within 1%. (Bit-exact unlikely given the multipole expansion truncation, but structural shape match would be enough.)

**Sprint M-vS-4 (4-6 weeks): theorem statement and proof.**
Synthesize §1-3 sprint results into the formal theorem. Prove the H₂ case directly. Prove the multi-atom case by induction on bond count (the LiH and NaH scale-up sprints supply the induction step empirically; formalize as theorem with sketch + numerical witness panel).

**Sprint M-vS-5 (4-6 weeks): paper drafting.**
LaTeX draft following the §3 outline above. Use the chemistry-pivot scoping memo + the verdict memo + the construction notes as raw material. Target J. Geom. Phys. or J. Noncommut. Geom.

**Total estimated effort:** 14-22 weeks (~4-5 months). Comparable to Papers 45 (K⁺-weak-form Lorentzian propinquity, 1-month effort) or 49 (OSLPLS strong-form bridge, 1-month effort). Within standard GeoVac math.OA standalone scope.

## §5 Risks and dependencies

### Risk 1: LiH 3-vertex doesn't match
The H₂ result is 2-vertex, which is a special case. LiH is the smallest non-trivial test. If the multi-vertex Bratteli combinatorics don't produce a bit-exact match, the paper arc retracts to "Marcolli-vS gauge networks work for H₂; multi-atom is a structural gap." Honest negative but kills the paper.

Mitigation: run M-vS-1 (LiH scale-up) FIRST, before committing to the theorem statement. Lowest cost decisive test.

### Risk 2: Two-body ERI doesn't have Bratteli structure
The corollary fails. Paper becomes one-body only. Less impactful but still a real paper — Marcolli-vS at the one-body level is original on GeoVac.

Mitigation: Sprint M-vS-3 is a separate sub-sprint. Frame the paper as one-body main result with two-body as future work if it doesn't close.

### Risk 3: Unitarity gap precludes "Bratteli network" claim entirely
The H₂ pilot showed L_e is not unitary (‖L_e† L_e − I‖ = 1.0). Marcolli-vS 2014's definition allows non-unitary bimodules in the inner-fluctuation sense (Connes-Chamseddine-style), but Perez-Sanchez 2024a requires unitarity. If a referee or collaborator pushes back on the non-unitarity claim, we may need to either (a) find a polar-decomposition rewrite that makes L_e unitary + a Hermitian "Higgs" piece, or (b) cite a Marcolli-vS variant that explicitly allows Hermitian bimodules.

Mitigation: the construction notes already document the issue. A short investigation (~1 day) into whether L_e = U_e ⊗ Φ_e factorization exists would resolve this. Cite Bertozzini-Conti-Lewkeeratiyutkul 2014 morphism axioms if they admit Hermitian intertwiners.

### Dependency 1: WH1 PROVEN lineage continuity
The paper relies on Paper 32, Paper 38, and the WH1 PROVEN structural framing. If those papers are still in flux, the new paper needs to thread through the existing math.OA lineage carefully. As of 2026-06-07, all WH1 papers are in stable state.

### Dependency 2: H₂ Bratteli pilot's driver and data
The empirical seed for the paper. The driver is on disk; data file `debug/data/bratteli_h2_pilot.json`. Reproducible in 5 seconds on stock hardware. Solid foundation.

## §6 Strategic positioning

This paper would:
- Provide the structural NCG home for GeoVac's chemistry composition that the conversation today centered on
- Validate the user's bond-sphere / composed-NCG intuition with a concrete published-framework match
- Complement the chemistry-engineering negatives (B.1 explicit-core HF doesn't bind NaH; LiH kwarg sweep clean STOP) by giving the framework a positive structural story to tell
- Sit naturally alongside Paper 20 (resource estimation) and Paper 57 (proposed chemistry-pivot paper) as the "structural why" complement to the "engineering what" of those papers

**Honest comparison to existing math.OA standalones:**
- Paper 38 (WH1 PROVEN): GH-convergence of SU(2) truncation. Proves a convergence theorem at the spectral-triple level.
- Papers 45/46/47/48/49 (Lorentzian arc): Krein-space propinquity / OSLPLS bridge / etc. Proves convergence theorems at the operator-system / pre-length-space level.
- This paper (proposed): GeoVac chemistry as Marcolli-vS gauge network. Proves a *structural identification* theorem — different shape from the convergence papers, but firmly in the same math.OA / J. Geom. Phys. lineage.

The result is also citation-friendly for chemistry/QC consumers: when someone asks "what's the framework's structural foundation?", we point at this paper. It's the chemistry-side companion to Paper 32 (general spectral triple) and Paper 38 (propinquity convergence).

## §7 Recommended starting move

Skip outreach (per PI direction this turn). Run **Sprint M-vS-1 (LiH scale-up)** as the next concrete piece. Single-week scope, decisive verdict on whether the multi-vertex Bratteli reading holds at all. If it PASSES, we have a viable paper arc. If it FAILS, we know to retract to "H₂-only structural observation" framing and the paper arc shrinks accordingly.

The PI's previous decision-gate from "all 3" remains in force: if engineering is exhausted (it is, per B.1 today) AND the structural framework matches (it does at H₂), the paper arc is the clear next investment.

## §8 Files

### To be created
- `debug/sprint_marcolli_vs_chemistry_paper_arc_scoping_memo.md` (this memo)
- `debug/bratteli_lih_pilot_driver.py` (Sprint M-vS-1, future)
- `debug/bratteli_lih_pilot_memo.md` (Sprint M-vS-1, future)

### Existing references
- `debug/bratteli_h2_pilot_memo.md` — pilot verdict
- `debug/bratteli_h2_pilot_construction_notes.md` — pilot construction details
- `debug/bratteli_h2_pilot_driver.py` — reproducible driver
- `debug/data/bratteli_h2_pilot.json` — pilot numerical results
- `papers/group3_foundations/paper_32_spectral_triple.tex` — framework's spectral triple base
- `papers/group1_operator_algebras/paper_38_su2_propinquity_convergence.tex` — propinquity / WH1 PROVEN
- `MEMORY.md` entry `wh1_marcolli_vs_lineage.md` — pre-existing identification
- arXiv:1301.3480 (Marcolli-van Suijlekom 2014) — primary external reference

### Memory file to update (per Bratteli verdict memo §4.3)
`memory/wh1_marcolli_vs_lineage.md` — note the bit-exact H₂ pilot confirmation. Don't touch in this scoping turn; queue for a follow-on update.
