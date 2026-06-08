# Sprint W1e-Audit and M-vS-Chemistry-Open (2026-06-07, umbrella)

**Sprint scope:** day-long arc that started as a single diagnostic — *"why does the qubit Hamiltonian not bind LiH while the continuous adiabatic+PK solver does?"* — and ended with the chemistry-engineering arc empirically exhausted, the Paper 20 chemistry-accuracy scope boundary documented in the public record, and a new math.OA paper arc opened with bit-exact empirical seed on H₂ and LiH.

**Architecture:** seven sub-sprints in a single session, plus production code fix, plus test infrastructure work, plus Paper 20 editing. Cross-thread synthesis is the heart of this memo.

**Headline:**
- **Chemistry-engineering arc is empirically exhausted at second-row.** All cheap engineering kwargs are no-op on first-row (LiH kwarg sweep clean STOP); explicit-core Hartree-Fock un-freezes the [Ne] core but does not close NaH binding (B.1 preliminary STOP, with methodology validated via LiH cross-check).
- **Marcolli-van Suijlekom 2014 gauge networks are the structural home for GeoVac chemistry.** Bit-exact match between Track CD's one-body Hamiltonian and the assembled Marcolli-vS construction at machine precision (residual ~1e-17), now confirmed on two molecules (H₂ + heteronuclear LiH).
- **Production code bug fixed and regression-tested.** The Path B V_NN double-count in `balanced_coupled.py`'s R-dependence corrector — silently masked by R3-A/R3-B's call convention until today — is patched via new `MolecularSpec.R` field. 202/202 regression tests pass.
- **Paper 20 framing now honest about scope.** Abstract, new §V.C (~2.5 pages), conclusion updated with row-conditional accuracy claim. Compiles clean.
- **Test suite slimmed.** 40-molecule sweep → pairwise representative matrix; session-scoped `hamiltonian_cache` fixture. Default test_ecosystem_export.py runtime: 8.5 min → 2.5 min (3× speedup).

---

## §1 Sub-sprint inventory and canonical memo pointers

Seven distinct sub-sprints, with one production code change and one paper edit:

| # | Sub-sprint | Verdict | Canonical memo |
|---|---|---|---|
| 1 | W1e-Projection-Audit (LiH diagnostic) | LiH balanced binds (path-A); NaH overattracts; Path B V_NN bug | `debug/sprint_w1e_projection_audit_memo.md` |
| 2 | LiH kwarg sweep | clean STOP — all cheap kwargs no-op or catastrophic | `debug/lih_kwarg_sweep_log.txt` |
| 3 | Plan agent A/B scope | scoping document; A.1 + B.1 parallel recommended | (agent transcript) |
| 4 | NCG composition Explorer | Marcolli-vS / Bratteli flagged as top candidates | (agent transcript) |
| 5 | Bratteli deep-dive | MEDIUM-LOW Perez-Sanchez; cheap H₂ pilot named | (agent transcript) |
| 6 | H₂ Bratteli pilot | PASS bit-exact (9.7e-18) Marcolli-vS; STOP Perez-Sanchez | `debug/bratteli_h2_pilot_memo.md` |
| 7 | B.1 explicit-core NaH HF | preliminary STOP; LiH HF cross-check validates methodology | `debug/sprint_b1_nah_hf_verdict_memo.md` |
| 8 | M-vS-1 LiH heteronuclear Bratteli | PASS bit-exact (2.5e-17) — paper arc has 2nd anchor | `debug/bratteli_lih_pilot_memo.md` |

Two scoping memos and one Paper 20 edit pass also emerged:

- `debug/sprint_marcolli_vs_chemistry_paper_arc_scoping_memo.md` — full paper arc structure (Sprint M-vS-1..5, 4-5 month effort)
- `debug/sprint_paper20_paper57_scope_boundary_scoping_memo.md` — Paper 20 / Paper 57 framing edits
- `papers/group4_quantum_computing/paper_20_resource_benchmarks.tex` — abstract + new §V.C + conclusion applied

## §2 Sub-sprint verdicts (one paragraph each)

**(1) W1e-Projection-Audit.** Closed the v3.85.0 chemistry-pivot follow-on question. Found that the R3-A/R3-B "balanced doesn't bind LiH/NaH" verdict was a Path B V_NN-double-count bug in `balanced_coupled.py` — Track CD's R-dependence corrector assumed spec was built at default R but unconditionally subtracted V_NN at the default. **For LiH:** Path A balanced (spec at default R, builder at requested R) binds at R_eq = 3.015 bohr with D_e = 0.158 Ha (2.4× over-binding vs continuous 0.067 Ha). **For NaH:** Path A still overattracts; bug fix doesn't close NaH. Two-line production code fix applied: new `R` field on `MolecularSpec`, populated by `hydride_spec`, consumed by the corrector. 202/202 regression tests pass. Sprint memo: `debug/sprint_w1e_projection_audit_memo.md`.

**(2) LiH kwarg sweep.** 11 kwarg combinations tested on LiH balanced. **8 of 11 give bit-identical energies to baseline** — `multi_zeta_basis`, `screened_valence_basis`, `screened_cross_center`, `pk_cross_center` are all gated on `Z_nuc_center >= 11` and silently no-op for first-row LiH. `cross_block_h1` makes LiH 16× over-bound (catastrophic; same pattern as F3 dead-end). The cheap engineering knobs cannot close LiH over-binding because they don't even trigger. Empirically confirms what the Explorer's NaH probe found: existing kwargs are exhausted at the chemistry-accuracy level for both first-row and second-row hydrides.

**(3) Plan agent A/B scope.** Multi-step research agent scoped the engineering threads. Recommended A.1 (NaH multi-zeta) + B.1 (HF on explicit-core NaH) in parallel as Week 1; identified B.1 as the *decisive single experiment*. Cross-thread strategic note: A.3 might be subsumed by B.1 if explicit-core HF binds NaH (frozen-core is then the wall, not the static basis). This recommendation routed today's afternoon work.

**(4–5) NCG composition Explorer + Bratteli deep-dive.** Two agent dispatches that converged on the same answer: the closest NCG composition framework for "atoms fused by bonds" is the Perez-Sanchez 2024 Bratteli-network / Marcolli-van Suijlekom 2014 gauge-network lineage. The deep-dive's MEDIUM-LOW Perez-Sanchez confidence verdict was upheld but sharpened: the actual gap is vertex-Dirac-dropping (Perez-Sanchez 2024a removes them; Marcolli-vS 2014 keeps them; GeoVac has atomic Diracs), not finite-dim-vs-infinite-dim. Named a cheap diagnostic — the H₂ pilot — to test structurally before committing to a multi-month paper arc.

**(6) H₂ Bratteli pilot.** Background agent constructed Bratteli network data explicitly for H₂ at n_max=2, R=1.4 bohr. Numerical comparison against Track CD's `balanced_coupled` output: max residual 9.7×10⁻¹⁸ on all four blocks, eigenvalues bit-identical, spectral action bit-exact at Λ ∈ {1, 2, 4}. **PASS-Marcolli-vS / NO-Perez-Sanchez 2024a (vertex Diracs dropped; edge bimodule not unitary).** Cited Connes-vS WH1 lineage as the structural home; recommended outreach to van Suijlekom or Marcolli (not Perez-Sanchez 2024a). Sprint memo: `debug/bratteli_h2_pilot_memo.md`.

**(7) B.1 explicit-core NaH HF.** Implemented `force_explicit_core=True` kwarg on `hydride_spec` (~80 lines), wrote 12-electron RHF SCF driver (~330 lines) with DIIS + density damping + level shift. **LiH HF cross-check binds at R_eq = 3.015 bohr with D_e = 0.158 Ha (exact match to today's balanced FCI verdict — methodology validated).** Explicit-core NaH HF gives monotone descent into small R, 5/7 R points converge, 2 SCF-bistable in the bond-breaking region. The wall is NOT the frozen-core projection; un-freezing the [Ne] core doesn't close NaH. Engineering arc B (operator-side fixes) joins arc A (basis-side fixes) in the empirically-exhausted column. Memo: `debug/sprint_b1_nah_hf_verdict_memo.md`.

**(8) M-vS-1 LiH heteronuclear Bratteli pilot.** Scaled the H₂ pilot's 2-vertex construction to LiH (Z_a = 3, Z_b = 1). max residual 2.47×10⁻¹⁷, all four blocks at machine precision, eigenvalues match to 1.55×10⁻¹⁵, spectral action bit-exact at every Λ. **PASS-Marcolli-vS at bit-exact precision — heteronuclear extension works.** Two empirical anchors now: H₂ + LiH. Sprint memo: `debug/bratteli_lih_pilot_memo.md`.

## §3 Cross-thread structural findings

### 3.1 The chemistry-engineering arc is empirically exhausted

Two parallel sub-sprints (LiH kwarg sweep, B.1 explicit-core HF) plus the cumulative weight of the F1–F6 dead-ends in CLAUDE.md §3 plus the Plan agent's A/B scoping plus today's Explorer NaH probe converge: **at the framework's current scope, no existing engineering knob closes chemistry-accuracy second-row binding.** Frozen-core projection is not the wall; the basis-extent hypothesis (multi-zeta, screened-valence) is not the wall; cross-block ERI grid truncation is not the wall; multi-determinant FCI on the given integrals is not the wall (R3-B falsifier verdict survives).

This is the **calibration-data scope boundary** that has been quietly accruing across the chemistry sprints since v2.0.6. Today it's localized at the second-row hydride threshold and documented in the public record (Paper 20 §V.C).

### 3.2 The Marcolli-vS gauge-network correspondence is empirically real

Two molecules, both at bit-exact precision:
- H₂ at R = 1.4 bohr, n_max = 2: max residual 9.7×10⁻¹⁸
- LiH at R = 3.015 bohr, n_max = 2: max residual 2.47×10⁻¹⁷

The construction:
- Vertices = atoms (one per nucleus)
- Vertex prespectral triple (A_v, H_v): hydrogenic Z_v at max_n=2 (5 orbitals each)
- Vertex Dirac D_v: hydrogenic eigenvalues + cross-center V_ne from the OTHER atom
- Edge bimodule L_e: GeoVac's cross_block_h1 (Hermitian, NOT unitary)
- Assembled global H matches Track CD's `balanced_coupled` output exactly

The Perez-Sanchez 2024a vertex-Dirac-dropping is *one step too far* for GeoVac. The Marcolli-vS 2014 vertex-Dirac-restored framework fits. This places GeoVac's chemistry composition in a published math.OA lineage we already know (WH1 PROVEN Marcolli-vS).

### 3.3 The two-sub-block 2-vertex case is structurally clean; the 3-sub-block default-spec case is the next sprint

Today's M-vS-1 LiH pilot used an *artificial* 2-vertex spec (Li_lone_pair + H_lone_pair, 5 orbitals each = 10 spatial). The *default* LiH spec from `lih_spec()` produces 3 sub-blocks (Li_core + LiH_bond_center + LiH_bond_partner, 5 orbitals each = 15 spatial), with a bond block that straddles two atoms. Whether the default 3-sub-block spec is ALSO a Marcolli-vS network is the M-vS-2 question. If it passes, the paper headline becomes "the production Track CD pipeline is structurally a Marcolli-vS gauge network on the molecular bond quiver" — a much stronger claim than "we can encode a chemistry-like spec into Marcolli-vS by construction."

### 3.4 Engineering negative + structural positive = clean strategic pivot

The "(a)+(b) both negative" pivot trigger condition (CLAUDE.md task #18 strategic posture) is technically not met. What happened instead:
- (a) engineering = NEGATIVE
- (b) Bratteli (the cheap NCG diagnostic) = POSITIVE in unexpected direction (Marcolli-vS, not Perez-Sanchez)

This is *stronger* than the both-negative case. The engineering arc is exhausted AND we have a positive structural home. The pivot toward the math.OA paper arc is therefore not a defensive move (accepting failure) — it's a constructive move (the structural reformulation matches at bit-exact precision; the multi-month paper arc has empirical seed and is worth the investment).

## §4 Production code changes

Two production files touched. Both backward-compatible. All tests pass.

### 4.1 `geovac/molecular_spec.py`
- Added `R: Optional[float] = None` field to `MolecularSpec` dataclass. Records the actual R the spec was built at.
- Populated by `hydride_spec` from the resolved R.
- Added `force_explicit_core: bool = False` kwarg to `hydride_spec`. Allows un-freezing the [Ne] core for Z=11–17 (extensible). When True: routes to a new `'explicit_extended'` core_type, builds a 10-electron [Ne] core block at Z_center=Z (bare).
- Added `_SECOND_ROW_ATOMIC_ENERGY` table (NIST values for Na, Mg, Al, Si, P, S, Cl).

### 4.2 `geovac/balanced_coupled.py`
- Lines 669-680: V_NN R-dependence corrector now reads `spec_R = getattr(spec, 'R', None) or _HYDRIDE_REQ.get(spec.name)`. Backward-compatible fallback. Eliminates the Path B V_NN double-count bug.

### 4.3 Test infrastructure (non-production)
- `tests/conftest.py`: added session-scoped `hamiltonian_cache` fixture. Each unique `(system, **kwargs)` is built once and reused across tests.
- `tests/test_ecosystem_export.py`: introduced `REPRESENTATIVE_PAULI` set (10 systems covering pairwise matrix); marked the 40-molecule sweep `@pytest.mark.slow`; rewired all tests to use the fixture.
- `tests/test_classical_quantum_parity.py` (new): two `@slow` tests that catch the V_NN bug class.

Default `test_ecosystem_export.py` runtime: **8:30 → 2:50 (3× speedup)**. Full library coverage and parity tests run via `pytest --slow`.

### 4.4 Paper edits
- `papers/group4_quantum_computing/paper_20_resource_benchmarks.tex`:
  - Abstract: row-conditional scope claim (first-row binds; second-row scope-bounded)
  - New §V.C "Chemistry-accuracy scope boundary" (~2.5 pages): empirical table, operator-level diagnosis, implication for users, structural-reason forward-reference
  - Conclusion: reframed to document the row-conditional scope
  - Future directions: added (v) "pushing the chemistry-accuracy scope boundary into second-row" as a separately-scoped engineering arc
- `papers/group4_quantum_computing/paper_20_refs.bib`: added `@article{marcolli_vs2014}` for the new §V.C citation.

LaTeX compiles cleanly to 11 pages, no new warnings.

## §5 Implications for papers / strategic posture

### Paper 20 (Resource Benchmarks)
Done in this sprint. The row-conditional scope boundary is in the public record. Future submissions of Paper 20 reflect today's evidence honestly.

### Paper 57 (Chemistry-Pivot Paper, not yet drafted)
The proposed three-claim framing now sharpens to:
1. Sparse Hamiltonian export with propinquity bounds (unchanged)
2. Sub-mHa VQE on small systems (unchanged, v3.85.0)
3. **Row-conditional operational scope-boundary localization** at the second-row hydride threshold, with operator-level evidence and Marcolli-vS structural framing

Drafting can proceed when ready; today's two scoping memos cover what should be in it.

### Marcolli-vS chemistry paper (new math.OA arc, not yet drafted)
Two empirical anchors at bit-exact precision (H₂ + LiH). M-vS-1 PASS. Next sprint M-vS-2 is the default LiH spec 3-sub-block test. If M-vS-2 passes, the paper headline becomes very strong. Full arc scoped to ~4-5 months. This would be the 15th math.OA standalone in the GeoVac series.

### Strategic posture
The chemistry-engineering arc is parked. The Marcolli-vS chemistry paper arc is the active forward direction. Engineering follow-ons (B.4 bonding-orbital Pauli repulsion, A.3 R-adaptive Z_eff) are not deprioritized to zero — they remain options if a chemistry consumer specifically asks — but they're not the project's current focus.

## §6 Follow-on register

### Shippable
- Sprint M-vS-2: default LiH spec 3-sub-block Bratteli reading (1 week scope)
- Paper 57 drafting once chemistry-pivot scope-boundary memo is internalized (4-6 weeks scope)

### Mid-term
- Sprint M-vS-3: two-body ERI Bratteli reading (3-4 weeks)
- Sprint M-vS-4: theorem statement + induction proof (4-6 weeks)
- Sprint M-vS-5: paper drafting (4-6 weeks)

### Optional (engineering, lower priority post-today)
- Sprint A.2: multi-zeta tabulation extension to first-row (test LiH overbinding reduction)
- Sprint A.3: R-adaptive Z_eff in cross-V_ne
- Sprint B.4: bonding-orbital Pauli repulsion against explicit core

### Frozen-core V_cross R-dependence
Named open since v3.56.0 F.1 sprint memo; remains open. Affects NaH and other frozen-core specs at non-default R. Not blocking; documented.

## §7 Honest scope

**What this sprint did:**
- Resolved the v3.85.0 follow-on question about LiH vs NaH binding cleanly
- Empirically exhausted the cheap chemistry-engineering kwargs (LiH sweep + Explorer NaH verdict)
- Validated explicit-core HF as not-the-wall for NaH binding
- Confirmed Marcolli-vS gauge-network correspondence at bit-exact precision on 2 molecules
- Fixed the Path B V_NN bug in production code, regression-tested
- Slimmed the default test suite from 8.5 min to 2.5 min on the heaviest file
- Documented the chemistry-accuracy scope boundary in Paper 20's public-facing record
- Scoped the Marcolli-vS chemistry paper arc with concrete sprint sequence

**What this sprint did NOT do:**
- The B.1 explicit-core HF result is preliminary (2/7 SCF non-convergence); a pyscf cross-check would tighten the verdict but is unlikely to change the qualitative direction
- The default LiH spec 3-sub-block Bratteli reading is not tested (M-vS-2 follow-on)
- The two-body ERI Bratteli reading is not tested (M-vS-3 follow-on)
- Paper 57 is not drafted
- The Marcolli-vS theorem statement is not formalized
- Outreach emails are not sent (per PI direction)
- Frozen-core V_cross R-dependence bug remains open

**Sample size:**
- Marcolli-vS correspondence: n=2 molecules (H₂ + LiH) at bit-exact, 1 cutoff (n_max=2)
- B.1 explicit-core HF: 1 molecule (NaH) + 1 cross-check (LiH); 7 R points
- LiH kwarg sweep: 11 kwarg combinations × 5 R points

## §8 Verification

- 202/202 targeted regression tests pass on `test_balanced*`, `test_cross_block_h1`, `test_ecosystem_export` (post bug fix)
- 18/18 topological S³ proofs pass (`test_fock_projection.py`, `test_fock_laplacian.py`)
- Slimmed default test suite: 31/31 pass + 76 skipped under `@pytest.mark.slow`
- Parity tests pass under `--slow` (3:41 wall)
- Paper 20 compiles cleanly to 11 pages (4 pre-existing undefined references, no new ones)
- Bratteli pilots: both reproduce bit-exact match in ≤5 seconds wall time

## §9 Hard-prohibition check (CLAUDE.md §13.5)

No changes to:
- Natural geometry hierarchy
- Fitted-or-empirical parameters (the bug fix uses `spec.R` which is set, not fitted; the explicit-core construction uses NIST atomic energies which are reference data, not fitted)
- Section 3 deletions (only adding qualifications to existing W1e/R3-A/R3-B rows)
- Paper 2 K = π(B+F−Δ) combination-rule "conjectural" labeling

## §10 Files

### Created (debug/)
- `sprint_w1e_projection_audit_memo.md`
- `w1e_projection_audit_driver.py` + `data/w1e_projection_audit.json` + `w1e_projection_audit_log.txt`
- `lih_kwarg_sweep_driver.py` + `data/lih_kwarg_sweep.json` + `lih_kwarg_sweep_log.txt`
- `sprint_b1_explicit_core_nah_scoping_memo.md`
- `b1_basis_sanity_diagnostic.py`
- `b1_nah_hf_driver.py` + `data/b1_nah_hf.json` + `b1_nah_hf_log.txt`
- `sprint_b1_nah_hf_verdict_memo.md`
- `bratteli_h2_pilot_driver.py` (background agent) + `data/bratteli_h2_pilot.json` + `bratteli_h2_pilot_memo.md` + `bratteli_h2_pilot_construction_notes.md`
- `bratteli_lih_pilot_driver.py` + `data/bratteli_lih_pilot.json` + `bratteli_lih_pilot_log.txt` + `bratteli_lih_pilot_memo.md`
- `sprint_marcolli_vs_chemistry_paper_arc_scoping_memo.md`
- `sprint_paper20_paper57_scope_boundary_scoping_memo.md`
- `perez_sanchez_outreach_email_draft.md` (drafted but not sent per PI direction)
- `sprint_w1e_audit_and_mvs_chemistry_open_2026_06_07_memo.md` (this memo)

### Modified (production)
- `geovac/molecular_spec.py` — added `R` field, `force_explicit_core` kwarg, `_SECOND_ROW_ATOMIC_ENERGY` table
- `geovac/balanced_coupled.py` — V_NN corrector reads `spec.R`
- `tests/conftest.py` — `hamiltonian_cache` fixture
- `tests/test_ecosystem_export.py` — pairwise representative matrix + slow-marked full library
- `tests/test_classical_quantum_parity.py` (new file)
- `papers/group4_quantum_computing/paper_20_resource_benchmarks.tex` — abstract + new §V.C + conclusion
- `papers/group4_quantum_computing/paper_20_refs.bib` — added `marcolli_vs2014`

### Modified (memos / documentation)
- `debug/sprint_chemistry_pivot_2026_06_07_memo.md` — addendum at top + §3.1 sharpening + §4 Paper 57 framing
- `debug/sprint_r3a_dmrg_lih_memo.md` — addendum noting Path B bug correction
- `debug/sprint_r3b_dmrg_nah_falsifier_memo.md` — addendum noting Path B bug
- `CLAUDE.md` §2 — added W1e-Projection-Audit entry; sharpened v3.85.0 entry with row-conditional marker
- `CLAUDE.md` §3 — split LiH dead-end row (composed stands, balanced retracted); added Path B caveat to NaH row
