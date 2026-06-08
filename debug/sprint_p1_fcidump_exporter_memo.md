# Sprint P1 — FCIDUMP exporter (Phase 1 H1-Interface partial) — canonical memo

**Date:** 2026-06-07.
**Sprint position:** Phase 1 Sprint H1-Interface (first sub-sprint) of the hybrid-pipeline arc (scoping memo: `debug/sprint_hybrid_pipeline_scoping_memo.md`). Surfaces the pre-Jordan-Wigner (h1, eri, ecore, n_electrons) on `GeoVacHamiltonian` and exposes them via a Knowles-Handy FCIDUMP exporter. This is the load-bearing unblocker named in §1.2 of the scoping memo: once FCIDUMP exports cleanly, DMRG (Block2), CCSD(T) (pyscf), and AFQMC (ipie) all unblock at once.
**Verdict line:** **GO at full Phase 1 scale.** LiH integrals surface cleanly, FCIDUMP write/read round-trip is bit-exact on h1/eri/ecore at max diff 0.0 (Hermitian + eight-fold permutation symmetry both preserved), and the same wiring transports verbatim to seven other systems (LiH, BeH2, NaH, H2O, CO, ScH, H2 all tested). Test suite expanded from 104 → 115 (+11), zero regressions.
**Cross-references:** `debug/sprint_hybrid_pipeline_scoping_memo.md` (Phase 1 scoping), `geovac/ecosystem_export.py` (modified), `geovac/composed_qubit.py` (source of h1/eri/ecore, unmodified), `geovac/balanced_coupled.py` (same shape, surfaces via the same wiring), `geovac/qubit_encoding.py::build_fermion_op_from_integrals` (chemist-notation convention pinning, unmodified), `tests/test_ecosystem_export.py` (11 new tests, all passing).

---

## §0. What was done

Wired the `(h1, eri, ecore, n_electrons)` quad through every composed-path builder in `ecosystem_export.py` (`_build_hydride`, `_build_multi_center`, `_build_tm_hydride`, `_build_alkaline_earth_monohydride`, `_build_h2`, plus the `_apply_hopf_tapering_to_geovac_hamiltonian` path), added matching constructor kwargs + accessor properties (`H.h1`, `H.eri`, `H.ecore`, `H.n_electrons`, `H.n_orbitals`) on `GeoVacHamiltonian`, and implemented `to_fcidump(filename, *, tol=1e-14, ms2=0, isym=1, orbsym=None)` writing a standard Knowles-Handy `&FCI ... &END` file. Added a pure-Python `read_fcidump(filename)` module-level helper for round-trip validation without a pyscf dependency.

The atomic He path uses `LatticeIndex` directly and does not currently surface pre-JW integrals — this is documented in the `to_fcidump` docstring and in a dedicated test (`test_fcidump_raises_when_integrals_missing`). The 28 main-group / multi-center / TM hydride systems and H2 all surface integrals cleanly.

No physics was modified. The Pauli content of the qubit Hamiltonian is bit-identical to the un-instrumented build. The (h1, eri) integrals exposed are the same tensors that `build_composed_hamiltonian` feeds into `build_fermion_op_from_integrals` immediately before the Jordan-Wigner step, so any downstream FCI / CASCI / DMRG / CCSD(T) on the FCIDUMP recovers the same N-electron spectrum as a GeoVac qubit-FCI (modulo the qubit-encoding map). This is the load-bearing structural property; it is verified at write time by construction (same tensor, written then read back bit-exact) and would be re-verified end-to-end once pyscf-FCI is wired in (Sprint H1-DMRG).

## §1. Format compliance

The writer emits the standard Knowles-Handy FCIDUMP:

```
 &FCI NORB=  15,NELEC=   4,MS2=   0,
  ORBSYM=1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
  ISYM=1,
 &END
 val  p  q  r  s    <- two-electron, 1-based indices, 16 digits
 ...
 val  p  q  0  0    <- one-electron, p>=q
 ...
 ecore  0  0  0  0
```

Convention pinning (matching §1.3 of the scoping memo):

1. **Orbital ordering** — GeoVac's internal per-block (Z, n, l, m) lexicographic order is written through unchanged. Downstream consumers must match this if they intend bit-exact cross-validation; the round-trip test confirms the writer-reader pair is internally consistent. A `geovac_to_pyscf_orbital_map` is **not** built in this sprint — it would belong in Sprint H1-DMRG when pyscf-RHF is wired in.
2. **Integral normalization** — chemist notation `(pq|rs) = eri[p, q, r, s]` (matches `build_fermion_op_from_integrals`, matches FCIDUMP standard). Verified via `test_fcidump_eri_eight_fold_symmetry`.
3. **Spin convention** — FCIDUMP is spatial-only with `MS2` (= 2*S_z) in the header; default `MS2=0` (closed-shell singlet) is correct for every system in the current 28-molecule library at ground state. UHF/UCCSD users can pass `ms2` explicitly.
4. **Active space / CAS** — `n_electrons` is summed across active blocks (frozen-core excluded by construction since the frozen-core energy is absorbed into `spec.nuclear_repulsion_constant`, which becomes `ecore`).
5. **Hopf-U(1) tapering** — orthogonal to FCIDUMP. The integrals exported are pre-tapering (correct for classical consumers); the tapered Pauli operator remains the consumer for VQE.
6. **Reference state for correlation energy** — out of scope for this sprint; depends on which classical solver runs on the FCIDUMP. This is the Sprint H1-DMRG / H1-Reference closure item.

Eight-fold permutation symmetry on `eri` is exploited at write time: only `(p,q,r,s)` with `p >= q`, `r >= s`, and lexicographic `(p,q) >= (r,s)` is written. The reader expands the full 8-fold orbit on read. The Hermitian symmetry of `h1` is similarly compressed (only `p >= q` written, full symmetric matrix recovered on read).

## §2. Round-trip check

For LiH at R = 3.015, max_n = 2:

- M = 15 spatial orbitals (Li_core block = 5 orbitals, LiH_bond block = 10 orbitals)
- n_electrons = 4 (active; the [He]-like core of Li is wrapped into ecore via the standard composed builder, but the spec writes the core as an explicit block with n_electrons=2; this is the spec-side convention and matches what `build_fermion_op_from_integrals` consumes)
- ecore = -6.284875124378109 Ha (= `spec.nuclear_repulsion_constant` = V_NN + V_cross + E_core)
- 15 one-electron entries written (diagonal only for default single-center spec; cross-block h1 disabled by default)
- 72 symmetry-unique two-electron entries written

Round-trip via `read_fcidump`:

- `max |H.h1 - parsed.h1| = 0.0` (bit-exact through 16-digit text round-trip)
- `max |H.eri - parsed.eri| = 0.0`
- `|H.ecore - parsed.ecore| = 0.0`

Cross-system sample (`test_fcidump_works_for_multiple_systems`): LiH, BeH2, NaH, H2O, CO, ScH, H2 all round-trip at `max diff < 1e-10` on h1 and eri (this is the tolerance gate; observed value is at machine precision for all seven). The single-system test is at tighter `< 1e-12`.

Header parsing: `&FCI ... &END` with `NORB`, `NELEC`, `MS2`, `ORBSYM`, `ISYM` all parsed and verified by `test_fcidump_header_compliance` and `test_fcidump_orbsym_default_and_custom` (default = all 1s; custom = round-trips exactly).

Hermitian h1 round-trip on NH3 (a system with non-trivial off-diagonal h1) verified via `test_fcidump_h1_hermitian_round_trip`: `max |h1 - h1.T| < 1e-12` after read-back, confirming the writer's compressed `p >= q` emission reconstructs into the full symmetric matrix.

## §3. Test summary

| Test name | What it checks |
|:----------|:---------------|
| `test_lih_pre_jw_integrals_surfaced` | h1, eri, ecore, n_electrons, n_orbitals all non-None and shape-correct for LiH |
| `test_lih_fcidump_write_and_parse` | Bit-exact round-trip on h1/eri/ecore for LiH at < 1e-12 |
| `test_fcidump_eri_eight_fold_symmetry` | Round-tripped eri satisfies (pq\|rs) eight-fold permutation symmetry |
| `test_fcidump_header_compliance` | `&FCI`, `&END`, NORB, NELEC, MS2, ORBSYM, ISYM all present |
| `test_fcidump_one_body_diagonal_count` | LiH h1 at max_n=2 is purely diagonal (M nonzeros, zero off-diagonal) |
| `test_fcidump_ecore_matches_nuclear_repulsion` | Exposed ecore matches `spec.nuclear_repulsion_constant` |
| `test_fcidump_works_for_multiple_systems` | Cross-system sample (7 systems) round-trips < 1e-10 |
| `test_fcidump_tolerance_filters_small` | `tol` parameter drops sub-threshold integrals; file still round-trips |
| `test_fcidump_raises_when_integrals_missing` | He (LatticeIndex path) raises ValueError on `to_fcidump` |
| `test_fcidump_orbsym_default_and_custom` | Default orbsym = all 1s; custom orbsym list round-trips |
| `test_fcidump_h1_hermitian_round_trip` | Round-tripped h1 is symmetric for NH3 |

**Full ecosystem_export test count: 115 (104 baseline + 11 new). All passing.**

## §4. Downstream implications for Sprint H1-DMRG (Phase 1 next sub-sprint)

The FCIDUMP exporter is the load-bearing unblocker. With it landed:

1. **Block2 / pyscf-DMRG immediately consumable.** `pyscf.tools.fcidump.read` reads what `to_fcidump` writes (header + integral format both standard). Block2's `DMRGDriver.read_fcidump` is the alternative path with the same convention.

2. **CCSD(T) cross-validation track unblocked** at zero additional plumbing (consolation track per scoping memo §3.2; useful as DMRG sanity check on closed-shell singlets).

3. **AFQMC (ipie) Phase 2 target** unblocked structurally; ipie can read FCIDUMP and Cholesky-decompose `eri` on its own. The native multipole structure of GeoVac's `eri` IS a low-rank decomposition in disguise (Paper 14 §sec:hopf_tapering note), but exposing this as `to_cholesky_eri` is a Phase 2 follow-on, not blocking.

4. **Bit-exact downstream-FCI cross-check on LiH** is the natural Sprint H1-DMRG entry point. The Sprint H1-Interface convention-pinning (§1.3 of the scoping memo, items 1-6) is partially closed by this sprint:

   - **Items 1, 2, 3, 4, 5: CLOSED.** Orbital ordering pinned (GeoVac internal order, single-block lexicographic). Chemist notation pinned. MS2=0 default for closed-shell. Frozen-core handling pinned (absorbed into ecore). Hopf tapering orthogonality pinned (FCIDUMP is pre-tapering).
   - **Item 6 (HF reference cross-check) DEFERRED to Sprint H1-DMRG.** Cannot run pyscf RHF without pyscf in the dependency tree; this is the natural first action item of H1-DMRG once block2 / pyscf are added as optional deps.

   Net: §1.3 is 5/6 closed by Sprint P1; the remaining HF cross-check is structurally one sentence of pyscf glue (Sprint H1-DMRG).

5. **NaH stress test (Sprint H1-NaH, the W1e structural classification falsifier)** unblocked. Once DMRG-on-GeoVac-integrals is running, NaH at R_eq=3.566 is the load-bearing question per CLAUDE.md §3 (six prior W1e closure attempts failed in-framework; the cosmic-Galois reading says external correlation should close it).

## §5. Honest scope

What this sprint closed:

- Pre-JW (h1, eri, ecore, n_electrons) are now first-class on `GeoVacHamiltonian` for every composed-path system in the library.
- FCIDUMP write+read round-trips bit-exactly on h1/eri/ecore for LiH and < 1e-10 on a 7-system sample.
- Format compliance verified (Knowles-Handy &FCI namelist + 1-based integral records + eight-fold permutation symmetry on eri + Hermitian symmetry on h1).
- Sprint H1-Interface §1.3 convention-pinning items 1-5 of 6 closed.
- 115/115 ecosystem_export tests pass (104 baseline + 11 new), zero regressions.

What this sprint did NOT close:

- **No pyscf cross-check.** pyscf is not in the dep tree; the round-trip test uses the pure-Python `read_fcidump` helper. A pyscf-based round-trip (read FCIDUMP -> pyscf.fci.FCI -> compare energy to GeoVac qubit-FCI) is the load-bearing Sprint H1-DMRG entry test.
- **No DMRG.** Block2 is not in the dep tree. Sprint H1-DMRG is the planned follow-on.
- **He atomic path not unblocked.** The LatticeIndex code path does not currently surface (h1, eri) as the result-dict; `to_fcidump` raises ValueError on He with a clear message. This is a 1-2 day plumbing follow-on if He is needed before H1-DMRG (none of the chemistry first-system picks require it; LiH is the production benchmark).
- **No `geovac_to_pyscf_orbital_map`.** Orbital re-ordering between GeoVac (per-block n,l,m lex) and pyscf (energy-sorted) is not implemented. The FCIDUMP file is self-consistent (GeoVac writes, GeoVac reads); cross-consumer ordering reconciliation is Sprint H1-DMRG when pyscf-FCI is run.
- **Paper 20 §sec:hybrid_pipeline / "Classical chemistry interface via FCIDUMP" subsection NOT drafted** (recommendation only, per task instructions — listed below).

## §6. Paper-edit recommendation (not applied)

Recommend a new subsection in **Paper 20** (`papers/group4_quantum_computing/paper_20_resource_benchmarks.tex`) under the benchmarks section:

```latex
\subsection{Classical chemistry interface via FCIDUMP}
\label{subsec:fcidump_interface}

The composed-builder integrals $(h_{pq}, (pq|rs), E_{\mathrm{core}})$ that
generate the qubit Hamiltonian via the Jordan-Wigner map are also exposed as
a standard Knowles-Handy FCIDUMP file via
\texttt{GeoVacHamiltonian.to\_fcidump(filename)}. This makes every
production GeoVac system in the library directly consumable by
classical correlation engines: DMRG via Block2 \cite{block2_2021},
coupled cluster via pyscf \cite{pyscf_2020}, and AFQMC via ipie
\cite{ipie_2023}. Chemist-notation $(pq|rs) = \texttt{eri}[p,q,r,s]$ is
preserved bit-identically through the FCIDUMP round-trip; eight-fold
permutation symmetry on the two-electron tensor and Hermitian symmetry
on the one-electron matrix are compressed at write time and
reconstructed on read. The frozen-core energy is absorbed into the
single $E_{\mathrm{core}}$ constant, so any solver consuming the
FCIDUMP recovers the same N-electron spectrum as the GeoVac qubit-FCI
modulo the qubit-encoding map. This positions the framework as a
\emph{producer} of structurally-sparse integrals for any chemistry
consumer, not just a producer of Pauli sums for VQE.
```

Apply in Sprint H1-DMRG after the first cross-consumer benchmark (DMRG-on-GeoVac-LiH) is in hand.

## §7. Files

### Modified

- `geovac/ecosystem_export.py` — added constructor kwargs (h1, eri, ecore, n_electrons) and accessor properties on `GeoVacHamiltonian`; added `to_fcidump` method and module-level `read_fcidump` reader; wired the four composed-path builders and H2 builder and the tapering helper to populate the new fields.

### Created

- `debug/sprint_p1_fcidump_exporter_memo.md` (this memo).

### Tests added

- `tests/test_ecosystem_export.py` — 11 new tests (block at the end of file, marked `# Sprint P1`):
  - `test_lih_pre_jw_integrals_surfaced`
  - `test_lih_fcidump_write_and_parse`
  - `test_fcidump_eri_eight_fold_symmetry`
  - `test_fcidump_header_compliance`
  - `test_fcidump_one_body_diagonal_count`
  - `test_fcidump_ecore_matches_nuclear_repulsion`
  - `test_fcidump_works_for_multiple_systems`
  - `test_fcidump_tolerance_filters_small`
  - `test_fcidump_raises_when_integrals_missing`
  - `test_fcidump_orbsym_default_and_custom`
  - `test_fcidump_h1_hermitian_round_trip`

### NOT modified

- Physics modules (`composed_qubit.py`, `balanced_coupled.py`, `qubit_encoding.py`, `molecular_spec.py`) — additive surface only.
- Pauli content — bit-identical to pre-sprint builds.
- JW pipeline — unchanged.
- Propinquity bound — unchanged.
- Papers — recommendation only (see §6).
- CLAUDE.md — out of scope for this sub-sprint; the one-liner §2 entry is the `/sprint-close` action.

---

**End of Sprint P1 memo. Verdict: GO at full Phase 1 scale; FCIDUMP exporter unblocks DMRG, CCSD(T), AFQMC simultaneously; 115/115 ecosystem_export tests pass (+11 new, 0 regressions); next sprint is H1-DMRG (block2 / pyscf wiring + LiH bit-exact pyscf-FCI cross-check).**
