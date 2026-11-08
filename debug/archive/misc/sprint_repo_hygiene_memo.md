# Sprint repo-hygiene + WH7 registration (2026-06-10, v3.110.0)

**Trigger:** PI concern that the repository "has become a graveyard of dead ends that
threaten to bury the model for anyone attempting to get into it," plus PI direction to
register WH7.

**Diagnosis:** the dead ends as *content* are fine — protected institutional memory
(§13.5) with an existing archive pattern (`tests/_archive/`, `geovac/_archive/`,
`papers/archive/`). The actual burial risk was structural: (a) `debug/` had 1,824 files
at one flat level (677 memos) with no active/closed distinction; (b) the file tree gave
newcomers no way to tell the ~12 keystone papers from the ~50 that are scaffolding,
observations, or descoped; (c) ~150 LaTeX build artifacts cluttered the paper folders
and repo root.

## Actions

1. **debug/ archive sweep** (`archive_sweep_2026_06_10.py`): 1,367 files moved into
   `debug/archive/<arc>/` (12 arcs). Keep rules — mtime ≥ 2026-06-01 (June frontier,
   407 files), referenced by any file in `tests/` (48 files, derived mechanically by
   regex scan so frozen-falsifier provenance can't break), pinned (3). Top level:
   1,824 → 458. Manifest `debug/archive/sweep_manifest_2026_06_10.json` maps every
   old path → new path. `pytest` already ignores `debug/` (pytest.ini), so no
   collection impact. Verified: frozen falsifiers + topological proofs green post-move.
   The `misc/` bucket (536 files) is classification-by-filename fallout; the manifest
   is ground truth. Script is re-runnable with a later CUTOFF as arcs close.

2. **papers/INDEX.md**: status map (KEYSTONE / ACTIVE / OBSERVATION / GUARDRAIL /
   PARTIAL / DESCOPED / DRAFT / HISTORICAL) for all ~60 papers + five-document reading
   path. Decision: descoped Lorentzian papers (45–49) are NOT physically moved — they
   are corrected in place per the de-versioning directive, cross-referenced from other
   .tex, and their Status notes are part of the honest record. The INDEX carries the
   status signal instead.

3. **README**: version badge → 3.110.0; "New here? Five documents, in order" block
   after the install line, pointing at field guide → claims register → INDEX →
   Paper 7 → Paper 14, with an explicit sentence that `debug/archive/` and
   `papers/archive/` are the record, not the product.

4. **Build-junk purge**: 147 untracked LaTeX intermediates (.aux/.log/.out/...)
   deleted from `papers/` and root. All were already gitignored; zero git impact.

5. **WH7 registered** (CLAUDE.md §1.7, explicit PI direction per §1.7 governance):
   *time-discreteness is observer-compactification.* Full entry in §1.7 with primary
   falsifier (Toeplitz temporal-compression program) and secondary falsifier
   (inherited Paper 35 π-prediction). Memory note
   `memory/wh7_time_observer_compactification.md`.

## Verification

- `pytest tests/test_p38_action_seminorm.py tests/test_p45_kplus_degeneracy.py
  tests/test_fock_projection.py tests/test_fock_laplacian.py` — green post-sweep
  (see session log).
- No production code in `geovac/` touched; no paper content touched (INDEX/README are
  documentation).

## Follow-ons

- The sweep's `misc/` bucket could be re-classified opportunistically when arcs are
  touched; not worth a dedicated pass (manifest already resolves lookups).
- Re-run the sweep with CUTOFF bumped when the June frontier closes (next natural
  point: after the outreach sends or the WH7 falsifier sprint).
