# paper_62 — debt (SMALL findings, owed guards, completeness gaps)

Owner of this list: the PM. SMALLs here are fixed in one scheduled sweep outside a cert cycle
(edit bar, `docs/qa/criteria.md`). Guards are written as their own activity and fire-tested.

## SMALL sweep (edit bar) -- paper_62 claim-search DELTA 2026-11-11 -- OPEN

- **Conclusion L1091-1093.** "and more than the equal-size count" is unscoped; true only for the bases tested
  (the ANO crossing between 24 and 35 functions is not located). Add "among the bases tested". (L786 is fine in context.)
- **Mechanism sentence L277-281.** "the exponential shape buys the mean-field part" is contradicted by ANO-R1, which
  recovers 21.0 of the 28.8 mHa RHF gap. Scope it "relative to cc-pVDZ".
- **Abstract L56-59.** "every twenty-four-function basis" needs "tested" / "here".
- **Prior art L986-989 (Lesiuk).** Scope the contrast to the water-tuned DZP.
- **zobel2020 (L741-743).** Cited for "ANO-R1 is a column truncation of ANO-R2"; that is our own numerical check
  (`tests/test_paper62_ano_matched_accuracy.py`). Reattribute. (Stored ANO-R2 file: O [4s4p2d1f], H [3s2p1d].)
- **HF record.** `docs/predictions/2026-10-31_hf_slater_second_molecule.md` L118-121 "the Gaussian basis needs more
  than aug-cc-pVDZ" is unscoped (cc-pV only).
- **Claim matrix row (SZ), `docs/claim_test_matrix.md` ~L898.** Says "AM ... not yet scored" (stale) and "every basis
  107,382 Pauli terms ... recomputed" (only ANO-R1's count is recomputed).
- **group2 synthesis ~L967.** Drops the "ANO-R1 falls short / crossing not located" clauses.
- **Paper 14 L1239-1243.** The equal-Pauli-count restatement lacks the canonical-orbital qualifier.
- **`docs/predictions/2026-11-01_hydrogen_ring.md` L230.** Omits HF from the molecules the ring contrasts with.

NITs: abstract L60-61 the 405-706 Ha 1-norm range includes cc-pVDZ's 674 without saying so; L32/L288
"cc-pVDZ is built for systematic convergence" uncited under a MEASURED tag (cite dunning1989); test docstrings
`tests/test_paper62_samesize_gaussian.py` L91-92 ("best same-size Gaussian" -> "tested"),
`tests/test_paper62_matched_accuracy.py` L102-103 ("the cc-pV matched-accuracy basis").

## Guards owed (separate activity; fire-test each against the wrong answer it excludes)

- `tests/test_paper62_samesize_gaussian.py` does not assert the pcseg-1 / def2-SVP Pauli counts, qubits, 1-norms or
  RHF energies (fire-test plants on those were asleep); the 21.0 mHa ANO-R1 RHF gap is unpinned.
- The "equal cc-pV Pauli counts" clause for the A35/A41 rungs is untested; A41 is pinned by C21 only (disclosed).

## OpenFermion accumulation drop (found 2026-11-11, v6.7.16)

- `jordan_wigner(build_fermion_op_from_integrals(...))` loses coefficients assembled from pieces below OpenFermion's
  `EQ_TOLERANCE` (1e-8): water Slater DZP 107,378 vs the true 107,382. Paper 62's counts use the closed form and
  stand; the library's composed LiH/H2O are unaffected. Audit any other registered Pauli count produced through
  OpenFermion on a large Gaussian-basis Hamiltonian (rebuild at a 2^20 scale and compare). PI call on scope.

## Surface staleness outside the paper (corpus-wide; PI call on scheduling)

- `viz/public/papers/*.html`: Paper 62's page was rebuilt 2026-11-11; the other 65 were last built 2026-11-09 and
  may carry pre-change abstracts. Five more (Papers 12, 14, 20, 23, 60) render raw `\gvq` macro text -- the builder
  fix (`debug/zenodo_manifest_build.py` `unwrap_registry_macros`) lands with v6.7.15, but their pages and manifest
  entries are not yet regenerated. Regenerating also refreshes the Zenodo descriptions the PI uploads.
- Reviewer-suggested C25 widening not yet done: scan docs/predictions Outcome sections and viz pages;
  `needs?\s+more\s+than\s+aug-cc-pVDZ\s+to\s+(match|reach)`; `brackets\s+the\s+Gaussian\s+cost` without a family name.
