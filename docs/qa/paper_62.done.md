# Paper 62 (Water at Equal Size: Slater versus Gaussian) — `/qa` profile

> **STATUS: NOT CERTIFIED -- never run (profile pre-registered 2026-11-01).**

> **Pre-registered 2026-11-01, before any review** (PI-delegated: "you can do the paper 62 dod"). Drafted by the
> PM; the PI may tighten or relax any line before the first run, but not during it (the /qa protocol, step 1). No
> run has happened. Not certified.

This file supplies only the Paper 62-specific scope and deltas. Everything else (C1-C25, LARGE/SMALL grading,
the edit bar, incremental stamps, the three-way verdict) is inherited from `docs/qa/criteria.md` and
`docs/qa/group2.done.md`, whose branch-defining criterion (benchmarking rule + guardrail-negative honesty) applies
here in full.

## Scope

**Paper:** `papers/group2_quantum_chemistry/paper_62_water_equal_size_slater.tex` (973 lines, 12 sections, 23
bibitems). Also in the `group2` scope. Deterministic gate: `--gate paper_62` (resolves to the `.tex` and its PDF).

**Summary surfaces (C9, GATING; measured non-empty 2026-11-01):**
- `papers/synthesis/group2_quantum_chemistry_synthesis.tex`, the subsection "Water at equal basis size" (from L974).
- `papers/INDEX.md` row 62 (L106) and the group2 row (L40).
- `CLAUDE.md` §1.5 composed-builder caveat and the §2 one-liners for v5.21.25 to v5.22.0.
- `docs/claim_test_matrix.md` rows `| 62 |` (L847 to L856 on 2026-11-01).
Do not record C9 as N/A.

**Records the paper rests on (C7):** `docs/predictions/2026-10-30_scale_ladder.md` (SL), its R3 addendum,
`2026-10-31_matched_accuracy_gaussian.md` (MA), `2026-10-31_hf_slater_second_molecule.md` (HS). Each is
pre-registered; the paper's "predicted before computing" sentences are checked against the record's commit order.

## Branch-defining criterion: **equal size is not equal accuracy, and canonical orbitals are not all orbitals**

The paper's central comparison is at equal basis size. Its cost comparison is at matched accuracy. Its Pauli-count
claim holds in canonical orbitals. A run PASSes this criterion only if all four hold:

1. **No sentence anywhere in scope** (paper, synthesis, INDEX, CLAUDE.md, matrix) presents an equal-size result as
   a matched-accuracy cost result, or the reverse.
2. **"Same Pauli count" carries "in canonical orbitals"** (or an equivalent scope) wherever it is stated as a
   finding. Koridon et al. 2021 show orbital rotations change the count; the abstract carries the qualifier, the
   synthesis and INDEX row 62 did not on 2026-11-01. Check each locus.
3. **"Nothing tuned on water" is true** for every number it is attached to. The Slater DZP's polarization
   exponents WERE tuned on water (fc-CCSD); only the ladder rungs are atom-set. Any sentence that blurs the two is
   a LARGE.
4. **Not a lattice result.** The paper says the bases are tuned exponential functions and the angular sparsity does
   not enter. Every summary of it must keep that (§1.5 rhetoric rule; the composed caveat of 2026-10-31).

## Dimensions to exercise (FULL run)

- **Code / test-backing (C1-C2)**, `code-reviewer` ×1, Opus. **Run with `--slow`** and report passed-vs-skipped per
  test. Files: `tests/test_water_ccsd.py` (solver), `tests/test_paper62_dzp_headline.py`,
  `tests/test_paper62_native_kits.py`, `tests/test_paper62_scale_ladder.py`, `tests/test_paper62_matched_accuracy.py`,
  `tests/test_paper62_hf_second_molecule.py`. The slow layer is about 2 h at 4 threads; budget it. Mandates: say
  for each headline whether its number is RECOMPUTED or READ as a literal; fire-test every guard.
- **Paper claims / prose (C3, C5, C6, C8)**, `claims-reviewer` ×1, Opus, enumeration-forced over every tier tag
  (9 bold tags on 2026-11-01: 4 MEASURED, 3 MEASURED pre-registered, 2 MEASURED negative) and every untagged
  numeric sentence in sec:native, sec:qubit, sec:prior and sec:limits.
- **External citations (C4)**, `citation-reviewer` ×1, Opus. Highest-risk items: `trenev2025` (the 107,382
  cc-pVDZ Pauli count is cited to its Table 5), `koridon2021`, `lesiuk2015`, `nemec2010`, `stark2025` (the
  -0.2558 Ha value is Table 3; geometry and basis are SI Tables S1-S2), `ema2003`, `forster2021`, `rico2004`,
  `dunning1989`, `crawford_projects`. Three-way verdicts; UNVERIFIABLE is not a pass.
- **Synthesis (C9)**, `claims-reviewer` ×1. GATING; see Scope.
- **Arithmetic audit**: every derived number (gaps, ratios, percentages, "covers about 94% / 45%") recomputed from
  the table it comes from. C21 covers the `\gvq` ones; the rest are by hand.
- **Deterministic (C10-C25)**, `--gate paper_62`. State grounded counts; record zero-surface criteria as
  UNMEASURED, not PASS (measure the surface at run start; do not copy it from here).
- **Completeness-critic** ×1, briefed from the reviewers' own reports.

## C8 headline claims (frozen goalposts, with backing as of 2026-11-01)

| # | Claim | Tier | Backing |
|---|---|---|---|
| 1 | Slater DZP 40.7 mHa below cc-pVDZ at fc-CCSD(T), 28.8 at SCF, 24 functions | MEASURED | RECOMPUTED: `test_paper62_dzp_headline.py` (slow), cc-pVDZ also in `test_paper62_matched_accuracy.py` |
| 2 | Fit converged: 12 -> 16 Gaussians moves every energy 0.07 mHa | MEASURED | RECOMPUTED: `test_paper62_dzp_headline.py::test_twelve_gaussians_are_converged_to_a_tenth_of_a_millihartree` (slow) |
| 3 | Hydrogenic Z/n kit 196 mHa above cc-pVDZ | MEASURED, negative | RECOMPUTED: `test_paper62_native_kits.py` (slow) |
| 4 | Single-scale Sturmian more than 1 Ha above; two-scale kit 130 mHa above | MEASURED, negative | two-scale = ladder R2 (`test_paper62_scale_ladder.py`); single-scale RECOMPUTED: `test_paper62_single_scale_sturmian.py` (atom fast; water slow, 878 s; added 2026-11-02) |
| 5 | Ladder: R5 3.1 mHa below cc-pVDZ, R5H 7.0; SL1 0.900, SL2 0.516 | MEASURED, pre-registered | `test_paper62_scale_ladder.py` |
| 6 | R3 diagnosis: correlation loss transfers one to one; 3d passenger; HF pieces cancel | MEASURED | `::test_r3_correlation_loss_is_a_one_centre_effect_carried_by_the_3d` (slow); the three water substitution runs are record-only |
| 7 | Pauli count 107,382 on 46 fc qubits for every 24-function basis, in canonical orbitals; 1-norm 405 to 706 Ha | MEASURED | counter vs openfermion (`test_water_ccsd.py`), cc-pVDZ count recomputed (MA guard); **the other kits' counts and the 1-norm range are READ** |
| 8 | MA: aug-cc-pVDZ 8.0 mHa short at 9.2x terms, 2.7x 1-norm, 80 vs 46 qubits; crossing before cc-pVTZ | MEASURED, pre-registered | `test_paper62_matched_accuracy.py` (aug recomputed; DZP target now recomputed by #1) |
| 9 | HF: Slater DZP 59.5 mHa below cc-pVDZ at 19 functions; aug-cc-pVDZ short at 9.3x terms | MEASURED, pre-registered | `test_paper62_hf_second_molecule.py` (Gaussian rungs recomputed; **Slater DZP READ**, 23 min) |
| 10 | External check: E_corr in Stark et al.'s TZ set agrees to 0.3 mHa | MEASURED, reported comparison | **record-only** (2.5 h): `debug/data/stark_zimmerman_hf_tz.jsonl` |
| 11 | Water-tuned DZP stays 34 mHa below the best atom-set rung; 21.7 with the atom-set 3d | MEASURED (cap) | derived from #1 and #5 (C21); the 21.7 from #12 |
| 12 | SL-3d: the 3d set on the atom's fc-CCSD(T) (Dunning's rule); transfer within 15% on four rungs; R5Hd 19.0 mHa below cc-pVDZ | MEASURED, pre-registered | `test_paper62_scale_ladder_3d.py` (ratios, band, SL3d-M, R5Hd position derived; R3 atom optimum and R3d water recomputed, slow; other rungs READ) |
| 13 | H6 ring (HR): Slater DZP 3.5 mHa below cc-pVDZ* at 30 functions, short of the 6.2 bar (STOP); cc-pVDZ+TZpol and aug-cc-pVDZ beat it, cheapest at 6.3x terms; canonical counts spread 1.4%, fixed labels equal; ring labels 5 vs 3 Z2, 4.37x fewer terms at 1e-10, site fewer at 1e-2 | MEASURED, pre-registered (STOP) | `test_paper62_hydrogen_ring.py` (all derived fast; cc-pVDZ* recomputed by the original solver and both labellings rebuilt, slow); **Slater DZP row READ** (45 min), 48-84 function rows READ |
| 14 | Torus (TR): toroidal winding exact; Z2 2/4/5 even N, 2/3/4 odd N; poloidal leak slopes rise toward 1 and 2; full/tor fraction fixed (observation) | MEASURED, pre-registered (GO) | `test_paper62_hydrogen_torus.py` (6x3, 7x3, slow 8x3 at the control geometry; slow: the scored step N = 10 -> 11 at s = 1.4 recomputed against the stored rows and TR3a's band); the other scored values and the full/tor ratio READ |
| 15 | Squashed ring (SQ): ring momentum dies as eps^p(c), charge-2 counting, wrap-around orders at odd N | MEASURED, pre-registered (GO) | `test_paper62_squashed_ring.py` (N = 7, 8 rebuilt, ten slopes); scored N = 9-13 values READ |

## Watch-notes (things the PM already suspects, ranked)

1. **Canonical-orbital qualifier.** The synthesis (water paragraph, first Pauli mention) and INDEX row 62 gained it
   2026-11-01, before any run. Check the remaining "same Pauli count" loci (synthesis HF sentence, CLAUDE.md §2).
2. **READ literals.** Rows 4, 7, 9 and 10 still rest partly on literals. A reviewer should say for each whether the
   prose claims more certainty than a READ number supports.
3. **sec:limits** says the solver is validated against an external reference only on STO-3G water and that the one
   larger-basis check is correlation-only on HF (edited 2026-11-01). Check that the solver section and the limits
   section say the same thing.
4. **Prior art (sec:prior)** was narrowed 2026-10-31 after a scan (Lesiuk, Nemec, Koridon credited). The claim of no
   prior equal-size Slater-vs-Gaussian qubit-cost comparison is a search result, not a proof; the prose must say so.
5. **R3 3d exponent rule** is an open PI decision (v5.21.28). The paper must present the tie to the 3p as a design
   flaw it diagnoses, not as a choice it defends.
6. **One geometry, one method per molecule.** Every summary keeps this cap.

## Change log

- 2026-11-02: rows 4 and 13-15 **approved by the PI** ("Ok, go for it", on the PM's recommendation to approve them
  as drafted). Row 14's backing gained the scored-step test the same day.
- 2026-11-02: rows 13-15 added and row 4's backing updated (PM draft, PI-directed "We should do all of these";
  **PI to approve before any run**). Sec. hring (HR, TR, SQ) postdates the 2026-11-01 draft; its abstract and
  conclusion now carry the H6 ring result, so rows 13-15 are summary-surface loci too.
- 2026-11-01: drafted (PM, PI-delegated), before any run. Backing column measured the same day, after three new
  guards (`test_paper62_dzp_headline.py`, `test_paper62_native_kits.py`) and the Stark-Zimmerman matrix row.
