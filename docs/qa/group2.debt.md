# group2 — debt (SMALL findings, owed guards, completeness gaps)

Owner of this list: the PM. SMALLs are fixed in one scheduled sweep outside a cert cycle (edit
bar, `docs/qa/criteria.md`). Guards are written as their own activity and fire-tested.
Source run: `/qa group2` FULL (incremental, first whole-group cert under stamps), 2026-11-03/04,
snapshot 706c5317 (`refs/qa/group2-full-2026-11-03`).

## Owed guards -- CLOSED 2026-11-04 (v6.2.1 guard pass)

Written as their own activity, fire-tested against a named wrong answer, then independently reviewed (fresh
code-reviewer). Files: `test_paper11_h2plus_pins.py`, `test_paper17_lih_headline_pins.py`,
`test_paper19_scheme_r_dependence.py`, `test_paper58_eri_census.py` (exact 2976),
`test_fcim_physical_h1_pes.py`, `test_fci_atoms_nnz_table3.py`, `test_paper8_sw_model_vs_genuine.py`,
`test_paper8_so4_algebra.py` (renamed; now collected), `test_paper13_schwartz_coefficient.py`,
`test_paper13_wall_and_li_nu0.py`, `test_h2_rovib_morse.py`, `test_algebraic_curve.py` (spectral un-skipped),
`test_paper15_wall_and_hf_coupling.py`, `test_paper12_bl_no_quadrature.py`, `test_paper59_pf_dictionary_guards.py`,
`test_paper60_owed_guards.py`, `test_paper62_dzp_headline.py` + `test_registry_crosscheck_p62.py`; C25
`he-bracket` fire-tested against the v6.1.5 text. Rows in `docs/claim_test_matrix.md`.

Review escapes fixed and re-fired: the Paper 12 quadrature spy missed a name bound before the monkeypatch (a
profiler spy added); the Paper 17 windows (+-0.01 bohr / +-10 cm^-1) passed a 2 % change of the PK width
(tightened to the printed precision). C22 caught a debug/ import in the Paper 62 paper-backing test (moved).

WEAK, open: `test_paper8_so4_algebra.py` is tautological (Pauli-matrix construction, no geovac import);
`test_paper8_sw_model_vs_genuine.py` hard-codes the model side; the FCI-M monotone / E(8) subtests hold for the
hop-on build too (physics checks, not discriminators); Paper 15's spectral 92.66 is a code value, not printed;
Paper 60 conditioning guards exercise test-file helpers; C25 `he-bracket` row-6 qualifier satisfied by the
raw-figure clause, owner window +-1 line.

## Pre-existing test rot (found, not caused, by this run)

- [RESOLVED v6.2.1: archived to tests/_archive/superseded/, feature removed in v2.7.0] `tests/test_composed_cusp.py` (slow, 9 errors): passes `cusp_correction=` to
  `ComposedDiatomicSolver`, which no longer takes it (both files last touched 2026-06-05).

## Routed / PI

- **Paper 16 vs Paper 13:** Paper 13 now has nu0 = 2 for Li 2S^e (nu = 1 harmonics are L=1, odd);
  Paper 16's nu_min = 1 is the Pauli-filling count without the L/parity constraint. Consistency
  check owed at the Paper 16 owner (group3).
- **Paper 15 H2 cusp model** (hand-set beta = 0.60, gamma = 0.3; disclosed): keep (tail 0.02 pp) or
  drop; recommendation keep, disclosed. The "2D+cusp" extended table's provenance (old coefficient?)
  is unverified.
- **CLAUDE.md §1.5** "zero-parameter construction" may need scoping for composed LiH after the
  Paper 17 disclosure (PI section).
- **Prospective re-derivation** of the PK B formula and the H2 cusp coefficient under the
  prediction-before-results protocol (PI raised 2026-11-04).

## SMALL sweep (edit bar) -- OPEN
- From the guard pass (2026-11-04): Paper 11 R_eq 1.9973 is a quadratic-fit bias (Brent minimizer 1.9971933, 2e-6
  relative — the prose understates); Paper 58 "29.8% measured vs 29.4% counted" -> exact 2976 = 2976 in matched
  real-Cartesian bases (also the test docstring); Paper 59 position-space fibre 1.8e-10 -> ~1e-15; Paper 60 "~4x the SW
  prefactor" is a condition-number ratio (prefactor ratio 9.2) and SW N^1.80 vs prose 1.85; Paper 12 Table I E_num
  N=27 -1.138836 vs recomputed -1.1388188; Paper 19 `build_balanced_hamiltonian` nuclear_repulsion carries the Li+ core
  energy (-7.2799) while all four electrons are explicit (disclosed offset; double-count-shaped convention, code
  review owed).

- Paper 8: body L418/L466 s-p mixing vs thm:selection; two beta conventions (thm:uniform vs
  corollary); thm:lih_beta is a measurement labelled Theorem, values drift (1.0306/1.9480/1.9732);
  history notes L79/L380/L1059; Avery2004IJQC pages.
- Paper 11: table L841 three-digit comparison; Lambda-selector tag.
- Paper 12: history notes ~L222/L934/L1307/L1737; 99.81% / 0.053 mHa bracketed not pinned; the
  headline r12 engine lives in `debug/prolate_r12_mpf.py`; appendix "0.01 Ha per element" vs 10%.
- Paper 13: DBOC "Berry phase" (L733), "non-trivial bundle" (L1415), L1462 gauge-invariant phase,
  L1413 cusp signature, mu(R) algebraic scope, Paper 15 attribution (L1687), sec ref (L1986),
  Gamma ratio provenance (L989); figure in `debug/plots/`; R50 test near-tautological.
- Paper 15: l_max=2 jump stated three ways; 20x/269x only >10x; figure scripts in `debug/`.
- Paper 17: L1322 "validates", L1691/L1773 "structurally optimal"; W1c/W1d/W1e block jargon + debug
  refs; L526 D_e grid claim; natbib "Note1" undefined footnote warning; rows not reproduced
  (R-dependent / self-consistent PK, 2D-solver rows of tab:solver_pk).
- Paper 19: L95 "convenience not necessity"; Discussion L2462 stale prospective text; L2394 "second
  packing axiom"; W1c jargon; n_max = 4 decider cites the engine test only; Route C cite ->
  `test_lih_golden.py`; "variational bound respected" (L633/L2486).
- Paper 58: K58-1 "at the level the basis convention permits" (exact in matched bases, 2,976);
  0.80 / 0.52 / 0.011 Ha unpinned; Eq. (1) "every census entry"; Tab II LIVE rows; 5x7 robustness
  driver-only; secondary numbers.
- Paper 59: integrated-class guard on frozen data; sec:scope cites records not tests; routeC
  docstring stale; zstar tautology; Matuzas 2026 multicentre? (abstract only read).
- Paper 60: preconditioner guard weak; docstring 0.78.
- Paper 62: lesiuk2015 published title from memory (verify); intro motivated by Paper 59 integrals
  not used; L493 "one-centre structure" undefined.
- FCI-A / FCI-M: remaining hop-era history rows (FCI-M v0.9.18/22; FCI-A v0.9.11/16-22/30; "7.1 ->
  10.7 eV"); pre-existing correction annotations (FCI-A L182/L400/L629-641/Note, FCI-M L78 header
  comment "Topological"); Fig. 2 image stale (NNZ).
- Synthesis: abstract "only the axial m-rule survives"; "zero-parameter" framing; historical
  "Withdrawn." paragraphs.
- Citations: P8 Avery2004 pages; P13 Kereselidze2016 -> Mol. Phys. 113, 3471 (2015); P13
  Aquilanti2001 -> IJQC 92, 212 (2003); P19 herbst2018 orphan; P59 avery 2018 (done); Gershtein-
  Krivchenkov title.
- Front door: README / INDEX remaining LiH -8.0490 history mentions; LiH energies unregistered in the
  numeric registry (G1).

## Completeness-critic gaps (2026-11-04)

- C25 covers only Papers 13 and 59 in group2: register headlines of Papers 11, 12, 15, 17, 19, 58, 62
  (the class that let the synthesis/front-door drift through).
- Coverage gaps on abstract claims: Paper 13 rovib (benchmark only), Paper 19 8,370-Pauli prolate
  count and VMC twin (debug only), Paper 15 HeH+ (no matrix row), FCI-A Li/Be (pinned by
  `test_direct_ci.py` but no matrix row).
