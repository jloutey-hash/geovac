# group2 — debt (SMALL findings, owed guards, completeness gaps)

Owner of this list: the PM. SMALLs are fixed in one scheduled sweep outside a cert cycle (edit
bar, `docs/qa/criteria.md`). Guards are written as their own activity and fire-tested.
Source run: `/qa group2` FULL (incremental, first whole-group cert under stamps), 2026-11-03/04,
snapshot 706c5317 (`refs/qa/group2-full-2026-11-03`).

## Owed guards (separate activity; fire-test each against the wrong answer it excludes)

OWNER-LARGE "test proves less than the prose / no test" findings and the guards the remediation
needs. The owner units stay unstamped until the guard exists.

- **Paper 8:** pin <1s|1/r_B|1s> R-dependence (0.896/0.610/0.328/0.125) and <1s|1/r_B|2s> != 0, so
  "exact" cannot return. `tests/geovac_paper8_tests.py` is never collected (name) — rename.
- **Paper 11:** pin |E(R=2) - (-0.6026342144949)| < 1e-12 at n_basis=20 and R_eq vs 1.99719 at 1e-4
  (current test: 0.01% / 1%; fire_test +1.21e-6 did not fire).
- **Paper 12:** quad-spy guard on `neumann_vee.compute_Bl_table`; update `test_Bl_against_quad`
  docstring (now an independent check); the 79.6/80.1% numerical comparator (open since
  2026-09-19); the r12^{0-1} grid table; HeH+ R_eq.
- **Paper 13:** He Schwartz coefficient vs Bromley-Mitroy; wall-converged coupled channel
  (R_min = 0.02); Li nu0 = 2; `test_h2_rovib_morse` bands/docstring to the recomputed values; the
  spectral radial test skipped although `solve_radial_spectral` exists (un-skip).
- **Paper 15:** Level-4 R_e,min wall dependence; Hellmann-Feynman P formula; 96.0% (no test since
  2026-09-14).
- **Paper 17:** `test_ab_initio_pk_v2` band 2.5-4.0 -> pin 3.222 +- 0.01 and omega 1465 (fire_test
  B inv2 -> median did not fire); l-dependent 3.176; Outlook l_max 3-4 (slow); tab:B candidates;
  4N 63.5% unbound (asserts has_minimum only).
- **Paper 19:** R-independence of the three non-balanced schemes; n_max = 3 energy on the corrected
  phase (-8.054515).
- **Paper 58:** `test_paper58_eri_census` density band 25-35 -> exact nz == 2976 against a real-basis
  symmetry count (fire_test -400 did not fire).
- **Paper 59:** eq:pf exact identity (sympy: integrand = -d/dx[sqrt(Q) e^{-Dx}]); euclidean dictionary
  at five points; the 1.8e-10 position-space fibre.
- **Paper 60:** 1-norm exponents 0.8193 / 0.7261 on the converged box (bands only; docstring cites
  the retired 0.842 / 0.78); conditioning numbers (1.47e6, N^1.57, 4.2x); preconditioner lemma guard.
- **Paper 62:** JW Pauli count 107,382 for an exponential basis (registry `p62_dzp_pauli`); HF
  canonical count equality; kit 1-norms 405/511/620.
- **FCI-M / FCI-A:** physical-h1 n_max = 2 LiH PES pin (monotonic, E(8) < E_sep - 0.2, BSSE -0.0884
  R-independent), fire-tested against a hop-on build; lambda = 0.1 interior minimum; zeta_B
  relaxation; Hypothesis-2 shifts; FCI-A NNZ vs Table III.
- **C25 `he-bracket`** was edited by a fixer in the remediation pass (anchor 0.020%, overclaim
  0.004|marginally below|bracket): review it as a guard and fire-test against the v6.1.5 text.

## Pre-existing test rot (found, not caused, by this run)

- `tests/test_composed_cusp.py` (slow, 9 errors): passes `cusp_correction=` to
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
