# group4 — debt (SMALL findings, owed guards, completeness gaps)

Owner of this list: the PM. SMALLs are fixed in one scheduled sweep outside a cert cycle (edit bar,
`docs/qa/criteria.md`). Guards are written as their own activity and fire-tested.
Source run: `/qa group4` FULL (incremental), 2026-11-04, snapshot c613a229 (`refs/qa/group4-full-2026-11-04`).

## Owed guards (separate activity; fire-test each against the wrong answer it excludes)

- **Paper 23 (encoder fix, v6.3.0):** fire-test the new spectrum-vs-FCI tests (deuteron, He-4 with and without
  Coulomb), the composed Pauli-vs-matrix test and the OpenFermion round trip against the retired real-coefficient
  XY/YX encoding; `_finalize_hermitian_pauli`'s raise; the updated pins (904 / 641.39 / 632.19 / 383.87); the
  Minnesota V_S = -91.85 pin against -91.4.
- **Paper 14:** regenerate `test_paper14_mpo_rank` JSONs (pair-diagonal vintage: 334 / 878, {9,16}) and pin the live
  profile {16,15,8,15,20,19,10,5}, the balanced lift 15 and DF rank 15; a test for tab:accuracy (LatticeIndex exact,
  n_max 2-5: 2.43/2.10/1.99/1.93%); the per-block tapering counts (838 -> 424, ...); re-pin or retire
  `test_paper14_rel_lambda` (retired 40.59/143.96/18.68); `test_ecosystem_export` ratio / isostructural tests (35 slow
  failures pinning 11.10 / 9.23 / 222...); `test_gaussian_reference` fixture (`fci_method=` no longer accepted); the
  TC "standard converges 5.3 -> 2.0%" sentence; the lambda band [1.72,1.90] cannot tell 1.774 from 1.784.
- **STO-3G baseline:** an artifact + test for LiH STO-3G raw JW 631 / QWC 179 / lambda 16.48 (12.34 non-identity)
  (drivers in the run scratchpad: lih_sto3g.py, sto_ni.py).
- **Paper 20:** `test_paper20_library::test_beh2_h2o_qpe_regime_one_norms_pinned` (slow, red, pins retired
  2652/306.4/5798/1511.1/373.4 at the wrong BeH2 geometry); tab:pes_nmax2 live energies; the 2.816 two-point exponent;
  SrH/BaH/multicenter systems; tapering 54%.
- **Paper 16:** `test_paper16_dirac_metric` cannot fail (recomputes delta = x/(1-x) and a Z-independent nu);
  `test_structure_type` is a lookup of the hardcoded table (tautological).
- **Paper 23:** the Fock rigidity theorem has no test (now [PRIOR ART]).
- **Paper 27 C17 entropy families** (`p27-ep2b-entropy-nmax*`): the exemption vocabulary includes "corrected",
  which tab:ep2b's caption carries, so a drift inside that table can be exempted; switch to same-line exemption and
  fire-test. Also fire-test the new `p27-ep2b-vs914-retired` family and the inverted `p27-minnesota-contrast`.

## Cross-group consequences

- **Paper 34 (group6):** its §V He 2^3P precision entry rests on Paper 14's relativistic rationals, which are a fit
  to NIST (selected by rational search, `geovac/breit_integrals.py`, Track BF-D). Route to the group6 run.
- **Papers 24 (group3) and 27 (group6):** the Minnesota V_S correction moves the HO two-fermion entropy values
  (followed through in v6.3.0).

## SMALL sweep (edit bar) -- OPEN

- Paper 14: abstract L85 ref; basis-axis exponents called "system size"; "Q^3.17 floor"; "~50x-320x" without the
  matched-accuracy caveat (ft_resources); d-block s+p row vs prose; "18 ... extended to 28"; BeH2 11.7 vs 32%;
  balanced 0.20% "composed accurate"; nested "correct binding"; M = n_max^2; LiH 837 vs 838 identity mixing in prose;
  figure title "1/M^2 decay" vs caption, generator hard-codes data (provenance in debug/); jargon; debug/ cites;
  citations: trenev bibitem order; "JW with 2-qubit reduction" attributed to Trenev (only JW stated); zeta 1.6875
  He "STO-3G" label; qubitization O(lambda t/eps); ChildsBerry lambda; Chawla S1-S3 contents, author list, journal
  ref; CaH 0.47% is a dipole; BLISS-THC 130.9.
- Paper 16: eq:metric_correction undefined source; acknowledgments credit; classifier `.nu` lowest-spin S_N-only vs
  tab:mu_values Hund's-rule spin; claim-test matrix row 65 stale; per-pair-4 NO-TEST row covers a retracted claim.
- Paper 20: L651 DF/THC claim; L1100 "five orders tighter"; abstract GH framing; Childs2021 undefined citation;
  Wecker shots citation; jargon; history prose in captions.
- Paper 23: "(Corrected 2026-...)" annotations at L258 and the D-HFS paragraph; "framework-native Bohr-Fermi";
  "Connes axioms" claim; "structurally sparse" with no nuclear baseline; the "universal" list extends Coulomb
  measurements to HO/WS/NN; L362 asymptotic from finite data; Friar-Sick page number from memory; Sarma-Stevenson
  scope; PPY2025 published ref; jargon; test docstring history.
- Synthesis: retired-history parenthetical "(334 Pauli and 13x QWC)".
- claim_test_matrix Paper 20 section retired-era (L464 334, L467 878, L468 O(Q^2.5), L469 11.10/9.23, L513/521, L803
  34.5 vs 34.0).

## Self-check round 1 SMALLs (2026-11-04; LARGEs fixed in v6.3.0)

- Paper 14:
  - The conclusion's "320x ... structural and persistent" lacks the matched-accuracy caveat, and so do the
    sec:ft_gaussian cc-pVDZ clause and the sec:composed_onenorm 862-vs-369 Trotter comparison.
  - The STO-3G provenance says "published molecular integrals", but they were recomputed at R = 3.015;
    FCI -7.881 vs -7.882 is unregistered.
  - The v2.0.24 note puts 37 Ha next to 838 (it should be 38.6 with PK) and says "the quantum computing approach".
  - Dated retired-value parentheticals remain in the abstract and conclusion (edit bar).
  - The abstract omits the unfavourable He Q = 10 count.
  - The abstract's sec:accuracy cross-ref is wrong.
  - The v2.5.0 note's nested "correct binding energies" and 17.1%.
  - The conclusion lead-in says "three cost metrics".
- Paper 16:
  - L548-553 attributes noble-gas IE to the irrep (D shares it).
  - L205/208-210 and the tab:nu caption carry withdrawn-rule history.
  - L664-667 says "B/D visible only at Level 2" but should cite occupancy.
  - Eq. ionization_surgery does not balance dimensionally (S^{3N-1} vs S^{3N-4} x S^2 x hyperangle).
- Paper 20:
  - The L1153-1162 "Basis completeness" paragraph mixes composed and balanced and repeats 0.20%.
  - The intro says 51-317x "raw" but the cc-pVDZ counts are reduced.
  - The conclusion's "universal 3.17 ... driven by Gaunt" drops the small-basis-fit qualifier (the exponent
    identity is forced by N = c(n_max) Q).
  - Dated history parentheticals remain throughout, plus the section 955-983 correction block (PI call:
    point it to Paper 14).
- Paper 23:
  - The deuteron "1-norm dominated by the diagonal" is 58% Z-only (He-4 72%, "about 70%" holds); both
    shares are unpinned.
  - The [PRIOR ART] tag's sources (Fock, Bander-Itzykson) may not carry the uniqueness clause; consider
    Wybourne/Iachello as Paper 24 does.
  - L789 "reproduces the leading-order recoil at 2.86%" lacks the calibration qualifier.
  - The L1045 residual is attributed to a_0 rounding.
  - History annotations remain at L274, L1117.
- Paper 27:
  - The tab:ep2b caption narrates the 2026-08-22 retraction ("corrected", "Compare the retracted values").
  - The N_max = 4 row is measured but untested (the JSON and tests cover N_max = 2, 3).
  - "Convention-independent commutator (0.75)" is N_max-dependent.
  - Same in Paper 24: L1144 "Magnitude caveat, 2026-08-29".
- group6 synthesis L325-327: the 21x ratio is restated without the owners' "load-bearing evidence is the
  commutator" caveat.
- group4 synthesis:
  - L263-267 calls the He cc-pVDZ/cc-pVTZ rows synthetic (only the dagger rows are).
  - The dated history parentheticals at L300-305, L474-482 and in the abstract remain.
  - "Lower edge of, not below" Q^3.9 (3.79 is below).
- Front door:
  - README L30/192/231 drop "not matched accuracy"; L41 "cheaper at the same qubit count" is unscoped (the
    atomic path is unfavourable at small Q).
  - The positioning page "Figures withdrawn" block and the one-pager's retired-figures parenthetical
    (11.11 vs owner 11.10).
- Field guide: the 37-system library is cited to Paper 14 (Paper 20 owns it); L276/L803 lack the
  matched-accuracy qualifier.

## Completeness-critic gaps (2026-11-04)

- C25 covers no group4 claim: register the matched-qubit caveat (load-bearing), the composed "zero-parameter"
  provenance and the grouping scope (front-door drift on all three went ungated).
- claims_register has no Paper 16 or Paper 23 row.
- Paper 23 D-HFS values cite only debug/ provenance; the backing test `tests/test_paper34_autopsy_baselines.py` is
  uncited.
- INDEX P14 row "at equal qubits" lacks "not matched accuracy".
- `docs/qa/stamps/group4.plan.json` unit keys use the old Paper 23 section title (retitled "A Two-Species").

## PI items

- CLAUDE.md §1.5: the STO-3G market-test sentence and "zero-parameter construction" (drafts given 2026-11-04).
