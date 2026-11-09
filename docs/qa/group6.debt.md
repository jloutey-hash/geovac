# group6 — QA debt (SMALL sweep list and owed items)

Opened by the seeded FULL run of 2026-11-08 (FAIL; record `debug/qa/group6_full_run_2026_11_08_notes.md`).
Edit bar: these SMALLs ride along only when a LARGE fix opens their unit; the rest are fixed in one scheduled
sweep outside the cycle. Items marked OWED need a computation, a guard, or a source.

## Owed guards (write and fire-test separately from the fixes)

- P27 `test_paper27_entropy_locus.py`: the control edge is parity-forbidden and cannot fail; rebuild it on the
  largest in-sector off-vertex edge (expect |dS/S| ~ 0.14 at n_max=3) and fix the runner-up comparison.
- P27: pin the physical-Hamiltonian commutator (10.6% / 9.3%) and block fractions (0.945 / 0.926); pin the physical
  gamma_inf extrapolation and the large-Z return toward 2 (local slope minimum ~1.905).
- P26: hub migration (Sec. IV) has no test; B and C values are member-dependent.
- P34 `test_paper34_l2_presence_bound.py`: rebuild to the Observation (D row recorded as a violation; drop the H 21cm
  one-input row or reclassify it; remove the fitted L2_CLASS_SCALE envelope).
- P34 `test_paper34_autopsy_baselines.py`: read the printed reduced-mass factors; tighten the chain tolerance (a
  1.9 ppm row error passes today).
- P35 `test_paper35_predictions.py`: the aggregate-tally test hard-codes its booleans; the Heisenberg-Euler test asserts
  a constructed pi; the Stefan-Boltzmann guard pins the wrong sign of F.
- P36 `test_paper36_lamb_budget.py`: Table XX closure is the only source-anchored check (add the 0.805 / 0.891 fm rows
  and per-formula checks); `_pct` helper is sign-blind; the 0.12 kHz G_SE,7 2P term is unguarded.

## SMALL sweep list

### Paper 26
- Docstrings in `tests/test_paper26_entanglement.py` and `tests/test_paper26_molecular_entanglement.py` carry retired
  numbers (0.303, ~50x, 4/6/4, 2e-3, Z>=5, S~Z^-2.56 as the ion's).
- eq:icv defines pairwise I_cv while Table II and the test use bipartite 2 S(1s); He 0.751 reproduces in neither.
- 40x core/bond ratio needs pk_in_hamiltonian=False (default gives ~57); the "different electron counts" explanation is
  false; "consistent with S~Z^-2.52" (predicts 16).
- R-independence of S_bond is by construction (tier MEASURED overstates).
- L493 "holds as a theorem" vs the provenance tier; L47-50 / L668-670 "<0.2%" vs L255 7%; L663-665 unqualified
  "separable operator cannot entangle"; L723-724 "diverges at dissociation"; L396-407 cutoff dependence; title / INDEX
  row 133 still headline the decoupling.

### Paper 27
- Title lacks the non-degeneracy qualifier; L1046 / L1105-1107 "bare-graph ground state" entropy; L1101-1103 wrong
  concavity bound; L803-805 vs L990-996 composed vs balanced; L583-585 "irrational" tagging; L36 / L432-433 "rigorously
  recovered" (circular); "saturates, not vanishes" from two points; Be S_kin=1.23 guard pins a LAPACK gauge.

### Paper 34
- Remaining exposure-class vocabulary at L4586-era rows, L5947, L5979, L6477 ("New V.D-class exposure").
- He 2^3P autopsy L5447 "same architecture that closes hydrogen Lamb", L5344-5345 "cross-register multi-focal".
- He-3 structural reading L6094-6100 (stale wall attributions; "50-500x").
- D/T overclaims (L5692 "reproduces at the depth of the budget"; L5884 "full mass x spin matrix").
- D autopsy ppm table mixes a_e conventions (L5565-5586).
- 2S Bethe log listed as a match without "not converged" (L3784); off-precision rows print Drake values as computed.
- Row "Fock + kappa" for hydrogen E_n (L3651); Schwinger row "33% overshoot" stale; HF-4 error code S.
- Mu HFS itemization (+4.703 MHz Karshenboim) does not add and the source is unverified; misplaced r_Z(t) caveat.
- Cs autopsy F_R vs Bohr-Weisskopf conflation; sum 1233.6 not 1233.8.
- mud "-0.50% vs Krauth" route; muHe "Diepold 1666.305" unverified; Li-7 "confirmed Walls 2003" unverified.
- Many sentence-level overstatements logged by the claims reviewers (status words "only/every/all four", "bit-exactly"
  with nonzero residuals, Lorentzian "no structural cost", alkali cliff reason, Ps closure by construction, muHe
  "close the gap", conclusion universality vs Zeeman/Stark).
- Bibitems: drake_swainson1990 title, strohmaier2006 title, Fee 1993 author list, Bonilla "Tab. IV" -> Table 3,
  Nieuviarts years.
- Abstract "depth-3 at -211 ppm" vs Remark "-456 ppm".
- `geovac/magnetization_density.py` docstrings (0.10124, -6657 vs -6665.2, eq 95/96); operator-level Zemach is a
  definitional identity (matrix rows 505, 932 overstate "verifies"); r_Z(p) uncertainty omits the kernel systematic;
  "remaining +7.2" is a plug.

### Paper 35
- L254-256 "iff m^2 in the same ring" is false (m^2 = sqrt 2); L755-760 / L923 / L966 / L1450 Hurwitz-3/2 driver
  claimed for the scalar; L956-966 tally "ten"; L610-617 "Paper 34 is not amended" stale; L1389-1396 refined
  Prediction 1 overreach; L429 / L752 / L1446 "newly computed" priority unverified; module docstring
  thermal_tensor_triple.py:394-404 false sign remark; L1425 / L1432 "bare graph" for the continuum spectrum;
  phantom section refs XIV.A/B; per-Weyl 7 pi^2/720; L1161 "theorem" heading.

### Paper 36
- Conclusion "alone ... +7.1 MHz" vs abstract "+6.9"; title "Convention Erratum" vs double count; the
  "An earlier decomposition ... is withdrawn" history sentence; Barker-Glover not in Table VIII; FNS scales rounded
  145 kHz; "every term the review tabulates" omits ~5 Hz terms; `two_loop_self_energy.py` docstrings; uncited
  CODATA constants.

### Synthesis
- N_Pauli = 27.90 x Q unannotated; L376-379 "one-body operator ... cannot entangle orbitals" unqualified;
  "Six independent instances" of the multi-focal wall vs seven.

### Criteria (docs/qa/group6.done.md)
- The C8 headline lists for P26 and P27 below the FREEZE carry retired values (gamma 1.93-1.96, edge 42.5/40.0,
  rho 0.94, S_HO values, 0.3033/49.6, 4/7/6); criterion 4.4 audits the retracted HO rigidity. Currency fix before the
  next run.

### Self-check round 2 SMALLs (2026-11-08, unswept)
- P26: step title; hub equation and abstract/conclusion still read the hub value as member-independent; L710 hub
  "robust" and the conclusion restatement.
- P27: abstract "large-Z edge" phrasing; L793 "every tabulated slope".
- P34: Lyman-alpha residual quoted as 35 ppm (recomputation gives 32, within rounding of the NIST inputs); abstract
  -211 vs body -456 ppm figure; Uehling rows sit under the projection-depth headings; L1849; "Refined Prediction 1"
  label at L567; residual "cleanly" wording at L4988 / L6201.
- P35: L1356 "Refined Prediction"; L1400 "uniformly confirm".
- Synthesis: L426 "focal-length-tagged"; hub value stated member-independent.

### Third-look SMALLs (2026-11-08, unswept)
- P34 obs:single_origin vs presence bound: the muH Lamb row (textbook assembly) still sits under "L2-count = 2+"
  while the carve-out says textbook-assembly rows lie outside the bound's chains; move it to an outside-the-bound
  note as the H Lamb row already is.
- P34 PK passage: "Paper 19 reads" -> "Paper 19 suggests" (owner hedges twice).
- P34 entropy restatement drops "on the physical Hamiltonian" (synthesis keeps it).
- Synthesis L509-512: split the presence-bound parenthetical into its own sentence.
- Pre-existing: r_Z is both a projection variable (P34 sec:layer2_boundary) and a presence-bound Layer-2 input;
  one sentence should state the bound counts r_Z's value as an input.
- P34 l.2834-2836: S_full(He) 0.0408 (hop model) set beside Paper 27's physical S_B 3.77e-2; confirm the two are the
  same quantity for He.
