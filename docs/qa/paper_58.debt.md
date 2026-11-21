# Paper 58 — SMALL debt from /qa paper_58 FULL #3 (2026-11-20, seeded)

Edit bar: fixed in one scheduled sweep, not inside the cert cycle. SMALLs in units a LARGE fix opened were fixed
with it (CHANGELOG v6.8.1). Guard items are guard-writing work (separate pass, fire-tested).

## SMALL sweep — paper text
- l.~46 abstract "We settle the question of whether that trade was necessary": say the group-theoretic necessity.
- l.~233-237 W1e "replacement rather than a correction": add "on this system" (abstract qualifier).
- l.~183 "GeoVac adds nothing to that continuous apparatus": the paper does add the R-operator orders and the
  master-integral weight; restrict to "the integrals themselves".
- l.~215, 223-225 the l-selection "still holds / exact by construction" for generation 2: scope to the two-body tensor
  (the builder's h block has no l-zeros beyond the m rule); l.~324-325 "its zeros there are threshold decisions".
- l.~228-229 "The sparsity was never a discovered property": say "the l-generated sparsity".
- l.~176-178 "standard Coulomb--Sturmian practice" for two-centre evaluation: general two-centre Slater practice.
- l.~328-331 "no two-center electron-repulsion engine exists ... no ERI was numerically evaluated": stale against
  l.~347-356 and l.~784 (the census predates the engine; say so).
- l.~329, 854 "per-side Gaunt feasibility removes 5.6%": the 176 removals are the L-overlap condition (l.~842).
- l.~385-386 "order = number of distinct decay rates": scope to the exponential-only classes.
- l.~404-405, 484-485 "every finite nonzero singular point is apparent": "reads as apparent" (mod-p reading).
- l.~425-429 the nine-function independence lemma is a hand proof (reviewed): say so inline, as l.~545 does.
- l.~528-529 s-type pair elementary: pre-registered P3b missed on case H9; disclose next to "pre-registered panel".
- l.~598-599 "the constants cancel, because f0(0)=0": justify f0(0)=0 or cite the check.
- l.~672 "none of the 625 is Gaunt-forbidden": counted over the 195 m-passing quartets.
- l.~687 "node-bearing n_max=3 sample ... 195/195": block size differs at n_max=3.
- l.~749 "exactly a complete elliptic integral" vs l.~763 numerical identification.
- l.~751-753 Paper 59 restated without its tier.
- l.~756-759 "no pi, no dilogarithm" as two-centre properties: scope per Neumann term.
- l.~770-776 "(AA|BB) carry no transcendental seed": say no special-function seed.
- l.~801, 818 "missing Neumann engine / awaits the engine": stale against l.~356, 783.
- l.~851 "corroborated below" -> "above".
- l.~1015-1023 Corollary point_group: scope to cross-centre label-generated sparsity (same-centre Gaunt survives,
  l.~912-920); add a tier tag and a tab:backing row.
- Thm (ii) l.~902-904 "every admissible l couples": no l-selection rule; nonzero at generic R.
- l.~1058 "factors of 3.4 and 5.4 for the m-rule": these are dense/genuine; the m-rule factors are 3.2 and 5.0.
- l.~1129 "predicts exponent 1/2" mixes R*(n) and n*(R).
- l.~1204-1208 DLMF top point c^2 = 300 at the double-precision floor (true local slope -1.87 vs -1.89; prefactor
  c^2.12 vs 2.07): drop the point or recompute at 40 digits; the verdict holds (asymptotic constant ~31).
- l.~1267-1268 "Screening moves the front; nothing else does": tested on H2/HeH+ only; scope it.
- l.~1325-1366 Lowdin and cost-conservation inference untagged; "currency / invariant" law language in an Observation.
- l.~1446 "correlation cannot have been what was missing": add "on this system".
- l.~1454-1459 NaH exponents "optimized on the isolated Na atom": E(Na) not stationary at the production shapes
  (dE ~ 0.5-2.7 mHa for 0.02 shifts); say "set on the isolated atom".
- l.~1508, 1526 LiH 53% and 15.0% have no test here; l.~1541 -8.887049 is an electronic energy (label it).
- l.~1537-1538 the twelve-electron leg is pinned by test (the paper says driver); tab:backing row Obs truncation.
- l.~1587 "variational zeta = 1.197" is the HF optimum; the FCI optimum is 1.2005.
- l.~1658 hybrid domain: "two independent quadratures <= 7e-12 on twelve cases" backed on four cases, one quadrature.
- Table tab:backing: rows missing for Thm (iii), the Corollary, the H2 closed-form PES, the DLMF dictionary and rate
  bound, the many-electron and coherence fronts, the Ei-domain engine finding; caption "see Sec. scope" dangles.
- NaH table "-38.0 %" for an unbound row without its R.
- Bibliography: lesiuk2014a cited as source of the (AA|BB) closed forms (it is a pointer); NOCI "fragment-native" gloss.

## SMALL sweep — derived documents
- group2 synthesis: "conservation law" stated flat (l.~762, 935, 1208); [OBSERVATION] dropped (l.~757);
  13.8x density tier is COUNTED (l.~737); abelian-residue scope statements dropped (l.~87-91, 749); NaH R_eq
  qualifiers (l.~769); currency gaps: QFD, exact-order results, one-centre root proof, Obs truncation; open questions
  omit the a > b weight, the form of f0, fixed-R values.
- field guide l.~469, 737; synthesis l.~781: weight-two restatements drop "seven integrals" / "0<a<b".
- docs/chemistry_period_inventory.md: row 3a Gaps cell stale (understates); row 3b omissions; row 8 EQ entry.
- docs/claims_register.md row 22 tier understates (abelian residue is INTERNAL THEOREM); row 5 qualifiers.
- docs/certified_reference_values.md l.~125-127 weight-two clause unscoped.
- docs/claim_test_matrix.md: rows 912/913 still say exactly 9 is not proved; row 607 BACKED-SOUND and OWED at once;
  rows 613/629 point to scripts now in debug/archive/misc/; row 614 "no pi / no dilogarithm two-centre-only" stale;
  row 600 "Exact decidability needs the general-m Neumann build" stale.
- Paper 32 l.~5354 drops "on this system".
- geovac/two_center_eri.py l.~1254-1256 "sigma != 0 ... ARGUED, not verified" vs l.~1279 and the sigma-pole test.
- Test docstrings: test_paper58_pachucki_exchange.py l.27-28, test_paper58_pachucki_extension.py l.34-35, 1005 say
  exactly 9 is not proved; test_paper58_noci_engine.py l.20 stale; test_paper58_nah_ladder.py l.25 "against experiment";
  test_paper58_census.py l.104-107 mis-cites 195/195 for the S/h decision.

## From the self-check (round 2), not fixed
- Weight-two restatements dropping "for 0<a<b" / the seven-integral scope: Paper 34 l.~1319, Paper 18 l.~363, Paper 57
  l.~158, field guide l.~469 and l.~737, group2 synthesis l.~782. One-centre restatements dropping "at rational exponents".
- Paper 57 l.~161 "decided pi-free by Lindemann-Weierstrass" (LW decides vanishing).
- Paper 59 l.~99-101 "the Neumann expansion terminates" (not at unequal exponents); l.~109 per-term dilogarithm claim.
- Paper 32 l.~5354 "(the neglected inter-center overlap)"; Paper 58 l.~1684 "whose neglected overlap is the subject".
- Paper 58 l.~459-461 the exchange-point value written before "at the exchange point"; l.~144-149 Paper 14 exponents
  without "matched qubits, not matched accuracy".

## Guard pass (separate, fire-tested)
- Many-electron front: the four front residuals, the zeta_eff trajectory, basis doubling and the coherence shifts are
  asserted nowhere; the slow test runs quick mode at tolerance 0.01 against the paper's 0.006, and depends on debug/
  (debug/aha_t1_core.py restored 2026-11-20; the 2026-11-08 sweep had archived it and the test errored).
- m != 0 exchange: 47 digits at the unused point (guard 1e-29 / 1e-19); a > b 34 digits (guard ~32.3).
- R^9 ln R coefficient: the power of a is untested (a = 1 at both pairs); the reviewer derived the onset coefficients
  for all seven integrals symbolically (upgrade candidate).
- Operator precisions above their guards: 10^-81, 6e-79, 3e-75; "121 orders" vs N_SLOW = 64; "degree 60" vs 40;
  "recovered on a further prime" unguarded.
- exact_order C1: R values not asserted distinct; strengthen to gcd degree 0 after stripping t, t-1, t+1.
- master_weight D2 tolerance 1e-3 (160x slack); l.106 tautological half; self-pinned literals l.288/295.
- 3ddel Q(t) frozen and checked at five ratios only; the 576-index runs pass outside pytest only.
- hybrid_every_ratio main file l.~186 checks Z positivity at u = 1 only (panel file is right); V tie point-sampled;
  24 of 30 panel region forms not tied to the evaluator (verified by the reviewer to 90-100 digits).
- bra/ket certificate: I_A never computed on the census configuration; blind to a uniform scale.
- (AA|BB) 5x7 grid and n_max = 3 sample: driver only; four-entry test.
- QFD: tau termination asserted numerically (fires at 1e-45, not 1e-70); certified R = 1.6, 2.0, zeta = 1.197 rows
  absent from benchmarks; "dps 60 vs 90 agree to 6e-85" untested.
- 38/38 panel order = rates checked on ~10 cases; swap Pauli ratios untested; LiH truncation values unpinned.
- NaH zeta provenance test accepts re-stored molecule-tuned values.
