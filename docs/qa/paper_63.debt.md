# paper_63 — debt (SMALL findings, owed guards, completeness gaps)

Owner of this list: the PM. SMALLs are fixed in one scheduled sweep outside a cert cycle (edit bar,
`docs/qa/criteria.md`). Guards are written as their own activity and fire-tested.
Source run: `/qa paper_63` FULL (first cert, unseeded), 2026-11-06, snapshot 5a1a5b17
(worktree `../geovac-qa-p63`).

## Owed mathematics (PI call)

- **The dictionary as a bound on the rung (self-check round 1, F1).** prop:rung, correctly stated, puts a
  general-orbital integral in the D-module M of its graph's 1s integral and bounds its one-parameter factors by
  those of M's RESTRICTION to the line. The missing step that would restore "rung at most the 1s rung": the
  restriction of the tree family's module to the one-scale line has no composition factor above F_n (for trees
  plausibly because K_0^(n-1) is a cyclic vector of the symmetric-power Bessel module, so every x^k K_0^a K_1^b
  weight lies in the module it generates, and the K_0(tx) transform carries that to F_n's module). Not argued in
  the paper; the claim is re-tiered to module membership plus "every computed case sits on the 1s rung".

## Owed guards (separate activity; fire-test each against the wrong answer it excludes)

- **Headline 6, record-only clauses (code D4; raised to the PI as a C1 gap on a headline).** The correlated
  triangle's order-nine annihilators, composition series and K3 maps at gamma = 2 (exponential), gamma = 1/2 and 3
  (mechanism), and for the cusp-correct geminal at gamma = 3/2 and 2, are on the record only; the guards
  (`test_paper63_{ga1_geminal,ga3_geminal,mg_mechanism}.py`, `test_paper59_c2b_feynman_g2.py`) cover the order-two
  and algebraic pieces. The paper now says so inline. Owed: a K3-map guard per gamma, or a guard on the stored
  annihilator's composition series. The K3 rung at gamma = 1 and 3/2 is guarded (slow).
- **Dictionary three-density check (code D7).** `tests/test_paper63_dictionary.py` `_w` hard-codes the Delta^-3 and
  Delta^-5 weights instead of applying eq:euler_weights to 2 c^-1/2 K_0; the two-density test does apply the
  operators. Owed: compute the Euler image in the three-density test. The "contains only K and E" assertion is
  automatic (derivatives of K are K and E) and carries no evidence.

## SMALL sweep

- Paper 63 bibitem `monkhorst_jeziorski1979`: title is "No linear dependence and multi-center integral problems in
  momentum space quantum chemistry" (J. Chem. Phys. 71, 5268), not "... molecular Sturmian calculations".
- Paper 63 four-density section (~:606): "(both on the prediction record's computation, not a regression test)" is
  wrong for c_{4,0} and c_{6,0}, which `test_ising_moments_are_the_registered_values` guards to 1e-25; only the
  5e-261 precision is record-only (code D6, an upgrade).
- Paper 63 three-electron section (~:369-371): "opens no new transcendence class" — scope to the holomorphic period
  (claims S8).
- Paper 63 four-centre section (~:334-335): "the one this paper was written on" is a split leftover (Paper 59's).
- Paper 63 sec:scope (~:1496): "the P3e-2 basis" is an undefined insider label.
- Paper 63 lem:integral proof (~:220-227): "derivatives of every order in those parameters" dominated termwise fails
  for l >= 3 and for centre derivatives of the Yukawa (distributional terms); the lemma is classical (Weniger); the
  justification needs replacing with a correct one, not a hedge (code D1, SMALL part). Blocked on writing that
  argument carefully.
- Paper 63 conclusion: carries cap three only; caps two (3-D reduction through four densities), four (only the
  triangle among rings) and five are in sec:scope and the abstract (claims S7).
- Done-file headline 1 says "polynomial in the Euler operators"; the paper now says "in c, c^-1 and the Euler
  operators" (code D8). The done-file is frozen for this run; carry the wording at the next co-write.
- Six guards backing Paper 63 keep `test_paper59_*` names and docstrings pointing at "Paper 59 sec:density_ring"
  (`correlated_triangle_g32`, `c2b_elliptic`, `offslice_y_operator`, `gaussian_fit_geminal`, `c2b_feynman`,
  `c2b_feynman_g2`); rename or re-point, and cite them in the prose they back (code D5).
- group3 synthesis (derived): lines ~1537-1539 drop the gamma = 1 mechanism cap and the BBBG Eq. 30 credit (F3);
  the pricing sentences drop the one-machine / own-grid caps (F4); "elliptic squared" and the Calabi-Yau operator
  carry no Bailey/BBBG/AESZ credit (F5); the intro roster of reconvergence papers (~:205-214) omits Paper 63.
- group2 synthesis :827 "The rungs beyond the third center" — the owner indexes rungs by linked densities (F6).

- Completeness critic, unguarded or unrun body numbers: the kite's Lee-Pomeransky critical-point counts (22 generic,
  21 at z_9; ten drop points) and the "41 sets / 29 unknowns / 18 held-out" fit vs the guard's 53 rate sets; the
  count's relation to the eighteen masters needs one sentence (sector convention, record 2026-10-26 l.43); the
  ring's -0.161 / -2.33, 191 / 8.9 and the -9216 t^2 source; the two rejected j-families and the Frobenius-trace
  non-isogeny (~:987-991); the pricing figures (slow operator_evaluator and production_geminal_fit not run this
  cycle); the Gegelia-van Straten conductors (2145/2465, 945/960) were confirmed by the citation reviewer but are
  not in the guard set.
- C5 scorecard reconciliation not audited this run: confirm every registered gate that missed (MP1-2, PF1/3/4, GF1,
  GR1/4 partial; GY3 disclosed) is disclosed where its result is stated.
- Paper 57 J4 "(F, pre-registered and confirmed)" against the owner's "a calibration rather than a hard test": the
  same row says "nearly a corollary of Paper 59"; judged faithful, recorded.

## Citations left UNVERIFIABLE (not fetched; never upgraded on plausibility)

Bonisch et al. 2021 "Eq. (3.7)"; Fresan-Sabbah-Yu parity of the pairing; Avery & Avery 4-centre "App. B" and the
2018 book's specifics; Davydychev-Delbourgo Sec. 5.1 recovering Nickel (nickel1978 bibitem lacks a title);
Lim-Tu-Yu specifics; Kutzelnigg 1994 / Braess 1995 claims (Braess pages 93-103); Vanhove annihilation (the paper's
own check); hehre1969; chetyrkin_tkachov1981; remiddi1991 lacks a title (Sims-Hagstrom 2003 report misprints in
Remiddi's eq. 33; the paper uses the PPR2004 form).
