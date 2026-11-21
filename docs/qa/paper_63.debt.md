# paper_63 — debt (SMALL findings, owed guards, completeness gaps)

Owner of this list: the PM. SMALLs are fixed in one scheduled sweep outside a cert cycle (edit bar,
`docs/qa/criteria.md`). Guards are written as their own activity and fire-tested.
Source run: `/qa paper_63` FULL (first cert, unseeded), 2026-11-06, snapshot 5a1a5b17
(worktree `../geovac-qa-p63`).

## Owed mathematics (PI call)

- **The rung bound beyond s-type trees.** RESOLVED for trees of s-type densities on the one-scale slice at the
  holomorphic-period level: prop:tree_bound (record RB 2026-11-06, guard `tests/test_paper63_tree_bound.py`), with
  attainment for n <= 8 when the period is nonzero (PI-approved 2026-11-06 after FULL #2; Basset guard). Still open:
  rings, l > 0 (momentum factors carry directions), other slices (cyclicity of K_0(b_2 x)...K_0(b_n x) in the tensor
  product of Bessel modules), the fibre, the real 3-D integral.

## Owed guards (the planned guard pass; fire-test each against the wrong answer it excludes)

- **Headline 6, record-only clauses. PI disposition 2026-11-06: accepted as a disclosed coverage gap; the guard is
  written in the guard pass.** The order-nine annihilators, composition series and K3 maps at gamma = 2
  (exponential), 1/2 and 3 (mechanism), and for the cusp-correct geminal at 3/2 and 2 are on the record only (said
  inline); the K3 rung at gamma = 1 and 3/2 is guarded (slow). Owed: a K3-map guard per gamma.
- **Dictionary three-density check (code D7).** `_w` hard-codes the Delta^-3 / Delta^-5 weights; compute the Euler
  image instead. The "contains only K and E" assertion is automatic.
- **Fuchs bound equality** (`test_fn_operators_irreducible_by_fuchs`): a `<=` -> `<` plant does not fire; add a control
  at margin 0.
- **F_5 operator on quadrature** (`_annihilation_residual(5, ...)`), re-verified by FULL #2 to 1e-15 / 4e-20.
- **Six `tests/test_paper59_*` guards backing Paper 63** (`correlated_triangle_g32`, `c2b_elliptic`,
  `offslice_y_operator`, `gaussian_fit_geminal`, `c2b_feynman`, `c2b_feynman_g2`): rename or re-point their
  docstrings ("Paper 59 sec:density_ring") and cite them in the prose they back (code D5).
- The tree identity with non-1s weights is guarded only for the chain of three; prop:tree_bound uses it at every n.
- Pricing: the evaluator's "a thirteenth" and the ~500x speed test, the fit law, the plain/weighted fit numbers have no
  guard (tagged [MEASURED], one machine).

## SMALL sweep

**Swept 2026-11-06 (v6.4.5):** Monkhorst-Jeziorski title; Remiddi title; Ising L-values now cited as guarded;
"opens no new transcendence class" scoped to the holomorphic period; "the one this paper was written on";
[MEASURED / classical] -> [MEASURED + PRIOR ART]; quadric counts "one in the guard"; bubble check "two in the
guard"; the g32 sentence moved back to the gamma = 3/2 paragraph; pricing fit errors "relative"; GF1's miss
disclosed; conjecture (i) appositive scoped to n <= 5; status "onto the Bessel-moment family"; placement "1s" and
[OBSERVATION] on the causal reading; "the P3e-2 basis"; prop:rung deg P <= m; abstract correlated-triangle gamma
scope, [MEASURED] on the attainment observation, "all s-type", the dictionary cap; intro "all s-type"; group2
synthesis "climb further" scoped and "rungs of more linked densities"; field guide reading path -> Paper 63; group3
synthesis "1s densities" and the gamma = 1 mechanism cap; Paper 18 "rungs of more linked densities"; Paper 59 :869
re-pointed to its sec:genus; Paper 57 J8 and J10 "1s" and the path reduction; seam md J8, J10 and the header;
matrix Hamiltonian row "over s-type orbitals".

**Still open:**
- lem:integral proof (~:220-227): termwise dominated convergence fails for l >= 3 and centre derivatives of the
  Yukawa; the lemma is classical (Weniger); replace the justification with a correct one (blocked on writing it).
- Conclusion carries caps one and three; caps two, four and five are in sec:scope and the abstract.
- :1460 calls eq:tree proved while its tag is [SYMBOLIC + MEASURED] (tag upgrade candidate); :1517 "The operators
  of the ladder are known" unscoped (true for F_n's operators).
- Abstract "on every slice tested the three-electron fibre is a K3 function integrated once" drops the proved /
  mod-p split.
- mg_mechanism guards more than its clause says (kite order, pair-mass singular points).
- The kite's critical-point counts (22 / 21 / ten drop points) and their relation to the eighteen masters (one
  sentence owed); the 41/29/18 fit vs the guard's 53 rate sets; the ring's -0.161 / -2.33, 191 / 8.9, -9216 t^2;
  the rejected j-families and the Frobenius non-isogeny.
- C5 scorecard: GY3 and GF1 now disclosed; MP1-2, PF1/3/4, GR1/4 not yet reconciled.
- group3 synthesis: pricing sentences drop the one-machine / own-grid caps; "elliptic squared" and the CY operator
  carry no Bailey/BBBG/AESZ credit; the intro roster of reconvergence papers omits Paper 63.
- `docs/claims_register.md` has no Paper 63 line.

## Citations left UNVERIFIABLE (not fetched; never upgraded on plausibility)

Bonisch et al. 2021 "Eq. (3.7)"; Fresan-Sabbah-Yu parity of the pairing; Avery & Avery 4-centre "App. B" and the
2018 book's specifics; Davydychev-Delbourgo Sec. 5.1 recovering Nickel (nickel1978 bibitem lacks a title);
Lim-Tu-Yu specifics; Kutzelnigg 1994 / Braess 1995 claims (Braess pages 93-103); Vanhove annihilation (the paper's
own check); hehre1969; chetyrkin_tkachov1981; remiddi1991 lacks a title (Sims-Hagstrom 2003 report misprints in
Remiddi's eq. 33; the paper uses the PPR2004 form).

## From /qa paper_63 FULL 2026-11-21 (unseeded) — SMALL debt
Paper text (units not opened by a LARGE fix):
- Abstract l.~48, 73-75 and Scope l.~2192: the order bound restated without the slice (t;1,...,1).
- Abstract l.~80-83 tags the genus-zero reading of the Averys' route [PRIOR ART]; the owner (l.~2172) says [OBSERVATION].
- Abstract l.~49-50 and Conclusion l.~2229 "any tree of n linked two-center densities has the same period": the owner says 1s densities.
- Abstract l.~63-64 "operators follow from integration by parts on the kite and the tetrahedron": Y at gamma = 1, 2 fitted from series.
- l.~1300 "sits on a rung" defined as an upper bound only, against prop:tree_bound's "at or below / on".
- l.~1902 "eighteen masters" vs l.~1909 "22 at generic z, 21 there" — reconcile or delete the "equals the number of masters" clause.
- l.~2115-2118 the integrated three-centre value cited as consistent with part (iv): Paper 59 shows only that it leaves genus zero.
- l.~1346 "one step above the weight-one two-center seed set" — per Neumann term (C16 two-centre-weight-one-whole-integral; add Paper 63 to cited_by).
- l.~1759 "equals the guessed one exactly" at four gamma: guard at gamma = 3/2 modulo a prime; the rest on the record.
- l.~720-724 path-of-four integrand "checked to 20 digits": no guard (reviewer re-derived it, 1.6e-9); l.~583-585 second-slice fibre agreement unasserted; l.~643-644 F_n quadrature 1e-28 "each": guard n = 4 (+ 6-8 slow).
- [MEASURED] side clauses quoting more cases / tighter precision than their guards: prop:tree_bound (4e-38, nine products), prop:tree_rung (m <= 4, 1e-27), prop:odd_banana (lambda = 3 unguarded, 2e-31), prop:star_attain (2e-16 at 40 digits; C to 1e-11 for nine trees), prop:tree_attain_hom (eleven mod-p trees; 1e-11 for nine), prop:tree_rung_l (6e-17 traces to debug/).
- l.~1068-1074 history prose ("An earlier version of this proof"); l.~847 [MEASURED] on exact identities; l.~1061 exponent tables for n = 6-8 not printed; l.~692 "(for l>0 this is not shown)" partly stale (n = 2 bounded); l.~2209 "open even for s-type trees" — name the Slater weights.
- B3: prop:tree_attain_hom step (ii)'s local rule is now written into the proof (2026-11-21); prop:star_attain (ii) still cites "the local rule stated after the conjecture".
- l.~640 and prop:tree_rung's Sym^m Bessel module construction: no citation at the locus (BBBG / Fresan-Sabbah-Yu).
- Bibliography: monkhorst_jeziorski1979 title is "No linear dependence or many-center integral problems ..."; l.~2180 claim unverified; weniger2005 venue; nickel1978 title.
Derived documents and infrastructure:
- viz/public/papers/paper_63_density_graph_ladder.html (generated from the Zenodo manifest) is stale against prop:tree_reduction and the attainment results: rebuild with debug/build_paper_pages.py at the next Release (the manifest is PI/release territory).
- .zenodo.json description ("a tree ... is annihilated") — refresh at the next Release.
- papers/INDEX.md row 63 omits the proved fibre propositions; group3 synthesis l.~1509 "every computed case" drops "s-type"; Avery map l.~346 same; period_ladder l.~296/390 "sits on" definition; inventory row 40 omits tree_attain_hom; docs/topic_to_paper_lookup.md has no Paper 63 row; CLAUDE.md section 6 tiers omit Paper 63 (PI-editable section? no: section 6 is PM-editable).
- docs/claim_test_matrix.md rows 175, 906 partly stale; row 906 cites a debug/ memo; mg/ib rows 195/214 contradict each other.
- check_retracted_terms.py notes at l.~9387 and l.~9479 describe pre-proof states (history; patterns unaffected). p63-tree-bound-level may trip on a correct future restatement of prop:tree_rung (mirror direction).
- docs/qa/paper_63.done.md caps 2 and the fibre half of 3, criterion 5 and W1 contradict the current owner: re-freeze with the PI at the next run.
- Unsent Avery draft (docs/outreach/avery_message_draft_2026_11_06.md, uncommitted): stale on attainment and "are the banana integrals" — fix before sending.
Guard pass (separate, fire-tested):
- The Pi != 0 extension (wrong answer to reject: a root nabla^2 evaluated as 1, giving pi instead of 5 pi/2); mixed-degree step (i) at n >= 5.
- Record-only mod-p legs (fork5_root0, path5_centre, seven second-prime results); five-star amplitude guard at 15% (1.12 C5 not caught).
- mg kite guard compares stored data only (a perturbed operator passes); ga1 distinctness accepts a rescaled operator; geminal_shapes a < n vs "most primes"; pricing figures (timings, 2e-69 two routes, M13 five t, the plain and weighted fit figures).
- test_paper63_four_centre_matched.py and test_paper63_fxk_per_key.py import debug/ modules and skip silently if debug/data files vanish.
- Tautological asserts: tree_rung:122, odd_banana:409-410, tree_attainment:406/411; literal pin odd_banana:401.
