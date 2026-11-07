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
