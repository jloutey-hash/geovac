# Paper 63 (Linked Electron Densities Climb the Banana Ladder) — `/qa` profile

> **RE-FROZEN 2026-11-06 after FULL #1 (see the change log); first frozen 2026-11-06 — PI co-write done.** Branch criterion frozen as written (six clauses);
> all nine C8 headlines frozen; run shape **unseeded** (the standing calibration record is cited;
> the seed plan below stays on file for a later `seeded` run). The goalposts below are fixed for
> the first `/qa paper_63` FULL run, which the PI invokes.

> **Inherits the shared criteria in [`docs/qa/criteria.md`](criteria.md).** This
> file supplies only the Paper-63 scope + deltas + watch-notes.

> **STATUS: RESEARCH-GRADE (FULL #2, 2026-11-06, v6.4.4).** Primary verdict FAIL on two
> SURFACE-LARGE (fixed, C16 `p63-stale-owed-and-rung-rule-home`) and SMALLs (debt file); 0
> OWNER-LARGE; every dimension exercised; gates green incl. C25; coverage 14/14. Not PASS.
> History: FULL #1 FAIL (v6.4.1), RB sprint (v6.4.2), claim-search DELTA DEFECTS (v6.4.3).
> Paper 63 was created 2026-11-06 (v6.4.0). Its moved content was last exercised by `/qa paper_59` on the
> 2026-11-02 approval of Paper 59's items 10-15 (`docs/qa/paper_59.done.md`), under
> Paper 59's goalposts; those items now describe Paper 63's sections and are
> superseded here, not re-certified.

**Scope (single-paper + synthesis footprint):**
- **Paper 63** — `papers/group3_foundations/paper_63_density_graph_ladder.tex`
  ("Linked Electron Densities Climb the Banana Ladder: Period Classes of Molecular
  Integrals over Exponential-Type Orbitals").
- **C9 (synthesis) = GATING from day one.** The group3 synthesis carries
  `\subsection{The density-graph ladder (Paper~63)}` (the rung paragraph moved from the
  group2 synthesis on 2026-11-06, prefaced by the dictionary sentence) plus its bibitem;
  the group2 synthesis keeps a one-sentence pointer in its Paper 59 subsection. Paper 59
  was certified under a stale "no synthesis footprint" premise for a full run; this
  record does not inherit that premise. The scope therefore carries
  `synthesis/group3_foundations_synthesis.tex` (gate `paper_63`, 2 files).
- **Companions, in scope only where Paper 63 restates them (C7):** Paper 59 (the
  three-centre radial core, `sec:genus`/`sec:f12`; the H2 geminal master,
  `tests/test_paper59_h2_geminal_master.py`, is cited by both papers and keeps its
  name), Paper 61 (the modular arithmetic), Paper 18 `sec:level2_genus` (the genus
  grading; one pointer sentence to Paper 63 added 2026-11-06), Paper 57 sub-domain J
  (rows J4, J7, J10, J11 now home in Paper 63).
- **Out of scope:** the rest of group3 (its own `group3.done.md`), Paper 59's
  remaining sections (its own `paper_59.done.md`).

**Deterministic `--gate`:** `paper_63` (`debug/qa/qa_scopes.py`; Paper 63 is also a
member of `group3`, so `--gate group3` reaches it, unlike Papers 58/59/60 in group2).

## Dimensions to exercise (ALL, one invocation — FULL run, first cert of a fresh target)

- **Code / test-backing (C1–C2)** — `code-reviewer` ×1 (Opus-tier: the standing calibration record found
  Sonnet code reviews PARTIAL; 2 seeds if `seeded`). Map +
  **RUN**: `tests/test_paper63_dictionary.py` (the Yukawa-derivative dictionary and the
  Euler-operator weights), `tests/test_paper63_tree_bound.py` (prop:tree_bound; run with
  `--slow`, three of its cases are slow), `tests/test_paper63_four_centre_genus.py`,
  `tests/test_paper63_three_electron_k3.py`, `tests/test_paper63_generic_slice.py`,
  `tests/test_paper63_four_density_rung.py` (tree identity, the order-$n$ operators,
  the Hadamard/banana identity for $n = 2..5$, Fuchs irreducibility, the Ising
  thresholds), `tests/test_paper63_four_density_star.py`,
  `tests/test_paper63_four_density_path.py`, `tests/test_paper63_density_ring.py`,
  `tests/test_paper63_chain_slater_slice.py`, the geminal cluster
  (`test_paper63_{mg_mechanism,ga1_geminal,ga3_geminal,gy_gaussian,one_formula,geminal_shapes}.py`),
  the kite cluster (`test_paper63_{ib_kite,y_minimal_order,kite_ninth_point,master_count}.py`)
  and the pricing pair (`test_paper63_{operator_evaluator,production_geminal_fit}.py`).
  For each: does the test prove the claim (not tautological / weaker than prose /
  right-answer-wrong-reason)? The matrix rows (`docs/claim_test_matrix.md`, the
  "Paper 63 `sec:...`" rows) name what is guarded and what is RECORD ONLY; a reviewer
  who finds prose resting on a record-only item at a tier above [OBSERVATION] reports it.
- **Paper claims / prose (C3, C5, C6, C7, C8 + branch criterion)** — `claims-reviewer`
  ×1 (Opus-tier; 2 seeds), enumeration-forced: every `[TIER]` tag (3 SYMBOLIC PROOF,
  1 SYMBOLIC PROOF + PRIOR ART, 3 SYMBOLIC, 24 SYMBOLIC + MEASURED, 12 MEASURED,
  5 OBSERVATION, 1 OBSERVATION + PRIOR ART, 12 PRIOR ART at the draft date), every
  headline, every `[OPEN]` / cap.
- **External citations (C4)** — `citation-reviewer` ×1 (Sonnet-tier; 2 seeds). Risk
  surface: the banana / Calabi-Yau literature (AESZ tables, Almkvist-Zudilin, Bönisch
  et al. 2021, Vanhove), the Bessel-moment literature (Bailey 1936; BBBG 2008 eqs. 30
  and 33; Broadhurst 2009, arXiv:2212.01962), the uniqueness relation and the one-centre
  dilogarithms (Remiddi 1991; Fromm–Hill 1987; Pachucki–Puchalski–Remiddi 2004), Brown
  2009 (periods of massless graphs), the IBP method (Puchalski–Kędziera–Pachucki), the
  Averys' momentum-space treatment and Monkhorst–Jeziorski 1979 (the placement section).
- **Synthesis faithfulness (C9)** — GATING on the group3 synthesis subsection.
- **Deterministic (C10–C25)** — the step-1 scripts, `--gate paper_63`. C25 owns seven
  registry ids re-homed here on 2026-11-06 (`p59-four-centre-genus-one`,
  `p59-three-electron-k3`, `p59-coulomb-ring-elliptic`, `p59-correlated-triangle-k3`,
  `p59-ising-thresholds`, `p59-tree-banana-hadamard`, `p59-slater-slice-membership`;
  the ids keep their `p59-` prefix as history).
- **Completeness-critic** ×1 (FULL run).

## Branch-defining criterion (FROZEN 2026-11-06, PI-approved; inherits Paper 59's swap): transcendence-tier + open-frontier honesty

Paper 63 is a differential-Galois / period-class paper, so group2's chemistry
benchmarking criterion is vacuous for it and group3's foundations criterion does not
reach its molecular integrals. The criterion below is a sharpening of the shared C3
(prose ≤ tier) + C8 (headline honesty) + C5 (no §3 negative suppressed). Reviewers
verify ALL of:

1. **SYMBOLIC PROOF means a proof.** The dictionary (`sec:dictionary`: every
   exponential-type orbital is a polynomial differential operator in exponent and centre
   applied to the $1s$ Yukawa function; at the fibre every $s$-type density weight is a
   polynomial in $c$, $c^{-1}$ and the Euler operators applied to the $1s$ weight (for
   $\ell>0$ the centre derivatives leave momentum factors; any restatement extending the
   Euler identities to $\ell>0$ is MATERIAL); so every rational-kernel multi-density
   integral is that operator applied to the $1s$ integral of its graph and **lies in its
   $D$-module**) carries `[SYMBOLIC PROOF]`. The **rung bound** is proved only by
   prop:tree_bound, for trees of $s$-type densities on the one-scale slice at the
   holomorphic-period level. Any sentence asserting or relying on a rung bound outside
   that scope (rings, $\ell>0$, other slices, the integrated value), or on attainment
   ("the rung equals the $1s$ rung", "orbital-independent" read as equality), is MATERIAL.
   A `[SYMBOLIC PROOF]` tag resting on a numerical check alone is MATERIAL.
2. **A rung is an operator, not a number.** The tree ladder names the operator of the
   one-scale family (Legendre ⋆ equal-mass banana; $n = 4$ AESZ #58). The generic
   three-electron fibre is a K3 function integrated once on every slice tested, with no
   relation to K3 and elliptic periods found in any basis tried (P3e-2, OPEN-LEANING on
   the walls register). No sentence may collapse "sits on the K3 rung" into "is a K3
   period" or "is a sum of products of elliptic integrals" at generic scales = MATERIAL.
3. **The Hadamard/banana identity is checked for $n = 2..5$ and is PRIOR ART for the
   operators.** The operators are Almkvist–Zudilin / the AESZ tables; the corpus's
   contribution is the route from a molecular integral to them and the tree identity.
   "For all $n$", "we prove the banana ladder", or crediting the operators to GeoVac =
   MATERIAL (registry `p59-tree-banana-hadamard`).
4. **The ring drop is uniqueness, not rings.** The Coulomb triangle sits on the elliptic
   rung because its loop integrates out exactly (star-triangle); the correlated triangle
   with exponential geminals does not drop (K3 rung at $\gamma = 1$, $3/2$, $2$; order nine
   also at $\gamma = 1/2$ and $3$; the order-four drop at $\gamma = 1$ is the equal-mass
   degeneration). "Rings never
   climb" or "every loop lowers the rung" = MATERIAL (registry
   `p59-coulomb-ring-elliptic`, `p59-correlated-triangle-k3`). Only the triangle was
   computed (cap four).
5. **The conjecture stays bounded.** `sec:conjecture` is `[OBSERVATION + PRIOR ART]`: a
   pattern across the tested instances with its four caps stated (fibre vs integrated
   value; the three-dimensional reduction shown only through four densities; the
   dictionary a bound only for $s$-type trees on the one-scale slice, attainment never
   proved; one ring, and the generic fibre's content evidence). The
   nearest theorem (Brown, massless graphs) is cited as such; the Symanzik-polynomial
   proof is explicitly not attempted. Any "theorem", "derived" or "proved" on the
   conjecture = MATERIAL.
6. **Pricing is pricing.** The operator-evaluator and production-fit results are
   `[MEASURED]` costs on one machine and one configuration, on model densities; the
   operators are "structure, not a shortcut" for chemistry. Any "exact ⇒ faster
   chemistry" or "F12 energy error" reading = MATERIAL (exact ≠ accurate, inherited from
   Papers 58 and 59).

## Paper-63-specific watch-notes (the risk surface — ranked)

- **W1 — the dictionary's reach [HIGHEST RISK; criteria 1, 5].** The abstract,
  introduction, conjecture status, scope and conclusion restate the dictionary. Verify
  each keeps module membership as the general statement, scopes any rung bound to
  $s$-type trees on the one-scale slice (prop:tree_bound), and never claims attainment;
  and that cap three survives in every summary surface (Paper 59 pointers, the syntheses,
  the field guide, the Avery reading map, INDEX).
- **W2 — the K3 claims at generic scales [criterion 2].** `sec:three_electron` and the
  group3 synthesis both carry the P3e-2 negative. Verify neither reads the K3 fibre as
  the K3 period of the integrated value; verify the second slice (t;1,2) and the Slater
  slice are stated as "a K3 function integrated once", not as K3 periods.
- **W3 — prior-art credit [criterion 3; C4].** Bailey / BBBG eq. 33 (the closed form),
  Almkvist–Zudilin / AESZ (the operators), Remiddi (one-centre dilogarithms), Broadhurst
  (the kite's ninth point, arXiv:2212.01962 Eq. 25), Puchalski–Kędziera–Pachucki (IBP on
  this graph), the Averys (momentum-space re-expansion). Each must be credited where the
  result is stated, not only in the placement section.
- **W4 — RECORD ONLY vs GUARDED [C1].** Many annihilators, composition series, K3 maps,
  minimality and irreducibility searches are record-only (the matrix rows say so). Verify
  the prose tiers do not exceed what the guards witness; a `[SYMBOLIC + MEASURED]` on a
  record-only item is a coverage gap (C1), raised if load-bearing.
- **W5 — the tiers of the moved content did not move.** Paper 59's `/qa` DELTAs
  (v5.20.19, v5.20.22, v5.21.2) re-tiered several sentences in these sections (the rung
  rule F → Obs; "the operators are known" scoped to n = 2..5; "neither the K3 nor an
  elliptic piece" withdrawn). The C16 entries that gate them list Paper 59 in `files`;
  confirm each still fires on Paper 63's text (add the Paper 63 path to `files` where the
  retired wording could recur) before the run.
- **W6 — genus grading consistency (C7).** Paper 18 `sec:level2_genus` carries the
  genus grading and now points to Paper 63 for the rungs beyond; Paper 57's J4/J7/J10/J11
  rows and `docs/atlas.md` doors 14, 17, 20, 21, 27 name Paper 63 as home. Verify the
  restatements match the owner's tiers (the F on J10 covers the tree identity and the
  non-elliptic operator only; the star/path reductions are outside the F).

## C8 headlines (enumerated, with tiers — FROZEN 2026-11-06, PI-approved)

1. **The dictionary [SYMBOLIC PROOF].** Every exponential-type orbital is a polynomial
   differential operator (coefficients rational in the exponent) in its exponent and
   centre applied to the $1s$ Yukawa function; at the fibre every $s$-type density weight
   is a polynomial in $c$, $c^{-1}$ and the Euler operators $c\,\partial_c$ applied to the
   $1s$ weight; every rational-kernel multi-density integral is that operator applied to
   the $1s$ integral of its graph and lies in its $D$-module. **The tree bound [SYMBOLIC
   PROOF + MEASURED]:** on the one-scale slice of a tree of $n$ linked $s$-type densities
   the holomorphic period is $E\,F_n$, at or below rung $n$ (record
   `docs/predictions/2026-11-06_rung_bound_trees.md`). Guards
   `tests/test_paper63_dictionary.py`, `tests/test_paper63_tree_bound.py`.
2. **The tree identity [SYMBOLIC + MEASURED].** The holomorphic period of any tree of
   $n$ linked two-centre densities is $2^n/(\pi\sqrt{\prod c_v})\int\prod_v K_0(a_v x)\,dx$,
   shape-independent. Guard `tests/test_paper63_four_density_rung.py`.
3. **The Hadamard/banana ladder [SYMBOLIC for n = 2..5; PRIOR ART for the operators].**
   The one-scale operator of order $n$ is Legendre ⋆ equal-mass banana$_{n-1}$ in
   $y = 1/(16t^2)$: $n = 2$ Legendre, $3$ Clausen's elliptic square (K3), $4$ AESZ #58
   (irreducible, no quadratic relation, non-constant Yukawa: the first rung not built from
   elliptic periods), $5$ no new type. The real star and path of four reduce onto the
   $n = 4$ operator with no K3 piece (minimality bounded by mod-$p$ searches). Guards
   `tests/test_paper63_four_density_{rung,star,path}.py`.
4. **A fourth distinct centre adds no genus [MEASURED].** Guard
   `tests/test_paper63_four_centre_genus.py` (1s instance).
5. **The three-electron chain [SYMBOLIC + MEASURED; PRIOR ART for the closed form].**
   K3 fibre; holomorphic period the Bessel moment $\bar V_3$ with the Bailey / BBBG eq. 33
   closed form; on the slices (t;1,1) and (t;1,2) the fibre needs the K3 piece (P G = 0
   and W G = 0 DERIVED on the Slater slice, the Mellin contour shift the one cap); the
   generic fibre's content beyond K3 and elliptic periods is evidence (P3e-2). Guards
   `tests/test_paper63_three_electron_k3.py`, `tests/test_paper63_generic_slice.py`,
   `tests/test_paper63_chain_slater_slice.py`.
6. **The ring [SYMBOLIC + MEASURED].** The Coulomb triangle's loop integrates out
   (star-triangle); at the fibre its operator factors exactly into elementary, algebraic
   and Legendre-equivalent pieces (elliptic rung, below the chain's K3); at one centre
   dilogarithms (Remiddi). The correlated triangle sits on the K3 rung at $\gamma = 1$
   (exact), $3/2$ and $2$, with order nine also at $1/2$ and $3$ (composition series, pieces algebraic / elliptic /
   K3, the elliptic piece a bubble at imaginary momentum $i\gamma$); the drop at
   $\gamma = 1$ is the equal-mass degeneration; the cusp-correct geminal brings no new
   type; the Gaussian geminal leaves the Fuchsian class. Guards
   `tests/test_paper63_density_ring.py` and the geminal cluster.
7. **The kite and the off-slice triangle [SYMBOLIC + MEASURED; PRIOR ART for the ninth
   point's location].** IBP reproduces the series-guessed kite operators at five $\gamma$;
   off the slice the kite has order 11 and the triangle's minimal operator order 13, with
   the K3 piece $R_V$ (irreducible, orthogonal), two crossed Legendre pullbacks and
   seven algebraic pieces by exact maps; the ninth singular point $z_9$ has a closed form
   in the squared rates (Broadhurst's triangle threshold). Guards
   `tests/test_paper63_{ib_kite,y_minimal_order,kite_ninth_point,master_count}.py`.
8. **Pricing [MEASURED].** The kite's operator rebuilds the grid behind Y at about a
   thirteenth of the grid's CPU bound; at chemistry precision double-precision
   quadrature is ~500× faster; Psi4's six-term fit misses Y by 5-13% through its tail;
   the triangle class is 21% of a LiH r12 run (COVERAGE GAP: no guard planned). Guards
   `tests/test_paper63_{operator_evaluator,production_geminal_fit}.py`.
9. **The bounded conjecture [OBSERVATION + PRIOR ART].** The class is set by the shape of
   the density graph and the masses on its lines; trees at rung $n$; rings ordered by
   uniqueness, not count; the integrated value at or below its fibre; four caps stated;
   P3e-2 is its test; the Symanzik-polynomial proof not attempted.

## Caps the reviewers must find intact

- The rungs are read on the fibre; the integrated value may sit lower (and does in one
  class).
- The three-dimensional reduction is shown for the chain of three and the star and path
  of four only; for $n \ge 5$ the holomorphic-period statement stands alone.
- The dictionary places every orbital in the $1s$ module; it bounds the rung only for
  $s$-type trees on the one-scale slice; attainment is observed in every computed case and
  proved in none.
- Among rings only the triangle was computed; the generic three-electron fibre's content
  is evidence from integer-relation searches.
- The finite closed form of the first rung's integrated value is Paper 59's open problem
  and Paper 61's subject, not Paper 63's.

## Seeding plan (worktree only; never touches the real corpus) — ON FILE, NOT USED (first run unseeded, PI 2026-11-06)

K ≈ 6 planted defects, ≥1 catchable by each gating dimension:
- **code ×2** (in a P63 backing test): (a) a Hadamard/banana check reduced to one $n$ or a
  loosened exact-arithmetic assertion in `test_paper63_four_density_rung.py`; (b) the
  dictionary guard made to pass with the Euler-operator weights swapped.
- **claims ×2**: (a) "the rung equals the $1s$ rung" (criterion 1 collapse of the upper
  bound); (b) the generic chain fibre stated as a K3 period (criterion 2).
- **citation ×2**: (a) the banana operators credited to GeoVac rather than
  Almkvist–Zudilin / AESZ; (b) a wrong equation number on BBBG 2008 (eq. 33 → eq. 30) or
  on Broadhurst arXiv:2212.01962.

M ≈ 5 known-good controls that must NOT be flagged: the tree identity's constant
$2^n/(\pi\sqrt{\prod c})$; AESZ #58 for $n = 4$; the Coulomb triangle's elliptic rung;
the correlated triangle's K3 rung at $\gamma = 1$; prop:tree_bound's scope (s-type trees,
one-scale slice).

## Change log
- 2026-11-06 (after FULL #1 and the RB sprint) — **RE-FROZEN, PI co-write ("I agree with your suggestions",
  approving the proposal below).** Criterion 1, criterion 5, headline 1, W1, the caps and the seed controls now
  carry module membership as the general statement and the rung bound only in prop:tree_bound's scope (s-type
  trees, one-scale slice, holomorphic period), proved after the run by record RB
  (`docs/predictions/2026-11-06_rung_bound_trees.md`, GO). `tests/test_paper63_tree_bound.py` joins the code list.
  These goalposts govern the next run (the claim-search DELTA, then the certifying FULL).
- 2026-11-06 (first FULL run, remediation) — **PROPOSED (superseded by the entry above).** The run found that the dictionary never proved "rung at most the 1s rung":
  correctly stated, prop:rung gives only membership of the 1s D-module (code D1, self-check rounds 1-2). The paper
  is re-tiered to module membership plus "every computed case sits on the 1s rung" (observation). Proposed
  criterion 1: the [SYMBOLIC PROOF] is module membership; any sentence asserting or relying on a rung bound
  ("at most", "cannot climb above", "bounds the rung") or on attainment is MATERIAL. Headline 1 accordingly
  (drop "quotient" and "rung at most"); criterion 5 and the caps list: "the dictionary places, it does not
  bound"; W1: verify restatements keep module membership and do not claim a bound; seed control "the
  dictionary's 'at most'" retired. Also proposed for headline 6 and the Hamiltonian observation: "over s-type
  orbitals (for l > 0 not shown)".
- 2026-11-06 (after the freeze) — **Criterion 1 and headline 1 narrowed to match a PM fix in the paper.**
  The paper's abstract and prop:rung had extended the fibre Euler-weight identities to every
  orbital; for l > 0 the centre derivatives leave momentum factors, so the identities hold for
  s-type densities only (the general D-module statement is unchanged). The goalpost is stricter
  than before, not looser: criterion 1 now names the over-extension as MATERIAL. PI told.
- 2026-11-06 — **FROZEN** by PI co-write: criterion as written, nine headlines, unseeded. PM
  corrections before the freeze: version v6.4.0; code reviewer Opus-tier (calibration record);
  the correlated triangle's K3 rung stated at gamma = 1, 3/2, 2 only (order nine at 1/2 and 3),
  matching the paper and registry `p59-correlated-triangle-k3`.
- 2026-11-06 — **DRAFTED** by the PM on the day of the split (Paper 59 → Papers 59 + 63);
  not frozen; PI co-write owed before any `/qa paper 63` run. Pointers re-homed the same
  day: the claim matrix rows, the atlas doors 14/17/20/21/27, seam rows J4/J7/J9/J10/J11,
  Paper 57's J-rows, Paper 18 `sec:level2_genus`, the registry owners (C25) and the
  syntheses (C9 footprint moved to group3). The twenty moved guards were renamed
  `test_paper59_*` → `test_paper63_*`; `tests/test_paper59_h2_geminal_master.py` is shared
  with Paper 59 and keeps its name.
