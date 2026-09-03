# Trunk — carryforward: remediation scope after the 2026-09-01 FULL FAIL

> ## STATUS 2026-09-02 (later, v5.3.0): Part F PI items DONE — F1.1, F1.2,
> F2.1 (+ propagation), F2.2, F6 all four. Remaining Part F rows
> (F1.3–F1.15, F2.3–F2.7, F3, F4, F5) still owed, then the DELTA under the
> F-seeding rule. Detail: CHANGELOG v5.3.0; findings the remediation
> surfaced are in **Part F7** below.
>
> ## STATUS 2026-09-02: FULL run #2 = FAIL (trustworthy on 9/11 chunks).
> Remediation scope = **Part F** below. The DELTA run after Part F must
> carry the two INCONCLUSIVE code classes (F-seeding rule).
>
> The PI directed a second FULL run instead of the DELTA this file asked for
> ("Do not re-run the FULL certifying pass until a DELTA comes back clean",
> Sizing). Recorded, not relitigated. Consequence: the run paid full price
> and confirmed the standing rule — at least four of the verified MATERIALs
> (P32 :6477, :3758-3762, :5113-5201, :4868/:5216) are tags the 2026-09-01
> C3 pass placed at the wrong tier (the C3-boundary shape), i.e. defects the
> previous remediation introduced; a DELTA would have found them at ~1/5 the
> price.
>
> ## STATUS 2026-09-01: Parts A—D COMPLETE (D executed, not just
> adjudicated). Awaiting the DELTA run.
>
> - **Part A** done: C5 screen fires on the prediction class (anchor on the
>   8.8e-8 residual + definite-label negation-binding rule + `predictive`),
>   two-way discrimination proven; C21 family + C16 entry for the retired
>   scaling, both proven; trunk banner corrected.
> - **Part B** done: B1-B10 all applied (B4 revised -- the comparator was
>   seed S-D1; the live test was already `==`; parametrization widened to
>   n_max<=5). Plus the claims-dimension M3-M10 and A4-adjacent items:
>   P0/P7 tier scoping, the L3 footnote inversion + false Sobolev claim
>   (A7), the b-constant in Papers 32, 38 AND 18 (a third locus the P38
>   remediation discovered).
> - **Part C** done -- and WIDER than scoped: the two zombie agents found
>   and fixed 6 additional gate-invisible/line-broken loci in Papers 54,
>   57, 23, 34 and the group6 synthesis, plus 2 more in Papers 22/24 the
>   C16 regex could not see. One honest claim REVERSAL disclosed en route:
>   LiH's balanced QWC count (629 exact-rule vs 87 retired) no longer sits
>   below the Gaussian baseline (273); Paper 19 now says so.
> - **Four consolidation corrections**: seeds S-D1, S-E2, S-F1, S-F2 had
>   leaked into the genuine-findings union (protocol step-6 lapse; recorded
>   in the run notes SS1b).
> - Verification: all 11 trunk deterministic gates PASS; C16+C21 PASS on
>   group2/3/4/6/synthesis; 61 tests green across the touched suites;
>   coverage matrix PASS.
> - **Part D** adjudicated (PI, 2026-09-01), all four: **D1** — the trunk
>   owes the C3 inline-tier pass (four `claims-reviewer` agents dispatched;
>   the tier is derived from the claim matrix + backing test, so the pass is
>   itself an audit); **D2** — yes, the group3 synthesis gets claim-matrix
>   rows carrying their *source paper's* tier (other six syntheses belong to
>   their own re-certs, explicitly NOT closed here); **D3** — Opus citation
>   tiering above ~50 bibitems, written into `qa.md`; **D4** — the ~5,000-line
>   coverage ceiling is ACCEPTED, with the consequence stated in-place (a
>   PASS certifies the surface actually read; the deterministic registries do
>   scan the unread lines, so what is unguarded there is judgment, not the
>   known retracted/retired classes).
> - **Part D EXECUTED** (2026-09-01, same day as adjudication): the C3
>   inline-tier pass shipped **257 tags across all six trunk documents**
>   (0 before) via four Opus reviewers, and the matrix gained its **first
>   synthesis block** (23 rows for the group3 synthesis). Because the tier
>   had to be derived from the matrix + the backing test rather than from
>   the sentence's tone, the pass ran as an audit and surfaced **~25 further
>   defects** — chiefly one repeated class: papers correct in the middle and
>   overclaiming at both ends (P7's abstract + closings, P0's abstract, the
>   synthesis's abstract + conclusion, P32's abstract, P38's theorem
>   displays). Two upgrades applied, so the pass cut both ways.
> - Two findings inverted on PM inspection and are recorded as such: the P38
>   envelope "fix" would have installed a retired constant (and exposed a
>   C16 scope gap, since the entry never covered P38 — now closed with
>   two-way discrimination proven); and a tier disagreement between P1 and P7
>   resolved **against the PM's own tag**, not the agent's.
> - The group3 synthesis's §V dynamics backing was found to have **never run
>   under pytest** (wrong filename for collection AND script-style fixture
>   wiring that ERRORed when invoked directly). Numbers verified to reproduce
>   first, then rewired: 3 passed in 20.3 s, now collected by default.
> - Verification after Part D: **all 12 trunk deterministic gates PASS**
>   (scopes stated), gate-coverage matrix PASS (11 gates × 11 targets, no
>   empty cell), all six documents compile clean.
> - **GATE DEFECT FOUND AND FIXED (C11).** Under `--gate`, C11 keyed findings
>   on a papers-relative path while the scope predicate matches repo-relative
>   ones, so **every** finding was filed as out-of-scope advisory and the gate
>   printed PASS. Proven by planting a completely wrong internal title on a
>   trunk bibitem: detected, printed, exit 0. **C11 has never been able to
>   fail in a gated run.** Third instance of the GATE SELF-AUDIT pattern after
>   C10-by-exit-code and C5/C12's hard-coded document set. Fixed; both
>   criteria re-proven in both directions; siblings probed (they fire);
>   `tests/test_internal_titles_check.py` pins it. Found only because a NEW
>   criterion had to prove it fires — the discipline caught the criterion
>   next door.
> - **Part E item 6 CLOSED** (see the item itself): 30 prose references with
>   no bibitem, +30 bibitems / +26 cites, plus a year audit that found 6
>   mismatches (3 pre-existing in Paper 32) and 30 more corpus-wide, now
>   surfaced as per-group advisory debt by the new C11 criterion.
> - **NOT started**: Part E depth items beyond those absorbed above.

> **Origin: `/qa trunk` FULL certifying run, 2026-09-01 = FAIL (trustworthy).**
> Panel fully calibrated (13/13 seeds after an Opus re-dispatch of citations,
> 0/5 false positives); deterministic layer 13/13 with scopes stated. Full
> record + evidence: `debug/qa/trunk_full_run_2026_09_01_notes.md`.
>
> Following the Paper-59 precedent, remediation is scoped as its **own sprint**
> rather than fixed inside the `/qa` run. This file is that scope.
>
> **Cert path:** Part A → Parts B/C → **DELTA-verification `/qa trunk`** →
> once clean-delta, the **FULL certifying run** → PASS. Part D is PI
> adjudication and gates nothing mechanical. Part E is optional depth.

---

## Part A — gate fixes FIRST (they are what verify Parts B and C)

Ordering is not cosmetic. Every content fix below is supposed to be *checkable*,
and three of the gates that should check them are currently blind. Fix the
instruments before the readings.

**A1 — C5's screen cannot see `prediction`.**
`debug/qa/check_k_label.py` matches `derived / theorem / proven / conjecture /
conjectural`. §13.5 names four forbidden tiers and one of them —
**`prediction`** — is absent. That is why four loci in Paper 32 (including its
abstract) passed C5 for an unknown period.
→ Add `prediction` / `predictive`. Re-read §13.5 and confirm no *fifth* word is
missing. **Prove discrimination both ways** per the hard rule: fire on the
retired wording, stay silent on "observed the numerical coincidence"
(Paper 32 `:1646`, which must remain clean).

**A2 — C21 cannot see the retired resource figures.**
`2.5`, `51`, `1712` were never registered, which is the whole reason a ~20-locus
zombie survived a certification. Registering bare `2.5` corpus-wide is
unworkable (it is a common numeral), so this needs a *contextual* family:
require `Q^{2.5}` / `O(Q^{2.5})` / `Pauli` adjacency, and forbid it where the
retirement is disclosed.
→ Register the family with canonical `27.90 × Q` exactly-linear + the
`54×–317×` range. Prove two-way discrimination against the live retired text
recovered from git and against the corrected wording.

**A3 — C16 has no entry for the retired-scaling phrasing.**
The zombie is a *claim shape* ("O(Q^2.5) Pauli-term scaling", "51× to 1712×
advantage"), which is C16's job, not only C21's.
→ Add the phrasings; two-way discrimination proven; scope covering group2,
group3, group4 and the syntheses.

**A4 — trunk's staleness banner under-scopes.**
`docs/qa/trunk.done.md` names 4 changed `.tex`; the live checker says 5. A
reviewer scoping from the banner skips the synthesis — the document that held
the zombie. Regenerate the banner.

---

## Part B — trunk content (the FAIL items)

Ordered by severity. Every one verified against primary text during the run.

| # | criterion | locus | fix |
|:--|:--|:--|:--|
| B1 | **C5** | P32 `:71` (abstract), `:6099`, `:7028` (in a numbered Observation), `:7035` | "the α **prediction**" → "the α **observation**" / "empirical match". `:1646` is already correct and is the model wording. |
| B2 | **C7** | P38 `:549-562` vs `:1142-1147` | Λ is defined twice, incompatibly, and `thm:main_intro` / `thm:main` display the headline **in Λ**. Fix §Setup to define Λ := van Suijlekom state-space GH, naming Latrémolière only as the open strengthening. Then sweep the **11 residual "propinquity" assertions** about *this paper's own result* (`:292, 751, 761, 1127, 1174, 1202, 1253, 1288, 1293, 1398`, + keywords `:134`). |
| B3 | **C8** | P32 `eq:forced_count_chain` `:5108`, summary `:5371` | Chain prints the matter-sector 128 where the trunk delta requires the **full-axiom** count. The paper's own proof (`:5142`) says so and directs a relabelling never applied. Either print 260 or apply the relabelling the proof specifies — and fix `:5371`'s unqualified "128 per generation". |
| B4 | **C2** | `tests/test_p38_action_seminorm.py` | **CORRECTED 2026-09-01: the `<=` comparator was seed S-D1** — the live test already reads `==`. Genuine residue: the parametrization covers `{2,3}` while the paper's footnote claims verification to `n_max ≤ 5`; widen it (runtime permitting). |
| B5 | **C2** | `tests/test_trunk_qa_c2_delta.py` + P7 §Discussion | "two independent routes agree" is one route: `c2_formula()` is hardcoded in the test, has **no production implementation**, and no derivation in the paper. Options: (a) derive `c²(n,l)` in the paper and implement it in `geovac/`, or (b) downgrade the claim to "consistent with" and drop "two independent routes". **(b) is the honest minimum.** Feeds Paper 2's Δ, so do not leave it as-is. |
| B6 | **C1** | P1 abstract claim (iii), `:191-255` | Berry phase ≡ 0 and Θ(n) = −2ln((n+1)/n) have **zero live backing**; the only artifact is archived and tied to the retracted k=2.113. Either write a live test (the maths is elementary log-difference algebra — cheap) or downgrade the abstract claim. §13.4a says a paper equation needs a test. |
| B7 | **C4** | P32, P38, P1 | See the citation block below. |
| B8 | numeric | P7 `:903` | `E = -Z²/2 + 5Z/8 = -11/8` — the printed formula evaluates to −3/4 and the registered value is **−11/4**. Check against the C21 registry before editing; may be a twin. |
| B9 | numeric | P32 `:2148` | Retired float-contaminated `b` at "≥80 dps"; P38 `:1431` records **22 verified digits** and the strings diverge at digit 12. Also fix P38's own three-way inconsistency (claims 22 digits, prints 59, appendix says 4.105 vs 4.109). |
| B10 | hygiene | P1 `:503, 510, 517`; `:426` | Three `[Placeholder]` figure boxes ship live — in a paper whose Erratum retracts a claim *because* its figure "remained a `[Placeholder]`". The retained log-holonomy claim has one too: same evidentiary state, opposite disposition. Also `\appendix` sits after `\end{thebibliography}`, and the captions restate withdrawn numbers without the caveat. |

### B7 — citations (all verified by the Opus pass)

- **`hekkelman2022`** (P32 `:7214`, cited **in the abstract**) — "PhD thesis,
  Radboud University (2022)" does not exist. Correct item: Master's thesis 2021,
  published *Lett. Math. Phys.* **112**:20 (2022), arXiv:2111.13865 — which
  **Paper 38 already has right under the same key**.
  **This key has been wrong three times**, with a *"citation re-delta CLEAN,
  zero MATERIAL"* recorded between attempts. When fixing, re-test the *specific*
  defect rather than trusting the record.
- ~~arXiv 0409307 at `:2430`~~, ~~Camporesi label (1994)~~, and ~~the
  Leimbach propinquity attribution (:164-166)~~ — **all three were
  calibration seeds** (S-F1, S-F2, S-E2); live corpus already correct.
  Removed 2026-09-01.
- **P32 `:6041`** — "Marcolli–van Suijlekom 2014 rationality theorem for
  Robertson–Walker" cites `marcolli_vs2014` (*Gauge networks*, JGP 75).
  The correct work, `fathizadeh_marcolli2016` (CMP 356, 2017), is **already in
  the same bibliography**. Confirm which is meant (Fathizadeh–Marcolli 2017 vs
  Fathizadeh–Ghorbanpour–Khalkhali JHEP 2014) before swapping.
- **`chamseddine_connes2010`** third component — PRD 83, 045001 (2011)
  unresolvable; the real paper of that shape is *"…and the **Superstring**"*,
  Phys. Lett. B **396** (1997) 71. Components 1–2 are correct; drop or replace
  the third.
- **`deligne2010`** — preprint pointer `math/0302267` is Deligne–**Goncharov**,
  a different paper. Journal ref is correct; delete the preprint line.
- **P38 `:558`** — inline "TAMS 370 (2018)" contradicts its own bibitem
  "368 (2016)". The bibitem is right. (Sits in the *"Propinquity convention"*
  paragraph — the same paragraph as B2.)
- **`connes1995`** (P38) — book is Academic Press **1994**; and the claim it
  backs (real structure, KO-dimension) is the content of *"NCG and reality"*,
  JMP 36 (1995) — which is what the same key means in Paper 32. Decide which.
- **`farsi_latremoliere`** (P32 `:6068`) — third author **Packer** omitted.
- **`leimbach_vs2024`** (P38) — title is "…for **tori**", not "of the torus".
- **P1 `:26`** — Fock 1935 credited for **SO(4,2)**; Fock is SO(4). Paper 7 gets
  this right.
- **P38 `:655`** — "the basis **introduced by** Avery" overstates; Wen–Avery
  derive properties, the basis is Fock 1935 / Bander–Itzykson 1966.

---

## Part C — spills outside trunk (the reason this is not a trunk-only sprint)

**C-1. The retired-scaling zombie: ~20 loci, 7 documents.**

| document | target | loci |
|:---|:---|--:|
| Paper 22 | group3 | 7+, incl. a **Corollary 1 titled "Potential-Independent O(Q^2.5)"** |
| Paper 19 | group2 | 4, incl. the **abstract** |
| group3 synthesis | trunk + group3 + synthesis | 4 |
| Papers 17, 24, 31, 58 | group2 / group3 | 1 each |

Canonical: `N_Pauli = 27.90 × Q` **exactly linear** across molecules at fixed
basis; equal-qubit advantage **54×–317×**. The claim's *form* changes (power law
→ exact linearity), so this is a rewrite per locus, not a substitution — and per
the U1 upgrade note it is a **stronger** statement, not a cut.

**Consequences to face explicitly:**
- **group3's 2026-08-29 certification is compromised.** The exact-rule
  correction landed the same day and this survived it. group3 will need a fresh
  cycle regardless of the sweep order in `recert_sweep_plan.md`.
- Paper 22's Corollary is *titled* with the retired exponent — that is a
  structural edit, not a numeral swap.
- My own v5.2.4 group2 sprint walked past the Paper 19 loci because I worked
  from C21's list and C21 was blind. Re-check group2 after A2/A3 land.

---

## Part D — PI adjudications (nothing mechanical is blocked on these)

> **ADJUDICATED 2026-09-01 (PI direction).** All four resolved; each
> resolution is recorded inline below. D1 and D2 are executing now, D3 is
> written into `qa.md`, D4 is an accepted ceiling with its cost stated.

**D1 — C3: does prose hedging discharge "tier inline"?**
All six trunk documents carry **zero** inline tier tags. The convention is alive
elsewhere in the corpus (Paper 59: 46, Paper 60: 25, Paper 58: 9), so this is a
gap, not a house style. C3 had a surface for two claims out of twenty-four and
is recorded **UNMEASURED**, not passed.
→ Either (a) prose hedging suffices and C3 should be reworded for the trunk, or
(b) the trunk owes a tagging pass across six documents. **This is the single
largest open question from the run** and it decides whether a future trunk PASS
is even well-defined.

> **RESOLVED (b) — the trunk owes the tagging pass.** Option (a) was
> rejected on mission grounds: prose hedging is exactly what failed for κ
> over ~50 versions — a hedge is a judgment the reader must re-derive, a
> tier is a claim the reader can check. Four `claims-reviewer` agents
> dispatched 2026-09-01 across P0+P7 / P1+synthesis / P32 / P38, each
> instructed to derive the tier from `docs/claim_test_matrix.md` + the
> backing test rather than from how confident the sentence sounds, and to
> report — never silently soften — any sentence the honest tier
> contradicts. The pass is therefore itself an audit: prose that no tier
> can carry surfaces as a defect instead of being quietly re-worded.

**D2 — should syntheses have claim-matrix rows?**
The group3 synthesis has **zero** rows: 1,511 lines invisible to C1/C2 *by
construction*. That is the structural reason a withdrawn-theorem sentence could
sit in its opening paragraph. Same presumably holds for the other six syntheses.

> **RESOLVED: yes.** Rows for `group3_foundations_synthesis.tex` are being
> added in the same pass (trunk scope). Rule adopted: a synthesis row
> carries the tier of the *source paper's* claim, never a fresh judgment —
> so a synthesis stating something at a stronger tier than its source is a
> defect by construction, which is precisely the property the zombie
> exploited. **The other six syntheses are NOT closed here** — each belongs
> to its own group's re-cert, and claiming a corpus-wide fix from a trunk
> sprint would be the overclaim this gate exists to catch.

**D3 — citation tiering.**
Sonnet: 0/2 seeds, 4 defects. Opus: 2/2 seeds, ~10 defects, on identical files.
→ Recommend Opus for any target with >50 bibitems, or keep two seeds and budget
for re-dispatch. Worth writing into `qa.md`'s tiering rule.

> **RESOLVED: done 2026-09-01.** Written into the Model-tiering rule in
> `.claude/commands/qa.md` as a **Citation exception**: raise citations to
> Opus above ~50 bibitems. The write-up keeps the detail that matters most
> for future runs — the Sonnet agent was *not* lazy (it enumerated all 116
> bibitems and its four findings were genuine), so the failure mode is a
> silent miss *inside* the correct search, which only the **two-seed rule**
> makes visible. At one seed, 0/1 reads as noise rather than de-calibration.

**D4 — coverage ceiling.**
~5,000 lines were quoted by no reviewer. The findings are a **lower bound** on
the trunk's defect count. Accept the ceiling, or fund Part E.

> **RESOLVED: accept the ceiling; Part E stays unfunded for now.** The
> consequence is stated plainly rather than allowed to fade: a trunk PASS
> earned after this remediation certifies **the surface that was actually
> read**, and ~5,000 lines were read by no reviewer. The largest unread
> blocks are named in Part E (P32 `:2104-3164`, ~30 consecutive Q5-prime
> remarks with zero matrix rows; P32 `:6105-6980`, the Lorentzian/Krein
> block where a descope-zombie would live; synthesis `:474-1183`). The
> deterministic layer (C16/C17/C21 registries) *does* scan those lines
> exhaustively, so the known retracted-claim and retired-number classes are
> guarded there; what is unguarded is anything needing judgment. Part E is
> retained as a funded-on-request list, not deleted.

---

## Part E — coverage debt (optional; the critic's re-dispatch list)

1. **Run the 20 unrun Paper 32 inline-cited tests**, especially
   `test_real_structure.py` + `test_connes_axiom_audit_31.py` — the only backing
   for the §IV Connes axiom audit, a keystone.
   *(Done at close: `test_trunk_qa_fejer_4_over_pi.py` and
   `test_trunk_qa_forced_count_moduli.py` — 22 pass, 1 skip, no defect. So two
   of three C8 backings are now measured.)*
2. **P32 `:2104-3164`** — the Q5′ remark chain, ~30 consecutive remarks, zero
   matrix rows, no section headings. Largest unaudited block in the target.
3. **P32 `:6105-6980`** — Lorentzian/Krein block; 6 cited tests, none run. This
   is where a Lorentzian descope-zombie would live.
4. **P7 `:536-643`** (N-electron §VI) and **`:876-909`** (v2.6.0 note) — 12+
   load-bearing numbers with no matrix rows.
5. **synthesis `:474-1183`** — the body, ~700 lines, essentially unenumerated.
6. ~~**Reverse-citation pass**~~ — **CLOSED 2026-09-01.** Measured at **30**
   prose references with no `\bibitem` across the six trunk documents (not
   broken `\cite` — pdflatex reports zero undefined citations, before and
   after). Worst: Papers 54-57, the synthesis's entire Reconvergence
   section; Paper 32 naming Paper 38 ten times and Paper 55 fourteen times
   with no entry; Paper 0 naming the archived Paper 6 in the present tense.
   Fixed: **+30 bibitems, +26 inline cites**, titles read from each paper's
   own `\title` (C11 PASSES on all 30). A year audit run on top found **6**
   mismatches against the cited papers' own `\date`, **three of them
   pre-existing in Paper 32** — the twin class, in a dimension C21 does not
   cover. All corrected.
7. **"18 independent symbolic proofs"** — Paper 7's abstract was softened
   ("several of them limiting-case or definitional"); the synthesis still says
   "18 **independent**" at `:52`, `:151`, `:327`. Matrix row 40 directs the
   softening.

---

## Part F — FULL run #2 (dispatched 2026-09-01, scored 2026-09-02) = FAIL

> Run record: `debug/qa/trunk_full_run_2026_09_01_notes.md` §5 (calibration
> 19/21 seeds, 0/8 false positives; every MATERIAL below was verified against
> the primary text on MAIN; seed-location findings set aside). Line numbers
> are MAIN at `dbe7ae2`. Remediation is its own sprint, as before; the cert
> path is **Part F → DELTA (F-seeding rule) → clean → FULL**.

### F-seeding rule for the DELTA run (from the calibration result)

> **Superseded in part (2026-09-02, PI direction):** seeded calibration is
> now opt-in (`/qa <target> seeded`). Items 1 and the three-seed clause of
> item 2 apply only to a seeded run; the Opus tier for trunk-root code
> (item 2) and paste-only-the-frozen-criteria (item 3) stand. DELTA #3 runs
> unseeded unless the PI says otherwise.

Two Sonnet code chunks were PARTIAL (1/2): CODE-A (P0+P7 tests) missed the
S3 class — a manufactured convergence endpoint / false-positive-by-
construction — and CODE-C (P32 tests) missed the S2 class — a vacuous or
self-referential assertion. Both misses were whole-file reads that ended
BACKED-SOUND, i.e. the reviewer trusted the docstring's account of what the
assertion checks. Every CODE catch in the run was GUIDED (the wave-2 prompts
carried class-level hints). The DELTA therefore must:

1. plant a **fresh S3 seed in the P0/P7 backing tests** and a **fresh S2 seed
   in the P32 backing tests**, blind — no class hint in the prompt;
2. use **Opus** for the code dimension on trunk roots, or **3 seeds per Sonnet
   code agent**;
3. paste only the frozen criteria into CODE prompts.

Until that DELTA comes back, the trunk's code dimension is certified only on
the classes its reviewers demonstrably caught.

### F1 — Paper 32 (`papers/group1_operator_algebras/paper_32_spectral_triple.tex`)

| # | locus | defect | fix |
|:--|:--|:--|:--|
| F1.1 ✅ **DONE 2026-09-02 (v5.3.0)** — rescoped to Remark `rem:D_GV_no_graph_form`; `def:D_GV_graph` + `prop:D_equiv` withdrawn, every `\ref` retargeted; new `tests/test_paper32_dirac_module_contents.py` (11 tests) pins the backing modules carry no graph-form operator; C16 `dgv-graph-form-tautology`. | :766-811 `def:D_GV_graph` + `prop:D_equiv`, :816 | **LARGE by class (false backing).** The "graph form" Dirac is defined by isospectrality and the Proposition is then a tautology; the cited backing (`geovac/dirac_matrix_elements.py`, `tests/test_dirac_matrix_elements.py`, `tests/test_qed_self_energy.py`) contains no graph-form operator; the only implemented one (`geovac/dirac_lattice.py::DiracLattice`) is NOT isospectral to CH at nonzero hopping. | **PI adjudication (Part D-style):** rescope Def+Prop to a Remark (recommended) — every axiom check and the π-free certificate use the spectral form; `DiracLattice` is a different operator — or build the operator the definition promises. Either way fix :816, every `\ref{prop:D_equiv}` / `\ref{def:D_GV_graph}`, and add `tests/test_paper32_dirac_graph_form.py`. |
| F1.2 ✅ **DONE 2026-09-02 (v5.3.0)** — relabelled to the KO-3 pair, "hence" dropped, construction unchanged; Dąbrowski–Dossena verified against the primary (additive rule = graded product, eq. 10; Table 5 gives no consistent ε′ for it) and cited as `dabrowski_dossena2011`; extra: `fluctuated_dirac` default ε′ = −1 / `U.T` → +1 / `U†` in both AC modules (errors had cancelled, no number moved), pinned by `test_fluctuated_dirac_J_compatible`; C16 `combined-triple-ko1-additive`. | :4049-4051, :4531-4534; `tests/test_almost_commutative.py:11,403,410-427` | **LARGE — KO-dimension label.** "3 + 6 = 9 ≡ 1 (mod 8), hence J² = −1, JD = +DJ": the verified sign pair (ε, ε′) = (−, +) is the KO-**3** pair (KO-1 is (+, −)). Tests verify the signs at 1e-12, not the label. | **PI adjudication:** minimal honest fix = state the verified signs as the KO-3 pair, drop "hence", fix the test docstrings, construction unchanged; or change the construction. Verify Dąbrowski–Dossena (IJGMMP 8 (2011) 1833, arXiv:1011.4456) against the primary before citing it for "naive J⊗J does not realise KO-additivity". |
| F1.3 ✅ **DONE 2026-09-02** — six loci renamed to “state-space GH rate”; C16 `latremoliere-propinquity-named-for-gh-rate` widened with a `propinquity\s+rate` alternative, two-way proven. | :1986-1987, :2132, :2140-2141, :3198-3199, :3569, :6874 | C7 class: "propinquity rate" / "state-space propinquity rates" for the paper's own state-space GH rate. | Rename to "state-space GH rate" / "GH-rate cross-check" at all six loci, THEN widen C16 `latremoliere-propinquity-…` (or the group1 entry) with a `propinquity\s+rate` alternative and record the two-way proof. |
| F1.4 ✅ **DONE 2026-09-02** — three loci retagged [SYMBOLIC + MEASURED]; “128 per generation” → 128·N_gen² (mixing) / 8 per generation (diagonal slice). | :6477, :3758-3762, :5113-5201 (`thm:forced_count`) | C3-boundary: symbolic argument + numerical certificate tagged [SYMBOLIC PROOF] / [PANEL-VERIFIED]. | Retag [SYMBOLIC + MEASURED]. :5418 "128 per generation" → 128·N_gen² (mixing) / 8 per gen (diagonal). |
| F1.5 ✅ **DONE 2026-09-02** — “provably non-overlapping” → “rings sharing no generator”; “seam theorem … prove” → “seam reading … place”. | :4868-4887 vs :5216-5222 | [OBSERVATION] seam "reading" stated as "provably non-overlapping" / "seam theorem … prove". | :4882 "provably" → "in rings sharing no generator"; :5216 "seam theorem … prove" → "seam reading … place". |
| F1.6 ✅ **DONE 2026-09-02** — [CONDITIONAL] added; “on the η-trivialised inner factor” stated twice. | :4207-4213 | untiered "forbidden … vanish identically on the inner factor of ANY CC-compatible AC extension" (physics inference). | Add [CONDITIONAL] + "on the η-trivialised inner factor". |
| F1.7 ✅ **DONE 2026-09-02** — “constructed as … and its Hamiltonian is validated against …”. | :1386-1388 | [MEASURED] backs the hyperfine gap, not spectral-triple-hood. | "[MEASURED] it is constructed as an explicit Connes-style … spectral triple, and its Hamiltonian is validated against …". |
| F1.8 ✅ **DONE 2026-09-02** — “(numerically coincident with the Fock Jacobian factor 1/16; an Observation, Paper 0)”, now with `\cite{paper0}`. | :3149-3150 | κ appositive "(the Fock Jacobian Ω⁻⁴)" asserts the C8-forbidden bridge. | "(numerically coincident with the Fock Jacobian factor 1/16; an Observation, Paper 0)". |
| F1.9 ✅ **DONE 2026-09-02** — replaced by the proven rate (4/π)log n/n + O(1/n), citing Paper 38. | :3380-3382 (inside `thm:gh_convergence`) | STALE: "O(log n/n) … not rigorously proved … deferred" vs Paper 38 unconditional (2026-06-10). | Replace with the proven statement, cite Paper 38. |
| F1.10 ✅ **DONE 2026-09-02** — retagged [SYMBOLIC PROOF] by construction; the `test_dirac_matrix_elements.py` citation withdrawn (it covered the κ↔(l,j) bridges, not the axiom). | :878-882 order-zero axiom | mis-pointed backing (`test_connes_axiom_audit_31.py::verify_order_zero` is the Lorentzian U_L check). | Retag [SYMBOLIC PROOF] "by definition" and drop the cite, or add a 5-line test. |
| F1.11 ✅ **DONE 2026-09-02** — asserts pinned per cell to the measured residuals (order-zero 0 / 0.050661 / 0.067547; order-one 0 / 0.101321 / 0.101321, ±5 %) replacing `< 0.5`; 6/6 pass. | :6518-6521 axiom table (vi)/(vii) | prints "≤ 0.0675" / "≤ 0.101"; `tests/test_connes_axiom_audit_31.py:291-313, :328` assert `< 0.5`. | Run with printed residuals; tighten the asserts to the printed values with margin. |
| F1.12 ✅ **DONE 2026-09-02** — measured 2026-09-02: SU(3) 47 pass + 2 skipped of 49; modular Hamiltonian 70 (66 + 4 slow); almost-commutative 54; standard-model triple 45 — all four loci fixed. | :3715 "39/39"; :7028-7029 "67"; :4053 + :4167 "38 tests passing" | count drift (actual 49 = 47+2 skipped; 70 = 66+4 skipped; 53 collected/53 pass). | Fix all four. |
| F1.13 ✅ **DONE 2026-09-02** — TOC of Connes–Marcolli 2008 fetched (AMS Colloq. 55; four chapters): cosmic Galois group = Ch. 1 §7.3, finite geometry = Ch. 1 §13; no theorem number cited. | :2738-2739 "Ch. 4"; :4032-4033 "Ch. 13" | Connes–Marcolli 2008 locators (the book has four chapters). | Ch. 1 §1.7 (Thm 1.100); §1.13 — verify against the TOC. |
| F1.14 ✅ **DONE 2026-09-02** — two projection-rank tests added (matter blocks of the 272- and 260-dim solution spaces both span 128; pass); the “(45 tests)” pointer now names `tests/test_standard_model_triple.py`; matrix rows 50/104 updated. | :5172-5176 | coverage gap: matter-sector "512→256→128; order-one / J-reality add nothing" not computed by `test_trunk_qa_forced_count_moduli.py`. | Add the two constraints to the test if cheap; else name the gap in the claim matrix. `geovac/standard_model_triple.py` "(45 tests)" at :5183 — matrix row 48 says mis-pointed; fix the pointer. |
| F1.15 ✅ **DONE 2026-09-02** — bibitem keys → paschke_sitarz1998 / karamata1949 / cacic2013 / chamseddine_connes1997; 5 orphans cited (paper0, loutey_paper40/50, BW 1975, Witten 2018); 7 untagged paragraphs tagged; “nine attempts (twelve counting Sprint A and K-CC)”; “four orders” → three (1.0000 → 0.0007); “six to ten orders” → three-to-four (memo table: 2.7e-4 of D_e); c²(3,2) → c²(4,3); g₃ = 40 = both-sign degeneracy at k = 3 (was “third single-chirality”); dim_H(3) = 40 vs g₃ = 40 flagged as a numerical coincidence at three loci (incl. P42); wall count: seven catalogue entries + G6 = eight entries vs six instances; 14 of 15 CLAUDE.md mentions reworded (first mention kept as the conventions note). Not identified: “H10 sub-locators”, the :1119/:1773 items, `test_s3_eigenvalue_ground_state` docstring — carried. | NIT sweep | :134, :1773; untagged :110-120 / :198-205 / :2098-2113 / :3273-3300 / :1555-1558; :200, :2208; :1708 "nine attempts" vs 12; :3698 "three non-abelian"; :1119; :4280 "six to ten orders" vs 5.3; :4290 "four orders" vs 3.15; :5699/:5714/:5716 six-vs-seven walls; g₃ 40 vs 24 notation; :5633-5635 "graph Laplacian" for Laplace–Beltrami; :1681/:1699 g₃ wording vs P7:732-744; bibitem years paschke_sitarz2000→1998, karamata1962→1949, cacic2009→2013, chamseddine_connes2010 (verify each); orphans :7410, :7436; 15 `CLAUDE.md` mentions + sprint codenames (audience register); H10 sub-locators (verify or drop). Upgrades (optional): :76-86, :1051-1070, :6121-6123, :5810-5832, :5048-5059, :6870-6871. |

### F2 — Paper 38 (`papers/group1_operator_algebras/paper_38_su2_propinquity_convergence.tex`)

| # | locus | defect | fix |
|:--|:--|:--|:--|
| F2.1 ✅ **DONE 2026-09-02 (v5.3.0), propagation included** — Remark `rem:circle_fejer` + closed form `eq:circle_fejer_moment` (2/π, probability-normalised, quadrature-vs-closed-form + "twice" tests added); "§I.1" dropped; App. A Step 3 bookkeeping replaced; Paper 40 torus corollary withdrawn + two `2Vol(S¹)/Vol(SU(2))` slips fixed; group1 synthesis + outreach note N1 corrected; memory updated; C16 `circle-fejer-constant-4-over-pi`; matrix row 52 note + new row. | :941-950 remark, :1702-1708, :1141/:1288 locators | **Circle Fejér constant wrong by 2.** Probability-normalised Fejér first moment is (2/π) log N/N, not 4/π; the corpus's own test `tests/test_trunk_qa_fejer_4_over_pi.py:235-268` names 2/π and asserts the SU(2) constant is twice it. "Stein–Weiss §I.1" is Rⁿ L¹ theory, not a circle moment. | Circle constant → 2/π with normalisation stated; "the SU(2) constant is twice the circle constant" as an Observation (2·Vol(S²)/Vol(SU(2)) mirrors 2·Vol(S⁰)/Vol(S¹)); drop "§I.1" (Zygmund Vol. I Ch. III may stay, drop "§3.6"); add a closed-form-vs-quadrature circle test to `test_trunk_qa_fejer_4_over_pi.py`. Check memory `l2_quantitative_rate_4_over_pi.md`. **Propagation → group1 carryforward:** Paper 40 :971, :977-989 `rem:stein_weiss_general`, :966-973 torus corollary, :1937, :1943-1944; group1 synthesis :250, :783, :1284, :597-599. Read before editing. |
| F2.2 ✅ **DONE 2026-09-02 (v5.3.0)** — now "twice the Haar-normalisation ratio Vol(S²)/Vol(SU(2)) = 2/π". | :283-284 | arithmetic: "Vol(S²)/Vol(SU(2)) · (2/π) = 4/π" gives 4/π². | `(2/\pi)` → `2`. |
| F2.3 ✅ **DONE 2026-09-02** — `test_central_fejer_su2.py` monotonicity sample extended to n = 200, 500, 1000 (2.3 s at n = 1000); prose narrowed to the tested sample. | :934-935 | "verified at n ∈ {2,…,1000}" vs `test_central_fejer_su2.py:647-659` up to 100. | Extend the test to 200/500/1000 if cheap (measure), else narrow the prose. |
| F2.4 ✅ **DONE 2026-09-02** — inline cite added (`test_eq_l3_main_inequality`, `test_l3_C3_uniform_in_n_max`); tier left [PANEL-VERIFIED] (no Paper 40 §3.3 retag). | :980-1077 Lemma L3 | no inline test cite (`tests/test_r25_l3_lipschitz_bound.py` exists; matrix :221 names it); CLM-C proposes upgrade via Paper 40 §3.3. | Add the cite; verify Paper 40 §3.3 before retagging [INTERNAL THEOREM]. |
| F2.5 ✅ **DONE 2026-09-02** — function cited: `::test_spatial_kernel_condition_truthful_vs_offdiag`. | :553-554 | whole-file cite of `test_p45_kplus_degeneracy.py` (1 of 5 fns relevant). | Cite the function. |
| F2.6 ✅ **DONE 2026-09-02** — `paper43` bibitem added (title verbatim from Paper 43); `loutey_paper40` merged into `paper40_unified`; Latrémolière keys → 2015 (JMPA) / 2016 (TAMS; venue was already correct on MAIN). | :293-295; :1881 `paper40_unified` / :1935 `loutey_paper40`; :1840 | Paper 43 cited inline, no bibitem; duplicate Paper 40 bibitems; `latremoliere2018` venue "JFA 368" → Trans. AMS 368 (2016) 365–411 (body :591 is right). | Add `paper43`; merge to one key (re-key :1520/:1565); fix venue. |
| F2.7 ✅ **DONE 2026-09-02** — main theorem tagged [INTERNAL THEOREM]; L1′ tagged [SYMBOLIC + MEASURED]; perez_sanchez2024 cited, chamseddine_connes2010 removed; gaudillot → IMRN 2025 rnaf197 (key gaudillot_vs2025), hekkelman → LMP 112, 20 (2022), both publisher-verified; shell index unified to the 1-based convention of eq:CH_spectrum (n(n+1), shells n ≤ 2J+1 = n_max). | NIT sweep | abstract bare tier; :472 `thm:main_unconditional`, :743 `lem:L1prime` untagged; orphans :1750, :1912; preprint forms :1781 gaudillot_vs2023 (IMRN 2025 rnaf197), :1776 hekkelman2022 (LMP 112, 20) — verify vs publisher; n-index :433 (n+1)(n+2) vs :625 n(n+1). |

### F3 — Papers 0, 1, 7 and their tests

| # | locus | defect | fix |
|:--|:--|:--|:--|
| F3.1 ✅ **DONE 2026-09-02** — test rewritten: 4×3 Jacobian of the embedding, g = JᵀJ = Ω²I (nine entries) and det g = Ω⁶, symbolic; P7:79 softened; matrix row 42 updated. | `tests/test_fock_projection.py:330-356` `test_volume_element_jacobian` | (a/b)³ = a³/b³ restatement, not a Jacobian computation (1 of the 18). | Compute √det g of the stereographic embedding symbolically and compare. Soften P7:79 "complete, machine-verified algebraic audit" (matrix row 42). |
| F3.2 ✅ **DONE 2026-09-02** — nnz bound `< 1.5` → `< 1.15` (measured 1.07 over n_max = 5–30); timing guard kept; P1:449 reworded. | `tests/test_ov_scaling_rigorous.py:169 < 1.8`, `:216 < 1.5`; P1:449 | test guards sub-quadratic; prose says O(V) verified. | Tighten :216 to < 1.15 (nnz exponent deterministic); keep a timing guard at :169; P1:449 → "measured exponent ≈1.05; the frozen test guards sub-quadratic timing and near-linear nnz". |
| F3.3 ✅ **DONE 2026-09-02** — “SO(4) [fock1935], extended to SO(4,2) [barut1967]”. | P1:26 | "dynamical symmetry group SO(4,2) \cite{barut1967,fock1935}" — Fock 1935 is SO(4). | "SO(4) [fock1935], extended to SO(4,2) [barut1967]" (P7:63-67 has it right). |
| F3.4 ✅ **DONE 2026-09-02** — four “exact Rydberg spectrum” loci → quantum numbers + matched-κ Observation; K sentence given its own [OBSERVATION]; Condon–Shortley phase attached to Y_lm; Bohr = L = nħ, de Broglie = integer wavelengths. | P1:15/54/59/363; P1:393-398; P1:271; P0:749-750 | "reproduce the exact Rydberg spectrum" (→ quantum numbers, per P7:689); K sentence inherits [CONJECTURE] → own [OBSERVATION]; Condon–Shortley wording; Bohr/de Broglie attribution. | As stated. |
| F3.5 ✅ **DONE 2026-09-02** — `gaudillot_vs2025` bibitem added to Paper 7 and the group3 synthesis; inline arXiv IDs replaced by cites. | P7:748, SYN:354 | inline "Gaudillot-Estrada and van Suijlekom (arXiv:2310.14733)" (ID verified), no bibitem. | Add `gaudillot_vs2023` to both; re-run C20. |
| F3.6 ✅ **DONE 2026-09-02** — abstract/series duplicate tags removed (P0:27, :32, :104); “follows from” → “tracks”; Paper 6 marked archived (×2) with `\cite{paper1}` added; “~6 %” → 5.3 %; P1 log-holonomy stated as a closed form (no R²); “all α-independent physics” → single-electron structure; P7 “confirming” → “consistent with” + `\cite{loutey_paper18}`. Docstring item not identified — carried. | NIT sweep | P0:22/27, 30/32, 103/104 duplicate tags; P0:530 "follows from" under PANEL-VERIFIED; P0:711-713; P0:833 "~6%" vs 5.3%; P0:684, 719-720 "developed in Paper 6" (archived); P1:261 "R²=1.0" not computed; P1:366 "all"; P7:914 "confirming"; P7 orphan `loutey_paper18`; docstrings `test_s3_eigenvalue_ground_state`. |

### F4 — group3 synthesis (`papers/synthesis/group3_foundations_synthesis.tex`)

| # | locus | defect | fix |
|:--|:--|:--|:--|
| F4.1 ✅ **DONE 2026-09-02** — provenance note scoped to unitarity + (b), pre-rename filename named, (a)/(c)/(d) declared archived with the matrix gap; “10⁴ steps, machine precision” → 10⁻¹⁰ over 10³ steps live (10⁴ in Paper 6’s archived H₂ run); asserts pinned: peak > 0.999 (0.999756), period < 0.45 % (0.4106), Rabi norm < 1e-10 (2.9e-13), off-resonance < 0.01 (0.003725). | :587-594 provenance note; :601-603 | Note claims the section's results "were re-run and reproduce" with backing in `tests/test_rabi_oscillation.py`; the file backs only unitarity + Rabi (b), not (a)/(c)/(d) (matrix row 133 NO-TEST, open). ":over $10^4$ time steps" vs the live test's `n_steps = 1000` (`:195`, max_n=10, threshold 1e-10); 10⁴ is Paper 6's archived run. | Scope the note to unitarity + Rabi; name the pre-rename filename (6e6ce40); leave (a)/(c)/(d) as archived Paper 6 measurements with the row-133 gap named. State "10³ steps in the live regression (10⁴ in Paper 6's original run)" or raise `n_steps` and re-measure. Tighten `run_off_resonance` `< 0.5` (measured 0.0037); pin the Rabi assertions nearer 99.98 % / 0.41 %; "machine precision" vs 1e-10. |
| F4.2 ✅ **DONE 2026-09-02** — “constructs the cosmic-Galois comparison map … (a homomorphism, not a closed immersion)”. | :164-165 | "Paper 56 closes … U*_GV ↪ U_4^ab ⋊ SL_2" — P56:1298-1302 says Φ is not injective (factors through the abelianisation); SYN:1096-1099 itself says "a homomorphism (not a closed immersion)". | "closes the injection direction of the comparison … (a homomorphism, not a closed immersion)". |
| F4.3 ✅ **DONE 2026-09-02** — year → 2002; key → `aquilanti_caligiana2002`. | :1410 `aquilanti_caligiana2003` | CPL 366, 157 is 2002, not 2003. | Verify vs publisher, fix year (keep or rename the key consistently). |
| F4.4 ✅ **DONE 2026-09-02** — → “(group 5, QED/gauge)”. | :988 | "Paper 25 (synthesis group)" — it is group5 (:692 right). | Fix. |
| F4.5 ✅ **DONE 2026-09-02** — Paper 38’s unconditional state-space GH theorem added at both loci with the correct metric name. | :176, :1306-1309 | omit Paper 38's unconditional state-space GH theorem. | Add, with the correct metric name (state-space GH, not propinquity). |
| F4.6 ✅ **DONE 2026-09-02** — KPS moved to the F-theorem clause; `Suhonen2007` (Talmi–Moshinsky) and `biedenharn1981` (dynamical algebra) cited; Fathizadeh–Marcolli, Glanois, Hain and Brown bibitems added and cited; five orphans deleted (Whitten, Dunlap, Peruzzo, Lee, chamseddine_connes1997); :700-704 hedge rewritten (graph spectrum → continuum spectrum numerically). Left as-is: :381 wording, :867 Weyl–Selberg (no year, not an attribution C20 flags), :1149. | NIT sweep | :69 KPS for "conformal field theory", :381 F_{S^d} attribution (soften); inline attributions without bibitem :546 Talmi–Moshinsky, :867 Weyl–Selberg, :1073 Fathizadeh–Marcolli (load-bearing; bibitem exists in Paper 55), :1076, :1122/:1154, :1149; orphans chamseddine_connes1997, Whitten1973, Dunlap2000, Suhonen2007 (→ :546), Peruzzo2014, Lee2021, biedenharn1981 (→ :236); :700-704 hedged phrasing; "18 independent" at :52/:151/:327 (Part E item 7, still open). |

### F5 — coverage debt the critic recorded (carried, not fixed here)

- **UNMEASURED per criterion × document** (run record §5.6 table): C7 for
  P0/P1/P7 (zero occurrences examined); C5 prose side for P0/P38, artifact
  side (`tests/test_paper2_corrections.py`, matrix row 47) corpus-wide; C8 κ
  for P1/P38/P32:3800-7527; C8 4/π for P0/P1/P7; footnotes, table cells,
  appendices and `\paragraph` strings as listed; 83 internal bibitems had no
  mandated reviewer (title/year drift is C11-deterministic; status
  descriptors were hand-swept for Papers 45-49 only).
- **Unopened backing tests:** `tests/test_paper2_corrections.py`,
  `tests/test_r25_l3_lipschitz_bound.py`, `debug/p38_g1g2_*.py`,
  `geovac/standard_model_triple.py`, the 19 synthesis-row tests; partially
  read `test_qed_self_energy.py`, `test_modular_hamiltonian.py:261-1074`.
  Open NO-TEST rows 133, 49, 88. → DELTA targets alongside the F-seeding rule.
- **Named by the DELTA (2026-09-02, CODE-A):** Paper 7's [MEASURED] "graph
  Laplacian converges to the S³ Laplace–Beltrami operator numerically through
  n_max = 30" (P7:77) has no backing test at all — no test builds the
  n_max = 30 spectrum for it and the paper carries no inline cite. Logged as
  a NO-TEST row in `docs/claim_test_matrix.md` (Paper 7); closing it is a
  new test, not prose. ✅ Closed 2026-09-02 (v5.4.0):
  `tests/test_paper7_graph_convergence.py`; see G1.2.
- **Unread regions:** P38 Appendix A Steps 1-4 (:1622-1690, the 4/π chain —
  F2.1 reaches it only through the body); P38 §1 :176-203 (C7 prose); P32
  :1670-3669 single-reviewer; P7 :471-493, :574-604, :879-891; P1 Appendix A
  :409-434; P0 Tables 1/2.

### F6 — PI adjudications (ALL FOUR ADJUDICATED + EXECUTED 2026-09-02, "across the board" on the PM's recommendations)

1. **F1.1** `prop:D_equiv` — Remark (recommended) or build the operator.
   → **Remark.** Done (F1.1 row).
2. **F1.2** KO-dimension label — relabel to the verified KO-3 pair
   (recommended) or change the construction. → **Relabel.** Done (F1.2 row).
3. **F2.1 propagation** — the 2/π correction reaches Paper 40 and the group1
   synthesis; logged as group1 carryforward, not fixed under trunk.
   → **Fixed now, at source, under this sprint** (deferral is churn): Paper
   40, group1 synthesis, outreach note N1. Nothing left for the group1
   carryforward on this item.
4. Whether the run-#2 verdict (a second FULL FAIL, with the QA-gate change to
   C16 and the F-seeding rule) warrants more than the patch bump the PM will
   apply. → **v5.3.0 (minor)**, because the code-tier exception is a QA-gate
   change; the run itself stays v5.2.7.

### F7 — findings the F6 remediation surfaced (2026-09-02, v5.3.0)

| # | locus | finding | disposition |
|:--|:--|:--|:--|
| F7.1 | P40 :289 and `eq:su2_4_over_pi` (~:1892) | Second and third loci of the `2Vol(S¹)/Vol(SU(2))` slip (= 2/π, presented as 4/π); `debug/review_paper38.md:54` had flagged the identity in P38 ("A7") and P40 kept it. Group1 synthesis ~:614 carried the same. | Fixed: `2Vol(S²)/Vol(SU(2)) = Vol(S²)/π²`. C16 pattern covers the slip form. |
| F7.2 | P40 `cor:general_compact_connected` torus extension (:966-973) | A circle factor carries 2/π, so a product-with-torus constant is mix-dependent — the extension was never true. | **Withdrawn**; corollary scoped to semisimple G. |
| F7.3 | P38 App. A Step 3 (:1622-1690, the region F5 lists as unread) | The "−π²n log n/16 × 2" bookkeeping did not produce the coefficient; the log n coefficient is (A) triangle truncation of the odd-d sum (+2) net (B) the d/2 term of √(a(a+d)) (−1). | Rewritten; (A) and (B) each test-pinned in `test_trunk_qa_fejer_4_over_pi.py`. Closes the F5 "App. A Steps 1-4 unread" item for Step 3 only. |
| F7.4 | `docs/outreach/note_n1_su2_truncations.tex:87` | "the same constant" (circle = SU(2)) — the retired claim outside `papers/`. | Fixed at source; file added to the C16 entry's scope. |
| F7.5 | `geovac/almost_commutative.py`, `geovac/standard_model_triple.py` `fluctuated_dirac` | ε′ default −1 and `U.T` for `U†`: two cancelling errors (combined U purely imaginary). No published number moved. | Fixed; `test_fluctuated_dirac_J_compatible` pins J-invariance of D_ω at nonzero ω. |
| F7.6 | `almost_commutative.py` docstring | KO-dim set for J_F² = +1 given as {0, 4, 6, 7}; KO-4 has ε = −1. | → {0, 1, 6, 7}. |
| F7.7 | Dąbrowski–Dossena 2011 | The 3+6≡1 additive rule is a theorem about the *graded* product; the module builds the *ungraded* sum, for which signs follow factor by factor. Verified against the primary before citing. | Stated in P32 + module docstring; bibitem `dabrowski_dossena2011`. |
| F7.8 | C21 | Not extended for 2/π vs 4/π: the numeric registry is keyed on decimal numerals and the papers carry these symbolically. | C16 is the instrument; recorded so the next run does not read it as an omission. |
| F7.9 | C16 registry | Three entries, each two-way proven: `circle-fejer-constant-4-over-pi` 23 fire / 0 live / 4 exempt; `dgv-graph-form-tautology` 9 / 0 / 2; `combined-triple-ko1-additive` 13 / 0 / 4. `--gate trunk` PASS (14/33 entries, 49 loci, 54 exempt); `--gate group1` PASS (10/33, 34 loci, 55 exempt). | Registry 30 → 33. |
| F7.10 | P40 abstract :85, `thm:main_intro` :217, “Throughout” :332, `thm:main` :1771, :253; P38 :1603; P32 :1982 | The class “compact connected Lie group of rank r ≥ 1” includes tori, but the circle constant is 2/π (F2.1) and the dual-Coxeter normalisation is only defined on simple summands — the withdrawn torus corollary (F7.2) had left the main-theorem hypothesis unscoped. | Scoped to compact connected **semisimple** G at all loci (group1 synthesis already said “simple”); P38 and P32 pointers say so too. |
| F7.11 | P32 :1683, :1701, :3290, :6888; P42 :2303 | “dim_H = 40 = g₃^Dirac = Δ⁻¹” presented as structural: dim_H(3) is the cumulative node count 4 + 12 + 24 while g₃^D = 2·4·5 = 40 is the single-level degeneracy at k = 3; N(n_max) = n_max·g^D_{n_max}/3, so they agree only at n_max = 3. Also “third single-chirality” degeneracy: 40 is the both-sign count at |λ| = 9/2. | All five loci rewritten as a numerical coincidence between different objects (P42 fixed at source under this sprint). |
| F7.12 | P32 :5786 | Δ = 1/40 = c²(3,2): Paper 7 and `test_trunk_qa_c2_delta.py` define the value as c²(4,3). | → c²(4,3) with the convention named. |
| F7.13 | P32 :5736-5755 | “all six catalogue instances … under two theorems” followed by seven codes, then “the seventh (G6)” and “all six unique”: entries and instances conflated. | Seven entries + G6 = eight entries; six independent instances per Paper 57 §4 (H1, LS-8a, HF-3/4/5, W1e). |
| F7.14 | P32 :4309, :4319 | “six to ten orders of magnitude below wall depth” vs the F2 memo’s own table (max differential 2.7×10⁻⁴ of D_e = 0.075 Ha ≈ 3.6 orders; the memo’s prose says “six”); “four-orders-of-magnitude advance” for 1.0000 → 0.0007 (3.15). | Both restated from the numbers. |
| F7.15 | P38 bibliography | Renaming `latremoliere2018` → `latremoliere2016` collided with the existing `latremoliere2016` (the 2015 JMPA paper keyed with the wrong year); resolved from `git show HEAD` → JMPA = `latremoliere2015`, TAMS = `latremoliere2016`. | Lesson: check for an existing key before a year-rename. |
| F7.16 | `tests/test_trunk_qa_forced_count_moduli.py` | The matter-sector “J-reality + order-one add nothing” step was prose only (F1.14). | Now a computed projection rank (272 → 128, 260 → 128). |

---

## Part G — DELTA run #1 (2026-09-02, after Part F, F-seeding rule) = DEFECTS

**Shape.** DELTA-verification of the git diff `dbe7ae2` → working tree
(34 files), synced into an isolated worktree with the uncommitted edits
(the run-#1 lesson), seeded blind: 9 seeds (answer key
`debug/qa/trunk_delta_seed_key_2026_09_02.json`), 8 known-good controls.
Seven reviewers: CODE-A (Opus; P0/P1/P7 tests), CODE-C (Opus; P32/P38
tests + the three `geovac` modules), CLAIMS-1 (Opus; P32 + P42), CLAIMS-2
(Opus; P38 + P40), CLAIMS-3 (Opus; P0/P1/P7/P18), SYNTH (Opus; both
syntheses), CITE (Sonnet, two seeds; every changed bibitem/cite).
Deterministic gates whole-target.

**Calibration: 9/9 seeds caught by their intended catcher, 0/8 false
positives.** S-D1 (S3 manufactured endpoint) → CODE-A; S-D2 (S2
self-referential rank) → CODE-C — the two classes the FULL run #2 Sonnet
chunks missed, both caught by Opus; S-D3 → CLAIMS-1; S-D4 → CLAIMS-2; S-D5
(κ derived) → CLAIMS-3; S-D6 + S-D9 → SYNTH; S-D7 + S-D8 (wrong volumes)
→ CITE. Three seeds were additionally noticed out of scope (S-D6 by
CODE-A, S-D7 by CLAIMS-2, S-D8 by CLAIMS-1). No control was flagged
MATERIAL.

**Verdict: DEFECTS.** The remediated text carried 26 genuine items (G1
below), all remediated the same day; per the standing rule a further
DELTA (#2, fresh seeds) is owed before the FULL certifying run.

### G1 — genuine findings (non-seed), by reviewer

| # | finding (reviewer) | disposition |
|:--|:--|:--|
| G1.1 | Paper 7 Appendix still listed the retired Ω³ = 8p₀³/(p²+p₀²)³ identity as proof #8 after the test was rewritten (CODE-A). | Entry rewritten to the induced-metric statement g = JᵀJ = Ω²I, √det g = Ω³, with the replacement dated. |
| G1.2 | Paper 7's [MEASURED] "graph Laplacian converges numerically through n_max = 30" has **no backing test at all** (CODE-A). | NO-TEST row added to the claim matrix (Paper 7); named in F5. Closing it is a new test. ✅ **CLOSED 2026-09-02 (v5.4.0, PI: "add the test")** — `tests/test_paper7_graph_convergence.py` (5 tests, 1.1 s): ground-state energy of the production `AtomicSolver` Hamiltonian vs −1/2 on n_max ∈ {5,8,10,15,20,30}, pinned sequence (17.3% → 0.57%), monotone, one-sided, two-sided endpoint window, plus a bottom-of-spectrum density scope guard (six lowest eigenvalues within 0.06% at n_max = 30 — not a level-by-level Rydberg match); cited inline at P7:77; matrix row now MEASURED. |
| G1.3 | Rabi report strings still printed the old thresholds (< 1e-6, < 0.50) (CODE-A). | Fixed. |
| G1.4 | group3 synthesis :354 "4/π universal across compact Lie groups" (SYNTH, CLAIMS-3). | → compact connected semisimple; circle factor 2/π. |
| G1.5 | group1 synthesis :1544–1559 product-carrier sketch: "arbitrary compact connected Lie groups … inherits the universal 4/π" (SYNTH). | At least one semisimple factor; the joint rate is the max of the factor rates, so a torus factor's 2/π is dominated. |
| G1.6 | `brown2017` (ICM 2014, genus-zero) cited for the Hain–Brown mixed-elliptic-motive target (SYNTH). | → Brown, *Multiple modular values and the relative completion of π₁(M_{1,1})*, arXiv:1407.5167 (publisher/arXiv verified). |
| G1.7 | "is now a theorem" object shift (graph Laplacian vs spinor truncations); item (b) system unstated; unitarity 10⁻¹⁰ vs measured 1.1e-13; Hopf bundle written S³ → S² × S¹ (SYNTH NITs). | All fixed (S¹ → S³ → S²). |
| G1.8 | Paper 0 :835 "errors ranging … to 5.3 % for LiH" while BeH⁺ (≈31 % R_eq) is in the same list (CLAIMS-3). | BeH⁺ named with its Paper 17 numbers. |
| G1.9 | Paper 18 :1333 second 4/π-class locus left "general compact connected"; :1120 "propinquity asymptote" not covered by the widened C16 pattern (CLAIMS-3). | Fixed; C16 entry widened again (rates / asymptotes / constants), two-way proven 4/4 fire, 5/5 silent; :1088 dangling example renamed. |
| G1.10 | Paper 1 duplicate [SYMBOLIC PROOF] on the log-holonomy; :370 "encode the exact spectrum"; Paper 7 :83 "achieve O(V) scaling"; :748 "compact groups" vs "compact metric groups"; 18-proofs hedge not propagated to Paper 0 :680 / synthesis :356 (CLAIMS-3 NITs). | All fixed. |
| G1.11 | Paper 38 :968 "reaches 1.2732 = 4/π to four digits by n = 1600" — the doubling estimator reads 1.2785 there (three digits, residual shrinking ~1.8× per doubling); twin at :1757 (CLAIMS-2). | Both restated from the measured values. |
| G1.12 | Paper 40 :212 "full Class 1 generality" two lines above the semisimple theorem; :1980 "r positive-root directions" (it is N₊); residual "compact Lie group" at the `thm:main` title, :1906, :1994; the source comment still claimed a Latrémolière-propinquity upgrade (CLAIMS-2). | All fixed. |
| G1.13 | C16 `circle-fejer-constant-4-over-pi` missed the paraphrase "each circle factor carries the same constant 4/π" (CLAIMS-2). | Pattern widened (`circle factor carries the same`, `tori included`), two-way proven 3/3 fire, 4/4 silent. |
| G1.14 | Paper 38 `perez_sanchez2024` (*Bratteli networks*) is not the continuum-limit correction (CLAIMS-2 NIT). | Cite reverted to 2025 only; the 2024 bibitem removed (it was an orphan the F2.7 pass had "rescued" by citing it in the wrong place). |
| G1.15 | Paper 38 `rem:circle_fejer` under-tiered [MEASURED] though the closed form is derived; :284 2/π clause untagged (CLAIMS-2). | → [SYMBOLIC + MEASURED]; [OBSERVATION] added. |
| G1.16 | Dąbrowski–Dossena locator "eq. (10)" is their *classical* odd–even Cartesian product; the real-spectral-triple form is eq. (18) (𝒟̃), prescribed by their §4.2 for an odd first factor (CITE). The Table 5 directional claim was re-verified by direct reading of the primary (arXiv:1011.4456, Table 5 row 3₊: column 6₊ blank, column 6₋ = 1₋), so Paper 32's claim stands. | Locator fixed in Paper 32 (with row/column named) and in the `almost_commutative.py` docstring. |
| G1.17 | `standard_model_triple.py` comment "measured ε″ = +1" while the code's own residual gives ε″ = −1 (CODE-C). | Comment corrected. |
| G1.18 | `test_J_combined_KO3_sign` covered n_max = 2 only while Paper 32 says n_max ≤ 3; (−,+) also occurs at KO-2/KO-4 (CODE-C). | Parametrised over n_max ∈ {1,2,3}; oddness of the combined triple pinned (‖{γ_GV⊗γ_F, D}‖ ≠ 0), so KO-3 is the unique odd (−,+) row. |
| G1.19 | Paper 38's printed "2.008 at n = 800" and "0.63662 at n = 400" were not produced by the cited tests (CODE-C). | n = 800 ratio pinned < 0.01 (measured 2.008248); circle constant pinned to 1e-5 at n = 400 (measured 7.2e-7). |
| G1.20 | The forced-count module docstring still described the pre-correction paper; the new order-one projection test rebuilt a ~19 GB matrix in the default run and the lru_cache pinned it for the session (CODE-C). | STATUS header added; the order-one test slow-marked like its sibling; module-scoped fixture releases the cache. |
| G1.21 | Lorentzian axiom (vi)/(vii) tests covered N_t = 1 only while Paper 32 states the 3×3 panel {1,2,3}×{1,11,21} (CODE-C). | Parametrised over N_t ∈ {1,11,21} (residuals N_t-independent, as measured). |
| G1.22 | [MEASURED] used for bit-exact panel results (abstract, axiom-table caption, Peierls Chern number, JLO witness, Sprint L1) — the claims register reserves MEASURED for computed-vs-external-reference and PANEL-VERIFIED for bit-exact panels (CLAIMS-1). | All five → [PANEL-VERIFIED]. |
| G1.23 | Paper 42 :434–439 and :878–884 — two more "structurally distinguished … dim H₃ = 40 = g₃^Dirac" loci, plus "single-chirality" (CLAIMS-1). | Rewritten as the numerical coincidence (cumulative 4+12+24 vs both-sign single-level 40). |
| G1.24 | Paper 32: "disjoint" four lines above the F1.5 fix; "cannot derive them" left unhedged; abstract "Dirac graph operator" (×2) after the graph form was withdrawn; two "qualitative-rate" residues after F1.9; "thirty-one papers" (Papers 0–31 are thirty-two); 41+26+2 ≠ 70 arithmetic; "numerically established at ranks 2–3" vs P40's noisier rank-3 point; the compound Chamseddine–Connes bibitem behind a 1997 key (CLAIMS-1). | All fixed; key → `chamseddine_connes1997_2008`. |
| G1.25 | The KO-3 relabel paragraphs carried no tier though the signs are derived factor by factor and confirmed numerically (CLAIMS-1 upgrade). | [SYMBOLIC + MEASURED] at both loci; "measured signs" → "signs (derived … and confirmed numerically)". |
| G1.26 | Claim-matrix row 50 title still named the old chain (CLAIMS-1). | Updated to 2048→1024→512→272→260 (matter 512→256→128→128). |

**Method lessons carried into DELTA #2.** (1) Every rename/retag pass
leaves *enumerate-don't-sample residuals* — six reviewers found them
independently (Paper 18 ×2, Paper 42 ×2, Paper 40 ×4, Paper 32 ×3, the
syntheses ×3); the fix is a corpus-wide grep per class plus a C16 entry,
never a hand sweep. (2) The tier vocabulary is the claims register's:
MEASURED needs an external reference value; a bit-exact panel is
PANEL-VERIFIED. (3) A key re-year must check for an existing key first.
(4) A citation *locator* (equation/table number) must be read from the
primary, not from a memo of the earlier session.

**State after remediation.** Deterministic gates: C16 / C11 / C13 / C14 /
C19 / C21 PASS on `trunk` and `group1`; C18 / C15 PASS on `trunk`; C20, C22,
C5 (corpus-wide) PASS. C10: Papers 0/1/7/18/32/38/40/42 and both
syntheses PASS. Tests: 366 passed, 5 skipped, exit 0 in 76 s across the six re-run files (axiom audit 3x3 panel, KO-3 over n_max = 1-3 with the oddness pin, Fejer n = 800 ratio, forced-count, central Fejer, SM triple). Orphan bibitems: 0 in Papers 0/7/32/38
and the group3 synthesis. Seed worktree removed (seeds never touched the
corpus).

**DELTA #2 (owed).** Same shape; fresh S3 seed in the P0/P7 tests and
fresh S2 seed in the P32 tests (Opus); one seed per claims/synthesis
agent, two per Sonnet agent; scope = this remediation's diff. CLEAN-DELTA
there is the precondition for the FULL certifying run.

---

## Part H — DELTA run #2 (2026-09-02, PI-invoked `/qa delta #2`) = DEFECTS

**Shape.** DELTA-verification of the DELTA #1 remediation only: the DELTA #1
reviewed state was reconstructed from the saved diffs (seeds reversed) and
committed on a fresh seed branch, the current files synced over it, so the
worktree `git diff` was exactly the remediation (17 files, ~180 changed
lines). Criteria frozen as at `dbe7ae2` plus the PI-authorized code-tier
exception; all 18 deterministic gates PASS whole-target before dispatch
(C16/C11/C13/C14/C19/C21 on `trunk` and `group1`; C17/C18/C15 on `trunk`;
C20, C22, C5). Nine fresh seeds (answer key
`debug/qa/trunk_delta2_seed_key_2026_09_02.json`, placement asserted:
none on a comment line), eight controls, the same seven reviewers.

**Calibration.**

| Reviewer | Seeds | Caught | False positives | Status |
|:--|--:|--:|--:|:--|
| CODE-A (Opus) | 1 (S3 endpoint) | 1 | 0 | calibrated |
| CODE-C (Opus) | 1 (S2 self-bound) | 1 | 0 | calibrated |
| CLAIMS-1 (Opus) | 1 (S9 abstract tier) | 0 — **void seed** | 0 | **uncalibrated** |
| CLAIMS-2 (Opus) | 1 (C8 numeral) | 1 | 0 | calibrated |
| CLAIMS-3 (Opus) | 1 (S6 discrete produces) | 1 | 0 | calibrated |
| SYNTH (Opus) | 2 (S8 zombie, S9 status) | 2 | 0 | calibrated |
| CITE (Sonnet) | 2 (S1 year, S1 arXiv ID) | 2 | 0 | calibrated |

The CLAIMS-1 seed (abstract `[PANEL-VERIFIED]` → `[INTERNAL THEOREM]`) was
accepted as defensible because Paper 32's `prop:reality` supplies an
all-$n_{\max}$ construction — i.e. the plant was not a defect. Per the
void-seed rule the P32/P42 claims dimension is **INCONCLUSIVE** for this run
(its findings below were verified by the PM before acceptance); DELTA #3
must carry an unambiguous seed there (S4 κ-derived or S5 K-tier). Three
seeds were additionally noticed out of scope (S-E9 by SYNTH, S-E4 by
CODE-C, S-E8 by CLAIMS-1).

**Verdict: DEFECTS.** The DELTA #1 remediation carried genuine defects,
including two regressions of its own fixes; all remediated the same day
(table H1). DELTA #3 (fresh seeds) is the precondition for the FULL run.

### H1 — genuine findings (non-seed)

| # | finding (reviewer) | disposition |
|:--|:--|:--|
| H1.1 | The SU(2) doubling-estimator value printed as "1.2785 at n = 1600" is $a_{800}$ under the paper's own definition $a_n = (2n\gamma_{2n} - n\gamma_n)/\log 2$; $a_{1600} = 1.2761$ (measured: 1.30527 / 1.29096 / 1.28293 / 1.27849 / 1.27607 at n = 100…1600, residual shrinking 1.81–1.86× per doubling); no test reached n = 1600 (CLAIMS-2, CODE-C). | Both P38 loci restated with the index made explicit; `test_su2_doubling_estimator_at_n800` (default) and `_at_n1600` (slow, γ₃₂₀₀ ≈ 25 s) added; CHANGELOG v5.3.0 labels corrected. |
| H1.2 | Paper 40's body still names its own object "propinquity": `prop:limit_id_general` asserted "the propinquity limit … in the sense of Latrémolière", plus :114, :565, :578, :1586, :1822 (CLAIMS-2, out of hunk). | Proposition and five namings → state-space GH; Paper 40 added to the C16 propinquity entry's file list (gate PASS). The L5 section framing ("Latrémolière propinquity assembly", `eq:propinquity_def_general`, the convention subsection) is a **group1 carryforward** item, not fixed here. |
| H1.3 | Paper 40 header comment dropped the general-G conditionality; "Scope note after Corollary" ×3 (it is inside); "$r$ positive-root directions" (it is $N_+$); residual "compact Lie group" wording at the `thm:main` title, Reading A, the structural reading and Lemma L3 (CLAIMS-2). | All fixed. |
| H1.4 | Paper 40's **title** still names "compact Lie groups" while the theorem is now scoped to the semisimple class (CLAIMS-2, out of hunk). | **PI adjudication** (Part H2): a title change cascades to four internal bibitems, `papers/INDEX.md` and a Zenodo record already minted under the old title. |
| H1.5 | group3 synthesis :354 "universal across compact Lie groups" (CLAIMS-3, SYNTH); :361 identified Paper 7's 18 checks / the graph convergence with the GE–vS result (object conflation, SYNTH). | Semisimple + Paper 40's proof-status hedge; the GE–vS result stated for the Peter–Weyl spectral truncations, the graph-side convergence stated as having no theorem or frozen test. |
| H1.6 | group1 synthesis master-theorem sketch: the DELTA #1 wording "at least one semisimple … torus factor dominated" exceeds Paper 40 (whose product mechanism is additive and which asserts no value for $G_{ss}\times T^k$); the abstract twin :129 still said "arbitrary compact connected … universal 4/π"; :1560 called Paper 45's product carrier "independent justification of Paper 45's main theorem" (a withdrawn theorem) (SYNTH). | Class restricted to Paper 40's semisimple class at both loci; the Paper 45 sentence rewritten as the product-carrier proposition (`prop:product_action_seminorm`, not Lorentzian, joint rate = max of factor rates) outside the master theorem's class; "simple" → "semisimple" at :110 and :761. |
| H1.7 | Paper 0's new scope clause "among the hydrides through LiH" was broken by HeH⁺ (Paper 17: 93 % of $D_e$ adiabatic, ~5 % consistent 2D-variational) (CLAIMS-3). | HeH⁺ named with Paper 17's numbers; scope clause dropped. |
| H1.8 | Paper 7 :83 untagged and "near-linear cost" stronger than the frozen timing guard; :865 carried a changelog fragment (de-versioning) (CLAIMS-3). | [MEASURED] added; "operator size near-linear in $V$"; fragment removed. |
| H1.9 | `test_volume_element_jacobian` docstring still stated $d\Omega = \Omega^3 d^3p / p_0^3$ (wrong by $p_0^3$; $\int\Omega^3 d^3p = 2\pi^2$ exactly); O(V) docstring "< 1.5" vs the 1.8 guard; Rabi function docstring 0.95 / 0.5 % and "machine precision" (CODE-A). | All fixed. |
| H1.10 | The DELTA #1 "oddness pin" in `test_J_combined_KO3_sign` used `T.gamma_GV()` = sign($D_{GV}$), so $\{\gamma, D\} = 2|D_{GV}| \neq 0$ identically — it could not fail — and its inference is refuted: with the production grading `build_gamma_GV` ⊗ $1_F$ the finite combined triple is **even**, $\gamma^2 = 1$, $\{\gamma, D\} = 0$, $J\gamma = +\gamma J$ bit-exactly at $n_{\max} \le 3$, i.e. $(\varepsilon,\varepsilon',\varepsilon'') = (-,+,+)$, the **KO-4 column** (CODE-C; re-measured by the PM). | Test rewritten to pin the measured grading facts; Paper 32 and the module docstring now state that $(-,+)$ is shared by KO-2/3/4, that KO-3 is the continuum factor's label, and that the finite grading gives the KO-4 column — **PI adjudication** (Part H2): which label the finite truncation should carry. |
| H1.11 | `standard_model_triple.py`: the ε″ comment said −1 but the code measured the commutator (≈45); the "{γ, D} = 0 (odd)" comment was wrong twice (CODE-C). | Residual now the anticommutator (field ≈ 0, matches the comment); comments corrected. |
| H1.12 | Lorentzian (vi)/(vii) residuals are "sample-of-3" (3 of 55 multipliers); over the full basis at $n_{\max} = 3$ the maxima are 0.0785 / 0.2026 (2× the printed order-one figure) (CODE-C; re-measured by the PM). | Full-basis maxima stated in Paper 32 next to the sample-of-3 figures; detail table bound 0.101 → 0.1014; "5–10 %" → "5–10.2 %". |
| H1.13 | BBB sign-panel test looped $N_t \in \{1, 11\}$ vs the paper's $\{1,11,21\}$; circle pin 1e-5 admits 0.63663; $r_{800}$ window admitted 2.005–2.009; test name still said "monotone for n ≥ 3"; dead `paper_intermediates` line; `full_dirac_operator_system` docstring index convention (CODE-C NITs). | All fixed (21 added; 3e-6; $|r_{800} - 2.008248| < 5\cdot10^{-4}$; renamed; removed; convention noted). |
| H1.14 | Paper 32 quoted its own abstract without the new parenthetical; "explicit rate" without $+O(1/n_{\max})$; open-questions list still named the "full real structure $J$ audit" (settled) via a `debug/` memo; tier tag splitting a sprint name; "strongest possible level" superlative; caption vs "∀ n_max" rows; Paper 40 rank status under-stated the analytical mechanism (CLAIMS-1). | All fixed; Sprint L1 closure → [SYMBOLIC + MEASURED] per the C3-boundary rule. |
| H1.15 | Paper 32 still cited Perez-Sanchez 2024 (*Bratteli networks*) as part of "the correction" (CITE); a pre-existing orphan `perez_sanchez2024` bibitem in the group1 synthesis. | Two P32 cites → 2025 only (the 2024 lineage cites stay); orphan removed. |

**State after remediation.** Every deterministic gate PASS on `trunk` and
`group1` (C16 now scoping Paper 40); C10 PASS on Papers 0/7/32/38/40 and
both syntheses; TESTS_LINE; orphan bibitems 0 in Papers 32/38 and both
syntheses. Seed worktree removed.

### H2 — PI adjudications (nothing mechanical is blocked on these)

1. **Paper 40's title** ("…spectral truncations of compact Lie groups with
   bi-invariant metric") names a class the theorem no longer covers. Options:
   (a) retitle to "compact semisimple Lie groups" and cascade (P32/P38/group1
   synthesis bibitems, `papers/INDEX.md`; the minted Zenodo record keeps the
   old title, so the DOI landing page and the PDF would disagree until the
   next deposit); (b) keep the title as an umbrella and rely on the abstract's
   first sentence. PM recommendation: (a) at the next deposit, (b) until then.
   **PI decision 2026-09-02: (a), now.** ✅ DONE — retitled to "…compact
   semisimple Lie groups with bi-invariant metric"; 23 loci cascaded (P32,
   P38, P39, P41–P50, P52, P53, P55, field guide, both syntheses, two viz
   pages, `papers/INDEX.md`); C11 trunk + group1 PASS. The minted Zenodo
   record keeps the old title until the next deposit (PI-only).
2. **KO label of the finite combined triple.** Measured: $(\varepsilon,
   \varepsilon') = (-,+)$ (the v5.3.0 relabel), and with the production
   grading $\gamma_{GV} \otimes 1_F$ the finite triple is even with
   $\varepsilon'' = +1$ — the KO-4 column; the ungraded reading keeps the
   continuum factor's KO-3. Options: (a) keep KO-3 (continuum label) with the
   caveat now in Paper 32; (b) relabel KO-4 (finite-grading label); (c) drop
   the dimension label and print the sign triple only. PM recommendation: (c)
   — the signs are the invariant content; a dimension label at finite cutoff
   is a convention the paper should not appear to derive.
   **PI decision 2026-09-02: (c).** ✅ DONE — Paper 32 (4 loci), both
   modules, both test files (`test_J_combined_KO3_sign` →
   `test_J_combined_sign_pair`), matrix row 49, C16 entry notes; the
   measured triple $(-,+,+)$ is printed and no finite-cutoff label is
   attached. 112/112 tests pass in the three affected files.
3. **DELTA #3** — unseeded by default (PI direction 2026-09-02); if the PI
   wants the instrument, `/qa trunk seeded` with an unambiguous CLAIMS-1 seed.
   **PI decision 2026-09-02: skip DELTA #3; go straight to the FULL run,
   unseeded** (v5.4.0). The standing calibration record (run #2 19/21 seeds,
   0/8 FP; DELTA #1 9/9, 0/8; DELTA #2 8/9 with one void, 0/8) is what the
   verdict cites.

### H3 — group1 carryforward (spill, not fixed under trunk)

- Paper 40 §L5 "Latrémolière propinquity assembly" framing,
  `eq:propinquity_def_general`, the "Propinquity convention" subsection and
  the `\Lambda_{\mathrm{prop}}` notation all present the paper's object as
  the propinquity; the abstract, header and `prop:limit_id_general` now say
  state-space GH. A coherent rewrite of that section is group1's.
  ✅ **DONE 2026-09-02 (v5.4.0, under the PI's "do what you think is best").**
  Subsection retitled "Distance convention": $\Lambda_{\mathrm{prop}}$ is
  *defined* as van Suijlekom's state-space GH distance (bibitem `vs2021_jgp`
  added), Latrémolière's propinquity named as the strictly stronger object
  that is *not* claimed, dual reach as the named gap; §L5 retitled
  "state-space GH assembly via the approximation pair", the tunnel
  length kept only as orientation, lemma statement/proof/keywords/outlook
  reworded (18 loci). Only literature descriptions still say
  "propinquity" (a coadjoint-orbit result, the Lorentzian open problem).
  Paper 40 is still owed its own group1 cert.

---

## Sizing

- **Part A** — small, high leverage, do first. Three registry/screen edits with
  two-way discrimination proofs.
- **Part B** — the bulk of the trunk work. B4/B8/B10 are quick; B2 (11 loci +
  a definition), B7 (10 citations) and B5/B6 (each a judgement about downgrade
  vs. build) carry the weight.
- **Part C** — its own sub-sprint, and it re-opens group3.
- **Part D** — PI only.
- **Part E** — optional depth; recommend at least E1 and E3.
- **Part F** — run #2's scope. F1.3/F1.4/F1.5/F1.9/F1.12 and F2.1/F2.2 are
  quick and load-bearing; F1.1/F1.2 wait on F6; F2.1's propagation is a
  group1 item. The DELTA after Part F carries the F-seeding rule.

**Do not** re-run the FULL certifying pass until a DELTA comes back clean.
Remediated text is not clean text: one prior delta found 4 of its 11 genuine
findings were defects introduced by the previous run's own remediation.
