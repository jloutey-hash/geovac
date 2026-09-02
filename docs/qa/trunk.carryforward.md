# Trunk — carryforward: remediation scope after the 2026-09-01 FULL FAIL

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

## Sizing

- **Part A** — small, high leverage, do first. Three registry/screen edits with
  two-way discrimination proofs.
- **Part B** — the bulk of the trunk work. B4/B8/B10 are quick; B2 (11 loci +
  a definition), B7 (10 citations) and B5/B6 (each a judgement about downgrade
  vs. build) carry the weight.
- **Part C** — its own sub-sprint, and it re-opens group3.
- **Part D** — PI only.
- **Part E** — optional depth; recommend at least E1 and E3.

**Do not** re-run the FULL certifying pass until a DELTA comes back clean.
Remediated text is not clean text: one prior delta found 4 of its 11 genuine
findings were defects introduced by the previous run's own remediation.
