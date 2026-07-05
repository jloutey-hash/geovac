# Sprint memo — the layer-discipline arc (rules 11–13, C18, corpus sweep)

**Dates:** 2026-07-04 → 2026-07-05 · **Releases:** v4.67.8, v4.67.9, v4.68.0 · **PI-directed throughout.**

## 1. What happened

Reading the field guide, the PI surfaced three distinct register defects in the synthesis layer, each of which turned out to be a *class*, not an instance:

1. **Withdrawal-chronicle language** — process narrative ("was retracted in June 2026", "this synthesis has been brought into line with the retraction") living in the reader-facing synthesis instead of the papers of record.
2. **Fabricated/unreliable durations** — "Three years ago the project produced transcendentals…" in a project a fraction of that age; more broadly, wall-clock units attached to the project's course that LLM drafting cannot reliably produce.
3. **Insider terminology** — "Multi-year wall", "no sprint-scale handle", bare internal codenames presuming the reader holds the project chronicle.

Each became: a pass over the synthesis layer → a codified rule → (where grep-able) a deterministic gate → a corpus-wide sweep → calibrated seeded verification. The arc also absorbed two adjudications the verification loops surfaced (CvS-deferral attribution; Q1′/Q2′ label drift) and a working-tree tidy (PDF-tracking policy, Notes.bib artifacts, recovered clobbered data).

## 2. The three rules (canonical homes)

| Rule | Content | Homes |
|:--|:--|:--|
| **11 — withdrawal placement** | Chronicle (dates, process narrative, event-framing) lives with the object of record (paper Status notes, claims register, CHANGELOG); syntheses carry current state + the mathematical reason ("descoped by Paper 45's degeneracy theorem") + at most one-clause headline-grade flags. Claims are *withdrawn/descoped*; "retraction" reserved for whole-paper retraction (never occurred). Zombie rule unchanged: withdrawn-asserted-live = MATERIAL. | `authoring_conventions.md` r11; `criteria.md` NIT class |
| **12 — no project-course durations** | No wall-clock units (years/months/weeks/days) on the project's own course, historical or forward-estimate. Fine: sequence language, real dates, version anchors, unit-free effort vocabulary (sprint-scale / beyond sprint scale / bounded / substantial program / long-range / deep frontier). Exempt: external-world history, physics usage, compute-cost facts. | r12; **C18** gate; `memory/feedback_no_duration_language.md` |
| **13 — audience register** | Reader-facing prose readable without the internal chronicle: bare verdict-metaphors translated or defined at first use (group2's "a guardrail in the project's sense:" is the model); sprint/track/WH codenames kept **only** as parenthetical chronicle anchors under a per-document conventions note (field guide "Reading conventions" ¶; first-use footnotes elsewhere); formalized concepts defined at first use. | r13; `criteria.md` NIT class; `memory/feedback_audience_register.md` |

Synthesis DoD criterion 5 pre-registers all three for `/qa synthesis`; deterministic set for that target is now C10–C18.

## 3. The C18 gate (build → harden → promote)

`debug/qa/check_duration_language.py`, self-tested (16 positives / 11 negatives), scans all `papers/group*` + `papers/synthesis` `.tex`. Pattern classes: duration-ago, units-of-project-work, over-the-past-unit, across-n-units, compressed-into-units, project-took-units, multi-unit-estimate, unit-scale-estimate, unit-long, numeric-duration, wordnum-duration.

Two blind spots were found *by the verification loops themselves* and are now permanent patterns:

- **Line-join semantics** (v4.67.9): "Over several years / of investigation" spanned a LaTeX source line break and was invisible to line-based scanning. Fix validated against the pre-fix HEAD text (catches the exact line).
- **Word-number adjectivals** (v4.68.0 sweep delta): "a focused two-day session" — digit-based patterns miss "two-day". Found by a chunk reviewer as an out-of-delta observation.

Promotion history: built 2026-07-04 with historical classes FAIL + effort-estimate vocabulary ADVISORY (~190-instance standing debt, the C14-debug-refs precedent); promoted to **all-FAIL 2026-07-05** when the sweep retired the debt. Final state: **0 fail / 0 advisory over 59 papers.**

## 4. The corpus sweep (v4.68.0)

~230 instances retired across 27 documents: 190 advisory-tier + 38 numeric forms surfaced by the extended patterns + the delta-found wordnum residual. Uniform vocabulary mapping: multi-year → long-range; multi-month → substantial/extended; multi-week+ → beyond sprint scale; explicit ranges ("$\sim 6$--$12$ months") → unit-free scope statements preserving priority labels and technical parameters. First-use chronicle-anchor footnotes rode along where papers used sprint/track codenames unqualified.

**Execution note (process lesson):** 5 parallel sweep agents died at the monthly spend limit after clearing ~64 instances with no final reports. Recovery cost ≈ zero because the work-list was deterministic — remaining C18 hits = remaining work; the partial edits were verified by diff review (correct class only; no numbers/math/citations touched) and the remainder finished in main session. *Design principle: give fan-outs a deterministic, re-derivable work-list so interruption never requires transcript archaeology.*

## 5. Verification record (three seeded deltas, one /qa invocation)

| Run | Scope | Panel | Sensitivity | Specificity | Genuine MATERIAL |
|:--|:--|:--|:--|:--|:--|
| Retraction-pass micro-delta (7-04) | rule-11 diff, 2 files | 1 Opus claims + 1 Sonnet citation | 4/4 (S9 via sharpened recovery re-dispatch) | 4 controls clean | **2** — CvS-deferral misattribution (P44 pins the 2021 deferral as the *Riemannian* GH question, closed by P38; open item = its Lorentzian analog); 4/π universality prose>tier (abstract + field-guide twin) — both fixed |
| `/qa synthesis` post-cert delta (7-05, PI-invoked) | register-pass diff + criterion-5 conformance | 1 Opus claims | 2/2 (S9 first-pass — the sharpening baked in) | 4 controls clean | 0 — cert stands, rules 11–13 verified |
| Sweep delta (7-05, PI-invoked, pre-release) | all 267 hunks, 4 chunks | 4 Opus claims | 5/5 tailored seeds (false-closure, scope inversion, open-frontier zombie, content loss, scope downgrade) — each caught via cross-locus corroboration | 5 controls clean incl. Q1′/Q2′ verified vs Paper 49 | 0; 4 NITs fixed (doubled article, "continuum continuum", "two-day session", P55 vocabulary harmonized) |

Arc totals: **11/11 seeds, 13/13 controls, zero false positives.** Deterministic layer green throughout (C5-K/C10–C18; C16's 46 withdrawal flags survived the compression and the replace_all edits).

## 6. Honest scope

- **Theorem-grade:** nothing — this is a documentation-integrity / QA-infrastructure arc, not a physics sprint. No equation added, no claim tier changed upward.
- **Faithful-record corrections (verified against papers of record):** the CvS-deferral rewording (P38/P44 primary text); the Q1′/Q2′ adjudication (Paper 49 explicitly: closes Q1′, does *not* close Q2′ — "partially addresses" via the OSLPLS containment; cocycle-deficit closed form a separate open item). The 2026-05-31 memory ("Q2′ CLOSED") overstated the paper and was corrected — a verify-current-state incident spanning three documents and a memory.
- **Deterministically verified:** C18 0/0 over 59 papers; all touched papers three-pass clean; zero seed leakage in all three worktrees.
- **Judgment-verified (calibrated):** only the *changed loci* of each pass — unchanged certified surfaces rest on the prior certs + the C16/C17/C18 registries, standard delta semantics.
- **Named follow-ons:** (i) none structural; (ii) *optional standing item*: the q2_q2prime incident suggests an occasional memory-vs-corpus audit pass (memories are point-in-time; papers move); (iii) the next Zenodo DOI snapshot will carry compiled PDFs for the first time under the new tracking policy — worth a glance at the archive size when it lands; (iv) rule 12/13 conformance for *future* prose is carried by the memory rules + C18; the register rule (13) has no deterministic gate — it remains claims-reviewer judgment.

## 7. Artifacts

- Rules: `docs/authoring_conventions.md` r11–r13; `docs/qa/criteria.md` (3 NIT classes + C18 + change-log); `.claude/commands/qa.md` step-1 C18 entry; `docs/qa/synthesis.done.md` criterion 5.
- Gate: `debug/qa/check_duration_language.py` (self-tested).
- Memory: `feedback_no_duration_language.md`, `feedback_audience_register.md` (new); `q2_q2prime_closed.md` (corrected).
- Chronicle: CHANGELOG v4.67.8 / v4.67.9 / v4.68.0 (incl. all three delta records); seed keys `debug/qa/{synthdelta_retraction_pass,synthesis_delta3,sweepdelta}_seed_key.json` (gitignored).
- Repo tidy: PDF-tracking policy note in `.gitignore`; 15 Notes.bib untracked; 2 clobbered `debug/data` production analyses (H₂O/LiH Pauli) restored from HEAD.
