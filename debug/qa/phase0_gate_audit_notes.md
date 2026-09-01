# Phase 0 — gate scope audit + C21 widening (2026-08-31)

Pre-flight for the re-certification sweep (`docs/qa/recert_sweep_plan.md`,
Phase 0 items P0.1 and P0.2). No `/qa` was invoked; no LLM reviewer was
dispatched. No paper was edited.

---

## 1. What prompted it

The 2026-08-31 trunk run found C19 reporting

    RESULT: PASS (no eaten-escape corruption in 0 paper(s) in scope 'trunk')

`--gate` was a substring filter on the file path and no trunk paper's path
contains "trunk", so the gate had never examined a trunk document. Its own
notes recorded the open item: **seven of ten gates print no file count**, so
their coverage could not be confirmed from output for any target.

That makes "the deterministic layer is green" an unverified claim across the
whole sweep, not a fixed problem in one gate. Phase 0 tested it.

---

## 2. Findings

### F1 — C19 was one instance of a class; three more gates were mis-scoped

| Gate | Defect | Effect |
|:-----|:-------|:-------|
| **C10** compiles | substring `--gate` | `--gate trunk` compiled **0** documents; now 6 |
| **C18** duration | directory-based `BRANCH_DIRS` | `--gate trunk` scanned **27** papers (all of group1+group3) instead of 6; group runs silently omitted their synthesis; single-paper targets widened to the whole corpus while still printing the target's name |
| **C5** K-label | `--gate` narrowed a corpus-wide rule | `--gate group5` downgraded a §13.5 K-tier violation in *every other paper* to advisory and printed PASS. Same defect the 2026-08-22 run fixed on the default path, reintroduced through the flag. |
| **C15/C13/C14/C20** | substring `--gate` + a local duplicate TRUNK list | worked by luck for group names, examined nothing for `trunk` |

C5's fix was proven two-way: a planted K-rule violation in a group5 paper now
FAILs under `--gate trunk` (it did not before) and the gate returns to PASS
when the plant is removed.

### F2 — C16, the zombie gate, selected ZERO entries for five of eleven targets

C16 selects **registry entries** by a hand-maintained `scope` tag. No entry
carried the tag `trunk`, `synthesis`, `paper_58`, `paper_59` or `paper_60`, so:

| target | entries selected, tag-only (before) | after (locus-derived) |
|:-------|---:|---:|
| trunk | **0** / 27 | 7 |
| synthesis | **0** / 27 | 23 |
| paper_58 | **0** / 27 | 1 |
| paper_59 | **0** / 27 | 0 (legitimately — nothing retracted there) |
| paper_60 | **0** / 27 | 0 (same) |
| group3 | 6 | 7 |
| group4 | 3 | 4 |

So `/qa trunk` on 2026-08-31 recorded "C16 green" having checked nothing, and
the three single-paper certs (P58/59/60) and the synthesis cert were all issued
with C16 examining zero entries.

The sharpest instance: the `propinquity-as-achieved-metric` entry lists
`paper_38_*.tex` and `paper_32_*.tex` — both **trunk** papers, added by the
2026-08-31 trunk run itself — under the tag `group1`. **The trunk gate could
not see its own fix.**

**Root cause and fix.** The `scope` tag is a hand-maintained summary of the
`files` loci and had drifted from them. Selection is now **locus-derived**: an
entry runs whenever any of its declared loci lies in the gated scope, with the
tag kept as a widening fallback. Proven two-way with a planted zombie in
Paper 38 (fires under `--gate trunk`; silent when removed).

The same fix was applied to **C17** (identical `selected()` shape): families
selected rose trunk 2→3, group3 3→8, synthesis 2→19, paper_58 2→5.

### F3 — a severity gap C16 cannot close on its own

The registry entry covering Papers 38 and 32 is **advisory** severity, so even
now C16 cannot FAIL the trunk gate on the propinquity-as-achieved-metric class
— which is the exact overclaim trunk's C7 criterion names. Currently 0 live
hits and 2 correctly-flagged exempt ones, so promoting it to `fail` would be
safe today. **Left for PI decision:** it changes what a cert can fail on.

### F4 — C21's required-context guarantee was vacuous for one family

`65` (ERI block count) required context `r"ERI|nonzero quartet"`, searched
case-**insensitively** with no word boundaries — so it matched inside
`charactERIzes` and `inhERIt`, admitting any `65` in the corpus. Four false
positives (Papers 32, 51×2, 34), every one a percentage or ppm figure.
Fixed to `\bERI\b`; still fires on the genuine "65 at n_max = 2" ERI count in
Paper 8, silent on all four.

Also fixed: `1.69`/`1.690` admitted the deuteron-polarizability Layer-2 line in
meV (the required context `alpha` is satisfied by any nearby `\alpha` in a
precision passage) — forbidden context `meV|polarizability` added, proven to
stay silent there and still fire on the genuine `O(Q^{1.69})` in the field
guide.

And a **duplicate key**: `334` was defined twice in `RETIRED`; Python keeps the
last, so the earlier entry's forbidden-context guard was silently discarded — a
guard present in the source and absent at runtime.

### F5 — widening C21 found 30 genuine live retired values

C21 previously accepted only three targets (group3/4/6) and hard-exited on
every other name, so group1/2/5/synthesis/trunk and the single-paper targets
had **no numeric guarding at all**. Scope now resolves through `qa_scopes`.

After removing every false positive above, what remains is genuine:

| target | live retired | what |
|:-------|---:|:-----|
| **group2** | **24** | the pair-diagonal vintage, un-swept: LiH 334/878, BeH₂ 354.9, H₂O 778, He 120/2659, block 65, `O(Q^{3.15})`, `O(Q^{4.60})` — including a **whole resource table** in Paper 19 (§ lines 863–865) and a present-tense "854 Pauli terms (2.56× the composed value of 334)" |
| **synthesis** | 2 | the **field guide** twice asserts `O(Q^{1.69})` 1-norm (retired; 1.774) alongside `O(Q^{2.5})` |
| **paper_58** | 2 | "334 Pauli terms for LiH at Q = 30, O(Q^{2.5}) scaling" ×2 |
| **group3** | 1 | Paper 57 §forced_chem: "forces 1.44% ERI density", "N_Pauli = 11.10 × Q", "across 38 molecules" — three retired values in one sentence |
| **group4** | 1 | Paper 23: "measured O(Q^{3.15}) for Coulomb atoms" + "O(Q^{2.5}) composed" |

**This is the answer to the sweep plan's §2 question.** The Layer-2 exact-rule
re-pricing swept the loci C21 was watching and left the rest of the corpus
carrying the retired convention. group2 is therefore **not** a low-churn
"cause B" target as the plan classified it — it holds the largest single pocket
of Layer-2 debt in the corpus.

`tests/test_numeric_registry.py::test_gate_passes_on_current_corpus[group4]`
now FAILS. That failure is correct and is left standing as the honest record
until the value is remediated.

---

## 3. What was built

- **`debug/qa/qa_scopes.py`** — one declaration of all 11 target scopes, keyed
  to `docs/qa/<target>.done.md`. Scopes are declared by **paper number**, so a
  rename fails loudly instead of shrinking a scope silently. `describe()`
  always carries the file count. Self-test asserts every scope resolves
  completely, that none resolves to zero, that a missing declared paper warns,
  and that the predicate agrees with the resolver both directions.
- **`debug/qa/gate_coverage_matrix.py`** — runs every gate against every
  target and FAILs on any zero-coverage cell. Distinguishes a *bug* (zero
  because the scope resolved to nothing) from an *unexercised criterion* (zero
  because the registry legitimately holds nothing), and requires the latter to
  be reported as UNMEASURED rather than green.
- All 12 gates wired to the shared resolver; dead scope tables retired.

Current matrix: **PASS**, every gate examines a non-empty named scope on all
11 targets, with C16 on paper_59/paper_60 flagged UNEXERCISED.

---

## 3b. Remediation done this pass (PI-directed)

The **6 non-group2 loci are FIXED** — Paper 57, Paper 23, the field guide (×2)
and Paper 58 (×2). Each was rewritten with its surrounding claim re-read
rather than its numeral swapped, and that immediately paid: three of the six
passages carried retired content C21 had **not** flagged.

- **A seventh retired value, unregistered.** The field guide (×2) and Paper 58
  asserted a `190×` reduction against cc-pVDZ; canonical is **76×**
  (63,519 / 837). No registry entry existed, so no gate could ever have caught
  it — it surfaced only because the §15 rule forces re-reading the prose.
  Now registered, two-way discrimination proven (fires on the retired wording;
  silent on an unrelated count and on a journal page number).
- **A structural gap in the registry.** Registering it exposed that `resolve()`
  reached MEASURED and CITED only — so a RETIRED entry could not point at a
  DERIVED replacement, and **a ratio is inherently derived**. The whole class
  of retired ratios / exponents / per-qubit figures was therefore
  unregisterable, and C21 could never say what such a locus should read.
  `resolve()` now evaluates DERIVED, with cycle protection.
- **Two claims got weaker and now say so.** Paper 23's "reduction below the
  naive O(Q^4)" is 3.77 rather than 3.15 — still a reduction, by a narrow
  margin, now stated explicitly. Paper 58's composed-vs-Gaussian exponent
  comparison (3.17 vs 3.9–4.3) survives by much less than the retired figures
  implied, now stated.
- Papers 57/23/58 and the field guide compile clean with zero undefined
  references; C19 clean; 17/17 registry tests pass; C21 is **0 live** on every
  target except group2.

**group2's 24 are scoped as their own sprint:**
`docs/qa/group2_retired_value_sprint.md`.

---

## 3c. Severity promotion (PI direction, 2026-08-31)

`propinquity-as-achieved-metric` promoted **advisory → fail**. The trunk
criteria name that exact overclaim, so at advisory severity the gate could only
print a note about the one claim it most specifically guards. The original
"too noisy to gate" rationale was stale: the exemption list now covers the
legitimate framework / named-gap / descope / state-space mentions, and the
corpus showed 0 live hits against 2 correctly-exempted ones.

Verified: all 11 targets still PASS (no-op today); with the retracted wording
planted in Paper 38 the trunk gate now **exits 1 and blocks**, and returns to
exit 0 when removed. The module docstring, which used this entry as its
example of an advisory class, was corrected — with the general lesson recorded:
tighten the exemption, then gate; do not leave a named class permanently
advisory.

## 3d. P0.3 reassessed — the "142 numerals" task was built on a bad metric

The sweep plan inherited "burn down the 142 unregistered multi-document
numerals" from the Stage-4 memo. Widening C21 to the whole corpus turned 142
into **713**, and inspecting the list shows the number was never a worklist:
check D's only filters were `value >= 10` and `appears in more than one
document`. Across 62 papers that admits every year (**2026 leads the list with
~1,460 occurrences**), every paper and section number, every qubit count.

Registering those would be *worse than leaving them* — registry rule 3 says
never register a value you have not measured or cited, because a registry of
guesses launders them into authority.

Two fixes, so the number means something:

- **Measurement-shaped subset.** Decimals, thousands-separated counts, and
  integers large enough not to be structural labels; years excluded. Reported
  alongside the raw count, never instead of it.
- **arXiv/DOI stripping.** Spot-checking found `2401.04` listed as a quantity
  in three documents. It is `arXiv:2401.03705` — a modern arXiv ID
  (YYMM.NNNNN) parses as a decimal.
  **Measured effect: 27 of 728 entries, 3.7%.** The corpus holds 388 arXiv
  IDs, but an ID only reaches this list if it appears in **two or more**
  documents, so the great majority never enter it; what leaked was the small
  set of IDs cited across several papers.
  > *Correction (2026-08-31).* This bullet first read "every bibliography was
  > feeding the worklist." That was wrong — one verified instance generalised
  > into a corpus-wide magnitude without measuring it, which is the same error
  > this file documents in §F2's wrap-blindness retraction. The filter is
  > still worth keeping (it costs nothing and removes a known false-positive
  > class), but it is a 3.7% cleanup, not a structural finding. Note also that
  > the inline-arXiv citation gate is a *different axis*: it checks whether IDs
  > resolve and match their bibitems, and does not feed or filter this numeric
  > scan — so it neither caused nor could have prevented this leak.

Result: 713 → **686 raw / 520 measurement-shaped** corpus-wide, and per target
**16–78** rather than 92–186. The bulk of that reduction is the
measurement-shaped filter, not the arXiv strip. It is reviewable, but still an
uncurated list, not a defect list. **P0.3 is therefore NOT a Phase-0 blocker**
— it is per-target maintenance, best done inside each target's own cycle where
someone already has that paper open.

---

## 4. Open

1. **group2's 24 retired numerals** — scoped, not started
   (`docs/qa/group2_retired_value_sprint.md`). Prerequisite for `/qa group2`.
2. **P0.4, item 1 — balanced QWC columns: CLOSED, and it surfaced a bigger
   defect one column over.** The QWC item was already discharged on 2026-08-30;
   re-measuring all six second-row rows from the live builder confirms every
   printed Q / Pauli / QWC cell reproduces exactly
   (`debug/qa/verify_balanced_qwc.py`, 30/30 cells).
   **But the λ column beside it is computed at LiH's bond length for every
   molecule** — `build_balanced_hamiltonian` defaults `R = 3.015` and does not
   read the spec's geometry, so all six published λ values are at the wrong R
   (0.6–5.2% off). The backing test held the correct values and could not fire:
   it is `@pytest.mark.slow`, so 12 of its 13 cases skip by default, and the
   2026-08-30 pass moved the papers *away* from it while recording the move as
   a fix. Pauli and QWC are geometry-independent, so the blast radius is that
   two columns (QWC moves too — greedy grouping, so it is heuristic rather
   than structural). **RESOLVED and LANDED** under PI direction to take the
   mission-reinforcing option: per-molecule experimental geometry, because a λ
   projected through a fictitious shared R describes no physical system and so
   can carry no falsifier. Builder fixed, both tables corrected with an
   explicit R column, registry inverted (it had held the wrong-geometry values
   as canonical), and the forced/free split pinned by
   `tests/test_paper20_geometry_independence.py`:
   `debug/qa/balanced_lambda_geometry_finding.md`.
3. **P0.4, item 2 — angular-gradient assembly: DIAGNOSED, NOT FIXED, and
   NOT a certification blocker.** The defect is now located precisely rather
   than known-by-symptom: (i) the assembly mixes a complex-convention
   `_gaunt_integral` with real-harmonic reasoning that drops conjugates on
   both electrons, so their m-shifts add instead of cancelling — the L_z
   violation in one line; (ii) structurally worse, the bra orbital is
   delta-selected against `(L_eff, m_eff)` from an assumption the same block
   later disproves in its own comments, which forces `Q ≡ 0` and silently
   truncates the Neumann sum to a single term; (iii) a second electron-2 block
   repeats the construction and must be fixed with it. An attempted
   conjugate-restoring patch is recorded as a **dead end** — it makes the
   m-algebra close on paper but measures 26 → 46 added entries, all still
   violating, because the delta-selection sits upstream of it. A real fix is a
   re-derivation validated against independent quadrature, not a patch.
   **Does not block `/qa group4`:** the only live claim that depended on this
   operator is withdrawn and registered, and the paper's verdict is
   over-determined by the radial-only result.
   `debug/qa/tc_angular_gradient_diagnosis.md`.
2. **F3 severity gap** — PI decision.
3. `tests/test_numeric_registry.py` parameterises only over C21's legacy
   three-target `SCOPES`, so group2's 24 defects are not regression-protected.
   Widening it to all 11 targets needs a ratchet/baseline first (the C22
   pattern), or it simply goes red.
4. The plan's P0.3 (142 unregistered multi-document numerals) and P0.4 (the two
   live group4 measurement items) are untouched.
