# Re-certification sweep plan

**Written 2026-08-31.** Companion to `check_cert_staleness.py`, which is the
authoritative *state*; this file is the proposed *order*. Every `/qa`
invocation named here is PI-fired — nothing in this plan authorises the PM to
self-trigger the gate.

---

## 1. The situation

All 11 targets read CERTIFIED in their `.done.md` and all 11 are stale. That
uniform verdict hides three different causes, and they want different
treatment:

| Cause | Targets | What actually happened |
|:------|:--------|:-----------------------|
| **A — Layer-2 re-pricing** | group3, group4, group6, trunk | The exact-rule ERI correction, the numeric registry, and the Stage-4 cross-document ledger moved published numbers. The changed surface is large and the changes are *coupled across group boundaries*. |
| **B — paper arcs outran the group** | paper_58, paper_59, paper_60, group2 | Single papers advanced hard (P59 alone: 1160 lines) after their group was certified. P59 carries an explicit open carryforward. |
| **C — slow drift** | group1, group5, synthesis | Small edits accumulated against a distant cert date. group5's whole diff is 131 lines. |

Cause A is the expensive one, and the reason a naive target-by-target sweep
would be wasteful: one defect class introduced by the exact-rule correction
lives in group3, group4 *and* group6, so three independent panels would each
pay full price to find the same thing.

### Changed surface, measured

| Target | Certified | .tex changed | Lines changed | Cause | Papers w/ tests |
|:-------|:----------|-------------:|--------------:|:-----:|----------------:|
| trunk | 2026-06-16 | 5 | 387 | A | 4 of 6 |
| group3 | 2026-08-29 | 11 | 953 | A | 7 of 11 |
| group4 | 2026-08-30 | 4 | 2251 | A | 4 of 4 |
| group6 | 2026-08-30 | 5 | 2490 | A | 4 of 4 |
| paper_58 | 2026-08-17 | 1 | 177 | B | 1 |
| paper_59 | 2026-08-21 | 1 | 1160 | B | 1 |
| paper_60 | 2026-08-18 | 1 | 700 | B | 1 |
| group2 | 2026-06-28 | 7 | 664 excl. P58/59/60 | B | 7 of 12 |
| group1 | 2026-06-24 | 13 | 721 | C | 10 of 16 |
| group5 | 2026-07-04 | 4 | 131 | C | 7 of 8 |
| synthesis | 2026-07-05 | 4 | 985 | C | n/a (code dim waived) |

"Lines changed" double-counts deliberately: group2's 321-line synthesis diff
is also the synthesis target's, and P58/59/60's 2037 lines are also group2's.
That overlap is the argument for the ordering in section 3 — sequence so each
defect is found once, by the cheapest detector that can see it.

---

> **UPDATE 2026-08-31, after Phase 0 ran** (`debug/qa/phase0_gate_audit_notes.md`).
> Two corrections to the table above:
>
> - **group2 is a cause-A target, not cause-B.** Widening C21 beyond
>   group3/4/6 found **24 live retired pair-diagonal values** in group2 —
>   the largest pocket of Layer-2 debt in the corpus, including a whole
>   resource table in Paper 19. It should move up the order, not sit in
>   Phase 4.
> - **The deterministic layer was not green anywhere.** Four gates were
>   mis-scoped and C16 selected *zero* registry entries for five of the
>   eleven targets, so several past certs recorded gates as green that had
>   examined nothing. Fixed; `debug/qa/gate_coverage_matrix.py` is now the
>   standing proof.

## 2. The stability precondition (read before firing anything)

A certification is a statement about a **frozen** surface. If Layer-2 work
continues landing paper edits on group3/4/6 while their panels run, those
certs are stale the day they are issued — which is exactly the state the
sweep is trying to exit.

So the first PI decision is not "which target first" but **"is the
paper-facing half of the Layer-2 arc done?"** Three items say not quite:

1. **The angular-gradient assembly bug is unfixed** (Stage-4 open follow-on
   2). It is disclosed and fenced, and Paper 14's resource multipliers were
   *withdrawn* rather than re-measured. Any re-measurement lands back on
   Paper 14 and Paper 20 — i.e. inside group4's scope.
2. **The balanced QWC columns** print a "pending re-measurement" note that is
   false; the builder returns the values directly. Currently being measured.
   Lands on group4.
3. **142 unregistered multi-document numerals** (101 group4, 24 group3).
   Advisory today, but each is an unguarded twin/derivation of exactly the
   class C21 exists to catch.

Recommendation: close 1–3 as a pre-flight sprint (Phase 0), then sweep. The
alternative — sweeping now and re-certifying group4 afterwards — pays for
group4 twice.

---

## 3. Ordering principle

Two rules, both already project doctrine:

- **Roots before branches** (section 9 of CLAUDE.md): a finding at a root
  re-prices everything above it. Trunk holds Papers 0, 1, 7, 32, 38 and the
  group3 synthesis.
- **Cheapest exhaustive detector first**: the deterministic layer runs
  whole-corpus in seconds and most cross-paper drift dies there for ~0
  tokens. It goes before every LLM dispatch, every time.

And one rule specific to this sweep:

- **Cross-cutting before per-group.** Cause-A defects do not respect group
  boundaries. Sweep the coupled numeric surface once, corpus-wide, via the
  registry — then let the per-group panels review prose against numbers that
  have already stopped moving.

---

## 4. Phase 0 — pre-flight (no LLM dispatch, no `/qa` invocation)

This is the highest-leverage phase and none of it costs panel tokens.

**P0.1 — Finish the gate scope audit.** The 2026-08-31 trunk run found C19
reporting PASS having examined **zero** trunk papers: `--gate` is a substring
filter on the path, and no trunk paper's path contains "trunk". Seven of ten
gates still do not print a file count, so their coverage cannot be confirmed
from output for *any* target. Until this is fixed, every "deterministic layer
green" in this sweep is an unverified claim, and a gate that scopes its
verdict away is indistinguishable from a working one.
→ Make the file count and the resolved scope mandatory in every gate's RESULT
line; add the `[scope] WARNING` for unmatched scope members to all ten. Then
re-run all gates against all 11 targets and record the counts.

**P0.2 — Extend C21 scope beyond group3/4/6.** Every other group's numerals
are unguarded by check C (Stage-4 open follow-on 3). group1, group2, group5,
synthesis and trunk all get panels in this sweep; sending them in without
numeric guarding means paying an LLM to find registry-class defects.

**P0.3 — Burn down the 142 unregistered multi-document numerals**, group4
first (101 of them). Registering a numeral is cheap; having a panel discover
a twin drift is not.

**P0.4 — Close the two live group4 measurement items** (angular-gradient
assembly, balanced QWC columns) so group4's surface stops moving.

**P0.5 — Re-run `check_cert_staleness.py --detail`** after P0.1–P0.4 and
re-baseline the table above. The trunk scope gap it had (omitting the group3
synthesis) is fixed; assume the other scopes have not been audited.

Exit criterion for Phase 0: all ten gates report a non-zero, *named* file
count for all 11 targets, and C21 covers every target that will be swept.

---

## 5. Phases 1–5 — the sweep

Each target follows the same three beats unless noted. All 11 are
**re**-certs, so none is a first cert:

> **DELTA run** (diff-scoped since last cert, seeded, one reviewer per
> affected dimension, changed loci pasted in) → **remediate** → repeat until
> **CLEAN-DELTA** → **FULL certifying run** (all dimensions, enumeration
> forced, completeness-critic) → PASS.

Only the FULL run can produce PASS. Budget roughly 0.5M tokens per delta and
~2.5M per full run, from the group4 arc's measured cost.

### Phase 1 — trunk

Already half-done: the deterministic layer ran 2026-08-31 (12/12 PASS with
scopes stated, two real gate defects fixed). **The four LLM dimensions were
never dispatched**, so the verdict is INCONCLUSIVE by construction.

Small surface (387 lines, 6 documents, 4 with tests) and it is the root.
Fire the LLM panel here first; a finding in Paper 0/7/32 re-prices group1 and
group3 downstream.

### Phase 2 — the Layer-2 cluster, in dependency order

**group3 → group4 → group6.**

group3 is the foundations root for the other two: Paper 22's
potential-independence density, Paper 18's taxonomy and Paper 31's partition
are what group4's resource numbers and group6's precision numbers are stated
against. Certify the definitions before the numbers that use them.

group4 and group6 carry the largest changed surfaces in the corpus (2251 and
2490 lines) and are the direct blast radius of the exact-rule correction. Both
were certified 2026-08-30 and both went stale within a day — expect FAIL on
the first delta and plan for two remediation cycles each.

Note that group4's `carryforward.md` still describes CF-1 as an open
decision. It is not: the PI resolved it 2026-08-29 (exact rule everywhere,
re-price honestly) and CF-1 is DISSOLVED — the pair-diagonal "convention" was
a wrong-sign-q bug. Mark the carryforward closed during this phase so no
reviewer reopens it.

### Phase 3 — the single-paper targets

**paper_59 → paper_58 → paper_60.**

P59 first, because its carryforward names the exact next action: a fresh
DELTA run to verify the weight-3/disc-8 fallback. That is a cheap, well-scoped
run. P58 (177 lines) and P60 (700 lines) are small.

Doing these before group2 matters: all three papers live in group2's folder,
so certifying them individually shrinks group2's remaining surface to 664
lines across four documents.

### Phase 4 — the low-churn targets

**group5 → group2 → group1.**

group5 is the cheapest target in the corpus (131 lines, 4 documents) — a good
panel warm-up between the expensive Phase-2 runs.

group2 comes after Phase 3 has removed P58/59/60. group1 is 13 changed
documents but mostly small edits, and its 359-line Paper 32 diff overlaps
trunk, so Phase 1 will already have reviewed that document.

### Phase 5 — synthesis, last

The synthesis target (field guide + 6 group syntheses) is the roll-up. Its C9
faithfulness is measured *against* the papers, so it can only be certified
after the six groups are. Four of its documents changed (985 lines), and three
of those diffs are edits the group runs will already have reviewed — which is
fine and expected: the synthesis question is not "is this number right" but
"does this summary still faithfully represent that paper".

---

## 6. What the PI fires, in order

Phase 0 is PM work and needs no `/qa`. After that:

```
/qa trunk          (LLM panel; deterministic already green)
/qa group3
/qa group4
/qa group6
/qa paper 59       (delta first — carryforward names the scope)
/qa paper 58
/qa paper 60
/qa group5
/qa group2
/qa group1
/qa synthesis
```

Each line is 2–3 invocations, not one (delta → remediate → full). Eleven
targets at ~3 runs is ~30 invocations; at measured cost that is the dominant
expense of the coalescing effort, which is the argument for spending Phase 0
generously — every defect the registry and the gates catch is one a panel
does not have to.

---

## 7. Expected verdicts, stated in advance

Pre-registering the expectation makes the result informative:

- **Likely FAIL on first delta:** group4, group6 (largest re-priced surface,
  certified one day before going stale), group3 (the M9-class cross-document
  zombies were found *in* this surface).
- **Likely clean-delta:** group5 (131 lines), paper_58 (177 lines).
- **Unknown:** trunk (never had its LLM panel run this cycle), synthesis
  (985 lines, but roll-up prose rather than numbers).
- **Known-open going in:** paper_59 (carryforward fallback awaiting its delta).

If group5 or paper_58 comes back FAIL, that is information about the panel or
about drift we have not modelled — worth pausing the sweep to understand.

---

## 8. Honest scope of this plan

This plan orders the sweep and names the pre-flight. It does not:

- decide whether the Layer-2 arc's paper-facing work is finished (section 2 —
  PI call);
- re-open any frozen DoD (`docs/qa/<target>.done.md` criteria stay frozen; if
  a target's criteria need to change, that is co-written with the PI *before*
  its run, per the gate protocol);
- estimate wall-clock anything.
