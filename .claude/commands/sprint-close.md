---
description: End-of-sprint protocol — canonical memo, CHANGELOG, CLAUDE.md §2 one-liner, paper edits, MEMORY index
---

Close the sprint we just completed. Walk the protocol below. If a step does not apply, say so explicitly and skip it (do not silently omit).

**1. Canonical memo.** One memo per sprint, ≤ 5000 words, in `debug/sprint_<name>_memo.md`. If a memo already exists, extend it; do not create a second memo for the same sprint. If multiple sub-sprints landed in parallel, one memo per sub-sprint OR one umbrella memo with subsections (PI judgment).

**2. CHANGELOG.md entry.** Append under a new `## [vX.Y.Z]` heading with `### Added` / `### Changed` / `### Closed` subsections as appropriate. Full prose is fine here; CHANGELOG.md is the chronicle home (per CLAUDE.md §13.11).

**3. CLAUDE.md §2 one-liner.** Format:\ `**Sprint NAME (YYYY-MM-DD):** Verdict in one sentence. See debug/MEMO.md.` Strict ≤ 30 words. Do NOT compose a multi-thousand-word §2 bullet — that's the technical debt §13.11 is fighting.

**4. CLAUDE.md §3 dead-end row (only if a negative result landed).** Strict ≤ 2 sentences + memo link, in the failed-approaches table.

**5. Paper edits (only if structural findings worth keeping landed).** Edit papers in-place per §13.8. If a paper grew/shrank pages, note the diff (e.g.\ "Paper 50:\ 15 → 16 pages, three-pass clean"). Compile to verify three-pass clean before committing.

**6. MEMORY.md index entry (only if a cross-session-worthy fact landed).** Strict ≤ 200 chars in MEMORY.md, full content in a new `memory/*.md` file. Do NOT auto-create memory files for sprint outcomes — those go in CHANGELOG.md. Reserve memory for facts not derivable from the corpus.

**7. Reference tables.** If the sprint added new benchmarks, add rows to `docs/validation_benchmarks.md` (NOT to CLAUDE.md §10). If the sprint added new topic→paper mappings, add rows to `docs/topic_to_paper_lookup.md` (NOT to CLAUDE.md §11). These files are the canonical homes for those tables.

**8. Verification.** Have all relevant tests pass? Has every new equation in a paper edit got a verification test (§13.4a)? Has the PM not modified anything in the hard-prohibitions list (§13.5)?

**9. Honest scope check.** What was closed at theorem grade? What is structural sketch? What is numerical observation? What are the named open follow-ons? Write these into the memo §6 ("Honest scope") section explicitly.

After completing the protocol, summarize what you did and what is ready for `/release`.
