# Sprint: CLAUDE.md Compaction (2026-05-31)

## 1. Goal

Reduce CLAUDE.md context-load cost by compacting §2 sprint chronicles to one-liners and extracting §10/§11 reference tables to standalone files.

## 2. What was done

### Confinement reframing archived
- §2 banner ("OPEN REFRAME — READ EVERY SESSION") removed per PI direction
- §1.7 "Framing reorientation" block (2 long paragraphs) compacted to 1 short "Organizing observation" paragraph
- Memory file `confinement_reframing_open.md` updated to ARCHIVED status with rationale
- MEMORY.md index entry updated

### §2 compaction
- "Best results" prose converted to a 12-row table
- All sprint chronicle entries (κ, α, supertrace, nuclear, RH, QED, spectral triple, Lorentzian, gauge, chemistry, precision, gravity arcs) compacted from multi-paragraph chronicles to grouped one-liner summaries
- §2 reduced from 346 lines to 122 lines (65% reduction)

### §10 and §11 extraction
- §10 (Validation Benchmarks, 330 lines) extracted to `docs/validation_benchmarks.md`
- §11 (Topic-to-Paper Lookup, 404 lines) extracted to `docs/topic_to_paper_lookup.md`
- Both sections replaced with 6-line pointers in CLAUDE.md

### Sprint-close skill updated
- `/sprint-close` gained step 7 (reference table updates) pointing to the new standalone files

## 3. Quantitative result

| Component | Before | After | Saved |
|:----------|:------:|:-----:|:-----:|
| §2 | 346 lines | 122 lines | 224 (65%) |
| §10 | 330 lines | 6 lines | 324 (98%) |
| §11 | 404 lines | 6 lines | 398 (99%) |
| **Total CLAUDE.md** | **2207 lines** | **1263 lines** | **944 (43%)** |

## 4. Files modified

- `CLAUDE.md` — §1.7, §2, §10, §11 compacted
- `.claude/commands/sprint-close.md` — step 7 added
- `memory/confinement_reframing_open.md` — ARCHIVED
- `memory/MEMORY.md` — index entries updated

## 5. Files created

- `docs/validation_benchmarks.md` — extracted §10 table
- `docs/topic_to_paper_lookup.md` — extracted §11 table

## 6. Honest scope

- **Closed:** CLAUDE.md compaction is mechanical and complete. No information lost — all detail preserved in CHANGELOG.md, sprint memos, and the new standalone files.
- **Closed:** Confinement reframing archived as organizing reading per PI assessment (species-II predictions falsified, no predictive content generated).
- **Not applicable:** No theorems, numerical observations, or structural sketches — this was a documentation/infrastructure sprint.
- **Open follow-on:** §3 (Approaches That Failed, 69 lines) and §6 (Paper Series, 190 lines) could be further compacted but are lower priority and more genuinely useful in context.
