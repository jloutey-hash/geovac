# GeoVac Research Agent System

## Architecture Overview

Four agents, one human decision-maker. You (the PI) remain the approval gate for all research directions. The agents expand your field of view — they don't replace your judgment.

```
                    ┌─────────────┐
                    │  PI (You)   │
                    │  Approves / │
                    │  Redirects  │
                    └──────┬──────┘
                           │
                    ┌──────▼──────┐
                    │   LEADER    │  ← Synthesizes project state,
                    │  Strategic  │    proposes ranked directions
                    │   Brief     │
                    └──────┬──────┘
                           │
              ┌────────────┼────────────┐
              ▼            ▼            ▼
        ┌──────────┐ ┌──────────┐ ┌──────────┐
        │ EXPLORER │ │DECOMPOSER│ │ REVIEWER  │
        │ Find new │ │ Algebraic│ │ Critique  │
        │geometries│ │separation│ │ papers    │
        └──────────┘ └──────────┘ └──────────┘
```

## The Research Loop

### Step 1: Strategic Brief (Leader Agent)

**When:** At the start of a research sprint, or when you feel stuck, or after completing a track.

**How:** Load `agents/LEADER.md` into context along with `CLAUDE.md` and `SCOPE_BOUNDARY.md`. Say:

> "Generate a Strategic Brief for GeoVac. Here are the papers and project docs."

The Leader will produce a brief with 3-5 ranked research directions, each with:
- What the direction is and why it matters
- Connection to GeoVac's existing framework
- What a positive result looks like vs a negative result
- Estimated effort (1 sprint vs multi-sprint)
- Risk assessment (novel vs incremental, likely to hit known dead ends?)

### Step 2: You Decide

Read the brief. Pick a direction, modify one, or reject all and ask the Leader to reconsider with additional constraints. The Leader should explain its reasoning clearly enough that you can evaluate it even if the underlying math is unfamiliar — that's its job.

### Step 3: Execute (Explorer or Decomposer)

Depending on the direction:

- **"Find a geometry for X"** → Explorer Agent. Searches literature for published geometric structures, filters through GeoVac's natural geometry principle, returns candidates with references.

- **"Make this computation algebraic"** or **"Where do transcendentals enter X?"** → Decomposer Agent. Takes a mathematical expression or code path and systematically catalogs continuous/transcendental content, proposes algebraic replacements, checks against the exchange constant taxonomy.

- **"Implement and test X"** → Standard Claude Code workflow with CLAUDE.md loaded.

### Step 4: Review (Reviewer Agent)

After a paper draft or significant result:

> "Review this paper draft against GeoVac standards."

The Reviewer checks internal consistency, rhetoric compliance, benchmark honesty, transcendental cataloging, and mathematical correctness. Returns a structured critique you can act on.

### Step 5: Return to Leader

After completing a track (positive or negative), feed the results back to the Leader for the next Strategic Brief. The Leader incorporates what was learned and adjusts its recommendations.

---

## Practical Usage in Claude Code

Each agent is a markdown file in `agents/`. To invoke one:

1. **In Claude Code:** Reference the agent file in your prompt:
   ```
   Read agents/LEADER.md and CLAUDE.md. Generate a Strategic Brief.
   ```

2. **In Claude.ai:** Upload the agent file + CLAUDE.md + relevant papers as context.

3. **Combining agents:** You can chain them in a single session:
   ```
   Read agents/LEADER.md. Generate a brief.
   [You pick direction #2]
   Now read agents/EXPLORER.md. Execute direction #2.
   ```

## What This System Cannot Do

- **Generate correct physics autonomously.** Every result must be validated against known benchmarks, symbolic proofs, or analytical limits. The agents help you explore — they don't certify truth.
- **Replace domain intuition for dead-end detection.** The Leader can flag directions that overlap with known failed approaches (Section 3 of CLAUDE.md), but novel dead ends require human judgment.
- **Guarantee literature completeness.** The Explorer searches what's findable; niche or unpublished results may be missed.

## What This System Does Well

- **Reduces cognitive load.** You don't have to hold 24 papers and 30 tracks in your head to decide what's next.
- **Prevents reinventing failures.** The Leader cross-references against documented dead ends before recommending directions.
- **Enforces consistency.** The Reviewer catches rhetoric violations, missing benchmarks, and unacknowledged limitations.
- **Surfaces connections.** The Leader and Decomposer can spot structural parallels between different parts of the framework that might not be obvious.
