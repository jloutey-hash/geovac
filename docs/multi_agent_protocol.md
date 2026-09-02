# Multi-Agent Protocol — reference sections

> Moved **verbatim** from CLAUDE.md §13 on 2026-09-01. Nothing rewritten.
>
> These four subsections are *reference* material that restates files already
> in the repo (`agents/*.md`, `agents/WORKFLOW.md`, `.claude/commands/*.md`),
> so they were costing bytes in a file loaded into every session and every
> sub-agent dispatch without being something an agent must hold in working
> memory.
>
> **What deliberately stayed in CLAUDE.md §13:** 13.5 hard prohibitions, 13.11
> content discipline, 13.9 the mandatory session-summary format, 13.8 paper
> update policy, and 13.2/13.3/13.4/13.6/13.7/13.9a — the things an agent must
> act on *without* following a pointer.

---

### 13.1 Architecture

Four layers:

1. **Research agents (strategic layer):** Four specialized agents in `agents/` that formalize the strategic planning process. The PI invokes these by loading the agent file + CLAUDE.md into a Claude session. See `agents/WORKFLOW.md` for the full protocol.

   | Agent | File | Role |
   |:------|:-----|:-----|
   | Leader | `agents/LEADER.md` | Synthesizes project state, proposes ranked research directions with reasoning. Reads CLAUDE.md + SCOPE_BOUNDARY.md + recent results. |
   | Explorer | `agents/EXPLORER.md` | Searches published literature for geometric structures that fit the natural geometry principle. Returns structured candidate reports. |
   | Decomposer | `agents/DECOMPOSER.md` | Takes mathematical expressions/code paths, catalogs transcendental content using the exchange constant taxonomy (Paper 18), proposes algebraic separations. |
   | Reviewer | `agents/REVIEWER.md` | Critiques paper drafts against GeoVac standards (rhetoric rule, benchmarking rule, transcendental cataloging) and general scientific rigor. |

   The research agents operate on the three guiding principles:
   - **Natural geometry search:** If in doubt, search for new geometries to incorporate (Explorer).
   - **Algebraic deconstruction:** When code runs long, approach algebraically. Diagonalize. Deconstruct the continuum (Decomposer).
   - **Transcendental cataloging:** Catalogue where transcendental quantities enter and their relationship to geometries (Decomposer).

   The PI remains the approval gate: agents propose, the PI decides.

2. **Plan mode (human + agents → sprint plan):** The PI reviews the Leader's Strategic Brief, selects a direction, and defines the sprint plan (track definitions, PM prompts, exit criteria). This layer produces the directive that the PM agent executes. CLAUDE.md changes originate here.

3. **PM agent (main Claude Code session, full context):** Reads CLAUDE.md, README, and all papers relevant to the current track. **Default: does the work directly in main session — coding, computation, paper edits, commits.** Sub-agents are reached for only when the PI explicitly directs dispatch OR when the work is genuinely parallelizable AND context-heavy enough that main-session would clog. Evaluates sub-agent results against verification checklists when they are used.

4. **Worker sub-agents (scoped context, opt-in):** Available for cases where parallel dispatch genuinely beats sequential main-session work (e.g. 4+ independent computational tracks, each context-heavy). Default is to NOT use them. Each sub-agent reloads CLAUDE.md and its task-relevant files, which is expensive; reserve for cases where the cost is justified by the parallelism or context-protection benefit.

**Sub-agent cost discipline (2026-05-26 policy update).** Earlier sessions defaulted to sub-agent dispatch for every sprint. The cost (CLAUDE.md re-loaded per dispatch, parallel multiplier on context budget) has been substantial, and many sprints have been done equally well or better in main session. Default is now flipped:\ main-session by default, sub-agent on PI direction. The research agents (Leader/Explorer/Decomposer/Reviewer) in `agents/` remain available for big strategic moments; they are demoted from default to opt-in.

**The research loop:**
```
Leader Brief (when invoked) → PI picks direction → PM works in main session →
PI directs sub-agent dispatch IF needed → Reviewer critiques papers (when invoked) →
Results feed back to next Leader Brief
```

### 13.4a Equation Verification Protocol

Every equation that appears in a GeoVac paper must have a corresponding numerical verification in the codebase. This protocol bridges the research agents (which propose mathematical structures) and the PM/worker pipeline (which implements and tests them).

**The rule:** No equation goes into a paper without a test that verifies it computationally. "Verified" means one of:

| Verification type | What it checks | Example |
|:------------------|:---------------|:--------|
| **Analytical limit** | Equation reduces to a known result in a limiting case | H₂ → 2×H as R→∞ |
| **Symbolic identity** | Equation matches an equivalent form term-by-term | Gaunt integral = product of 3j symbols (symbolic, exact) |
| **Numerical cross-check** | Equation's output matches an independent computation | GeoVac He energy vs NIST reference |
| **Dimensional consistency** | All terms have consistent units | Energy terms in Ha, not mixed Ha/eV |
| **Symmetry verification** | Claimed symmetries hold numerically | Hermiticity of Hamiltonian matrix (H = H†) |

**Implementation:**

1. **When a sub-agent derives a new equation:** The PM must dispatch a verification sub-agent that implements the equation in code and checks it against at least one of the verification types above. The test goes into `tests/` alongside the production code.

2. **When the Decomposer proposes an algebraic separation:** The PM must verify that the separated form reproduces the original to machine precision (< 1e-12 relative error). Both the original and separated forms must be implemented and compared.

3. **When the Reviewer flags an unverified equation:** The PM treats this as a blocking issue — the equation is removed from the paper or a verification track is opened before the paper is finalized.

4. **Test naming convention:** Equation verification tests are named `test_paper{N}_eq{M}` or `test_paper{N}_{description}` so they can be traced back to the specific paper claim they verify.

**What counts as sufficient verification depends on the claim:**

- **Exact identities** (e.g., "the graph eigenvalue is n²-1"): must be verified symbolically or to machine precision across all relevant quantum numbers up to n_max=5.
- **Scaling laws** (e.g., "Pauli count scales as O(Q^3.8)"): must be verified across at least 4 data points with a log-log fit. Report the fit residual. **The example is a cautionary one:** this exponent read O(Q^2.5) for months because the fit was clean and nobody re-derived the evaluator underneath it. A log-log fit certifies the *shape* of a claim, never the *values* it is fitted to -- price the inputs by an independent route as well ([[feedback_independent_route_crosscheck]]).
- **Accuracy claims** (e.g., "0.004% error"): must report the computed value, the reference value, the reference source, and the actual percentage to enough digits to confirm the claim.
- **Negative results** (e.g., "this approach diverges"): must show the divergence numerically (e.g., error increasing monotonically with basis size across at least 3 points).

### 13.9b Project slash commands

Eight slash commands defined in `.claude/commands/`:

| Command | Type | Purpose | When to fire |
|:--------|:-----|:--------|:-------------|
| `/aha` | Trigger | Two-phase generative pass: reach wild (cross-domain Outlier + Inversion, judgment off), then filter hard (ground/verify/hedge/defend). | PI-chosen moment after a substantive result lands. Forces genuine lateral reach with the guardrails quarantined to a second phase. |
| `/sprint-close` | Trigger | End-of-sprint protocol:\ canonical memo, CHANGELOG entry, CLAUDE.md §2 one-liner, optional §3 row and paper edits, optional MEMORY index entry, verification check, honest-scope check. | When a sprint completes and is ready for release. Standardises the close pattern. |
| `/checkpoint` | Trigger | Version bump + commit + tag, with precondition checks (version string bumped, §2 entry exists, papers compile clean, tests pass, hard-prohibition check, untracked-file callout). **Push is opt-in (`/checkpoint push`) and never targets `main`.** Does NOT create a GitHub Release — those are manual, mint a Zenodo DOI, and are PI-only. | After `/sprint-close` has staged the content. Renamed from `/release` 2026-08-26. |
| `/qa` | Trigger (**PI-invoked ONLY**) | The QA certification gate: pre-registered criteria (`docs/qa/<target>.done.md`) + per-run calibration controls (seeded defects from `docs/qa/seed_defects.md`) + independent fresh reviewers (`code-reviewer` / `claims-reviewer` / `citation-reviewer`) → a three-way **PASS / FAIL / INCONCLUSIVE** verdict that distinguishes "target not done" from "reviewer not trustworthy." Skill: `.claude/commands/qa.md`. | When the PI judges a target (trunk / a branch) ready to *certify* as done. **The PM must NEVER self-trigger it, run it proactively, or nudge toward it each sprint** — control of timing is the PI's (qa.md hard rule). |
| `/walls` | Trigger (**PI-invoked**) | The negative-results refinement gate — sibling of `/qa` (which refines positive claims). Re-verifies each wall against current state (STANDING/SOFTENED/BREACHED/MIS-SCOPED), classifies HARD/SOFT/OPEN-LEANING, clusters walls into families, and mines shared mechanisms into structural findings; earned crystallizations are surfaced for PI promotion. Register: `docs/walls/register.md`. Skill: `.claude/commands/walls.md`. | When the PI wants the accumulated corpus of walls (§3 dead-ends + `debug/`/`memory/` memos) curated. PM never self-triggers; `/sprint-close` *files* a new wall, `/walls` *refines* the set. |
| `/audit-claim` | **Force-fire backup** for [[feedback_audit_numerical_claims]] | Curve-fit-audit pattern from `docs/curve_fit_audit_memo.md`. The standing rule should already fire on any "X matches Y" claim; this command is a backup if it didn't. | When the PI notices the PM made a numerical-coincidence claim without running the audit. |
| `/diag` | **Force-fire backup** for [[feedback_diagnostic_before_engineering]] | Diagnostic-before-engineering pass. The standing rule should fire automatically when ≥ 2 honest negatives accumulate; this command is a backup. | When the PI notices the PM is about to launch an implementation sprint past the 2-negatives threshold without running the diagnostic. |
| `/transcendental-tag` | **Force-fire backup** for [[feedback_tag_transcendentals]] | Paper 18 + master Mellin engine + Paper 34 classification. The standing rule should fire automatically on any transcendental appearance; this command is a backup. | When the PI notices a transcendental was introduced into production code or a paper without classification. |

The first five are genuine triggers (no corresponding memory rule, fired by PI at a specific moment). The last three are force-fire backups (their primary mechanism is a memory rule that should fire automatically; the slash command exists to manually invoke the discipline when the rule failed to trigger).

**Adding more.** Memory rule:\ write `memory/feedback_<name>.md` with the standard frontmatter + body, and add an index line to `MEMORY.md` (≤ 200 chars). Slash command:\ write `.claude/commands/<name>.md` with a single `description:` frontmatter line and the prompt body. No registration, no schema, no config.

### 13.10 Research Agent Integration

**When to invoke research agents:**

| Situation | Agent | What to provide |
|:----------|:------|:----------------|
| Starting a new research phase or feeling stuck | Leader | CLAUDE.md + SCOPE_BOUNDARY.md + any recent results |
| "What geometric structure handles X?" | Explorer | The specific question + relevant papers |
| "Where do transcendentals enter this computation?" | Decomposer | The expression/code path + Paper 18 |
| Paper draft ready for critique | Reviewer | The draft + CLAUDE.md + prior papers it builds on |
| After completing a track (positive or negative) | Leader | Updated CLAUDE.md with new results |

**Handoff from research agents to PM:**

When a research agent produces a result (e.g., the Explorer identifies a candidate geometry, or the Decomposer proposes an algebraic separation), the PI translates it into a sprint plan for the PM using the standard track format:

```
Track [XX]: [Agent output → implementation goal]
  Agent source: [which agent, what it proposed]
  PM prompt: [specific implementation task]
  Exit criteria: [what constitutes success/failure]
  Papers affected: [which papers would be updated]
  Equation verification: [which new equations need tests]
```

**Research agents do not write code or modify files.** They produce analysis, proposals, and critiques that the PI evaluates and the PM/worker pipeline implements. This separation is intentional: the agents explore, the PM executes, the tests verify.
