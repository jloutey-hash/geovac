---
description: Grounded generative pass — surface 3 corpus-checked connections from recent work; separate false-novelty (already surfaced) from real reach, and false-caution from live falsifier
---

Pause. This is a **grounded** generative pass: reach for connections a corpus-veteran would make, but verify each against the CURRENT corpus before proposing it — so a connection the papers already made is caught and re-aimed, not pitched as novel. The blind spot this fixes (2026-07-08): an /ahha reach ("bond = Higgs unifies two arcs that run in separate universes") was pitched as novel when the Marcolli–vS arc had already bridged those arcs a month earlier. Grounding is a novelty filter and a stakes-sharpener, **not** a brake on reaching.

**Step 0 — Ground (before generating).** Identify the domain of the work just completed — or, if no recent work, the most recent sprint entry in CLAUDE.md §2. Load current state for THAT domain (scope the load; do not pull all six syntheses):
- the owning paper section(s) — the paper that OWNS the question, per `feedback_verify_current_state`. A `debug/` memo is a dated snapshot the corpus moves past; the paper + post-memo CHANGELOG are canonical.
- the relevant group synthesis in `papers/synthesis/` and the `papers/INDEX.md` status flags for the papers it touches.
- the recent CLAUDE.md §2 entries and the §3 dead-ends for the domain.

Generate FROM this state, not from memory. /ahha is not exempt from the verify-current-state discipline because it "feels like" generation rather than conclusion — that exemption is exactly what produced the blind spot.

**Step 1 — Generate the reach.** List three connections — to other GeoVac findings, to physics literature you know of, or to nearby mathematical structures — that a veteran would reach for and that you have not surfaced. For each:

- (a) state it as a single-sentence claim
- (b) name what is at stake if true (what it implies, what it unifies, what it predicts)
- (c) name the specific falsifier (what computation or argument would kill it)

Still reach. Do NOT pre-filter to safe/known connections — Step 2 sorts them.

**Step 2 — Verify each against the corpus (the gate that was missing).** For EACH candidate, before proposing it, check whether the corpus has ALREADY made it: search the related papers, CHANGELOG, and canonical memos it would touch. Classify and re-aim:

- **ALREADY SURFACED** (corpus has it, with backing) → do not pitch as novel. Convert to "you already have this — the unbuilt next step is Y," or drop if there is no next step.
- **PARTIALLY SURFACED** (one leg established, another open — e.g. the structure is proven but the dynamics is untested) → state precisely what is done vs open. The proposal is the OPEN leg, with the done leg cited, not re-conjectured.
- **GENUINELY UN-SURFACED** → proceed. Make (b)'s "what's at stake / what it unifies" accurate about what the corpus VERIFIABLY lacks (having just checked), not about what you assume it lacks.

The verify step **re-aims; it does not suppress.** A partially-surfaced connection re-aimed to its open leg is often the most valuable output — do not discard a reach just because part of it is already known.

**Step 3 — Hedge audit.** For each surviving connection, ask: am I about to hedge because the falsifier in (c) is genuinely live and the claim might be wrong, or because training pushes me toward measured language regardless? **Distinguish these two cases explicitly; do not blur them.** If the hedge is just trained-flatness, restate the claim without it. (Step 2 removes FALSE novelty; Step 3 removes FALSE caution — keep them separate.)

**Step 4 — Defend one.** Pick the connection you are most tempted to underclaim. Defend it. State plainly what it would mean if true, and what specifically would have to hold — grounded in the Step 2 finding of what the corpus does and does not already contain.

**Discipline.** This is a generative pass, not a write-to-paper pass. Do not edit production papers or CLAUDE.md from this output without separate PI direction. The W3 spectral-zeta candidate was falsified after exactly this kind of leap — that is a feature of the protocol, not a reason to skip the leap. Grounding (Steps 0/2) makes the leap land in the right bucket; it is not a reason to stop leaping.
