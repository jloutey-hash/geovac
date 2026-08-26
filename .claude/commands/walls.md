---
description: The /walls pass — curate and refine the corpus of negative results (the walls). Re-verify each against current state, classify hard vs soft, cluster into families, and mine shared mechanisms into structural findings. The negative-results sibling of /qa. PI-invoked. Usage: /walls <cluster|all|new>.
---

# /walls — the negative-results refinement gate

`/qa` refines the **positive** corpus (papers, claims): it hunts for a claim that asserts more than its backing supports. `/walls` refines the **negative** corpus — the walls: CLAUDE.md §3 (the failed-approaches table), the `debug/` dead-end memos, and the crystallized meta-findings in `memory/`.

**Why this exists.** In the skeleton picture, physics splits into a *forced* part (the discrete, exact skeleton) and a *free* part (the transcendental calibration data injected at projection). **The walls are the survey stakes on that border** — each negative result marks a place where the forced side stops. A pile of un-curated negatives is just a graveyard; *refined*, it is the map of the boundary. And this is not hypothetical: two of the project's best **positive** structural insights — `multi_focal_wall_pattern` ("discrete labels couple cleanly; spatial focal lengths don't compose") and `two_kinds_of_sparsity` ("symmetry sparsity dies at two centers; distance sparsity GeoVac never had") — were *born* from clustering negatives until a shared mechanism fell out. `/walls` is the machine that does that on purpose.

> **PI-invoked.** Like `/qa`, this is a deliberate, PI-timed curation pass — **not** a per-sprint reflex, and the PM does not self-trigger it. Division of labour: `/sprint-close` **files** a new wall (appends a §3 row); `/walls` **refines** the accumulated corpus. (Genuine-trigger command, no corresponding memory rule — CLAUDE.md §13.9a.)

**Home.** `docs/walls/register.md` — the living register. CLAUDE.md §3 is **append-only** (§13.5 hard prohibition: negatives are never deleted or modified). The register is where re-verification, hard/soft tags, and clustering live; it *points back* to §3 rows and memos, and never rewrites them.

## Inputs
- **Scope** = `$ARGUMENTS` — one of:
  - a **cluster name** (refine one family, e.g. `/walls chem-accuracy`),
  - `new` (only walls filed since the register's `Last run` date),
  - `all` (the whole register — large; prefer running by cluster).
  - If absent, ask the PI for a scope, and show the cluster index from the register.
- **The register** `docs/walls/register.md`. If a wall in scope has no entry yet, create one.
- **Current corpus state** — for every wall touched: the **owning paper section** + **CHANGELOG since the wall's date**. A §3 row is a *dated snapshot*, not current truth (feedback_verify_current_state).

## The four moves — per wall, then per cluster

### 1. RE-VERIFY — is the wall still standing?
For each wall in scope, read its owning paper section + the CHANGELOG since its date, and assign exactly one status **with the evidence that decides it** (no bare verdicts):

- **STANDING** — still holds. Cite the current-state evidence checked (paper §, CHANGELOG range). *A "STANDING" with no cited check is a rubber stamp, not a verdict — reject it.*
- **SOFTENED** — holds in general, but a scope was carved out. Name the carve-out + its artifact. *(Canonical example: "the cheap electron-cusp treatment is impossible without losing sparsity" → SOFTENED — breached cheaply for **atoms** by xTC, Paper 14 §tc_atomic_sparsity, v5.0.9.)*
- **BREACHED** — a later result crossed it. Cite the result. Per §3-append-only, a breach is **annotated** (register + a new §3 row that references the original), **never** a deletion of the original row.
- **MIS-SCOPED** — the wall is real but its stated scope is wrong (too broad or too narrow). State the corrected scope.
- **SUPERSEDED** — folded into a sharper statement (a crystallized meta-finding). Point to it.

**Discipline (feedback_validate_before_reducing).** Claiming SOFTENED / BREACHED / MIS-SCOPED requires **exhibiting** the specific breach — a computed result, a paper section — never a hunch. Downgrading a wall's authority on assertion is the same failure `/qa`'s calibration guards against on the positive side.

### 2. CLASSIFY — hard or soft?
Tag each standing wall by *what kind of barrier it is*, because that decides strategy:

- **HARD** — a proven structural impossibility: a theorem, a category / rigidity obstruction, a conservation / parity argument, a cost-conservation (dual-basis) result. **Keepers**, and often *positive* impossibility results worth writing up as-is. *(e.g. symmetry-sparsity-dies-at-2-centers — SO(4)→axial→point-group, N3b m-rule-only; composition-wall — ‖[P_A,P_B]‖=0.50, dual-basis theorem.)*
- **SOFT** — engineering / convention / precision / basis-limited: our attempt wasn't good enough, or it's a numerical-precision ceiling, or a convention artifact. **Revisit candidates.** *(e.g. a PSLQ negative at a given precision budget; the Möbius α>1 form that turned out to be a fixed-a substrate artifact.)*
- **OPEN-LEANING** — tested-negative, but a *named* better attempt exists. **Frontier.** *(e.g. cheap-TC γ-selection: flexible Jastrow tested → caps ~8 mHa, but imported general-TC machinery may do better.)*

### 3. CLUSTER + MINE — what does a family of walls share?
Group the scope's walls by structural cause. For each cluster, attempt the **crystallization move** (the one that produced `multi_focal_wall_pattern` and `two_kinds_of_sparsity`): state the *single unifying mechanism* in one sentence, list the **independent** negatives that support it, and give an **operational consequence** (an A/B/C dispatch rule). Track cluster maturity:

- **CRYSTALLIZED** — a clean shared-mechanism statement exists (points to its meta-finding memo / paper §). Re-verify it still covers its members and that no member has drifted out.
- **FORMING** — ≥2 independent negatives visibly share a cause, but the one-sentence statement isn't sharp yet. **This is the working queue** — the single highest-value output of a `/walls` run is moving a cluster FORMING → CRYSTALLIZED.
- **SINGLETON** — a lone wall with no family yet.

### 4. PROMOTE (PI-gated) — surface earned crystallizations
When a FORMING cluster reaches a clean statement **+** independent support **+** a falsifier, surface it as a structural-claim candidate **for the PI**. The empirical bar the `multi_focal_wall_pattern` used: *"five independent observables converging is past the point of coincidence"* — treat **≥3 independent negatives** as the floor. Promotion into a paper claim or a WH is **PI-gated** (same governance as WH promotion, §1.7): `/walls` proposes, the PI decides. `/walls` **never** autonomously turns a negative into a positive paper claim, and **never** presents a soft wall as hard (or vice versa) without the exhibited evidence from moves 1–2.

## Output to the PI
- **Register delta** — which walls changed status this run, each with its exhibited evidence.
- **Cluster map** — each cluster's maturity, members, and shared-mechanism statement where one exists.
- **Promotion candidates** — any FORMING cluster that reached the bar: the one-sentence mechanism, its independent support, its falsifier — flagged for PI decision.
- **Strategic read** — which SOFT / OPEN-LEANING walls are the best **revisit** candidates and why; which HARD walls are **paper-worthy** impossibility results not yet written up.

## Hard rules
- **PI-invoked.** Not a per-sprint reflex; the PM never self-triggers it. `/sprint-close` files new walls; `/walls` refines the corpus.
- **§3 is append-only (§13.5).** `/walls` never deletes or rewrites a §3 row. Breaches / re-scopes are annotated in the register and, if load-bearing, appended as a **new** §3 row referencing the original. The historical negative is preserved intact.
- **No bare verdicts.** Every status carries the current-state evidence that decides it (paper § + CHANGELOG range). Every downgrade exhibits the specific breach (feedback_verify_current_state + feedback_validate_before_reducing).
- **Promotion is PI-gated.** `/walls` proposes crystallizations; turning one into a paper claim or WH needs PI direction.
- **A fresh adversary sharpens re-verification.** For a high-stakes or long-standing wall, the PI may have a fresh agent re-check it **without** being told the prior verdict (the `/qa` fresh-adversary principle). The failure mode this guards: the corpus assuming a wall is dead *forever* and never re-checking — exactly how an atomic-xTC-style breach hides in plain sight.
