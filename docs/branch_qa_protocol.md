# Branch QA Review Protocol — full text

> Moved **verbatim** from CLAUDE.md §9 on 2026-09-01. Nothing rewritten.
>
> Safe to relocate because the section says so itself: its three QA principles
> are *"memorialized in the `code-reviewer` / `citation-reviewer` agent types
> and `docs/authoring_conventions.md`"*, and the cycle is executed by `/qa`
> against `docs/qa/criteria.md`.
>
> CLAUDE.md keeps the five-step shape, the **claim → artifact rule**, the
> dependency-order rule and the PM-vs-PI disposition split — the parts an agent
> must act on without following a pointer.

---

### Branch QA Review Protocol

The corpus is QA'd branch by branch (the dependency tree of §6 / the field guide). Each branch runs the same cycle, in order:

1. **Synthesis update** — an agent re-reads the branch's papers and proposes an update to the branch's group synthesis (bring current; add missing reconvergence; fix stale claims). PM applies.
2. **Adversarial paper review** — a parallel agent audits the branch's papers + synthesis for overclaim, §1.5 rhetoric, zombie citations (withdrawn/descoped results), claims-register consistency, status drift, cross-ref hygiene. PM reconciles and **verifies any agent conflict against primary text before acting** (see the Hodge-SL₂/layer-count precedents — agents disagree, primary text decides).
3. **Cycle if needed** — re-run 1–2 if the adversarial pass forces material synthesis changes.
4. **Adversarial code review, per paper** — one `code-reviewer` agent (`.claude/agents/code-reviewer.md`) per paper: map each load-bearing claim to its backing test + code, RUN the tests, and audit whether the test actually *proves* the claim (not tautological / false-positive / weaker than the prose); flag claims with NO test. Results populate `docs/claim_test_matrix.md`.
5. **Disposition** — PM **fixes small issues directly** (status drift, cross-ref hygiene, precision, missing caveats) and **raises large issues to the PI** (a load-bearing claim with no/weak/false-positive backing; a test that proves less than the prose; a suspected bug in a keystone result; anything touching a hard prohibition or a keystone's status). Checkpoint (`/checkpoint`, patch grade) when the branch is current.

**Claim → artifact rule.** Every load-bearing paper claim maps to a backing test, recorded in `docs/claim_test_matrix.md` (the granular companion to the reader-facing `docs/claims_register.md`). A claim with no backing test is a **coverage gap** — logged in the matrix, and raised to the PI if load-bearing — never a silent omission. New equations follow the §13.4a naming convention (`test_paper{N}_*`). Run the review in dependency order: the **trunk roots (Papers 0, 1, 7, 32, 38) before the branches**, so a finding at a root re-prices everything above it.

**Three QA principles** (apply to every branch; memorialized in the `code-reviewer` / `citation-reviewer` agent types and `docs/authoring_conventions.md`):

1. **Provenance visibility.** Every load-bearing claim wears its register tier *inline in the paper* (SYMBOLIC PROOF / INTERNAL THEOREM / MEASURED / PANEL-VERIFIED / CONDITIONAL / OBSERVATION / CONJECTURE), not only in `docs/claims_register.md`, and the prose may assert no more than the tier. κ is the canonical failure — "derived" survived ~50 versions because no tier was pinned to the prose; a matching / convergence test backs "matched / converges," never "derived."
2. **Fresh adversary.** The reviewer is a fresh agent that never saw the discovery-mode rationalization — prose confidence is not evidence; the claim is decided by the artifact + primary text. This is why the per-paper review is a dispatched sub-agent, not main-session self-review.
3. **Two-way verdict.** Verdicts move both directions: a claim *under*-stated by its backing is UPGRADED (4/π earned a genuine derivation, not only downgraded) — the PM reconciles against primary text, which guards against demolition bias. A pass that only ever cuts is miscalibrated.
