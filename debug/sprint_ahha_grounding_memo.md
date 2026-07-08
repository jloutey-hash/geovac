# Sprint: grounding the /ahha skill (process fix)

**Date:** 2026-07-08
**Type:** process / tooling. No physics findings, no production code.
**Trigger:** an /ahha generative pass produced a connection pitched as novel that
the corpus had already resolved a month earlier. The PI named the blind spot and
directed the fix.

## 1. The incident (what exposed the blind spot)

An /ahha pass on the v4.73.x composition-wall work surfaced three connections. The
PI picked **Idea 1** to pursue: *"the chemical bond is a Connes–Chamseddine inner
fluctuation over a 2-point finite geometry — the two atomic centers are the two
sheets, the off-diagonal overlap block is the Higgs-analog connection, binding is
the spectral-action double-well."* The pass framed this as unifying two arcs that
"run in totally separate universes" (chemistry composition ↔ the WH1 spectral
triple).

On investigation (mandatory current-state check before concluding), that framing
was **stale**. Both legs of Idea 1 were already settled on 2026-06-07:

- **Structural leg — already PROVEN.** Sprint M-vS-2 (`debug/sprint_mvs2_lih_default_plus_rsweep_memo.md`)
  showed the production LiH `h1` **is** a Marcolli–van Suijlekom gauge network on
  the 2-vertex Li↔H bond quiver, bit-exact (residual 0.0). The edge intertwiner
  L_e (‖L_e‖_F = 1.378) IS the bond-as-connection object. "Two sheets + a bridge"
  is the literal content of that result — not a new conjecture.
- **Dynamical leg — already DEAD, structurally.** Sprint M-vS-2 Q2 found the
  spectral action S(D)(R) = Tr exp(−D²/Λ²) **monotone** in R at every Λ — no
  interior minimum, no binding. The companion Spectral-Action-Expansion sprint
  (`debug/sprint_spectral_action_expansion_chemistry_diagnostic_memo.md`) explained
  *why*: for a finite-dim Dirac there is no UV divergence, so the Chamseddine–Connes
  Seeley–DeWitt hierarchy (the term that would build the Higgs double-well) **cannot
  arise**. The CC binding mechanism is a continuum phenomenon; the n_max cutoff
  removes it. Independently consistent with W1e (binding is lost at the projection
  step, not the evaluation step).

So the correct disposition of Idea 1 was reachable from the corpus alone:
**structural leg already true; dynamical leg structurally dead at finite cutoff;
the only live residue is the continuum (M×F, continuum atomic S³ factors) version,
which is a research program, not a diagnostic.** The /ahha pass reached this only
*after* the PI sent the PM to check — the skill itself did not ground.

**One genuine byproduct (kept):** Idea 1 and the pass's Idea 3 (compactness ⇒
discreteness) are the *same statement* — the bond is a Higgs-analog structurally,
but the binding-as-symmetry-breaking dynamics is a continuum property the finite
compact skeleton cannot host. This is a clean re-reading, not a new result; it is
recorded here, not promoted.

## 2. The blind spot (root cause)

`/ahha` had a **generate → hedge-audit → defend** flow. It filtered for *false
caution* (trained flatness) but had **no step that checked candidates against the
current corpus** — no filter for *false novelty*. Connections were generated from
the model's mental picture of the corpus, which is a dated snapshot, and pitched as
un-surfaced without verification.

The standing rule `feedback_verify_current_state` (2026-06-14) exists precisely for
this, but it was scoped to "before forming a **verdict**." A generative pass *feels*
like generation, not conclusion, so it slipped under the rule's trigger. That
exemption is the root cause.

## 3. The fix (three homes, one discipline)

1. **`.claude/commands/ahha.md`** — inserted two steps around the existing
   generate/hedge core (not a rewrite):
   - **Step 0 (Ground):** load the owning paper section + the relevant group
     synthesis + §2/§3 for the domain *before* generating. Scoped to the domain.
   - **Step 2 (Verify):** for each candidate, search the related papers/CHANGELOG/
     memos and classify ALREADY / PARTIALLY / UN-SURFACED, then **re-aim** (not
     drop). Explicit guard: the verify step re-aims, it does not suppress — a
     partially-surfaced connection re-aimed to its open leg is often the best
     output (Idea 1 is the worked example). Old Steps 2/3 → 3/4.
   - Division of labor made explicit: **Step 2 kills false novelty; Step 3 kills
     false caution.**
2. **`memory/feedback_verify_current_state.md`** — widened trigger from
   "concluding" to "concluding **or proposing a connection/direction from corpus
   state (a generative pass such as /ahha)**"; added "your own mental model of the
   corpus" alongside the debug/ memo as an untrustworthy snapshot; appended the
   2026-07-08 incident to **Why**.
3. **CLAUDE.md §9 Current-State Check** — same widening in the canonical workflow
   doc; MEMORY.md index line updated to match.

**Walk-through proving the fix catches the incident:** under the new flow, Idea 1 →
Step 0 loads Paper 32 + the M-vS synthesis + §3 → Step 2 searches the M-vS sprints
→ classifies PARTIALLY SURFACED (structural proven, dynamical dead) → re-aims to
the continuum residue, citing the done legs instead of re-conjecturing them. That
is the disposition reached this session — reached first, not after a manual nudge.

## 4. Files touched

- `.claude/commands/ahha.md` — ground-then-verify rewrite (Steps 0–4).
- `memory/feedback_verify_current_state.md` — trigger widened to generative passes.
- `memory/MEMORY.md` — index line updated.
- `CLAUDE.md` — §9 Current-State Check widened; §1 version; §2 one-liner.

No `geovac/` production code, no `tests/`, no paper edits.

## 5. Verification

- No production code changed → `/regression` not applicable (nothing in the diff
  maps to a consumer test; the topological-integrity baseline is unaffected by
  doc/skill/memory edits).
- Hard-prohibitions (§13.5) clean: no geometry-hierarchy change, no fitted
  parameter, no negative-result deletion, no Paper 2 combination-rule change.
  CLAUDE.md edits confined to §1 (version only) / §2 (one-liner) / §9 (Workflow —
  PM-editable per the §13.5 access table). No edits to any NO-tier section.

## 6. Honest scope

- **Theorem grade:** none. This is a process sprint.
- **Structural (process):** the /ahha blind spot is characterized precisely
  (generation was exempt from the verify-current-state rule) and closed in all
  three canonical homes.
- **Numerical / physics:** none new. The Idea 1 investigation **re-confirmed**
  existing §3 results (M-vS-2 structural leg; spectral-action-expansion dynamical
  leg) via the current-state check — it produced no new physics and no new §3 row.
- **Named open follow-on (NOT launched):** the continuum bond-as-Higgs path —
  treat each atom as its full continuum S³ spectral triple (WH1/Paper 38) and the
  bond as the finite factor, i.e. a genuine M×F almost-commutative geometry, and
  ask whether the spectral-action double-well (binding) appears there where the
  n_max=2 truncation forbids it. This is the only place the dynamical leg could
  still live. Research program, not a diagnostic; needs explicit PI direction to
  scope.
