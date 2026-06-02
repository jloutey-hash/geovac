# Sprint: Confidence-Review Infrastructure — Memo

**Date:** 2026-06-01
**Type:** Meta / tooling sprint (not a physics sprint).

Built a reusable pre-broadcast paper-audit capability, calibrated it, ran a
partial foundations audit, then **deliberately wiped the audit discoveries (PI
direction)** so a fresh model can re-audit cold. This memo documents the
*durable* outcomes — instruments, methodology lessons, verified paper fixes —
NOT the wiped findings.

## 1. What this sprint was

A conversational session that reframed the question from "is GeoVac a theory of
everything?" toward "build justified confidence before broadcasting to science
communities." It produced a reusable confidence-review capability and applied
two verified corrections to Paper 2. The substantive audit findings (a partial
foundations audit of Papers 7/0/1) were intentionally cleared at PI direction to
enable an unbiased replication by a different model.

## 2. Instruments built (DURABLE, kept) — `agents/`

- **`CONFIDENCE_AUDITOR.md`** — external skeptical referee for pre-broadcast
  confidence. Design: (i) *invert the prior* — internal confidence labels
  (PROVEN, first-in-literature, STRONG, the WH register) are claims under audit,
  never evidence; (ii) every verdict grounds in an external result, a computation
  the agent ran, or a derivation it did — never "Paper X says so" / "the test
  passes"; (iii) per-claim taxonomy A externally-verified / B
  internally-consistent-only / C overstated / D unverifiable-here / E wrong;
  (iv) the **circularity map** (EXTERNAL vs GEOVAC-ONLY) as central deliverable.
- **`CITATION_CHECKER.md`** — web-enabled neutral citation + novelty checker.
  Honest ceiling: can *confirm-false* a novelty claim (find prior art) but never
  *confirm* novelty ("first in literature" → "to our knowledge"). Verdicts: OK /
  wrong-metadata / misattributed / doesn't-support / can't-find. Motivated by the
  documented fabricated-arXiv-ID incident in the §3 dead-ends table.

Mission framing (both): *de-risk a broadcast* (find embarrassing errors, grade
each claim by what it rests on, say how loudly each can be stated) — NOT certify
the physics, which no internal process can do (per `agents/WORKFLOW.md`).

## 3. Calibration — blind, on Paper 2

Calibrated by a **blind** run on Paper 2 (α observation), whose honest profile
was independently fixed first: noteworthy statistical observation; three
structural homes for B/F/Δ; combination rule explicitly conjectural; not a
derivation. Both agents, told nothing of the expected verdict, independently
landed there: YELLOW, honest-conjecture, headline cubic match 8.8e-8 reproduced
as EXTERNAL, citations GREEN (no fabrications). They surfaced real defects (§5)
and were correct even when a PM hand-check of one finding was briefly wrong.
Verdict: instruments are honest — neither rubber-stamp nor false-fraud.

## 4. Two methodology lessons (DURABLE — baked into `CONFIDENCE_AUDITOR.md`)

The most valuable durable output. Both are now mandatory checks in the auditor:

1. **Verify-the-verifier, against EXACT text.** Re-derive every consequential
   finding from the paper's exact wording (never a paraphrase) before it drives
   an edit. Incidents: (a) a PM hand-check "refuted" the auditor's §VIII.D
   finding by summing from k=0, but the paper's definition sums from k=1 — the
   auditor was right, the hand-check misread the index; (b) two flagged items
   (Paper 7 App. B arithmetic, Paper 0 "S⁵ six-dimensional") could not be located
   in the source on grep, so they were NOT edited.
2. **Cross-corpus check.** Per-paper-in-isolation audits systematically
   *over-flag*, because the corpus self-corrects across papers. Incident: an
   isolated audit of Papers 0/7 flagged the "graph Laplacian → spectrum /
   −1/16 topological" framing as a live crack — but Paper 18 §κ already
   reclassifies κ as conformal/calibration and states the spectrum comes from
   H = κℒ + W, "not from κ alone." The correct finding is "early papers LAG
   Paper 18," not "framework broken." The auditor must search the corpus (esp.
   Paper 18 taxonomy) for a later/companion treatment before calling a defect.

## 5. Verified paper fixes applied (Paper 2) — both PM-verified before editing

- **Table II + abstract `137.035987` → `137.036011`.** The printed value sat on
  the wrong side of CODATA (below). Solving the paper's own cubic
  x³−Kx+1=0 with K=π(42+ζ(2)−1/40) gives physical root 1/α = 137.036011 (above
  CODATA), same 8.8e-8 magnitude. Verified at 30 dps (mpmath).
- **§VIII.D summation index `k=1` → `k=0`** in the B_formal / N definitions. The
  paper's printed closed form N(2)=(d+1)(d+4)/2 and its stated quadratic
  m²+dm−2(d+2)=0 are both k=0-consistent; only the definition's lower limit read
  k=1. Verified symbolically (sympy): under k=0 the identity B_formal/N=d
  collapses to numerator d(m−2)(m+d+2)=d(m²+dm−2(d+2)), root m=2 — exactly the
  paper's claim; under k=1 it fails. One-character fix; no other text needed (the
  k=0 term contributes 0 to B_formal). Locked by `tests/test_paper2_corrections.py`.

## 6. Honest scope

- **Theorem grade:** nothing new proven. The §VIII.D fix restores an internal
  consistency the paper already claimed.
- **Verified (external/symbolic):** the two Paper 2 corrections (mpmath + sympy);
  the calibration verdict on Paper 2 (instruments honest).
- **Numerical observation (unchanged):** Paper 2's α cubic match stays a
  conjectural numerical observation — combination-rule "conjectural" label
  untouched (§13.5 hard-prohibition respected).
- **Deliberately NOT kept (wiped, PI direction):** the partial foundations-audit
  discoveries (Papers 7/0/1), including the convergence-framing question on
  Papers 0/7. Wiped so a fresh model audits cold; it is expected to re-derive
  them. The only carry-forward is the §4 instrument lessons.
- **Held, not applied (named follow-ons for the fresh audit):**
  - Paper 7 App. B arithmetic and Paper 0 "S⁵ six-dimensional" — auditor-claimed
    but exact text not located on grep; deferred.
  - Convergence-framing reframe on Papers 0/7 and CLAUDE.md §1/§1.5 —
    PI-territory; a sync to Paper 18's already-honest version, not a teardown.
  - `perez_sanchez_2024` citation annotation — checker-flagged, not PM-web-
    verified.
  - Full review campaign (Tiers 0–3, ~50 papers) is **unstarted** — only 3
    foundation papers were partially audited before the wipe.

## 7. Files

- **Added:** `agents/CONFIDENCE_AUDITOR.md`, `agents/CITATION_CHECKER.md`,
  `debug/confidence_review/LEDGER.md` (clean scaffold),
  `tests/test_paper2_corrections.py`, this memo.
- **Changed:** `papers/group5_qed_gauge/paper_2_alpha.tex` (value + §VIII.D
  index); `agents/CONFIDENCE_AUDITOR.md` patched with §4 lessons.
- **Cleared (PI direction):** 8 per-paper dossiers
  (`debug/audit_paper{2,7,0,1}.md`, `debug/citecheck_paper{2,7,0,1}.md`).
- **Left in place:** `debug/alpha_numerology_audit.py` (reusable tool, not a
  discovery).
