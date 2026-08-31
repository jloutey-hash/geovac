# `/qa trunk` re-certification — deterministic layer (2026-08-31)

**Status: deterministic layer COMPLETE and green. LLM panel NOT dispatched.**
Verdict is therefore **INCONCLUSIVE by construction** — four gating
dimensions (code/test-backing, claims, citations, synthesis) are unexercised,
and §7 says an unexercised gating dimension forces INCONCLUSIVE, never a
"PASS on the exercised dimensions."

This file exists so a fresh session can resume without re-deriving the
groundwork. The panel should run **cold** — the PM who did this layer spent a
long session forming opinions about this corpus, which is exactly the
condition §9 principle 2 says a reviewer must not be in.

## Criteria (frozen, unchanged)

From `docs/qa/trunk.done.md`. Scope: Papers **0, 1, 7** (group3) + **32, 38**
(group1) + the **group3 foundations synthesis**.

Branch deltas:
- **C7** — Paper 38 / WH1 is PROVEN *scoped to the van Suijlekom state-space
  GH distance* (translation-seminorm metrization), with no residual
  "Latrémolière propinquity" overclaim.
- **C8** — κ = −1/16 is an **Observation** (no bridge), not a derivation;
  **4/π** is **derived (numerics-pinned)**, not full symbolic proof; the
  Forced-Count moduli chain is stated at its **full-axiom** count.
- No branch-specific C14+.

## Why this re-cert was triggered

`check_cert_staleness.py`: trunk certified **2026-06-16**, and **4 `.tex`
have changed since** — `paper_32_spectral_triple`,
`paper_38_su2_propinquity_convergence`, `Paper_0_Geometric_Packing`,
`Paper_7_Dimensionless_Vacuum`. All 11 QA targets are similarly owed; trunk
goes first because §9 says a finding at a root re-prices everything above it.

## Deterministic findings — 2 real defects, both fixed

### 1. C19 examined ZERO trunk papers while reporting PASS

`check_latex_escapes.py --gate` is a **substring filter on the path**, and no
trunk paper's path contains the string `trunk`. It reported
`PASS ... in 0 paper(s) in scope 'trunk'`. **It had never inspected a trunk
document.**

Exactly the gate self-audit failure the protocol names: *"a gate that fires
correctly and scopes its verdict away is indistinguishable, from the outside,
from a working one."* It was visible only because C19's RESULT line prints a
file count — the other nine mostly do not, which is why their coverage still
cannot be confirmed from output alone (see Open below).

**Fixed:** named `trunk` scope from the DoD; substring behaviour kept for
every other value; a `[scope] WARNING` prints any scope member not found.
Now `6 paper(s)`. Selftest still fires (10 positive / 5 negative).

That warning immediately caught **two errors in the fix itself** — the glob
yields absolute paths while the scope list is relative, and Paper 1's `.tex`
is lower-case `paper_1_spectrum.tex` while its PDF is `Paper_1_Spectrum.pdf`.
A silent 5-of-6 would have looked like success.

### 2. Paper 38 carried NO C16 entry at all

`propinquity-as-achieved-metric` listed papers **39 and 40** explicitly —
not 38, the paper whose trunk C7 criterion names that exact overclaim.
(Severity is *advisory*, so it would not have gated regardless; the gap is
still worth closing.)

**Fixed:** papers 38 and 32 added to the entry. Two-way discrimination proven
per the hard rule — FIRES on three retired wordings, SILENT on two corrected
wordings **and** on Paper 38's live text (0 hits).

## A finding I reported and then had to retract

I initially reported C16's **line-wrap blindness** as "the one that matters
most." The mechanism is real: patterns match line-by-line with `[^.\n]{0,45}`
spans, and Paper 38 line 534 ends `stated the main theorem in the` while 535
begins `Latr\'emoli\`ere propinquity`, so that alternative never fires there.

**Measured, the impact is ZERO** — joining adjacent line pairs across every
registry entry and every gated file surfaces **0** additional live hits, **0**
on fail-severity entries. The Paper 38 straddle sits inside a
`\begin{remark}[history]` describing an earlier draft: disclosed history,
which the exemption exists for.

I generalised from one mechanism to a corpus-wide claim without measuring it.
Recorded in the C16 docstring as a **latent limit with its measurement**, so
nobody re-derives it as an alarm. (Scope: single-wrap straddles; spans are
≤45 chars so that dominates, but a two-wrap straddle was not tested.)

## Deterministic scorecard — 12/12 PASS

| gate | result | scope actually examined |
|:--|:--|:--|
| C11 internal titles | PASS | trunk |
| C5 K-label | PASS | trunk |
| C13 paper→test refs | PASS | trunk |
| C14 paper→file refs | PASS | trunk |
| C16 retracted terms | PASS | trunk (0 trunk-scoped entries; papers covered via group1/group3 globs) |
| C17 headline numbers | PASS | trunk (15 exempt/historical) |
| C18 duration language | PASS | trunk |
| C19 eaten escapes | PASS | **6 papers** (was 0) |
| C15 inline arXiv | PASS | trunk |
| C20 inline attributions | PASS | baseline |
| C21 numeric consistency | PASS | all scopes |
| C22 test-claim backing | PASS | corpus |

## What the fresh panel must do

Dispatch **all four** LLM dimensions (§7: the verdict is the AND across them):

1. **code / test-backing** — `code-reviewer`, one per trunk paper with tests.
   Apply the **independent re-derivation** rule: recompute constants by an
   independent route, and fire-test every regression guard by constructing
   the wrong value.
2. **claims / prose** — `claims-reviewer`. **Enumerate every K-sentence**
   (the K-prohibition is semantic and a hedge does not cure a tripwire), every
   κ and 4/π tier statement, and every status-table row.
3. **citations** — `citation-reviewer`, CONFIRMED / WRONG / **UNVERIFIABLE**,
   with the no-upgrade-on-plausibility rule.
4. **synthesis** — separate `claims-reviewer` on the group3 foundations
   synthesis.
5. **completeness-critic** — FULL runs only; it produced the highest-value
   finding of FULL cert #2.

Seeds: ≥1 per dimension, planted in a worktree (`git worktree add
../geovac-qa-seed-trunk`), answer key to `debug/qa/trunk_seed_key.json`,
**never** in the real corpus. Path-pin every reviewer prompt to the worktree.

**Specific attention for the claims reviewer, unresolved here:** Paper 38
line 291 reads *"Lemma~\ref{lem:L5} the Latrémolière propinquity assembly"*,
and the file header comment reads *"SU(2)/S^3 Latremoliere propinquity
convergence"*. The paper defines "Latrémolière propinquity" as a term of art
in its **Propinquity convention** paragraph (§Setup) and records the descope
in `rem:history38`. **Whether those two loci are a residual C7 overclaim is a
judgment call I deliberately did not make** — the PM here is not a fresh
adversary on it, and C16 is silent because the phrasings do not match its
retired patterns.

## Open / not done

- **Seven of ten gates do not report a file count**, so their trunk coverage
  cannot be confirmed from output. C19's bug was visible *only* because it
  prints one. Worth making the count mandatory in the RESULT line.
- `check_cert_staleness.py`'s trunk scope **omits the group3 synthesis**,
  which the DoD includes — same class of scope gap, one level up.
- The LLM panel, entirely.
