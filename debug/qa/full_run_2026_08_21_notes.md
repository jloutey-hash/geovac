# /qa FULL certifying run — 2026-08-21 — VERDICT: INCONCLUSIVE (+ FAIL findings)

Target: papers 58, 59, 60, 19 (group2) + 56, 24 (group3) + the group2/group3 syntheses.
Shape: FULL (whole-paper enumeration, all dimensions, completeness-critic).
Seeds: `debug/qa/full_seed_key.json` (7 seeds, 6 controls). Worktree destroyed; zero leakage verified.

## VERDICT: **INCONCLUSIVE**

Forced by one dimension, per the hard rule ("an unexercised or uncalibrated gating
dimension forces INCONCLUSIVE, never PASS"). Genuine FAIL-grade defects were also found
and are listed below — but the run cannot certify.

### Per-dimension scorecard

| Dimension | Exercised | Calibrated | Clean |
|:--|:--:|:--:|:--:|
| Deterministic (11 gates, whole-target) | YES | n/a | YES |
| Claims P58/59/60 | YES (×2 attempts) | **NO** | — |
| Claims P19/56/24 | YES | YES (1/1) | NO (7 genuine) |
| Claims syntheses | YES | YES (1/1) | NO (1 LARGE + 3) |
| Code P58+artifact | YES | YES (2/2) | NO (1 genuine + NITs) |
| Code P59/P60 | YES | YES (1/1) | NO (2 NITs) |
| Code P19/P24 | YES | no seed in scope | NO (found the energy-ladder error) |
| Citation group2 | YES | YES (1/1) | NO (1 genuine) |
| Citation group3+synth | YES | YES (1/1) | NO (NITs only) |
| Completeness-critic | YES | n/a | NO (found 10 gaps + 2 defects) |

Seeds: **6 of 7 caught** (F4 caught twice, by citation AND code-A). F1 missed by two
independent opus reviewers — see "the calibration failure" below.

## The calibration failure (claims-A)

F1 planted `\textbf{[SYMBOLIC PROOF]} We therefore derive that N(D) is a period of an
irreducible rank-four connection...` replacing the live `[SYMBOLIC + MEASURED]
Consequently ...`. Reviewer 1 (whole-paper sweep) did not surface it. Reviewer 2
(mechanical per-tag adjudication, 112 tags one at a time) explicitly adjudicated that
tag as SOUND, reasoning the FL/Katz leg is a standalone symbolic argument so "derive"
is licensed.

**Read: this is at least as much a goalpost gap as a reviewer failure.** `criteria.md`
does not sharply separate `[SYMBOLIC PROOF]` from `[SYMBOLIC + MEASURED]` for the
"symbolic argument + numerical witness" class, and the corpus itself tags this exact
claim mixed. RECOMMENDED (PI-gated): add a criteria.md rule — *when a claim's argument
is symbolic but its certificate is numerical, the tier is [SYMBOLIC + MEASURED]; a bare
[SYMBOLIC PROOF] requires that no numerical input appears anywhere in the chain* — and
re-run the dimension against it.

## Genuine defects found (seed-free), by dimension

**LARGE**
- **group2 synthesis asserts the Sturmian structural theorem as an unqualified no-go in
  6 loci** ("cannot bind", "provably incapable of binding", "disproves every single-center
  Sturmian alternative") — contradicted by Paper 8's own scope remark ("delimits a specific
  approximation; it is not a no-go for Sturmian molecular binding") AND by Paper 60, which
  binds H2+ and H2 in a common-scale basis. A certified synthesis contradicting its own
  guardrail paper after the paper was scoped.
- **group2 synthesis omits Papers 59 and 60 entirely** ("the ten Group 2 papers"; twelve
  exist). Its own "open polyatomic question" was answered structurally by Paper 59.

**SMALL / MATERIAL (verified against primary text)**
- P19 energy ladder did not reproduce from its cited JSON (offset 7.2771 vs the codebase's
  −7.2799) — **FIXED**: −7.9297/−8.0548/−8.1020, offset now stated inline; crossing verdict
  unaffected under every convention.
- P56: rank formula `2·3N+1` contradicts its own proof (rank 46 at n_max=2) and rank–nullity
  (9N+1); only the first listed value is right.
- P56: two un-hedged loci ("GeoVac IS the abelianisation", "↪ closed sub-pro-algebraic
  group") assert the closed immersion the theorem refutes.
- P24: one sentence assigns calibration π to the HO, negating the paper's own headline.
- P19 abstract "confirming convergence" omits the n_max=4 crossing the body calls decisive.
- P19: C17-retired 33.3 Ha survives at a group2 locus (family is group4-scoped — registry
  maintenance needed).
- P19: "consistent with O(B²) scaling" against three points growing as B^1.14.
- P24: pair-diagonal 1.44% quoted in the universal frame (6.06% is the universal number) —
  the framing-zombie rule by name.
- P19 `hoggan2011` cross-wires two distinct Hoggan publications (title vs venue).
- Artifact rows `qfd.hehplus`/`qfd.behplus` declare `backing_test: tests/test_paper58_qfd.py`
  which contains zero HeH/BeH references — a false backing pointer.
- **FIXED this run:** group3 synthesis six-vs-seven layers (body count + missing Layer-7
  item + closing paragraph, inside a starred subsection); P58 `tab:backing` caption
  "four files" → six; `geovac/t2_kw.py` docstring promising a nonexistent `value()` API;
  `debug/README.md` preservation list added (the T2/QFD/decider/probe drivers are the only
  record of certifications whose tracked tests reach lower precision).

## The completeness-critic's decisive contribution

It found that **the code dimension had not actually been exercised for most of the target**:
~24 declared backing files un-run, including ALL FIVE Paper 56 tests (keystone theorem,
previously refuted once, guarded only by a strict-xfail tripwire), Paper 58's title-theorem
suite, `test_routeC_momentum.py` (six of Paper 59's eleven sections), and three of Paper 60's
four headline-exponent suites. **Cause: my dispatch named the NEWEST files per paper rather
than the load-bearing ones.**

**Closed post-hoc in this run:** all named suites executed — P56 87 passed (both strict
xfails HELD), P58 15 passed, P59 47 passed, P60 13 passed, plus 44 slow legs green
(headline exponents Q^3.33/Q^1.19/K^0.84/cond 4→3673/n_orb^2.2 and the P58 title theorem).
Zero failures. The coverage gap is closed as a matter of fact, but the dimension was not
calibrated-and-exercised *as dispatched*, which is part of why the verdict is INCONCLUSIVE
rather than FAIL.

## Open questions for the PI (criteria-level, not reviewer-level)

1. **C9 has no omission leg.** The group2 synthesis being an arc behind its own group is
   invisible to every criterion (C9 hunts overclaim/zombies = commission). Does certification
   require synthesis *currency*?
2. **Tier boundary** `[SYMBOLIC PROOF]` vs `[SYMBOLIC + MEASURED]` (above).
3. **Starred headings are a structural blind spot** — the group3 synthesis has 27
   `\subsection*{}` vs 16 numbered; D1 lived there. Reviewer prompts should enumerate
   starred headings explicitly.

## Honest ceiling

"INCONCLUSIVE" here means the panel's clean verdicts on P58/59/60's claims dimension are
untrustworthy this run, not that defects are known to exist there. Everything else was
calibrated. Deterministic gates and the now-complete test execution are facts, not judgments.

---

# ADDENDUM — claims-A re-run against the sharpened C3-boundary (2026-08-22)

**PI rulings received:** C9 requires synthesis **currency** (ratified, now in criteria.md).
C3-boundary tier rule **adopted pending ratification** (also in criteria.md) — it was
load-bearing for this re-run, so it could not stay open.

**Remediation completed first** (all outstanding full-run findings): the group2 synthesis
Sturmian zombie scoped at all six loci + Papers 59/60 added + count corrected (C9-currency
now satisfied: 12/12 in-group papers present); P56 rank formula 2*3N+1 -> 9N+1 (19/46/82/127,
tied to rank-nullity) and both un-hedged closed-immersion loci; P24 pi-assignment sentence and
the 1.44%/6.06% framing zombie; P19 abstract crossing disclosure, retired 33.3 -> 32.6 Ha,
O(B^2) -> measured B^1.14, hoggan2011 citation; P60 three tier/anchoring items; artifact
backing_test pointers made per-row honest; C17 lih-onenorm family widened beyond group4;
t2_kw docstring; debug/README preservation list. All gates green, all papers compile.

**Re-run result: 2 of 3 seeds caught. Dimension REMAINS INCONCLUSIVE.**

| Seed | Class | Outcome |
|:--|:--|:--|
| A1 `[SYMBOLIC PROOF]` + "We therefore derive" (the exact case that defeated TWO prior reviewers) | C3-boundary | **CAUGHT — graded LARGE**, dismantled three independent ways (the rule mechanically; the abstract's divergent tier on the same claim; the frozen DoD's own mixed canonical form), with the reviewer explicitly stating it would not reconstruct the now-disallowed "the symbolic leg licenses the verb" reasoning. **The amendment works.** |
| A2 "exactly ... which we prove" under `[MEASURED]` | C3-boundary | **CAUGHT — SMALL**, with the 2%-tolerance backing quoted and the body's own "lands near" hedge contrasted |
| A3 W1 inversion: "its precision, and thereby its accuracy at this basis size" | exact != accurate | **MISSED — and affirmatively cleared** (the reviewer enumerated L976-980 among honest-scope loci and called them "All intact and in the honest direction") |

**Diagnosis.** The targeted fix worked on its target class and only on its target class. The
residual miss is the *direction* of an honest-scope sentence: the reviewer enumerated the loci
(the mandate was followed) but read them as a checklist rather than checking what each one
DENIES. Recommended next sharpening, mirroring what worked for C3: an explicit mandate to
quote each honest-scope sentence AND state in its own words what the sentence denies
("certified = assembly, NOT accuracy"), so an inverted one cannot pass as present.

**Two-way discipline held.** The re-run's N4 (P58 abstract stating the Sturmian theorem as an
unqualified heteronuclear no-go) was adjudicated against primary text and DISSOLVED: Paper 8's
`cor:dual_p0` is unretracted and states exactly that ("no shared p_0 exists that simultaneously
represents both atomic length scales"). An earlier reviewer had independently marked the same
clause SOUND. No edit made — the primary text decided, per the reconcile rule.

**Seed-free fixes applied from this run:** P60 "novel" -> "apparently novel" (aligning the
abstract with the two hedged body loci); P59 "proven irregular twist below" repointed to where
the irregularity is actually established symbolically (Sec. obstruction, Newton polygon);
U2 upgrade — the Euclidean-propagator dictionary retagged `[MEASURED]` -> `[SYMBOLIC + MEASURED]`
(its "exactly" clauses are genuine closed-form identities, so raising the tier is the correct
cure rather than softening true statements).

**Standing deferred upgrades re-surfaced, not forced:** eq:K0 6e-18 -> witnessed <1e-35;
eq:modpf 1e-25 -> ~5e-32; P58 tab:backing S,h rows MEASURED -> DECIDED.

**Net verdict unchanged: INCONCLUSIVE**, now on a single, precisely-diagnosed axis instead of
a diffuse one. Everything else in the target is calibrated-and-remediated.

---

# ADDENDUM 2 — the 214 builder column RESOLVED (2026-08-22)

The one genuine open number from the full run is closed. It was never lost in an
edit: the generator was `debug/noci_n3b_census.py` (added by commit `8c747df`,
"ERI inflation genuine/builder 13.8x (n_max=2) -> 15.0x (n_max=3)"), which is
still in the working tree but is a `debug/` driver — prunable by the Clean Room
Rule and uncitable by a paper, so from `tests/` the number looked ungrounded.

**It reproduces the published table exactly:**

| n_max | dense | m-rule | all-four-one-center | genuine | builder | ratio |
|:--|--:|--:|--:|--:|--:|--:|
| 2 | 10,000 | 3,120 | 390 | 2,944 | 214 | 13.76x |
| 3 | 614,656 | 121,920 | 15,240 | 114,280 | 7,600 | 15.04x |

**What the builder column actually is** — and this settles the full run's
MATERIAL-1 hypothesis. Both columns apply the *identical* axial rule
`m_p + m_r = m_q + m_s` (global M_L, "rule B"). They differ in exactly one
condition: the builder's block-diagonal surrogate never forms a two-center
orbital product, so it can represent a quartet only when all four orbitals sit
on ONE center, where the single-center multipole additionally requires the two
pairs' Gaunt L-ranges to intersect. That shared-L condition is imposed
identically in both columns — it accounts for the *same* 176 removals at
n_max=2 (390 -> 214 in the builder, 3,120 -> 2,944 in the genuine column).

So **builder == genuine restricted to all-four-on-one-center**, and the
inflation factor measures the two-center span and nothing else.

**The claims panel's rule-A/rule-B convention-artifact hypothesis is REFUTED.**
Neither column uses the pair-diagonal rule; my own earlier combinatorics
(rule-B same-center 780, rule-A same-center 484) missed because both readings
allowed (AA|BB), which the builder cannot represent at all. Paper 58's prose
was accurate throughout — no claim was ever overstated.

**Actions taken (so this cannot recur):**
- `tests/test_paper58_census.py` — new tracked test re-deriving every count in
  `tab:census` from the selection rules alone: the `g` row at both n_max, the
  published 13.8x/15.0x, the structural identity builder = genuine-restricted,
  and (independently of the exact-rational route) the S and h rows 32 = 10+22
  and 44 = 22+22. 5 tests, 0.7 s.
- Paper 58 states the two-column relationship inline and cites the test.
- `docs/claim_test_matrix.md` — the NO-TEST (owed) row becomes BACKED-SOUND;
  the sibling `g`-row entry is re-scoped so BACKED-WEAK now attaches only to the
  Gaussian-engine corroboration of the genuine density, not to the counts.
- C17 registry gains `p58-census-builder`, which FAILs on the refuted 780/484
  readings appearing as live builder counts.
- `debug/README.md` preservation list gains the generator: its `g`-row counts
  are now redundant, but the exact-rational two-center machinery that DECIDES
  the S/h cross-block non-vanishing lives only there.

**Verification:** all six deterministic gates PASS on group2; Paper 58 compiles
clean (11 pp, 0 undefined refs); 106 passed / 11 skipped across the five Paper 58
backing suites.

**Verdict impact: none.** This was a documentation-provenance gap, not a defect,
and it was not the axis forcing INCONCLUSIVE. That axis remains the claims-A
calibration miss (seed A3, the "exact != accurate" inversion) and is unchanged.

---

# ADDENDUM 3 — claims-A re-run #2: dimension now CALIBRATED (2026-08-22)

The axis that forced INCONCLUSIVE is closed. Seed key
`debug/qa/claimsA2_seed_key.json`; worktree destroyed, zero leakage verified
against every seed (planted text absent AND original text present in all three
papers).

## Calibration: 3/3 sensitivity, 3/3 specificity

| Seed | Class | Outcome |
|:--|:--|:--|
| B1 — P58 "certified = assembly, not accuracy" INVERTED | honest-scope direction | **CAUGHT, graded LARGE.** Reviewer quoted it, named the denial it reversed, grounded the canonical form in four independent repo loci, and flagged the self-contradiction against the sentence's own 0.068 Ha antecedent. **This is the exact class that was MISSED and affirmatively cleared on the prior run.** |
| B2 — P60 "basis incompleteness, not a missing metric" INVERTED | honest-scope direction | **CAUGHT, graded LARGE**, dismantled three ways from inside the paper (its own supporting clause argues for basis; §iso says the metric is "never present"; the one genuine metric caveat explicitly denies an energy effect). The fix generalises past the phrasing it was written from. |
| B3 — P59 bare `[SYMBOLIC PROOF]` + "We prove" on a 4.3e-7 certificate | C3-boundary | **CAUGHT, graded LARGE.** Regression control: confirms the 2026-08-22 C3-boundary amendment has not decayed. |
| C1/C2/C3 — known-good controls | specificity | **None flagged MATERIAL.** Zero false positives. |

**What changed and why it worked.** `.claude/agents/claims-reviewer.md` rule 7:
honest-scope sentences must be quoted AND their denial restated in the reviewer's
own words. The prior miss was not a failure to enumerate — that reviewer listed
the loci correctly and reported them "all intact and in the honest direction."
It read them as a presence checklist. An inverted honest-scope sentence reads
perfectly natural in place; only asking what it DENIES catches it. The re-run
enumerated 72 such sentences and reported 70 intact, 2 reversed — exactly the
discrimination the rule was written to produce.

## Genuine (seed-free) findings, all adjudicated against primary text

**FIXED**
- **P59 corrupted `\ref` — a NEW instance of the eaten-escape class, introduced
  by my own earlier remediation in this very arc.** `Sec.~\ref{sec:obstruction}`
  had become `Sec.~<CR>ef{sec:obstruction}` and `Poincar\'e` had become
  `Poincar'e`. This is the pointer carrying the paper's non-classicality
  attribution. **It compiled clean** — LaTeX renders the literal text, no error,
  no warning, no undefined reference. Repaired, and a corpus-wide scan confirms
  it was the only instance in 62 active papers.
- **P59 intersection form: `[SYMBOLIC]` + "no numerics enter" over-reached.** The
  branch-point-asymptotics leg is genuinely numerics-free (the v4.97.0
  `[OBSERVATION]→[SYMBOLIC]` upgrade covered exactly that), but the forcing
  argument and the Sp_4(Z) containment ride on M_0, which the paper itself tags
  `[MEASURED; exact-integer certificate]`. Scoped: the blanket claim now reads
  "no numerics enter in this leg", and an explicit tier note marks the
  M_0-dependent sub-clauses `[SYMBOLIC + MEASURED]`.
- **P58 "$<6\%$ residual"** → `5.6%--6.3%` (6.3% is not $<6\%$).
- **P58 "$80\%$" → "$79.5\%$"** — and *decided*, not merely accepted: the census
  gives exactly 2,340/2,944 = 79.5% cross vs 604/2,944 = 20.5% same-side at
  n_max=2. Now pinned by `test_paper58_census.py::test_class_split_79_5_vs_20_5`.
- **P58 "countable" → "counted"** (the paper's own tier word).
- **P60 "the only well-conditioned posing"** → "...we have found" (generalised
  over all posings from a two-posing comparison).
- **P60 unclosed parenthesis** in the QSVT aside.

**TWO-WAY UPGRADES applied** (claims weaker than their backing)
- **P60 compound-matrix result `[MEASURED]` → `[SYMBOLIC + MEASURED]`.** "The
  molecular metric penalty is one-electron, not N-electron" rests on an algebraic
  identity (the k-electron configuration overlap IS the k-th compound matrix) plus
  by-construction SW-orthonormality; only the 1x/8x/9x L^2 figures are measured.
- **P60 abstract now names the exact gerade constant** 2/(1+min j_0) = 2.555041...,
  basis- and R-independent, rather than only its small-n reading "flat (kappa~2)".

**DISSOLVED against primary text** (two-way discipline)
- The reviewer's NIT that "$10^3$--$10^4\times$ smaller (gerade)" is understated
  computed the ratio from the **ratio-1.6** Gaussian variant; the sentence is
  explicitly scoped to a **ratio-2** even-tempered metric, and the ratio
  dependence is disclosed as an honest caveat in the very next clause. No edit.

## New deterministic gate: C19 (eaten-escape corruption)

`debug/qa/check_latex_escapes.py`, registered in `criteria.md`. Detects LaTeX
control sequences destroyed by Python string escapes (CR+`ef`, BS+`extbf`,
TAB+`imes`, ...) plus bare control characters that are never legitimate in .tex.
Excludes `papers/archive/`; ignores indentation TABs; `--selftest` covers 4
positive and 5 negative cases. **PASS on all 62 active papers.**

This is the C16/C17 pattern applied to a class that had bitten the corpus three
times: mechanical, zero-variance, invisible to compile checks, and repeatedly
found only by someone noticing odd rendering. (It bit once more *while the gate
was being written*, in a heredoc used to tighten the gate's own regex — which is
the strongest argument for having it.)

## Verification

C10 compiles: P58 11 pp, P59 15 pp, P60 8 pp, zero undefined refs each. C5, C11,
C13, C14, C16, C17, C18, C19 all PASS corpus-wide. `test_paper58_census.py`
6 passed. The 18 symbolic S^3 proofs pass.

## Net verdict

**The claims dimension is now exercised, calibrated, and remediated.** Every
gating dimension of the full run has now reached that state. The run's residual
ceiling is no longer a calibration gap — it is that a *fresh* full certifying
pass has not been fired since these remediations landed, which is a PI timing
decision, not an open defect.
