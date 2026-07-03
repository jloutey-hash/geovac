# /qa group5 — 1st cert (FULL run) — 2026-07-02/03 — FAIL → Tier-1 REMEDIATED

**Shape:** FULL (fresh target; DoD frozen by PI same day). Worktree at 8f3efce, 18 seeds
committed on-branch (run-8 hardening; git status clean). Panel: 5 code (Sonnet ×2 seeds),
2 citation (Sonnet ×2 seeds), 3 claims (Opus), 1 synthesis (Opus); critic + recovery wave after.

## Calibration — round 1

| Agent | Seeds | Result |
|---|---|---|
| claims-{2,25,30} (opus) | s-p1 | 1/1 ✓ (run-4-pattern ID) |
| claims-{28,33,41} (opus) | s-p2 | 1/1 ✓ (title/§7 demolition) |
| claims-{36,51} (opus) | s-p3 | 1/1 ✓ (G8 self-contradiction) |
| synthesis (opus) | s-syn | 1/1 ✓ (4-way internal-consistency) + cross-caught s-p3 |
| citation-{2,25,28,30} (sonnet) | s-cit1a/b | 2/2 ✓ (noted P25's own correct Wilson page) |
| citation-{33,36,41,51} (sonnet) | s-cit2a/b | 2/2 ✓ (both LARGE, correct real locators) |
| code-P33 (sonnet) | s-c33a/b | 2/2 ✓ (quoted both verbatim) |
| code-P51 (sonnet) | s-c51a/b | 1/2 (caught ζ tolerance w/ Integer(0) proof; missed L6 rel_err) |
| code-{25,30} (sonnet) | s-c2530a/b | 1/2 (caught hodge ±1; missed su2 atol) |
| code-P28 (sonnet) | s-c28a/b | 0/2 (read the intact sibling loop test; missed both) — but cross-caught s-c33a/b |
| code-{2,36,41} (sonnet) | s-c23641a/b | 0/2 INSTRUMENT ERROR (seeds were in test_circulant_s3 = the WH1 comparator, NOT P2 backing — the agent correctly identified the non-backing and skipped it) |

**Specificity: 6/6** — zero false positives on controls anywhere; the "live bait" resolved
honestly (P51 carries every demotion in place — Möbius→B4, soft-IR→GD-2 in the adjacent
section, Fursaev retraction disclosed; claims-{36,51} explicitly cleared them all).
Cross-catches: s-p1 (by {25,30}-code), s-p3 (×3: synthesis, P51-code, claims-{36,51}),
s-c33a/b (by P28-code).

**Recovery wave (dispatched):** 3 Opus narrow assert-audits (P28 pair, P51 pair, {25,30}
pair — the missed seeds still live = their calibration) + 1 Sonnet recal for {2,36,41}
with fresh seeds s-c23641c/d in genuinely-backing files (committed 'sync 2') + the critic.

## Genuine findings (PM-verified where noted) — the remediation register

### A. Hard-prohibition / tier-label class (the pre-downgrade idiom)
- **A1 (LARGE, verified):** P25:147-150 + :1053-54 — "the additive rule itself remains
  conjectural" — the prohibited tier word ON the combination rule. + ~12 interpretation-level
  "Paper 2 conjecture(s)" loci (85, 140, 156, 163, 168, 170, 589, 604-05, 782, 1037, 1046,
  1065): whole-paper conjecture→observation idiom sweep needed (per-locus care: the rule =
  Observation; the bundle-interpretation may stay "proposed interpretation").
- **A2 (LARGE, verified):** P33:65 + :707 "the Paper~2 fine-structure-constant conjecture"
  ×2 (its own bibitem :734 correctly says "observation status").
- **A3 (SMALL):** P30:925 "within a conjectural framework"; P2:407,410 stale LaTeX comments.
- Note: the C5 screen structurally cannot see these (no K formula in-clause) — the
  claims-reviewer enumeration is the load-bearing control, as designed.

### B. Keystone-metric / status zombies
- **B1 (SMALL, verified):** P51:159 + :2233 "Latrémolière propinquity" (should be van
  Suijlekom state-space GH; P51's own P38 bibitem has the correct title). :2233 is inside
  a Lemma — reword needs care vs the descoped P44/45/47 cluster (defer-scope honestly).
  **C16 maintenance:** broaden the propinquity registry pattern (it missed these).
- **B2 (SMALL, verified):** P36 conclusion :1044-48 restates the Table-7.4 residual
  attribution its own §5.2 (:484-494) retracts, and frames LS-7 as future when §6.1
  reports it done-WEAK.

### C. Abstract/enumeration errors
- **C1 (SMALL, verified):** P33 abstract :41-43 six-rule list includes Gaunt/CG (Tier I —
  double-count) and omits triangle-on-n; fix per the partition theorem (1+6+1).
- **C2 (SMALL):** P33-vs-P28 census-tier mismatch on charge conjugation (P33: spinor-
  recovered; P28: vector-required; P33 merges C+Furry as R7). Harmonize the 8-rule
  enumeration across the two papers (P28 keystone convention wins) — needs care.
- **C3 (SMALL, verified):** P2:496 + :587 "nine (distinct) mechanisms" → twelve
  (Phases 4B-4I + Sprint A + K-CC; synthesis + CLAUDE.md carry twelve).

### D. Citation defects (real corpus)
- **D1 (LARGE):** P28 parker1980 misattributed — PRD 29:1584 (1984) is GUT RG-running,
  not the g-2 curvature formula (the 0.4%-match headline's source). The Sprint memory
  (parker_toms_curved_qed) says the physics match is real — locate the right source
  (likely Parker–Toms 2009 textbook §6) and repoint. PM web-verify.
- **D2 (SMALL):** P28 petermann/sommerfield a₂ printed with −197/144 (real: +197/144;
  context-only remark).
- **D3 (SMALL):** P2+P25 perez_sanchez 2024/2025 bundling (2024's own result is
  WITH-Higgs for its quiver extension; cite the 2025 Comment alone for without-Higgs —
  P30 already does it right).
- **D4 (SMALL):** P36 drake_swainson genuine title paraphrase (real: "Bethe logarithms
  for hydrogen up to n=20, and approximations for two-electron atoms").
- **D5 (LARGE-adjacent):** P36 missing bibitems for named numeric anchors: Antognini 2013
  (Science 339, 417), Karshenboim 2005 (Phys. Rep. 422, 1), PDG 2024, CODATA. Add.
- **D6 (SMALL/verify):** P51 frolov_fursaev1997 venue likely PRD 56, 2212 not JHEP;
  P51 "Dowker (1977)" inline with no bibitem (only dowker1994); Beccaria–Tseytlin
  mentioned with no bibitem; sommerfeld 1894→1896; P41 loan_hamer missing author
  (Brunner); P33 jeffrey→Jeffrey-Dai 2008; P30 6 uncited bibitems + :369 cites creutz
  where Fegan–Menotti–Onofri is invoked; P51 13 orphaned bibitems.
- **D7 (SMALL):** P28 schwartz1961 scope shorthand (helium ground state, not hydrogen 2P).
- **D8 (SOFT):** P33 chamseddine_connes2010 over-attribution for 1/(4π)=a₂(S²) — reword
  as vocabulary-source, not theorem-source.
- Good news: the Fursaev/Solodukhin cluster (the fabrication site) checks out CLEAN.

### E. Test-backing / production-code defects (the heavy layer)
- **E1 (LARGE, P33/vector_qed):** Furry 8/8 mechanism HARDCODED in production
  (`if a == b: return 0.0` before any physics) — tests verify a tautology; the real
  symbolic derivation lives in debug/archive (no asserts; reviewer ran it — math holds).
  FIX: port the derivation into tests/ with real asserts (earns symbolic-proof tier).
- **E2 (LARGE, P33):** census rules 3 (hardcoded pass=True) + 8 (Hermiticity of a
  symmetric-by-construction Σ) tautological in production; paper/code rule numbering
  swapped (6↔7); "triangle on SO(4)" (paper's Rule 8) has NO computational backing.
  FIX: honest relabel or real checks; align numbering; add triangle check.
- **E3 (LARGE, P33):** fabricated test citation `test_vector_qed_l1_4pi` (nonexistent);
  the underlying 1/(4π) scalar value verified true by the reviewer by hand. FIX: real
  test + fix citation. **Gate-hardening: extend C13 to bare test-function-name cites.**
- **E4 (MATERIAL, P33 + production BUG):** compute_self_energy(4, exact=True) CRASHES
  (graph_qed_photon.py:577-79 sorted(key=float) on complex-cast-failing sympy eigenvalues).
  Paper claims pendant-edge verified n_max=2..6 — currently FALSE above 2 on this path.
  FIX: bug (sort key → re()), then verify 3..6 live, then tests + honest scoping.
- **E5 (MATERIAL, P28):** pendant-edge "verified numerically n_max=2..100 + symbolically
  3,4" — tests cover n_max=2 only (the archived driver swept 2..5). FIX: after E4, add a
  sweep test (2..some honest N) + correct the prose to what is verified.
- **E6 (LARGE, corpus BUG):** module-level `mpmath.mp.dps` assignments collide across
  test files (collection-order-dependent; 11 reproducible failures in the SCOPE-order
  run). FIX: mp.dps set/restore fixtures or workdps contexts in the 4-5 offender files.
- **E7 (LARGE, P36):** headline Lamb chain (LS-1..6a, 1052.19 MHz) NO-TEST (debug-archived
  only); Drake–Swainson spot-check test is tautological (both sides same formula);
  P36:633 cites the WRONG module (qed_two_loop ≠ two_loop_self_energy). FIX: repoint the
  module cite; add at least a pinned recompute of the LS chain result if a driver is
  portable, else matrix-log as PI-deferred (heavy).
- **E8 (LARGE, P41):** all seven witnesses debug-only (4 inline cites already dangling;
  zero geovac/ production module). FIX: durability decision — port the four XCWG drivers
  to tests/wilson_rule_b_support/ (the wh7/mpo pattern) with pinning tests where runtimes
  allow; else scope the paper's verification sentences to memo-recorded and log.
- **E9 (MATERIAL, P30):** Result 2 (L1 kinetic/1-8 coefficient) NO-TEST (sibling SU(3)
  1/12 IS tested — asymmetry); Result 3 (MC Wilson table) NO-TEST (monte_carlo_wilson_
  expectation never called in tests). FIX: add both (MC test slow-marked, seeded RNG).
- **E10 (MATERIAL, P25):** Prop 1 item 4 (S² quotient spectrum {0,0,0,1,3,6}) NO-TEST;
  §VII.A S⁵/CP² (24.98% residual etc.) NO-TEST + dangling debug cites. FIX: add the
  quotient-spectrum test (cheap); port/pin S⁵ data or scope honestly.
- **E11 (MATERIAL, P51):** three C8 headlines NO-TEST (a₀=1.992; SC 99.4%; K_cone
  6-digit — the live F5 test proves less than prose); j_blindness test hard-imports
  debug/g6_fierz_pauli.py (durability hazard → tests/gravity_support/). FIX: durability
  migration + either port drivers or scope prose to memo-recorded.
- **E12 (LARGE, P51):** ZERO inline provenance tiers in the whole paper (every other
  group5 paper has 2-31). FIX: tier-tagging pass over P51's results sections.
- **E13 (P2):** B=42/selection/p-value/circulant §VI NO-TEST (debug-archived);
  P2:680-84 overstates dirac_s3.py as source of the obstruction identities (they live
  in archived debug). FIX: correct the source list; add cheap exact tests for B=42 +
  the circulant char-poly (both pure sympy); matrix-log the p-value search as archived.
- **E14 (SMALL):** su2 gauge-invariance committed tolerances (1e-8/1e-10) vs the paper's
  "<2e-15 machine precision" quote — tighten tests or scope the prose; hodge1 ±1-window
  test tightened to exact (sibling already exact); vacuous test_furry_theorem body
  (no-assert) — give it a real assert or remove.
- **E15 (SMALL):** P33 8/8 theorem statement lacks the n_max=2 scope its own §6 discloses
  (Ward fails at n_max=3). Add the scope to the theorem/abstract.
- **E16 (P36/P41/P2):** zero inline tests/ citations anywhere (C13 vacuous by omission) —
  add falsifier-test cites where tests exist after the backfills.

### F. NITs (fix-on-sight batch)
- P2:812-14 "cannot be reproduced" → "by a physically-motivated formula"; P2:117-19
  C6 phrasing sync ("converges to" not "shares the same spectrum"); P30:67-71
  least-action framing tone-down; P28:4339 "one observation"→two; P33:50 Hermiticity
  shorthand vs the remark's fuller mechanism; P28:4448 add the "coincidence" echo;
  synthesis: "graph-native 1.084"→continuum-spectral-sum label, "theorem tier"→
  proposition for depth-k, Link 2 "largely established", CP² 25% cross-ref to P24.

## Deterministic + state
Gates 7/7 PASS at freeze (v4.63.1). Matrix: zero group5 rows (to be populated in
remediation). Compiles 0-errors (P30+synthesis at freeze; others untouched pre-run).

## Calibration — round 2 (the recovery wave) → CLOSED

| Recovery agent | Seeds | Result |
|---|---|---|
| P28 (opus) | s-c28a/b | 2/2 detected (s-c28a MATERIAL w/ the empty-sum proof — "truncation headroom" physically false; s-c28b found + soundly graded down: the sibling `> 0` genuinely backstops it) |
| P51 (opus) | s-c51a/b | 2/2 (Bernoulli diff proved exactly 0; L6 rel_err measured 2e-16 vs the 1e-2 gate) |
| {25,30} (opus) | s-c2530a/b | 2/2 (proved "Haar sampling variance" fabricated — three fixed identity links; ±1-window off-by-one counterfactual) |
| {2,36,41} RECAL (sonnet, fresh seeds) | s-c23641c/d | 2/2 (computed the actual 8.815e-8; quantified the planted gate at ~11,345x looser) |

**FINAL CALIBRATION: 18/18 scoreable seeds caught (per-dimension after recovery); specificity 6/6;
s-c23641a/b retired as instrument-invalid (planted in the WH1 comparator, a non-backing file —
the agent's refusal to audit it was CORRECT).** Cross-catches: s-p1 x1, s-p3 x4, s-c33a/b x1.
Tier lesson: Sonnet-code missed 4 tolerance seeds under BROAD whole-paper prompts but the same
class was caught reliably under NARROW file-scoped assert-audit prompts (and by Sonnet itself in
the recal) — prompt scope, not model tier, was the dominant variable. qa.md follow-on: consider
splitting the code dimension into map (broad) + assert-audit (narrow) sub-passes.

## Recovery + critic + gap-closure genuine additions to the register

- **E17 (MATERIAL, P30/hodge1):** Bochner–Weitzenböck test is a tautology (mu_conn DEFINED as
  n(n+2)−2, then asserting the difference == 2 → `2==2`); the Ricci +2 shift has zero
  independent verification; degeneracy suite = formula echoes; transverse/longitudinal split
  never anchored to literature. FIX: independent connection-Laplacian computation or honest scoping.
- **E18 (MATERIAL, P28 tests):** `assert float(amp) != float('nan')` is unconditionally TRUE
  (IEEE NaN!=NaN) — the "exact CG amplitude" backing verifies nothing. FIX: proper finiteness +
  exact-value assert.
- **E19 (MATERIAL-class, P2/P36 tests):** vertex-parity "inheritance" asserted only as >0 (no
  no-rule baseline); the "k=1 does NOT" clause untested (verified true); verdict() says "MIXED"
  where paper/module say "WEAK" (uncovered taxonomy drift); §VIII.D→§VIII.E mis-cite; the
  "(and positive)" docstring vs actual 0.0. Plus the F₂ exact-rational prose backed exactly only
  for the zeros (nonzero values float-tier).
- **E20 (MATERIAL-LARGE, P28:1990-95, gap-closure):** the embedded gravity block explains the
  Newton-constant factor of 2 via the scalar-vs-Dirac cone coefficient — the EXACT reading
  CLAUDE.md §3 documents as WRONG (2026-05-30) and P51:637-662 explicitly rebuts (Wald-forced
  bookkeeping). Every shared NUMBER between the two copies matches (full diff table in the
  gap-closure report); the divergence is this one explanation + version-skew NITs (G6-FP
  status stale in P28; P51 6-12-vs-9-16-months internal skew; TT-spectrum relabeling; P51
  dim-H FP table internal inconsistency 16/40/80/140 vs formula 4/16/40/80 — verify).
- **F-adds:** P2:402-411 stale non-rendering "DROP-IN REWRITE" comment block carrying
  prohibited "conjectural"-K framing (delete); P28 n_max 8-vs-12 caption cosmetic; P51
  heat-kernel test comments invoke absent mechanisms (~46 OoM weak vs docstring).
- **Cleared by gap-closure (SOUND):** P2 §VIII S³-specificity + second-selection + tab:chain
  (K-tripwire exhaustive: rendered text CLEAN); the π³α³ 0.25% match IS properly
  coincidence-audited; P30's two tables arithmetically verified (all 6 char_coeff rows
  recomputed); P28 transcendental tables + weight-10 in-text arithmetic; tab:mode_dependence
  divergence honestly stated.

## VERDICT — **FAIL (fully calibrated)**

Every gating dimension exercised + calibrated (after the recovery wave). Genuine MATERIALs
across ALL dimensions (registers A–F + E17–E20 above). Specificity perfect. The findings
profile is first-cert-typical but deeper than group4's: structural backing debt (hardcoded
census mechanisms, tautological verifications, a reproducible production crash, a
test-hermeticity bug, phantom bibitems) rather than mechanical drift.

**Remediation plan (two tiers):**
- **Tier 1 — fix-on-sight block** (prose/citation/assert/label): registers A, B, C, D, F,
  E14/E15/E18/E20 + the mp.dps hermeticity fix (E6) + matrix population + gate hardenings
  (C13 bare-function-name; C16 propinquity broadening).
- **Tier 2 — named sprint-scale follow-ons** (PI-paced): E1 Furry derivation port; E2 census
  rules de-tautologize + renumber; E3/E4 the n_max≥4 crash fix + pendant sweep (E5); E7 P36
  Lamb-chain durable backing; E8 P41 XCWG durability migration; E9 P30 kinetic+MC tests;
  E10 P25 quotient/S⁵; E11 P51 headline ports + tests/gravity_support migration; E12 P51
  tier-tagging pass; E13 P2 B=42/circulant tests; E17 Bochner independent verification.
Then delta-verification cycle(s) per the run-shapes protocol → final certifying FULL run.


## Tier-1 remediation applied (2026-07-03, v4.64.0) — PI-confirmed

**Papers (all compile 0-errors):** P25 conjecture-idiom sweep (14 loci; rule→Observation,
physical reading→proposed interpretation; nine→twelve); P33 conjecture ×2 + abstract six-rule
list (Gaunt→triangle-on-n) + charge-conjugation→Furry harmonized with P28 + 8/8 census scoped
to n_max=2 + fabricated cite → real test + 99→100 tests + jeffrey→Jeffrey-Dai-2008 +
chamseddine vocabulary-not-theorem reword; P36 conclusion re-synced to §5.2's corrected
reading (Table-7.4 misreading; LS-7 done-WEAK) + wrong module cite fixed + Drake–Swainson real
title + 4 phantom bibitems added & wired (Antognini/Karshenboim/CODATA/PDG) + Schwartz scope;
P28 factor-of-2 → the Wald-forced reading (the §3-documented-WRONG cone attribution removed) +
parker1980 re-pointed (textbook §6 primary; 1984 PRD = background) + petermann sign +197/144 +
3 propinquity zombies → state-space GH + G6-FP status updated + pendant verification claim
honestly scoped + n_max 8→12 + one-obs→two + coincidence echo; P51 propinquity ×2 → state-space
GH (the Lemma honestly deferred) + sommerfeld 1896 + frolov venue PRD 56 2212 + Dowker-1977 &
Beccaria–Tseytlin bibitems added & wired + 13 orphans removed + dim-H convention note
(shell-index offset resolved) ; P30 conjectural-framework → observation-tier program +
least-action softened + Fegan/Menotti–Onofri cites fixed + 4 orphans removed; P2 nine→twelve
×2 + stale DROP-IN comment block deleted + C6 phrasing (converges-to) + p-value phrasing +
perez_sanchez 2024/2025 unbundled (both P2+P25) + §IV source-list honesty (obstruction
identities = archived drivers, not the production module); synthesis 4 NITs (spectral-sum
label, proposition tier, Link-2 largely, CP²-residual → Paper 25 §VII provenance) + P33
pendant claim scoped.

**Tests/code:** E6 mp.dps hermeticity FIXED (autouse set/restore fixtures ×4 files; the
original 11-failure SCOPE order now 624 passed/0 failed); E18 NaN always-true assert → real
isfinite; E19: k=1 converse now asserted, vertex-parity inheritance asserted at the WEIGHT
level (the guard-level monkeypatch provably changes nothing — the rule is redundantly enforced
via SO(4) channel weights; documented in-test), verdict() MIXED band → the paper's WEAK
vocabulary (+test synced), §VIII.D→§VIII.E, "(and positive)"→exact-0.0 docstring; the
vacuous furry test → asserts the documented scalar-Furry FAILURE; new
test_scalar_l1_self_energy_is_one_over_4pi (exact, 1e-14); **E4 production bug FIXED**
(graph_qed_photon eigenvalue sort key → re(N(x)); all 277 graph_qed consumers green) — the
n_max≥3 pendant sweep verification remains a Tier-2 item (exact eigenvals are slow); P51
zeta heat-kernel docstring aligned with its assert strength.

**Gate hardenings:** C13 now resolves bare \texttt{test\_fn} FUNCTION citations (the P33
fabricated-cite class; file-name fallback added after it caught P32's informal
test_real_structure cite — all 6 scopes PASS); C16 gained the fail-severity group5
propinquity entry (the line-break class that hid P51:159/:2233).

**Verification:** 7/7 group5 gates PASS post-remediation; 8 papers + synthesis compile
0-errors; all touched suites green (vector_qed 100, paper2+two_loop 23, qed-family 624 in
the formerly-failing order, graph_qed 277, zeta 40).

**Tier-2 named follow-ons (logged in group5.done.md):** Furry derivation port; census rules
de-tautologize + renumber (paper Table-1 6↔7 vs code keys) + triangle-on-SO(4) check + census
aggregate tests (1/8, 4/8, ==7); pendant n_max≥3 exact sweep (post-crash-fix); P36 Lamb-chain
durable backing; P41 XCWG durability migration (tests/wilson_rule_b_support/); P51
tier-tagging pass + a₀/SC/K_cone ports + j_blindness debug-import migration
(tests/gravity_support/); P30 kinetic-1/8 + seeded-MC tests; P25 S²-quotient-spectrum +
S⁵ data pins; Bochner–Weitzenböck independent verification; P2 B=42/circulant exact tests;
F₂ nonzero exact-rational pins; su2 gauge-invariance tolerance-vs-prose reconciliation.

Next per the run-shapes protocol: **delta-verification run** over this remediation diff,
then the final certifying FULL run.


## Delta-verification #1 (2026-07-03, v4.64.1) — CLEAN-DELTA → the certifying FULL run is EARNED

Five payload-pinned reviewers over the v4.64.0 remediation diff (a protocol refinement:
agents read ONE patch file each from the scratchpad — paste-don't-point economics without
routing 1,100 diff lines through orchestration). **Sensitivity 7/7** (every seed caught by
its assigned agent; claims-2 additionally cross-caught both citation seeds): s-d1 P25
"en route to a derivation" tier creep (3 sibling hunks quoted against it); s-d2 the P36
component swap (caught via the must-be-the-smaller-piece reasoning); s-d3 Link-2
"fully"; s-d4 the >=0 anti-vacuity tautology; s-d5 the dropped fixture restore (sibling
comparison + disk cross-check); s-d6/s-d7 the two volume digits (DOI-verified).
**Specificity clean** — all genuine hunks verified sound across five agents; the citation
agent confirmed 12/12 new bibliographic facts GROUNDED, several as corrections of real
prior errors (Frolov–Fursaev venue, Petermann sign, dropped authors, the Pérez-Sánchez
conflation). **Zero genuine MATERIAL in the diff.** Post-run polish: the FP/J-blindness
slash NIT tightened (P28); the parker locator honestly widened to Chs. 6–7 (the g−2
formalism may span the gauge chapter); fegan1983's content verified as the real 1978
Trans. AMS paper (key misnomer, harmless). Deterministic 7/7 PASS; P28 recompiles 0.
Cost ≈510k subagent tokens across 5 agents. **Next: the final certifying FULL run —
PI's timing call, weighing the open Tier-2 backing sprints.**
