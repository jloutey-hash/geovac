# Sprint: /qa certifying gate — the compounded v4.85–v4.106 delta re-review + remediation — 2026-08-21

**Origin:** PI fired `/qa` on the standing compounded "Phase-4 re-review OWED" (certified papers 58, 59,
60, 19 [group2] + 56, 24 [group3] + the group3 synthesis delta). Run shape: **DELTA-verification**
(the owed accumulated diffs; a full certifying run remains PI-gated on this clean delta).
Canonical memo; seed keys `debug/qa/delta85_seed_key.json` + `delta85r_seed_key.json`; both seeded
worktrees destroyed, zero leakage verified.

## Round 1 — the delta panel (8 dispatches, all dimensions)

**Calibration: sensitivity 10/10 planted seeds caught; specificity: zero manufactured faults**
(one designated control — the "84-digit" claim — was flagged and on verification the REVIEWER was
right: the control designation was flawed; counted as a genuine catch, not a false positive).
Two run-infrastructure errors (mine) corrected mid-run: citation diffs pre-dated the seeds
(dimension re-exercised against the worktree bibliographies → 2/2 seeds caught); `debug/` +
CHANGELOG not synced into the worktree (three "absent backing" LARGEs deflated to genuine cores).

**Verdict: DEFECTS** — ~22 genuine MATERIALs after seed/artifact removal. Headliners:
- P59 abstract+conclusion stated FALSE mathematics ("no invariant subspace" for a unipotent
  monodromy — must be "no invariant *coordinate* subspace").
- P19's withheld n_max=4 energy CROSSES exact (−7.9325→−8.0576→−8.1048 vs exact −8.0705; extrap
  ≈−8.13): the energy also converges to a wrong limit — disclosed; "converges excellently" scoped
  to n_max≤3 (4 loci).
- P58's 84-digit claim exceeded the artifact's own capping convention (→60) and its honest-scope
  understated the H₂ gap 2.5× (0.027→0.068 Ha, 0.027 scoped to ζ=1.197).
- Synthesis elliptic-layer paragraphs = pre-deflation zombies (central-fibre τ=i; meet-only-at-CM;
  "disjoint") → fully rewritten to the current family-level/X₀(2)/T-2-negative state.
- Six unretired-superseded-clause instances across P59/P56/P24 (+ seven-layer propagation to
  INDEX/P31/synthesis); citation SMALLs (V. K. Murty initials; harris_michels page range).
- Tracked-backing gaps: T2's in-suite pin was ~11 digits → `geovac/t2_kw.py` promoted (the (KW)+
  Watson-tail evaluator; executable ~21-digit witness `test_kw_mpmath_witness` excluding the old
  anchor); σ-law tie test to the paper's measured 30.1/15.4/4.4; Sugiura-1927 + `exchange_hp` legs;
  the −4.26 mHa phase-magnitude pin; cross-dps guard; tolerances tightened; C17 family
  `t2-collinear-anchor` added (caught a stale matrix row within seconds); C18 same-day pattern.
- Artifact: H₂ (60-digit, convention-stated) + LiH (30-digit, tail-bounded) rows added → 66 entries.
- DoD paper_59.done.md refreshed (66-digit canon; the intersection-form goalpost UPGRADED to
  [SYMBOLIC] per its own v4.97.0 closure — a two-way goalpost fix).

## Round 2 — the remediation-delta (2 verifiers, 4 fresh seeds, 4 controls)

**Claims verifier: calibrated (2/2 seeds: the reverted six-layer footnote; the "proven to be"
tier regression), 4/4 controls clean.** 37/41 checklist items FIXED first pass; residuals fixed:
"provably lives" tier slip; TWO raw-TAB corruptions my own heredoc scripts had injected into P59
(\to and \textbf eaten — the JSON/heredoc escape trap, now a named lesson); two more "converges
excellently" loci; the 4th P31 six-layer locus; P58 atom-list harmonization; DoD residuals.
**Code verifier: calibrated (2/2 seeds: the defanged 1e-12 witness tolerance; the self-comparing
tie test), 8/8 genuine items DELIVERED.** Genuine residual fixed: the same-day pattern now has
selftest coverage (17 pos/12 neg).

**Bonus catches from the corruption scan:** a PRE-EXISTING backspace corruption in C17's
pauli-advantage-floor regex (a silently-dead guard since it was written — repaired); the repaired
C18 same-day pattern immediately surfaced 2 more live loci (Papers 32, 17 — fixed).

## Final state

**REMEDIATION-DELTA: CLEAN.** All deterministic gates PASS (incl. the two repaired registries);
74 fast + slow suites green; all nine touched papers compile; both worktrees destroyed with zero
seed leakage. **The full certifying run (the only shape that can emit PASS) is now unlocked —
PI-gated.**

## Process lessons (for qa.md / future runs)
1. **Sync debug/ + CHANGELOG into seed worktrees** (three false-LARGEs and a V1 this run).
2. **Generate diff files AFTER planting seeds**, or point citation reviewers at the worktree
   files, never the diffs (the round-1 citation de-calibration).
3. **Never edit corpus text via bash-heredoc python** — the JSON+escape double-decode injects
   control characters (\t, \b, \r); use Write-tool scripts with raw strings (two TAB corruptions
   + my C18 pattern's backspaces all came from this; C17's pre-existing corruption suggests an
   earlier session hit the same trap).
4. Pre-registered digit thresholds must be calibrated against ring DIMENSION (10^(D/n)).
