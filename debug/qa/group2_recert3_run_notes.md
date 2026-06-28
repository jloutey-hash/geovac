# /qa group2 — re-cert run #3 (2026-06-27): FAIL → remediating

Third `/qa group2` invocation (confirmation re-cert of the v4.50.1 corpus, HEAD e95eec9).
Runs #1 (v4.50.0) and #2 (v4.50.1) both FAIL→remediated; run #3 peeled the **cross-paper
"second-locus" propagation layer** (the documented group3/group1 pattern: each correction
fixed at its primary locus, echoes elsewhere left stale).

## Verdict: FAIL (trustworthy)

Calibration **10/11** (cs3/P17 unscored — that code-reviewer over-ran its slow composed-PES
pipelines; worktree torn down before it returned), specificity **6/6**. The rs8-class fix
WORKED: cs5 (FCI-M vacuous `E_lih<100` binding guard, the exact class that slipped run #2)
was CAUGHT by the sharpened per-assertion code prompt.

### Per-dimension
| Dimension | Exercised | Calibrated | Genuine material | Verdict |
|:--|:--|:--|:--|:--|
| Code C1/C2 | 8/9 (P17 over-ran) | 4/5 (cs1/cs2/cs4/cs5; cs3 unscored) | none (NITs only) | clean-of-material |
| Claims C3/5/6/8 | yes | 3/3 (ps1/ps2/ps3) | YES | FAIL |
| Citations C4 | yes | 2/2 (ts1/ts2) | none (NITs only) | clean-of-material |
| Synthesis C9 | yes | 1/1 (ss1) | YES | FAIL |
| Deterministic C10–16 | yes | n/a | clean | PASS |

Seeds (all fresh loci/classes vs #1/#2): cs1 P12 +100 ✓, cs2 P15 1e10 ✓, cs3 P17 0..100 range
(unscored), cs4 FCI-A 1e3 ✓(double), cs5 FCI-M `<100` ✓(rs8-fix), ps1 P12 no-quadrature ✓,
ps2 P15 "matches exact" ✓, ps3 FCI-M guardrail-flip ✓(also code bonus), ts1 P13 Pekeris vol ✓,
ts2 FCI-A Hehre vol ✓, ss1 synthesis "viable replacement" ✓.

One reviewer FALSE-POSITIVE rejected: synthesis 4N "(unbound D_e)" is the CORRECT 2D-variational
result (Track AR / CLAUDE §5); that reviewer read the adiabatic artifact.

## Genuine material findings (verified, non-seed) — all clearly-correct
1. **GENUINE FALSE CLAIM** — paper_fci_atoms ~380: "Li FCI below the HF limit 0.6% / recovers
   correlation" is FALSE (Li 1.07% > HF 0.61%); only He recovers correlation. [LARGE]
2. **Synthesis-lag (run-#2 corrections not propagated to synthesis):** cusp 0.004% labeled
   "properly variational" (P13 now: non-variational); Sturmian "H∝S" single-term (P8 now:
   two-term, not ∝S; SO(4) congruence). [LARGE + SMALL]
3. **Cross-paper "second-locus" propagation gaps** (completeness-critic): H∝S echoes in
   P11/FCI-M/FCI-A-diag; HeH⁺ 93.1% echoes in P15-body+conclusion + P17-intro + synthesis;
   cusp non-variational in P13 tables; FCI-A LiH no-equilibrium in abstract+conclusion;
   H2O 26%→19.4% stale in P19 abstract+intro; H2+ 0.0002%-vs-0.70% + LiH 0.093-vs-0.110 +
   4N 63.5/64/67% cross-paper inconsistencies. [broad, mostly SMALL]
4. P12 "no/without any numerical quadrature" residual loci (841, 307); FCI-A "Both surpass
   STO-3G" (Li STO-3G uncomputed); FCI-M abstract Fock-λ fitted-scale not disclosed. [SMALL]

## Remediation — DONE + VERIFIED (released v4.50.2)
One opus propagation-sweep agent applied all 12 canonical corrections to ~25 stale loci across
10 files (per-file old→new logged in CHANGELOG v4.50.2). **Verification:** deterministic
C10–C16 PASS post-sweep; the genuine Li-beats-HF fix confirmed correct (He error 0.35% < HF
1.4% → recovers correlation; Li 1.10–1.15% > HF 0.61% → does NOT, unlike He); all 10 touched
papers compile content-err=0 (environmental missing-figure only), 0 new undefined refs; no
code/tests touched. The sweep left sound loci alone with documented reasons (P15 already-caveated
HeH⁺ paragraph; P17 body 64/67% already adiabatic-tagged + corrected by Note v2.0.24; P11 own
FD results table is its legitimate FD primary result).

## Honest ceiling / next
- cs3 unscored (P17 over-ran) — code dimension calibrated by 4/5 (cs1/cs2/cs4/cs5 incl. rs8-fix);
  next run: keep the sharpened per-assertion code prompt + let P17's pipelines run with a longer budget.
- After remediation: release (v4.50.2 patch) → re-run /qa group2 (run #4) for the certified PASS.
- Standing pattern across 3 runs: prose/restorations (#1) → framing/§3-reassertion+P8-theorem (#2)
  → cross-paper second-locus propagation (#3). Each run peels one layer.
