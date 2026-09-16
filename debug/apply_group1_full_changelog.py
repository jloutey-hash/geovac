r"""CHANGELOG + version + group1.done.md record for the group1 FULL cert run.

Patch (v5.12.2): paper corrections + a code guard. Flagged in the session
summary as possibly minor (the P53 Theorem 5.6 rescope changes a published
claim's tier) for the PI to decide -- not bumped unilaterally, and the /qa
skill file itself was not touched.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

ENTRY = """## [v5.12.2] - 2026-09-15

**`/qa group1 full` (PI-invoked certifying run) = FAIL, remediated. NOT a PASS.** The first whole-group FULL run of group1 since the 2026-06-24 bite-wise cert; live scope was the 10 papers left after the 46-49 archival (29, 39, 40, 42, 43, 44, 45, 50, 52, 53) + synthesis. 14 deterministic gates + a 10-agent LLM panel + a completeness-critic + 3 critic-driven gap re-dispatches. Verdict FAIL on verified material defects, all SMALL; the one LARGE dissolved under verification.

### The LARGE dissolved into a decomposition bug

The claims panel flagged a LARGE: Paper 40's Cor L3_closure / Lem L3_interior claim the Dirac-triangle inequality "rigorous at all ranks, every summand", contradicted by a documented counterexample in the backing (matrix row 214: "G2 (1,0)v(0,4) ratio ~2.4"). The P40 **code** reviewer, re-deriving independently, found the 2.4 is an artifact of `dirac_triangle_extended_verify.py`'s `tensor_product`, which is **not dimension-conserving outside its validated panel** -- it reports a Schur-impossible trivial summand for that pair. **The PM reproduced it**: (1,0)x(0,4) sums to 5236 != 4662; the dimension-correct DT value is < 1 (holds). No counterexample exists. Two reviewers converged on one locus with opposite readings and the code re-derivation settled it -- the exact discrimination the panel exists for.

### Verified material defects (all SMALL), remediated

- **A false counterexample in the backing.** Matrix row 214 AND the slow-test docstring cited the bug's "ratio 2.4" as real evidence justifying the panel restriction. Withdrawn; corrected to "panel-bounded by the driver's validated range", and a dimension-conservation guard added to `run_panel` (fire-tested: raises on the buggy pair, silent on clean panels) so a panel-widening can never again pass a mis-decomposition off as a result.
- **Retired "propinquity" labels** for GeoVac's own Paper 38/39/40 results, which are state-space GH: P42 section heading (contradicted its own body), P44 literature bucket, two synthesis loci.
- **Un-caveated pre-descope framing.** P43's intro "the Lorentzian content ... literally satisfied at (3,1), not merely a structural correspondence" -- signature-blind caveat added to match the abstract and Paper 42 Sec.10. P53's Theorem 5.6 asserted the strictly stronger Latremoliere propinquity via max(reach, height)->0 while its own remark says the height leg (Lebesgue constant 2.01, Lambda-independent) does not vanish; scoped to state-space GH with the Latremoliere form marked conditional.
- **Two synthesis clauses** in the archived-Paper-47 subsection describing descoped/degenerate content in present-tense "established" language -- both caveated to the surviving spatial-rate-formula reading.
- **P40 group corollaries** now carry the inline tier (rate rigorous rank-1, numerically pinned rank>=2; convergence conditional for general G); the P39 bib annotation of P40 gets the same rate caveat.
- **Toyota M. -> R.** (P42, P43 bibitems; the only citation defect on ~290 bibitems enumerated). **P53 stale bib titles** for the archived 45/46 corrected. **P43 sub-epsilon residual** given its one-line provenance caveat (degenerate scalar wedge unitary -> eps^2).

### The headlines all held, by independent re-derivation

Every load-bearing result was re-derived by a reviewer's own route: P50's F-theorem (PSLQ-recovered from the framework spectrum; both KPS sides fire-tested), P45's K+ annihilation (a genuine Krein-forced identity, the 15.87-norm spatial Dirac projected out; not a restriction artifact), P53's height 2.01 (sharp L1 norm, re-derived by a second route), P40's rank-1 4/pi (to 1e-8), P44's prop=2 (SVD-rank, guards reject k=1/k=3/hardcode), P29's 84->80 arithmetic (bit-for-bit). The C4 citation surface -- the branch's known fabrication history -- is clean end to end, both load-bearing inline theorem numbers (van den Dungen Prop 4.1, Nieuviarts Def 2.2) confirmed against primary sources, and the highest-risk 2026-dated arXiv ID resolves exactly.

### The completeness-critic earned the run

It caught a PM scoping error: Paper 39 fell out of all three claims chunks and both citation chunks, and Paper 44's code went unaudited because its test is not named `test_paper44_*`. Three focused re-dispatches (P39 claims, P39 citations, P44 code) all returned CLEAN, closing the gaps -- without them three gating dimensions were unexercised (INCONCLUSIVE, not clean).

### Owed to the PI (two items)

1. **Paper 40's "rigorous at all ranks" prose is not edited.** No counterexample exists, but the interior-summand closure (Lem L3_interior) is an analytical Steinberg/Brauer-Klimyk argument the code cannot verify; the honest computational tier is panel-verified + the Kumar-PRV asymptotic bound. Whether the analytical lemma certifies the all-ranks claim is a primary-math call. If it holds, only the (now-corrected) backing note was wrong; if it has a gap, the paper prose should be scoped to panel-verified.
2. **Paper 50's S7 "DONE" catalogue row has no backing test** and no external check is possible (framework-generated, and an S7 erratum history is on record). A regression test is the only pin; logged as a coverage gap.

Group1 is **not certified** by this run -- a FULL run that FAILs and is remediated needs a clean DELTA before a certifying FULL can PASS. Re-cert owed.

"""

CL = "CHANGELOG.md"
with io.open(CL, encoding="utf-8") as fh:
    t = fh.read()
anchor = "## [v5.12.1] - 2026-09-14"
if anchor not in t:
    print("FAILED changelog anchor"); sys.exit(1)
with io.open(CL, "w", encoding="utf-8") as fh:
    fh.write(t.replace(anchor, ENTRY + anchor, 1))
print("  + CHANGELOG v5.12.2")

CM = "CLAUDE.md"
with io.open(CM, encoding="utf-8") as fh:
    c = fh.read()
c = c.replace("**Version:** v5.12.1 (September 14, 2026)",
              "**Version:** v5.12.2 (September 15, 2026)", 1)
bullet = "- **/qa DELTA on the group1 archive = DEFECTS, remediated (2026-09-14, v5.12.1):**"
new = ("- **/qa group1 FULL certifying run = FAIL, remediated (2026-09-15, "
       "v5.12.2):** the one LARGE dissolved into a decomposition-driver bug "
       "(the '2.4' DT counterexample was a Schur-impossible mis-decomposition); "
       "all real defects SMALL. Headlines re-derived sound; C4 clean. Critic "
       "caught a P39/P44 coverage hole. 2 PI items. NOT certified. See CHANGELOG "
       "v5.12.2.\n")
if bullet not in c:
    print("FAILED sec2 anchor"); sys.exit(1)
c = c.replace(bullet, new + bullet, 1)
with io.open(CM, "w", encoding="utf-8") as fh:
    fh.write(c)
print("  + CLAUDE.md version + Sec. 2 one-liner")

# group1.done.md: record the run outcome (do NOT stamp CERTIFIED)
DONE = "docs/qa/group1.done.md"
with io.open(DONE, encoding="utf-8") as fh:
    d = fh.read()
stamp = """- 2026-09-15 — **WHOLE-GROUP FULL certifying run = FAIL -> REMEDIATED**
  (v5.12.2), the first whole-group FULL since the 2026-06-24 bite-wise cert.
  Live scope = 10 papers after the 46-49 archival + synthesis. 10-agent panel +
  completeness-critic + 3 gap re-dispatches (all CLEAN). One LARGE dissolved
  under verification into a `dirac_triangle_extended_verify.py` decomposition
  bug (the "G2 (1,0)v(0,4) ratio 2.4" is a Schur-impossible mis-decomposition;
  dimension-correct DT < 1). All real defects SMALL, remediated: false-2.4
  backing note withdrawn + dimension-conservation guard (fire-tested); retired
  "propinquity" labels -> state-space GH; P43 signature-blind caveat; P53
  Thm 5.6 scoped to state-space GH (Latremoliere form conditional); 2 synthesis
  present-tense-descoped clauses; P40 corollary tiers; Toyota M.->R.; P53 stale
  bib titles. Headlines all re-derived SOUND; C4 clean. **NOT certified** -- a
  FAILed+remediated FULL needs a clean DELTA before a certifying FULL can PASS.
  PI items: (1) P40 Lem L3_interior "rigorous at all ranks" analytical tier
  (no counterexample, but code can't verify -- primary-math call); (2) P50 S7
  "DONE" untested (no external check possible). Memo: CHANGELOG v5.12.2.
"""
marker = "## Change log\n"
if marker in d:
    d = d.replace(marker, marker + stamp, 1)
    with io.open(DONE, "w", encoding="utf-8") as fh:
        fh.write(d)
    print("  + docs/qa/group1.done.md run record (NOT certified)")
else:
    print("  ! group1.done.md change-log marker not found; record NOT written")
