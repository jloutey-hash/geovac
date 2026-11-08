"""Record group2 baseline FULL run - Batch 4 (Paper 58 + synthesis C9 +
completeness-critic) + baseline completion. v5.11.17. Idempotent."""
from __future__ import annotations
import sys, pathlib

CL = "CHANGELOG.md"
ANCHOR = "## [v5.11.16] - 2026-09-13\n"
MARKER = "## [v5.11.17]"
ENTRY = """## [v5.11.17] - 2026-09-13

**`/qa group2` baseline FULL run - Batch 4 (Paper 58 + group2 synthesis C9 + completeness-critic) = FAIL, remediated. The group2 baseline (4 batches) is COMPLETE.** Five reviewers, tree frozen. The LARGE was a stale-echo the run's own completeness pass caught in *this session's Batch-3 work*.

### Paper 58 - CLEAN on code and citations, one SMALL on claims

- **Code (all backing tests RUN):** no LARGE, no false-positive, no tautology. 25 slow legs + the non-slow suite all pass; the certificate is mutation-guarded, the (AA|BB) decider fire-tested both directions, the decompactification front cross-checked by an independent prolate quadrature. The three DoD MATERIAL items (blind decider, QFD precision, C22 syspath) are verifiably remediated.
- **Citations:** exceptionally strong - 25+ external citations individually source-verified (Curie/Bethe/von Neumann-Wigner mechanism, Shibuya-Wulfman, Ruedenberg exchange, Avery Sturmians, DLMF S30.3, Herring, two-projections theory), no WRONG-ID, no orphan bibitems, no fabricated IDs.
- **Claims (F1, MATERIAL/SMALL, fixed):** the polyatomic *scope* paragraph (S I.A) labelled the FCI anchor figures (+10.0%/+14.9%, matched to composed FCI 11.7%/19.4%) as "RHF" and spliced them onto the RHF-ladder endpoints (+1.4%/+2.5%) as one basis-driven curve. Per the backing memo the RHF ladder is 7.38->1.40% (BeH2) and 8.77->2.46% (H2O), and minimal-basis FCI is *worse* than minimal-basis RHF. Reframed: the FCI-vs-FCI margin and the RHF basis ladder are now separate, correctly labelled, and the "basis-limited not method-limited" claim is softened (RHF is not at its correlation limit). The code reviewer independently corroborated the muddle (SMALL-4). Does not propagate to the synthesis.

### The LARGE - a Batch-3 correction that was locus-incomplete (C9 synthesis + completeness-critic)

Batch 3 (v5.11.16) re-measured the FCI-atoms energies but the correction reached only the **authoritative surfaces** (abstract, Table I, convergence-detail table) - **not** the conclusion or the several prose echoes. The C9 synthesis reviewer and the completeness-critic together caught the retired pre-ERI-fix errors (He 0.35% / Li 1.07% / Be 0.90%) still live at: the group2 **synthesis** graph-native-FCI paragraph; the FCI-atoms **conclusion**; the He and **Li convergence sequences** in the figure discussion (the Li sequence was in *no* reviewer's line list - caught by a systematic re-grep, the exact locus-incompleteness the Summary-Surface Reading Rule exists to prevent); the hybrid-vs-exact-h1 comparison; the "remaining basis error" paragraph; and the graph-native comparison. **This is the Summary-Surface Reading Rule failure inside my own prior batch.** Swept exhaustively this time (re-grep confirms zero retired FCI-atoms values remain):

- He hybrid 0.56/0.45/0.38/0.35 -> **0.50/0.37/0.29/0.26%**; Li exact-h1 5.04/1.15/1.10 -> **5.00/1.12/1.06%**, n>=3 range 1.10-1.15 -> **1.06-1.12%**; He-vs-exact-h1 at n=4 0.35/2.08 -> **0.29/1.99%** (He n4 exact-h1 re-measured); conclusion + remaining-error He/Li/Be -> **0.26/1.03/0.71%**; synthesis line likewise.
- **The graph-native comparison needed a REFRAME, not a number swap.** With the grid (hybrid) He now at 0.26% (Table I) and graph-native at 0.25%, they **agree to 0.01 percentage points** - the old "0.35% grid vs 0.25% analytical, analytical more accurate" 0.10-point gap was the wrong-sign-q **grid bug**, not a real accuracy difference. Post-fix the two routes confirm each other rather than one beating the other.

### Completeness-critic - GAPS, all now closed

- **GAP-3 (the LARGE above):** FCI-atoms propagation - closed by the exhaustive sweep.
- **GAP-2:** the group2 DoD (`group2.done.md`) still ratified the superseded C8 literals "Be 0.90% / Li 1.07%" (the ".done.md ratifying a retired value" class) - updated to 0.71% / 1.03%.
- **GAP-1** (Paper 58 Batch-4 coverage "not evidenced") was a **false gap**: it ran before this record existed; Paper 58 received three fresh Batch-4 reviewers (code + claims + citations).
- **Deterministic layer 8/8 PASS** (C10 compile 13 papers, C13/C14/C16/C17/C19/C21/C22) - confirmed by the critic and re-confirmed post-remediation.

### Also

claim-matrix: FCI-A He-hybrid row 0.35%->0.26% (NO-TEST but reproduces this session); Paper-58 QFD row "84-digit"->60-digit (matching the paper + the row's own note); two stale test docstrings (`test_paper58_qfd` "84-digit"->60, `test_paper58_decompactification_front` header t_c 2.7456->2.664 - asserts were already correct).

### Carried to cert (declared NITs / OWED, none silently dropped)

Paper 58: the (AA|BB) "0 Gaunt-zeros is basis-forced" disclosure (code SMALL-1, not a false positive); the S/h census DECIDED-vs-MEASURED tier upgrade (earned but defensibly conservative, logged since 2026-08-17); LiH 30-digit + polyatomic scope numbers + many-electron front residuals are driver-level OWED; Clementi inline reference; smith1972 geometric-mean framing. Corpus-wide owed: the He 0.19%@n_max=7 and P15 96.0% in-paper extrapolation footnotes (Batch 2/3); herbst2018 orphan + uncited -8.071 (P19); the earlier per-batch NIT ledgers (all tracked in v5.11.14-16).

### Verdict

**Batch 4: FAIL, remediated** - Paper 58 clean (code + citations), 1 SMALL (F1); the LARGE was a locus-incomplete Batch-3 correction, now swept exhaustively; the completeness GAPs closed. **group2 baseline COMPLETE across four batches (v5.11.14-17): all 13 documents covered, every batch found and fixed defects (guardrail scope; KW literal; H2O reconcile; a P17 self-contradiction; the FCI-atoms LARGE + pair-diagonal zombie; and this run's incomplete-propagation catch), deterministic layer green, all papers compile.** It is a baseline re-measure, not a certification (58/59/60 in by reference; purpose per the DoD is to re-measure the whole group's surface). Clean deltas on 58/59/60 and the owed footnotes remain before any group2 cert.

"""


def main() -> int:
    p = pathlib.Path(CL); t = p.read_text(encoding="utf-8")
    if MARKER in t:
        print("skip"); return 0
    if t.count(ANCHOR) != 1:
        print(f"MISS anchor count={t.count(ANCHOR)}"); return 3
    p.write_text(t.replace(ANCHOR, ENTRY + ANCHOR), encoding="utf-8")
    print("ok CHANGELOG v5.11.17")
    c = pathlib.Path("CLAUDE.md"); ct = c.read_text(encoding="utf-8")
    ct = ct.replace("**Version:** v5.11.16 (September 13, 2026)",
                    "**Version:** v5.11.17 (September 13, 2026)", 1)
    # add a single group2-baseline one-liner to S2 (baseline now complete)
    s2_anchor = "- **/qa paper_61 DELTA = DEFECTS, remediated (2026-09-13, v5.11.12):**"
    s2_bullet = ("- **/qa group2 baseline FULL run COMPLETE (2026-09-13, v5.11.14-17):** 4 batches, "
                 "13 docs; every batch found+fixed defects (guardrail scope; pair-diagonal zombie; "
                 "FCI-atoms stale energies re-measured; a locus-incomplete same-session fix caught by "
                 "the completeness pass). Baseline re-measure, not a cert. See CHANGELOG v5.11.14-17.\n")
    if s2_bullet.split(":**")[0] not in ct and s2_anchor in ct:
        ct = ct.replace(s2_anchor, s2_bullet + s2_anchor, 1)
        print("CLAUDE.md S2 group2-baseline bullet added")
    c.write_text(ct, encoding="utf-8"); print("CLAUDE.md v5.11.17")
    return 0


if __name__ == "__main__":
    sys.exit(main())
