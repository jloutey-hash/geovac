r"""CHANGELOG + version for the group1 archive DELTA and its remediation.

Patch (v5.12.1): this touched a gate's IMPLEMENTATION (guard coverage) and two
paper loci, but did not add, remove, or rescope any C-criterion and did not
touch .claude/commands/qa.md -- so the standing "qa-skill change => minor" rule
does not fire.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

CL = "CHANGELOG.md"
CM = "CLAUDE.md"

ENTRY = """## [v5.12.1] - 2026-09-14

**`/qa` DELTA on the group1 archive (v5.12.0) = DEFECTS, remediated.** First verification of the Lorentzian-tail archiving. Claim-impact scope: the four archived papers plus the new C24 gate. Two reviewers (both Opus), unseeded, against the standing calibration record (trunk DELTA #1/#2, 17/18 seeds, 0/8 false positives). All deterministic gates green whole-group at close.

### Per-dimension scorecard

| Dimension | Exercised | Result |
|---|---|---|
| Deterministic (C10-C24, group1) | yes, whole-group | CLEAN |
| Claims / status (C14) | yes, 7 live citers of the archived papers | 2 material, remediated |
| Code (the C24 gate + follow-ons) | yes | 2 material-small + 3 nit, remediated |

### Claims: two summary surfaces credited a descoped claim as achieved

Both are pre-existing staleness the archiving pass did not sweep, and both are the paraphrase class C16 is structurally blind to — a citer restating a descoped claim in its own words, which is exactly why the claim-impact reviewer and not the phrase gate found them.

- **Field guide** said "Paper 49 closes the strong-form Krein-MS bridge (Q1')" with no decomposition, while the group1 synthesis makes the identical claim correctly by leading with "the Lambda-inheritance ... claims ... are descoped ... the cocycle-deficit algebra and the OSLPLS category design survive." Rewritten to the decomposed form: a construction at the categorical/algebra level, metric-level Lambda-inheritance descoped.
- **Paper 50's** "Place in the series" retrospective listed the arc as having "established ... convergence theory at ... strong-form Lorentzian ... via OSLPLS", contradicting the same paper's own descope-aware introduction. The descoped levels are now flagged as descoped and not convergence results.

No C16 entry was added for either: "strong-form Lorentzian" appears in dozens of correct descope-aware sentences ("... is descoped"), so any pattern broad enough to catch the crediting form fires on the correct form too, which the discrimination rule forbids. This defect class is owned by the claim-impact reviewer by design.

**The clean half, reported with equal weight.** The reviewer enumerated every citation of the four archived papers across all seven live documents. Papers 43, 45, 52, 53 and the group1 synthesis body are descope-accurate throughout — Paper 45 in particular carries "the Paper 49 cocycle machinery, which survives the descope" correctly. Both archive-note edits from v5.12.0 were verified **faithful**: the four surviving pieces map one-to-one to the register. No live document describes 46-49 as current or forthcoming.

### Code: the gate's own regression net had the holes it was built to close, one level up

The gating check (A, register integrity) was independently confirmed sound — it anchors on the filesystem glob, not the register's self-report, so no archived paper can escape it, the exact opposite of the C11-could-not-fail failure. The gaps were in the guards around it:

- **The ratchet's known/fresh partition lived only in `main()`**, untested; a future edit inverting the set-difference would silently stop surfacing NEW re-derivation loci. Extracted to a pure `partition_probe_hits` helper with three unit tests, one fire-tested by inverting the diff.
- **The selftest advertised "all five probes fire" while exercising 3 of 5 check-A branches.** The nonexistent-file and empty-reason branches fired but were unguarded. Added both as selftest probes and mirror-test assertions, both fire-tested.
- The page builder would wipe every page on an empty-but-valid manifest; guarded on a non-empty manifest. A loop variable shadowing the `page()` helper was renamed.

Every new guard was fire-tested against the specific wrong answer it excludes, in a pass separate from the fix, per the guard-writing rule. Mirror test 10 -> 15 assertions.

### Verdict

DELTA = DEFECTS, remediated. Not a clean delta this run; a re-run would be the precondition for a group1 FULL certifying pass, which is separately owed (the group1 cert record is stale to 2026-06-24). The archiving itself is faithfully reflected across the corpus; the two claims defects were pre-existing summary-surface staleness the archive pass surfaced rather than caused.

"""

with io.open(CL, encoding="utf-8") as fh:
    t = fh.read()
anchor = "## [v5.12.0] - 2026-09-14"
if anchor not in t:
    print("FAILED: changelog anchor")
    sys.exit(1)
with io.open(CL, "w", encoding="utf-8") as fh:
    fh.write(t.replace(anchor, ENTRY + anchor, 1))
print("  + CHANGELOG v5.12.1")

with io.open(CM, encoding="utf-8") as fh:
    c = fh.read()
if "**Version:** v5.12.0 (September 14, 2026)" not in c:
    print("FAILED: version anchor")
    sys.exit(1)
c = c.replace("**Version:** v5.12.0 (September 14, 2026)",
              "**Version:** v5.12.1 (September 14, 2026)", 1)

bullet = "- **Paper-retirement policy + C24; Lorentzian tail archived (2026-09-14, v5.12.0):**"
new = ("- **/qa DELTA on the group1 archive = DEFECTS, remediated (2026-09-14, "
       "v5.12.1):** two summary surfaces still credited a descoped Lorentzian "
       "claim as achieved (the C16 paraphrase blind spot); the C24 gate's own "
       "guard coverage had gaps. Archive-note edits verified faithful. See "
       "CHANGELOG v5.12.1.\n")
if bullet not in c:
    print("FAILED: Sec 2 anchor")
    sys.exit(1)
c = c.replace(bullet, new + bullet, 1)
with io.open(CM, "w", encoding="utf-8") as fh:
    fh.write(c)
print("  + CLAUDE.md version + Sec. 2 one-liner")
