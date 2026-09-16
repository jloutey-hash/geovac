r"""CHANGELOG + version for the paper-retirement policy and the C24 gate.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

ENTRY = """## [v5.11.21] - 2026-09-14

**Paper retirement now has a policy and a gate. It did not have either, and that is why four supplanted papers sat in the live set indefinitely.** CLAUDE.md Sec. 9, `docs/retired_papers.md`, C24, and its mirror test.

### The gap, stated precisely

Papers 46 through 49 rest on a model the Paper 45 annihilation theorem withdrew. They are **honestly labelled** -- each states its own descope in its own abstract, which I checked before claiming otherwise -- and every claim-level mechanism worked correctly on them. What never happened is anyone deciding whether they should still be in the live set, because no step in the process was responsible for asking. **A decision gap, not an honesty gap.** The distinction matters: nothing here needs fixing in the papers.

### What the policy says

- **Retire claims, not papers.** The precedent is Paper 2, which was demoted Conjectures to Core to Observations, had its label downgraded, and is still live and cited by 16 others. Tier demotion is the default; retirement is the exception.
- **Neither trigger is sufficient alone.** A withdrawn foundation is not grounds (Paper 45 is descoped and load-bearing, with five healthy dependents). Zero dependents is not grounds (Papers 52 and 53 have none because they are DRAFT, which is unfinished rather than supplanted).
- **Archive, never delete** -- but not for the reason first given. Preservation is already handled: the Zenodo deposit holds the PDF bytes, so the archival copy exists independently of the repo. The actual reason is that the repo copy is what future work greps, and pruning has already cost this project real time. Deleting is not catastrophic; archiving is simply free.
- **The decision is the PI's.** C24 informs and never nominates a verdict.

### The metric took three attempts, and the first two are recorded so nobody rebuilds them

| Metric | Result |
|---|---|
| Plain citation count | Misses the cluster entirely. Papers 46-49 carry 6, 8, 3 and 2 citers while nothing outside the group depends on any of them. |
| Citations from outside the paper's own strongly connected component | Collapses. 49 of 54 papers form ONE component, so "external" is empty for almost everything and the check reported 28 papers including three keystones. |
| **Health of the dependents** | Separates with no false positives. |

A paper propped up mainly by other descoped papers is in a dying subtree. On this corpus the four Lorentzian papers take the top four slots by unhealthy-dependent count and **no other paper has a single one**. Paper 45 reads five healthy dependents, so the metric visibly distinguishes the one that stays from the ones worth examining.

Check B stops short of a verdict deliberately. Whether a citation is a real dependency or a see-also is the Sec. 13.8 test, and that is not automatable.

### Keeping the archive discoverable

An archived paper's approaches would otherwise be invisible to anyone starting new work, which is the Sec. 3 re-derivation problem one level up. The answer is the same one Sec. 3 uses: record the attempt with the phrase a future sprint would actually use, and let the gate report it. **No reminder**, because reminders rot -- this corpus has measured that repeatedly. `docs/retired_papers.md` carries trigger terms for every archived paper with no live successor, and C24 probes them. Reading the seven archived papers to write those terms surfaced a distinction worth keeping: two are covered by a live successor and carry no re-derivation risk at all, while Paper 6 is **valid and merely orphaned**, so the right response to a hit on it is reuse rather than rebuilding.

### Guards

C24's `--selftest` proves all five probes fire against synthetic defects rather than real debt. `tests/test_paper_retirement_check.py` adds ten assertions in the FIRE direction, including the discrimination that matters: a descoped-but-well-supported paper must not read as a dying subtree. That is the Paper 45 case, and a metric that could not tell it from Paper 46 would nominate the wrong paper.

**PI note.** Bumped as a patch, but this is arguably a minor: Sec. 9 gains a standing policy and the `/qa` deterministic list gains a gate, both of which Sec. 9 lists as corpus-significant. Flagging rather than bumping unilaterally, per the version rule. **I also wired C24 into `.claude/commands/qa.md` so it actually runs** -- gate changes have historically been PI-directed, and a gate that never runs is the failure mode this corpus has hit four times, so I chose running over inert. Revert that one line if you would rather it stayed out.

**Nothing was retired.** The Lorentzian decision is yours and remains open.

"""

CL = "CHANGELOG.md"
with io.open(CL, encoding="utf-8") as fh:
    t = fh.read()
anchor = "## [v5.11.20] - 2026-09-14"
if anchor not in t:
    print("FAILED: changelog anchor")
    sys.exit(1)
with io.open(CL, "w", encoding="utf-8") as fh:
    fh.write(t.replace(anchor, ENTRY + anchor, 1))
print("  + CHANGELOG v5.11.21")

CM = "CLAUDE.md"
with io.open(CM, encoding="utf-8") as fh:
    c = fh.read()
if "**Version:** v5.11.20 (September 14, 2026)" not in c:
    print("FAILED: version anchor")
    sys.exit(1)
c = c.replace("**Version:** v5.11.20 (September 14, 2026)",
              "**Version:** v5.11.21 (September 14, 2026)", 1)

bullet = "- **/qa DELTA #3 = DEFECTS, remediated (2026-09-14, v5.11.20):**"
new = ("- **Paper-retirement policy + C24 (2026-09-14, v5.11.21):** four "
       "supplanted papers sat live because nothing owned the decision; "
       "register + dependency-profile gate added. Two metrics failed first. "
       "See CHANGELOG v5.11.21.\n")
if bullet not in c:
    print("FAILED: Sec 2 anchor")
    sys.exit(1)
c = c.replace(bullet, new + bullet, 1)
with io.open(CM, "w", encoding="utf-8") as fh:
    fh.write(c)
print("  + CLAUDE.md version + Sec. 2 one-liner")
