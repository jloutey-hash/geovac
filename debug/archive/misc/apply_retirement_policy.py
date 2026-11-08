r"""Add the paper-retirement policy (CLAUDE.md Sec. 9) and the C24 criterion.

Sec. 9 is PM-editable (Sec. 13.5 table: "7-9 (Code/Coding/Workflow) | Yes").
The policy is deliberately short -- the PI asked for about a third of the first
draft, after correctly pointing out that most of that draft was built on DOI
permanence, which does not apply to a personal Zenodo-stamped repo.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

CM = "CLAUDE.md"
CRIT = "docs/qa/criteria.md"
QA = ".claude/commands/qa.md"

POLICY = r"""### Paper Retirement (added 2026-09-14, PI-directed)

**Retire claims, not papers.** The corpus's most-tested case is Paper 2: its
combination rule was demoted Conjectures -> Core -> Observations and then
relabelled conjecture -> observation, and the paper is still live and cited by
16 others. Tier demotion is the default move; retirement is the exception.

**Neither trigger is sufficient alone.**
- A withdrawn foundation is not grounds. Paper 45's main theorem was withdrawn
  and it remains load-bearing, with five healthy dependents.
- Zero dependents is not grounds. Papers 52 and 53 have none because they are
  DRAFT, which is unfinished, not supplanted, and needs finishing or dropping
  rather than archiving.

**The signal is the health of the dependents.** A paper propped up mainly by
other descoped papers is in a dying subtree. Two cheaper metrics were tried and
both failed, recorded so they are not rebuilt: plain citation count misses a
supplanted cluster entirely, because its members cite each other; and
"citations from outside the paper's own strongly connected component" collapses,
because 49 of 54 papers form one component. **C24** computes the surviving
metric and prints the profile.

**Archive, never delete.** Preservation is not the reason -- the Zenodo deposit
holds the PDF bytes, so the archival copy already exists independently of the
repo. The reason is that the repo copy is what future work greps, and pruning
has already cost this project real time
(`memory/feedback_resurrect_pruned_artifacts.md`). Deleting is not catastrophic;
archiving is simply free.

**Every retirement declares itself.** Move the `.tex` to `papers/archive/`, add
a row to `docs/retired_papers.md` (class, reason, and -- for anything with no
live successor -- trigger terms), keep the `papers/INDEX.md` row as a tombstone,
and stamp the documents that cited it, exactly as a retracted *claim* stamps its
dependents.

**Trigger terms are how the archive stays discoverable.** An archived paper's
approaches would otherwise be invisible to anyone starting new work, which is
the Sec. 3 re-derivation problem one level up. The answer is the same: record
the attempt with the phrase a future sprint would actually use, and let C24's
probe report it. No new habit, and no reminder to forget.

**The decision is the PI's.** C24 informs; it never nominates a verdict and
never blocks.

"""

CRITERION = r"""## C24 -- paper retirement: register integrity + dependency profile (added 2026-09-14, PI direction)

`debug/qa/check_paper_retirement.py`. Three checks, one of which gates.

- **A. Register integrity (FAILS).** Every `.tex` in `papers/archive/` has a row
  in `docs/retired_papers.md` declaring a class and a reason, and every row
  whose class has no live successor declares trigger terms. A file in the
  archive with no row fails; a row naming a file that is not there fails.
- **B. Dependency profile (advisory).** For every paper whose `INDEX.md` status
  is DESCOPED / PARTIAL / DRAFT, print how many of its citers are themselves
  in that set. Ranking, not verdict.
- **C. Re-derivation probe (advisory).** Grep each retired paper's trigger terms
  across the live corpus and report hits, so an approach attempted in an
  archived paper surfaces when it is attempted again. Over-broad terms are
  visible by their hit count.

**Why it exists, measured.** Four papers (46-49) sat in the live, DOI-stamped
set for months after the model they rest on was withdrawn. They were honestly
labelled -- each states its descope in its own abstract -- and no step in the
process was ever responsible for asking whether they should still be there.
That is a **decision** gap, not an honesty gap, and no existing gate covered it.

**Two metrics were tried and failed** before B settled, recorded so nobody
rebuilds them: plain citation count (misses a self-citing cluster: Papers 46-49
carry 6, 8, 3 and 2 citers while nothing outside the group depends on them), and
citations from outside the paper's own strongly connected component (collapses:
49 of 54 papers form one component, so the check reported 28 papers including
three keystones).

**B stops short of a verdict on purpose.** Whether a citation is a real
dependency or a see-also is the Sec. 13.8 test -- *if the cited claim were
withdrawn tomorrow, would this sentence have to change?* -- and that is not
automatable. Retirement is a PI decision.

Self-test: `--selftest` proves all five probes fire against synthetic defects,
not against real debt, so it keeps working once the debt is cleared.

"""

applied, failed = [], []

# ---- CLAUDE.md Sec. 9 ------------------------------------------------------
with io.open(CM, encoding="utf-8") as fh:
    t = fh.read()
anchor = "### Benchmarking Rule"
if "### Paper Retirement (added 2026-09-14" in t:
    applied.append("policy already present")
elif anchor in t:
    t = t.replace(anchor, POLICY + anchor, 1)
    with io.open(CM, "w", encoding="utf-8") as fh:
        fh.write(t)
    applied.append("CLAUDE.md Sec. 9: Paper Retirement policy")
else:
    failed.append("CLAUDE.md anchor")

# ---- criteria.md -----------------------------------------------------------
with io.open(CRIT, encoding="utf-8") as fh:
    c = fh.read()
if "## C24 --" in c or "## C24 —" in c:
    applied.append("C24 criterion already present")
else:
    c = c.rstrip() + "\n\n---\n\n" + CRITERION
    with io.open(CRIT, "w", encoding="utf-8") as fh:
        fh.write(c)
    applied.append("docs/qa/criteria.md: C24 entry")

# ---- wire C24 into the /qa deterministic list ------------------------------
with io.open(QA, encoding="utf-8") as fh:
    q = fh.read()
hook = "   - **C15 inline arXiv-ID resolvability**"
addition = """   - **C24 paper retirement** -- `debug/qa/check_paper_retirement.py`
     (FAILs if an archived paper has no row in `docs/retired_papers.md`, or a
     row lacks its class, reason, or -- where there is no live successor --
     trigger terms. Also prints, as **advisory**, the dependency profile of
     every DESCOPED/PARTIAL/DRAFT paper, and any retired paper's trigger term
     that is live in the corpus again). **The gate informs; retiring a paper is
     a PI decision and C24 never nominates a verdict.** Corpus-wide;
     `--selftest` built in.
"""
if "C24 paper retirement" in q:
    applied.append("C24 already wired into qa.md")
elif hook in q:
    q = q.replace(hook, addition + hook, 1)
    with io.open(QA, "w", encoding="utf-8") as fh:
        fh.write(q)
    applied.append("qa.md: C24 added to the deterministic check list")
else:
    failed.append("qa.md anchor")

print("applied %d" % len(applied))
for a in applied:
    print("  +", a)
if failed:
    print("UNMATCHED:")
    for f in failed:
        print("  -", f)
    sys.exit(1)
