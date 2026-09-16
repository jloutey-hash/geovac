r"""Finish C24: correct the docstring, drop the discarded SCC code, rewrite the
check-B display and its selftest probe.

The module docstring still claimed the SCC metric is "the principled definition
and nothing has to be hardcoded".  That was written before the metric was
measured against the real corpus, where 49 of 54 papers turned out to form a
single component.  Leaving that sentence would have been the exact defect this
session spent all day remediating: a replacement reading asserted before it was
checked.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import re
import sys

P = "debug/qa/check_paper_retirement.py"
EDITS = []


def edit(old, new, label):
    EDITS.append((old, new, label))


# ---- 1. docstring: the SCC route was tried and failed --------------------
edit(
    """A plain citation-orphan test cannot find them, and that is the whole difficulty:
they cite each other, so from the inside the cluster looks healthy.  The
measurement that works is citations from OUTSIDE the paper's own strongly
connected component of the citation graph.  A self-citing island IS an SCC, so
the SCC is the principled definition and nothing has to be hardcoded.

THREE CHECKS
  A  register integrity ............ FAILS the gate
  B  retirement candidates ......... advisory (nominates; never decides)
  C  re-derivation probe ........... advisory

Only A fails.  B must never block: a brand-new paper legitimately has zero
external citers, and so does a deliberately terminal one.  The gate nominates;
retirement is a PI decision.""",
    """A plain citation-orphan test cannot find them, and that is the whole difficulty:
they cite each other, so from the inside the cluster looks healthy.

TWO METRICS WERE TRIED AND FAILED, recorded here so nobody rebuilds them:

  1. Plain citation count.  Misses the cluster entirely -- Papers 46-49 carry
     6, 8, 3 and 2 citers apiece while nothing outside the group depends on any
     of them.
  2. Citations from outside the paper's own strongly connected component.  A
     self-citing island IS an SCC, so this looked principled;  measured against
     the real corpus it collapses, because 49 of the 54 numbered papers form
     ONE component.  "External" is then empty for almost everything, and the
     check reported 28 papers including three keystones.

What carries signal is the HEALTH OF THE DEPENDENTS.  A paper supported mainly
by other descoped papers is in a dying subtree.  On this corpus that separates
with no false positives: the four Lorentzian papers take the top four slots by
unhealthy-dependent count, and no other paper has a single one.

THREE CHECKS
  A  register integrity ............ FAILS the gate
  B  dependency profile ............ advisory (informs; never decides)
  C  re-derivation probe ........... advisory

Only A fails.  B must never block, and it deliberately stops short of a verdict:
whether a citation is a real DEPENDENCY or a see-also is the Sec. 13.8 judgment
("if the cited claim were withdrawn tomorrow, would this sentence have to
change?"), and that is not automatable.  Retirement is a PI decision.""",
    "docstring: the SCC route recorded as tried-and-failed")

# ---- 2. drop the now-unused Tarjan implementation ------------------------
edit(
    '''def sccs(nodes, edges):
    """Tarjan strongly connected components. edges: node -> set(node)."""''',
    '''def _sccs_REMOVED(nodes, edges):
    """REMOVED 2026-09-14 -- see the module docstring.

    Kept only as a named tombstone so the approach is not rebuilt: on this
    corpus 49 of 54 papers form one strongly connected component, so the
    "outside my own component" metric it was written for reports nothing
    useful.  Deleted rather than repaired.
    """''',
    "sccs: tombstoned")

# ---- 3. main(): new call + display ---------------------------------------
edit(
    """    print()
    print("B. retirement candidates (ADVISORY -- nominates, never decides)")
    cands = check_b(paths, cited_by, cites)
    if not cands:
        print("   [ok] every live paper has at least one citer outside its own component")
    for c in cands:
        mark = "<<" if c["flag"] else "  "
        print("   %s Paper %-3d %-46s status=%-9s %s"
              % (mark, c["num"], c["file"], c["status"] or "-",
                 "SELF-CITING ISLAND" if c["island"] else "isolated"))
    flagged = [c for c in cands if c["flag"]]
    if flagged:
        print()
        print("   %d candidate(s) marked << : no external dependents AND a"
              % len(flagged))
        print("   descoped/partial/draft status. Retirement is a PI decision.")""",
    """    print()
    print("B. dependency profile of every not-healthy paper (ADVISORY)")
    cands = check_b(paths, cited_by)
    if not cands:
        print("   [ok] no live paper carries a descoped/partial/draft status")
    else:
        print("   %-7s %-11s %-8s %-10s %s"
              % ("paper", "status", "healthy", "unhealthy", "unhealthy dependents"))
        for c in cands:
            print("   P%-6d %-11s %-8d %-10d %s"
                  % (c["num"], c["status"], len(c["healthy"]), len(c["unhealthy"]),
                     ", ".join("P%d" % x for x in c["unhealthy"]) or "-"))
        print()
        print("   Read it as a ranking, not a verdict.  A paper with NO healthy")
        print("   dependents, or whose dependents are mostly themselves descoped,")
        print("   is in a dying subtree and is worth a look.  Whether any given")
        print("   citation is a dependency or a see-also is a judgement call.")
    flagged = [c for c in cands if not c["healthy"]]""",
    "main: dependency-profile display")

# ---- 4. selftest probe for the new check B -------------------------------
edit(
    """    # B: a self-citing island must be nominated, a well-cited paper must not.
    #    91<->92 cite only each other; 93 is cited by 94 (outside its component).
    cites = {91: {92}, 92: {91}, 93: set(), 94: {93}}
    cited_by = {92: {"paper_91_a.tex"}, 91: {"paper_92_b.tex"},
                93: {"paper_94_d.tex"}}
    fake = ["x/paper_91_a.tex", "x/paper_92_b.tex",
            "x/paper_93_c.tex", "x/paper_94_d.tex"]
    got = {c["num"] for c in check_b(fake, cited_by, cites)}
    ok = {91, 92}.issubset(got) and 93 not in got
    bad += 0 if ok else 1
    print("  [%s] B: self-citing island nominated, externally-cited paper not"
          % ("OK" if ok else "DEAD"))""",
    """    # B: a descoped paper propped up by other descoped papers must surface with
    #    its unhealthy count;  a healthy paper must not appear at all; and a
    #    descoped paper with healthy support must appear with unhealthy == 0.
    status = {91: "DESCOPED", 92: "PARTIAL", 93: "ACTIVE", 94: "DESCOPED"}
    cited_by = {91: {"paper_92_b.tex", "paper_94_d.tex"},   # both unhealthy
                93: {"paper_91_a.tex"},                      # healthy paper
                94: {"paper_93_c.tex"}}                      # healthy support
    fake = ["x/paper_91_a.tex", "x/paper_92_b.tex",
            "x/paper_93_c.tex", "x/paper_94_d.tex"]
    rows = {c["num"]: c for c in
            check_b(fake, cited_by, status_of=lambda n: status.get(n, ""))}
    ok = (93 not in rows                                   # healthy: absent
          and len(rows.get(91, {}).get("unhealthy", [])) == 2
          and len(rows.get(91, {}).get("healthy", [])) == 0
          and len(rows.get(94, {}).get("unhealthy", [])) == 0
          and len(rows.get(94, {}).get("healthy", [])) == 1)
    bad += 0 if ok else 1
    print("  [%s] B: dying subtree surfaces with its unhealthy count; a healthy"
          % ("OK" if ok else "DEAD"))
    print("       paper is absent and a well-supported descoped paper reads 0")""",
    "selftest: probe rewritten for the dependency-profile metric")

with io.open(P, encoding="utf-8") as fh:
    t = fh.read()

applied, failed = [], []
for old, new, label in EDITS:
    if old in t:
        t = t.replace(old, new, 1)
        applied.append(label)
    else:
        failed.append(label)

with io.open(P, "w", encoding="utf-8") as fh:
    fh.write(t)

print("applied %d of %d" % (len(applied), len(EDITS)))
for a in applied:
    print("  +", a)
if failed:
    print("UNMATCHED:")
    for f in failed:
        print("  -", f)
    sys.exit(1)
