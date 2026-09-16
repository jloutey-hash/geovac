r"""C24 -- paper retirement: register integrity, candidates, re-derivation probe.

WHY THIS GATE EXISTS (measured 2026-09-14).
Four papers (46, 47, 48, 49) sat in the live, DOI-stamped set, unreviewed as a set, after
the model they rest on was withdrawn by the Paper 45 annihilation theorem.  They
were honestly labelled -- each states its descope in its own abstract -- and
nothing in the process was ever responsible for asking whether they should still
be in the live set.  That is a DECISION gap, not an honesty gap.

A plain citation-orphan test cannot find them, and that is the whole difficulty:
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
change?"), and that is not automatable.  Retirement is a PI decision.

GATE SELF-AUDIT RULE: `--selftest` proves every check FIRES, against synthetic
probes rather than against real debt, so it keeps working after the real debt is
cleared.
"""
from __future__ import annotations

import argparse
import collections
import json
import glob
import io
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
REGISTER = os.path.join(ROOT, "docs", "retired_papers.md")
ARCHIVE = os.path.join(ROOT, "papers", "archive")
INDEX = os.path.join(ROOT, "papers", "INDEX.md")
PROBE_BASELINE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              "retirement_probe_baseline.json")

SUMMARY_DOCS = {
    "geovac_field_guide.tex",
    "group1_operator_algebras_synthesis.tex",
    "group2_quantum_chemistry_synthesis.tex",
    "group3_foundations_synthesis.tex",
    "group4_quantum_computing_synthesis.tex",
    "group5_qed_gauge_synthesis.tex",
    "group6_precision_observations_synthesis.tex",
}

NEEDS_TRIGGERS = {"CLOSED", "CLOSED-VALID"}
CANDIDATE_STATUSES = ("DESCOPED", "PARTIAL", "DRAFT")

CITE_RE = re.compile(r"\\cite[a-zA-Z]*\{([^}]*)\}")
NUM_RE = re.compile(r"paper[_]?(\d+)", re.I)
ROW_RE = re.compile(r"^\|\s*`([^`]+)`\s*\|\s*([A-Z-]+)\s*\|([^|]*)\|([^|]*)\|(.*)\|\s*$")


# ---------------------------------------------------------------- helpers
def live_papers():
    return sorted(glob.glob(os.path.join(ROOT, "papers", "group*", "*.tex"))
                  + glob.glob(os.path.join(ROOT, "papers", "synthesis", "*.tex")))


def paper_num(path):
    m = NUM_RE.search(os.path.basename(path))
    return int(m.group(1)) if m else None


def read(path):
    return io.open(path, encoding="utf-8", errors="replace").read()


def parse_register(text):
    """Rows under '## Register'.  Returns list of dicts."""
    rows = []
    body = text.split("## Register", 1)[-1].split("\n## ", 1)[0]
    for line in body.splitlines():
        m = ROW_RE.match(line.strip())
        if not m:
            continue
        fname, cls, retired, reason, triggers = m.groups()
        trig = [t.strip() for t in triggers.split(";") if t.strip() and t.strip() != "—"]
        rows.append({"file": fname.strip(), "class": cls.strip(),
                     "retired": retired.strip(), "reason": reason.strip(),
                     "triggers": trig})
    return rows


def citation_graph(paths):
    """number -> set of basenames citing it, and number -> set it cites."""
    cited_by = collections.defaultdict(set)
    cites = collections.defaultdict(set)
    for f in paths:
        me = paper_num(f)
        text = read(f)
        for m in CITE_RE.finditer(text):
            for key in m.group(1).split(","):
                n = NUM_RE.search(key.strip())
                if not n:
                    continue
                t = int(n.group(1))
                if t == me:
                    continue
                cited_by[t].add(os.path.basename(f))
                if me is not None:
                    cites[me].add(t)
    return cited_by, cites


def _sccs_REMOVED(nodes, edges):
    """REMOVED 2026-09-14 -- see the module docstring.

    Kept only as a named tombstone so the approach is not rebuilt: on this
    corpus 49 of 54 papers form one strongly connected component, so the
    "outside my own component" metric it was written for reports nothing
    useful.  Deleted rather than repaired.
    """
    index = {}
    low = {}
    on = {}
    stack = []
    out = []
    counter = [0]

    def strong(v):
        work = [(v, iter(sorted(edges.get(v, ()))))]
        index[v] = low[v] = counter[0]
        counter[0] += 1
        stack.append(v)
        on[v] = True
        while work:
            node, it = work[-1]
            advanced = False
            for w in it:
                if w not in nodes:
                    continue
                if w not in index:
                    index[w] = low[w] = counter[0]
                    counter[0] += 1
                    stack.append(w)
                    on[w] = True
                    work.append((w, iter(sorted(edges.get(w, ())))))
                    advanced = True
                    break
                if on.get(w):
                    low[node] = min(low[node], index[w])
            if advanced:
                continue
            work.pop()
            if work:
                low[work[-1][0]] = min(low[work[-1][0]], low[node])
            if low[node] == index[node]:
                comp = set()
                while True:
                    w = stack.pop()
                    on[w] = False
                    comp.add(w)
                    if w == node:
                        break
                out.append(comp)

    for v in sorted(nodes):
        if v not in index:
            strong(v)
    return out


def index_status(text, num):
    m = re.search(r"^\|\s*%d\s+`[^`]+`\s*\|\s*\**([A-Za-z-]+)\**" % num, text, re.M)
    return m.group(1).upper() if m else ""


# ---------------------------------------------------------------- checks
def check_a(rows, archive_files):
    """Register integrity. FAILS."""
    problems = []
    listed = {r["file"] for r in rows}
    for f in sorted(archive_files):
        if f not in listed:
            problems.append("archived file has NO register row: %s" % f)
    for r in rows:
        if not os.path.exists(os.path.join(ARCHIVE, r["file"])):
            problems.append("register row names a file not in papers/archive/: %s"
                            % r["file"])
        if r["class"] in NEEDS_TRIGGERS and not r["triggers"]:
            problems.append("%s is %s but declares no trigger terms"
                            % (r["file"], r["class"]))
        if r["class"] not in NEEDS_TRIGGERS | {"SUCCESSOR-COVERED"}:
            problems.append("%s has unknown class %r" % (r["file"], r["class"]))
        if not r["reason"]:
            problems.append("%s declares no reason" % r["file"])
    return problems


def check_b(paths, cited_by, status_of=None):
    """Dependency profile of every paper whose own status is not healthy.

    ADVISORY, and deliberately NOT a verdict.

    Two metrics were tried and discarded before this one, both recorded so
    nobody rebuilds them:

      * plain citation count -- misses a supplanted cluster entirely, because
        the members cite each other (Papers 46-49 had 6, 8, 3, 2 citers apiece
        while nothing outside the group depended on any of them);

      * citations from outside the paper's own strongly connected component --
        collapses, because 49 of 54 papers form ONE component, so "external"
        is empty for almost everything and the metric reports 28 papers
        including three keystones.

    What carries signal is the HEALTH of the dependents: a paper supported
    mainly by other descoped papers is in a dying subtree.  On this corpus that
    separates cleanly -- the four Lorentzian papers hold the top four slots and
    no other paper has a single unhealthy dependent.

    Whether a citation is a real DEPENDENCY or just a see-also is the §13.8
    judgment ("if the cited claim were withdrawn tomorrow, would this sentence
    have to change?") and is not automatable.  So this prints the profile and
    stops.
    """
    if status_of is None:
        idx = read(INDEX) if os.path.exists(INDEX) else ""
        status_of = lambda n: index_status(idx, n)
    out = []
    for f in sorted(paths, key=lambda p: (paper_num(p) or 0)):
        n = paper_num(f)
        if n is None or os.path.basename(f) in SUMMARY_DOCS:
            continue
        status = status_of(n)
        if status not in CANDIDATE_STATUSES:
            continue
        healthy, unhealthy = [], []
        for src in cited_by.get(n, ()):
            if src in SUMMARY_DOCS:
                continue
            sn = paper_num("x/" + src)
            if sn is None or sn == n:
                continue
            (unhealthy if status_of(sn) in CANDIDATE_STATUSES
             else healthy).append(sn)
        out.append({"num": n, "file": os.path.basename(f), "status": status,
                    "healthy": sorted(set(healthy)),
                    "unhealthy": sorted(set(unhealthy))})
    out.sort(key=lambda r: (len(r["healthy"]), -len(r["unhealthy"])))
    return out


def partition_probe_hits(hits, baseline):
    """Split check_c hits into (known, fresh) against a recorded baseline.

    Pure and side-effect-free so the mirror test can drive it directly -- the
    ratchet discrimination this computes (a NEW locus taking up a retired
    paper's topic vs. a live paper that has always narrated the arc) is the
    C22-critical half of check C, and until this was extracted it lived only in
    main() with no test (DELTA code review D1, 2026-09-14).

    `hits` is a list of (fname, term, where); `baseline` maps "fname|term" ->
    list of previously-seen docs.  A hit is FRESH iff it names a doc not in its
    baseline entry.  Returns (known, fresh), each a list of
    (fname, term, where, new_docs).
    """
    known, fresh = [], []
    for fname, term, where in hits:
        key = "%s|%s" % (fname, term)
        was = set(baseline.get(key, []))
        new_docs = sorted(set(where) - was)
        (fresh if new_docs else known).append((fname, term, where, new_docs))
    return known, fresh


def check_c(rows, paths):
    """Re-derivation probe. ADVISORY."""
    blobs = {os.path.basename(f): read(f).lower() for f in paths}
    hits = []
    for r in rows:
        for t in r["triggers"]:
            key = t.lower()
            where = sorted(b for b, txt in blobs.items() if key in txt)
            if where:
                hits.append((r["file"], t, where))
    return hits


# ---------------------------------------------------------------- selftest
def selftest():
    print("C24 selftest -- every check must FIRE on a synthetic defect\n")
    bad = 0

    # A1: archived file with no row
    p = check_a([], {"ghost.tex"})
    ok = any("NO register row" in x for x in p)
    bad += 0 if ok else 1
    print("  [%s] A: archived file with no register row" % ("OK" if ok else "DEAD"))

    # A2: CLOSED row with no triggers
    rows = [{"file": "nope.tex", "class": "CLOSED", "retired": "x",
             "reason": "r", "triggers": []}]
    p = check_a(rows, set())
    ok = any("no trigger terms" in x for x in p)
    bad += 0 if ok else 1
    print("  [%s] A: CLOSED row declaring no trigger terms" % ("OK" if ok else "DEAD"))

    # A3: unknown class
    rows = [{"file": "nope.tex", "class": "WHATEVER", "retired": "x",
             "reason": "r", "triggers": ["t"]}]
    p = check_a(rows, set())
    ok = any("unknown class" in x for x in p)
    bad += 0 if ok else 1
    print("  [%s] A: row with an unknown class" % ("OK" if ok else "DEAD"))

    # A4: register row naming a file not in papers/archive/ (D2, 2026-09-14)
    rows = [{"file": "does_not_exist_in_archive.tex", "class": "CLOSED",
             "retired": "x", "reason": "r", "triggers": ["t"]}]
    p = check_a(rows, set())
    ok = any("not in papers/archive" in x for x in p)
    bad += 0 if ok else 1
    print("  [%s] A: row naming a file absent from the archive"
          % ("OK" if ok else "DEAD"))

    # A5: row declaring no reason (D2, 2026-09-14)
    rows = [{"file": "nope.tex", "class": "SUCCESSOR-COVERED", "retired": "x",
             "reason": "", "triggers": []}]
    p = check_a(rows, set())
    ok = any("no reason" in x for x in p)
    bad += 0 if ok else 1
    print("  [%s] A: row declaring no reason" % ("OK" if ok else "DEAD"))

    # B: a descoped paper propped up by other descoped papers must surface with
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
    print("       paper is absent and a well-supported descoped paper reads 0")

    # C: a trigger term present in the live corpus must be reported
    class Fake(dict):
        pass
    rows = [{"file": "old.tex", "class": "CLOSED", "retired": "x", "reason": "r",
             "triggers": ["zzz unique probe phrase"]}]
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, "live.tex")
        io.open(f, "w", encoding="utf-8").write("we revisit the ZZZ Unique Probe Phrase here")
        hits = check_c(rows, [f])
    ok = bool(hits)
    bad += 0 if ok else 1
    print("  [%s] C: trigger term live in the corpus is reported"
          % ("OK" if ok else "DEAD"))

    # B probe now also asserts the fresh/known partition directly (D1).
    known, fresh = partition_probe_hits(
        [("old.tex", "widget", ["a.tex", "b.tex"])],
        {"old.tex|widget": ["a.tex"]})
    ok = (len(fresh) == 1 and fresh[0][3] == ["b.tex"] and not known)
    bad += 0 if ok else 1
    print("  [%s] C: partition surfaces a NEW doc, suppresses the baselined one"
          % ("OK" if ok else "DEAD"))

    print()
    if bad:
        print("RESULT: SELFTEST FAIL -- %d probe(s) could not fire" % bad)
        return 1
    print("RESULT: SELFTEST PASS -- all probes fire")
    return 0


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--gate", default=None, help="accepted for symmetry; C24 is corpus-wide")
    ap.add_argument("--update-baseline", action="store_true",
                    help="record the current probe hits as the baseline")
    args = ap.parse_args()
    if args.selftest:
        return selftest()

    if not os.path.exists(REGISTER):
        print("RESULT: FAIL (no docs/retired_papers.md; C24 requires the register)")
        return 1

    rows = parse_register(read(REGISTER))
    archive_files = {os.path.basename(p)
                     for p in glob.glob(os.path.join(ARCHIVE, "*.tex"))}
    paths = live_papers()
    cited_by, cites = citation_graph(paths)

    print("A. register integrity (FAILS the gate)")
    problems = check_a(rows, archive_files)
    if problems:
        for p in problems:
            print("   [FAIL] %s" % p)
    else:
        print("   [ok] %d archived paper(s), all registered with the required fields"
              % len(archive_files))

    print()
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
    flagged = [c for c in cands if not c["healthy"]]

    print()
    print("C. re-derivation probe (ADVISORY, ratcheted)")
    hits = check_c(rows, paths)

    baseline = {}
    if os.path.exists(PROBE_BASELINE):
        try:
            baseline = json.load(io.open(PROBE_BASELINE, encoding="utf-8"))
        except ValueError:
            baseline = {}

    known, fresh = partition_probe_hits(hits, baseline)

    # The baseline size is ALWAYS printed. A ratchet that hides its own size is
    # how debt becomes permanent (the C22 rule).
    print("   [baseline] %d known hit(s) across %d recorded term(s) -- these are"
          % (len(known), len(baseline)))
    print("              the live papers that narrate an archived arc, not"
          " re-derivations")
    if not hits:
        print("   [ok] no retired-paper trigger term is live in the corpus")
    if not fresh:
        print("   [ok] no NEW locus has taken up a retired paper's topic")
    for fname, term, where, new_docs in fresh:
        print("   [NEW] %-38s %-34s now also in: %s"
              % (fname[:38], '"%s"' % term[:32], ", ".join(new_docs[:4])))
        print("         ^ that topic was attempted before. Read the archived"
              " paper before rebuilding it.")
    for fname, term, where, _ in known:
        if len(where) > 8:
            print("   [broad] %-38s %-34s %d docs -- term too generic to be useful"
                  % (fname[:38], '"%s"' % term[:32], len(where)))

    if args.update_baseline:
        snap = {"%s|%s" % (f, t): sorted(w) for f, t, w in hits}
        io.open(PROBE_BASELINE, "w", encoding="utf-8").write(
            json.dumps(snap, indent=2, sort_keys=True) + "\n")
        print("   [baseline updated] %d term(s) recorded" % len(snap))
    hits = fresh

    print()
    if problems:
        print("RESULT: FAIL (%d register defect(s); %d candidate(s), %d probe hit(s) advisory)"
              % (len(problems), len(flagged), len(hits)))
        return 1
    print("RESULT: PASS (register intact; %d retirement candidate(s) and %d "
          "NEW re-derivation hit(s) reported as advisory)" % (len(flagged), len(hits)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
