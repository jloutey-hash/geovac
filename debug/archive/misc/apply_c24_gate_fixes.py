r"""DELTA code-review remediation for the C24 gate (v5.12.0).

The reviewer confirmed the gating check (A) is sound and genuinely fires. The
findings are all in the guards AROUND the machinery, which is the "one level up"
version of the very holes C24 was built to close:

  D1 (MATERIAL-SMALL) the ratchet's known/fresh partition lives only in main(),
     so no test guards it -- a future edit inverting the set-difference would
     silently stop surfacing NEW re-derivation loci. Fix: extract it to a pure
     helper the mirror test drives directly.
  D2 (MATERIAL-SMALL) --selftest advertises "all five probes fire" while
     exercising 3 of 5 check-A branches; the nonexistent-file and empty-reason
     branches fire today but are unguarded. Fix: add both probes + mirror
     assertions.
  D3 (NIT) an empty-but-valid manifest ([]) would make the prune wipe every
     non-index page. Fix: guard the prune on a non-empty manifest.
  D5 (NIT) the prune loop variable `page` shadows the page() builder helper.

This script does the CODE fixes only. The new TESTS go in a separate applier and
are fire-tested there, per the guard-writing rule (a guard written in the same
pass as its subject inherits the subject's reasoning).

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

GATE = "debug/qa/check_paper_retirement.py"
PAGES = "debug/build_paper_pages.py"
EDITS = []


def edit(path, old, new, label):
    EDITS.append((path, old, new, label))


# ---- D1: extract the partition into a pure, testable helper ---------------
edit(GATE,
     '''def check_c(rows, paths):''',
     '''def partition_probe_hits(hits, baseline):
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


def check_c(rows, paths):''',
     "D1: partition_probe_hits helper extracted")

edit(GATE,
     '''    known, fresh = [], []
    for fname, term, where in hits:
        key = "%s|%s" % (fname, term)
        was = set(baseline.get(key, []))
        new_docs = sorted(set(where) - was)
        (fresh if new_docs else known).append((fname, term, where, new_docs))

    # The baseline size is ALWAYS printed.''',
     '''    known, fresh = partition_probe_hits(hits, baseline)

    # The baseline size is ALWAYS printed.''',
     "D1: main() calls the helper")

# ---- D2: two missing check-A probes in --selftest -------------------------
edit(GATE,
     '''    # A3: unknown class
    rows = [{"file": "nope.tex", "class": "WHATEVER", "retired": "x",
             "reason": "r", "triggers": ["t"]}]
    p = check_a(rows, set())
    ok = any("unknown class" in x for x in p)
    bad += 0 if ok else 1
    print("  [%s] A: row with an unknown class" % ("OK" if ok else "DEAD"))''',
     '''    # A3: unknown class
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
    print("  [%s] A: row declaring no reason" % ("OK" if ok else "DEAD"))''',
     "D2: A4 (nonexistent file) + A5 (empty reason) probes added")

edit(GATE,
     '''    print()
    if bad:
        print("RESULT: SELFTEST FAIL -- %d check(s) could not fire" % bad)
        return 1
    print("RESULT: SELFTEST PASS -- all five probes fire")
    return 0''',
     '''    # B probe now also asserts the fresh/known partition directly (D1).
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
    return 0''',
     "D1/D2: selftest gains the partition probe; count language de-hardcoded")

# ---- D3 + D5: guard the prune, de-shadow the loop variable ----------------
edit(PAGES,
     '''    live_ids = {e["id"] for e in entries}
    pruned = []
    for page in sorted(OUT_DIR.glob("*.html")):
        if page.stem == "index" or page.stem in live_ids:
            continue
        page.unlink()
        pruned.append(page.stem)
    if pruned:
        print(f"pruned {len(pruned)} orphan page(s): {', '.join(pruned)}")''',
     '''    # Guard: an empty-but-valid manifest ([]) must NOT be read as "every page
    # is an orphan" -- that would delete the whole directory (D3, 2026-09-14).
    # A real manifest always has entries; an empty one is a build error
    # upstream, not a signal to wipe.
    if not entries:
        print("WARNING: manifest is empty; skipping prune to avoid wiping pages")
    else:
        live_ids = {e["id"] for e in entries}
        pruned = []
        for html_page in sorted(OUT_DIR.glob("*.html")):
            if html_page.stem == "index" or html_page.stem in live_ids:
                continue
            html_page.unlink()
            pruned.append(html_page.stem)
        if pruned:
            print(f"pruned {len(pruned)} orphan page(s): {', '.join(pruned)}")''',
     "D3+D5: prune guarded on non-empty manifest; loop var de-shadowed")

by_path = {}
for path, old, new, label in EDITS:
    by_path.setdefault(path, []).append((old, new, label))

applied, failed = [], []
for path, items in by_path.items():
    with io.open(path, encoding="utf-8") as fh:
        t = fh.read()
    for old, new, label in items:
        if old in t:
            t = t.replace(old, new, 1)
            applied.append(label)
        else:
            failed.append("%s   [%s]" % (label, path))
    with io.open(path, "w", encoding="utf-8") as fh:
        fh.write(t)

print("applied %d of %d" % (len(applied), len(EDITS)))
for a in applied:
    print("  +", a)
if failed:
    print("UNMATCHED:")
    for f in failed:
        print("  -", f)
    sys.exit(1)
