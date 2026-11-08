r"""Close DELTA #2's one advisory SMALL: the empty-manifest prune guard in
build_paper_pages.py had no regression test.

Same move as the D1 fix: extract the prune DECISION into a pure helper
(`orphan_stems`) so it can be tested in isolation, leave the destructive
unlink() in main(), and add a `--selftest` that pins the wipe-guard. The
reviewer verified the guard works by full-run fire-test; this makes that
guarantee survive a future refactor.

The critical assertion: empty entries -> ZERO orphans (never "every page is an
orphan"). That is the wipe the guard prevents.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

P = "debug/build_paper_pages.py"
EDITS = []


def edit(old, new, label):
    EDITS.append((old, new, label))


# ---- 1. extract the decision to a pure helper -----------------------------
edit(
    '''def main() -> None:''',
    '''def orphan_stems(entries, existing_stems):
    """Page stems to prune: those with no manifest entry, never index, never a
    live page.  Pure and total so the wipe-guard is testable (DELTA #2 D3).

    THE LOAD-BEARING INVARIANT: an empty `entries` returns an EMPTY set, never
    "every stem is an orphan".  A real manifest always has entries; an empty one
    is an upstream build error, not a signal to wipe the crawlable pages.
    """
    if not entries:
        return set()
    live_ids = {e["id"] for e in entries}
    return {s for s in existing_stems if s != "index" and s not in live_ids}


def main() -> None:''',
    "orphan_stems helper extracted")

# ---- 2. main() uses the helper --------------------------------------------
edit(
    '''    if not entries:
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
    '''    if not entries:
        print("WARNING: manifest is empty; skipping prune to avoid wiping pages")
    existing = {p.stem for p in OUT_DIR.glob("*.html")}
    pruned = sorted(orphan_stems(entries, existing))
    for stem in pruned:
        (OUT_DIR / f"{stem}.html").unlink()
    if pruned:
        print(f"pruned {len(pruned)} orphan page(s): {', '.join(pruned)}")''',
    "main() prunes via the helper")

# ---- 3. --selftest that pins the guard ------------------------------------
edit(
    '''if __name__ == "__main__":
    main()''',
    '''def _selftest() -> int:
    """Pin the wipe-guard and the prune decision (DELTA #2 D3)."""
    bad = 0

    # THE guard: empty manifest must prune NOTHING, even with pages present.
    got = orphan_stems([], {"paper_1", "paper_2", "index"})
    ok = got == set()
    bad += 0 if ok else 1
    print("  [%s] empty manifest -> prune nothing (no wipe)"
          % ("OK" if ok else "DEAD"))

    # A normal manifest prunes only true orphans, never index, never live.
    entries = [{"id": "paper_1"}, {"id": "paper_2"}]
    got = orphan_stems(entries, {"paper_1", "paper_2", "index", "paper_old"})
    ok = got == {"paper_old"}
    bad += 0 if ok else 1
    print("  [%s] normal manifest -> only true orphans pruned, index+live kept"
          % ("OK" if ok else "DEAD"))

    # A single-entry manifest does not wipe the rest to zero (reviewer's case).
    got = orphan_stems([{"id": "paper_1"}], {"paper_1", "index"})
    ok = got == set()
    bad += 0 if ok else 1
    print("  [%s] single-entry manifest -> no spurious prune" % ("OK" if ok else "DEAD"))

    print()
    if bad:
        print("RESULT: SELFTEST FAIL -- %d probe(s) could not fire" % bad)
        return 1
    print("RESULT: SELFTEST PASS -- the wipe-guard holds")
    return 0


if __name__ == "__main__":
    import sys as _sys
    if "--selftest" in _sys.argv:
        _sys.exit(_selftest())
    main()''',
    "--selftest added")

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
