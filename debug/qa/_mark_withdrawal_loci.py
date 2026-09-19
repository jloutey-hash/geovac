r"""Append id-carrying withdrawal markers to the four records that QUOTE the
retired headlines in order to withdraw them.

WHY THESE ARE FALSE POSITIVES OF MY OWN MAKING.  Two turns ago I narrowed both
C17 families' `exempt_if_nearby` from a generic alternation
(`retired|superseded|withdrawn|...`) down to the per-entry marker plus specific
corrected phrasings -- because the generic version was excusing two LIVE loci
on incidental vocabulary in neighbouring table rows (the 2026-09-04 DELTA #5
class).  That narrowing was right.  Its consequence is that my own sweep text,
which says things like "(0.0002% retired 2026-09-19)", is now correctly seen as
a live occurrence: the gate cannot tell "asserted" from "explained", and the
registry's designed answer is the marker that names the entry it withdraws.

So each of these gets `[retracted YYYY-MM-DD: <entry-id>]`, which exempts THAT
entry and no other.

Run:  python debug/qa/_mark_withdrawal_loci.py
"""
from __future__ import annotations

import io
import sys

P11ID = "p11-h2plus-0002pct-retired"
P13ID = "p13-he-0019pct-nonexistent"
MARK = "[retracted 2026-09-19: {}]"

EDITS = [
    # docs/claim_test_matrix.md:268 -- quotes both the retired % and the
    # retired 5000x ratio in order to retire them.
    ("docs/claim_test_matrix.md",
     '| 11 | spectral Laguerre **machine precision** (retired: 0.0002%,',
     '| 11 | spectral Laguerre **machine precision** ' + MARK.format(P11ID)
     + ' (retired: 0.0002%,'),

    # docs/paper_notes_archive.md:108
    ("docs/paper_notes_archive.md",
     'via spectral Laguerre (0.0002% retired 2026-09-19) |',
     'via spectral Laguerre (0.0002% retired 2026-09-19 '
     + MARK.format(P11ID) + ') |'),

    # docs/paper_notes_archive.md:110
    ("docs/paper_notes_archive.md",
     'cusp-extrapolated (0.019% retired 2026-09-19), fiber bundle',
     'cusp-extrapolated (0.019% retired 2026-09-19 '
     + MARK.format(P13ID) + '), fiber bundle'),

    # docs/project_closeout_plan.md:26
    ("docs/project_closeout_plan.md",
     'the 0.0002% recorded here was retired 2026-09-19)',
     'the 0.0002% recorded here was retired 2026-09-19 '
     + MARK.format(P11ID) + ')'),
]


def main() -> None:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    by_file: dict[str, list[tuple[str, str]]] = {}
    for path, old, new in EDITS:
        by_file.setdefault(path, []).append((old, new))

    failures = []
    for path, pairs in by_file.items():
        s = io.open(path, encoding="utf-8").read()
        n = 0
        for old, new in pairs:
            if old in new and old not in s:
                failures.append(f"{path}: anchor NOT FOUND -> {old[:70]!r}")
                continue
            if s.count(old) != 1:
                failures.append(f"{path}: anchor count {s.count(old)} -> {old[:60]!r}")
                continue
            s = s.replace(old, new, 1)
            n += 1
        if n:
            io.open(path, "w", encoding="utf-8", newline="").write(s)
        print(f"  {path}: {n}/{len(pairs)} marked")

    if failures:
        print("\nFAILURES (nothing written for these):")
        for f in failures:
            print("   " + f)
        raise SystemExit(1)
    print("\nall withdrawal markers applied")


if __name__ == "__main__":
    main()
