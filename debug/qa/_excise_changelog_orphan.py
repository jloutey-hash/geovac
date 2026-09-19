r"""Excise the orphaned fragment of the old corrupted bullet 4.

WHAT WENT WRONG IN MY OWN REPAIR, which is worth recording because it is a
second-order instance of the same trap.

The original bullet 4 was corrupted to contain a bare form feed (`\f` eaten
from `\footnotemark`).  My repair iterated `s.splitlines(keepends=True)` and
replaced the line that started with the bullet prefix AND contained a control
character.  But **Python's line-splitting treats `\x0c` as a line boundary**,
so the corrupted bullet had ALREADY been split in two:

    line A: "4. **Collateral defects ...** `<TAB>ext{}` in a `dcolumn` cell, a `"
    line B: "ootnotemark` with no `<FF>ootnotetext`, and a `\cite{bates1953}` ..."

I replaced line A with the corrected bullet and left line B untouched.  Result:
the good bullet, immediately followed by a dangling fragment of the bad one.
The very same property (`\x0c` is a line terminator) is also why my earlier
`for ln, line in enumerate(s.splitlines())` scan printed nothing while the
whole-file count said one form feed remained.

LESSON: form-feed corruption cannot be repaired -- or even located -- with
line-based operations.  Use raw-offset surgery on the whole string.

Run:  python debug/qa/_excise_changelog_orphan.py
"""
from __future__ import annotations

import io
import sys

PATH = "CHANGELOG.md"
FF = chr(12)


def main() -> None:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

    s = io.open(PATH, encoding="utf-8", newline="").read()
    before = s.count(FF)
    if before == 0:
        print("no form feeds present; nothing to do")
        return

    i = s.find(FF)
    # The orphan begins at the start of its physical line (after the previous
    # "\n") and ends at the next "\n".  Delete exactly that span.
    start = s.rfind("\n", 0, i) + 1
    end = s.find("\n", i)
    if end < 0:
        raise SystemExit("orphan runs to EOF -- inspect manually")
    orphan = s[start:end]

    # Safety: the orphan must be a FRAGMENT, not the good bullet.  The good
    # bullet starts with "4. **Collateral"; the fragment must not.
    if orphan.lstrip().startswith("4. **Collateral"):
        raise SystemExit("refusing: that line is the GOOD bullet, not the orphan")
    if "ootnotemark" not in orphan and "ootnotetext" not in orphan:
        raise SystemExit(f"refusing: unexpected orphan content -> {orphan[:90]!r}")

    print(f"  excising orphan ({len(orphan)} chars): {orphan[:100]!r}")
    s = s[:start] + s[end + 1:]
    io.open(PATH, "w", encoding="utf-8", newline="").write(s)

    s2 = io.open(PATH, encoding="utf-8", newline="").read()
    print(f"  form feeds: {before} -> {s2.count(FF)}")
    for ch, name in ((chr(11), "VT"), (chr(8), "BS"), (chr(13), "CR")):
        if s2.count(ch):
            print(f"  NOTE residual {name}: {s2.count(ch)}")

    # The good bullet must have survived intact.
    good = ("4. **Collateral defects while fixing defects:** `" + chr(92)
            + "text{}` in a `dcolumn` cell")
    print(f"  good bullet 4 intact: {good in s2}")
    entry = s2.split("## [v5.14.3]", 1)[1].split("## [v5.14.2]", 1)[0]
    print(f"  control chars remaining in the v5.14.3 entry: "
          f"{[(repr(c), entry.count(c)) for c in (FF, chr(11), chr(8)) if entry.count(c)] or 'none'}")


if __name__ == "__main__":
    main()
