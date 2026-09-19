r"""Repair eaten escapes in the v5.14.3 CHANGELOG entry.

WHAT HAPPENED, and it is the funniest possible instance of it.  The paragraph
documenting that I keep mangling backslashes in bash heredocs was itself
written through a bash heredoc, and got mangled.  `cat -A` shows:

    `^Iext{}`          -- `\t` ate itself: TAB + "ext{}"   (wanted \text{})
    `^Lootnotemark`    -- `\f` ate itself: FF  + "ootnotemark"
    `^Lootnotetext`    -- same
    `0\.0002`          -- wanted the DOUBLED form `0\\.0002`, so the sentence
                          now displays the CORRECT pattern while claiming it
                          was the broken one, inverting its own point.

This is the C19 class exactly (`debug/qa/check_latex_escapes.py`): a LaTeX
control sequence destroyed by a Python string escape, where the result is a
bare control character that no compiler complains about.  C19's own docstring
carries the hard rule I broke: **never apply a LaTeX-bearing edit through a
bash-heredoc Python replacement string -- write the script to a file with raw
strings and execute it.**  Fourth violation this session.

Run:  python debug/qa/_fix_changelog_eaten_escapes.py
"""
from __future__ import annotations

import io
import sys

PATH = "CHANGELOG.md"

# Bullet 4: rebuild with the control characters replaced by real backslashes.
BAD4_PREFIX = "4. **Collateral defects while fixing defects:** `"
GOOD4 = (
    "4. **Collateral defects while fixing defects:** `" + chr(92) + "text{}` in a "
    "`dcolumn` cell, a `" + chr(92) + "footnotemark` with no `" + chr(92) +
    "footnotetext`, and a `" + chr(92) + "cite{bates1953}` in Paper 11's abstract "
    "when its own bibitem is `Bates1953` (natbib: \"Citation `bates1953' ... "
    "undefined\"). All three caught by compiling; all fixed."
)

# Bullet 1: show the DOUBLED escape, which is what actually shipped.
BAD1 = "so escapes doubled (`0" + chr(92) + ".0002`) and matched nothing"
GOOD1 = ("so escapes doubled (`0" + chr(92) * 2 + ".0002` -- two literal "
         "backslashes, matching nothing a corpus ever contains) and matched nothing")


def main() -> None:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

    s = io.open(PATH, encoding="utf-8").read()
    fixed = []

    # --- bullet 4: locate the corrupted line and replace it wholesale ---
    lines = s.splitlines(keepends=True)
    for i, ln in enumerate(lines):
        if ln.startswith(BAD4_PREFIX) and ("\t" in ln or "\f" in ln):
            lines[i] = GOOD4 + "\n"
            fixed.append("bullet 4 (TAB/FF control chars -> real backslashes)")
            break
    s = "".join(lines)

    # --- bullet 1: single -> doubled escape ---
    if BAD1 in s:
        s = s.replace(BAD1, GOOD1, 1)
        fixed.append("bullet 1 (single escape -> doubled, restoring its point)")

    io.open(PATH, "w", encoding="utf-8", newline="").write(s)

    for f in fixed:
        print("  fixed: " + f)
    if not fixed:
        print("  nothing matched -- inspect manually before retrying")
        raise SystemExit(1)

    # Verify: no bare TAB/FF anywhere in the v5.14.3 entry.
    entry = s.split("## [v5.14.3]", 1)[1].split("## [v5.14.2]", 1)[0]
    bad = [(c, entry.count(c)) for c in ("\t", "\f", "\v", "\r")
           if entry.count(c)]
    print(f"  residual control characters in the entry: {bad or 'none'}")
    for want in (chr(92) + "text{}", chr(92) + "footnotemark",
                 chr(92) + "footnotetext", chr(92) * 2 + ".0002"):
        print(f"    contains {want!r}: {want in entry}")


if __name__ == "__main__":
    main()
