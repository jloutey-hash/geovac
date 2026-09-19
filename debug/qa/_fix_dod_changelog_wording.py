r"""Match the group2 DoD Change-log wording to C17's declared exempt phrasing.

The Change-log entry I added documents the retraction and therefore QUOTES the
retired number.  C17 correctly flags that as a live occurrence -- it cannot tell
"asserted" from "documented".  Its family declares these exempt phrasings:

    \[retracted \d{4}-\d{2}-\d{2}: p11-h2plus-0002pct-retired\]
    | reproduced to machine precision
    | reference-limited residual
    | wrong by ~?7\.6 orders
    | retired 0\.0002

My sentence says "falsified by ~7.6 orders", a near-miss on the fourth.  A
fifth id-carrying marker would also clear it, but the sentence already IS the
withdrawal, so matching the declared vocabulary is cleaner than stacking
another marker -- and it keeps the exemption surface equal to what the family
actually documents.

Run:  python debug/qa/_fix_dod_changelog_wording.py
"""
from __future__ import annotations

import io
import re
import sys

PATH = "docs/qa/group2.done.md"
OLD = "**0.0002%** falsified by ~7.6 orders"
NEW = "**0.0002%** wrong by ~7.6 orders"


def main() -> None:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

    s = io.open(PATH, encoding="utf-8").read()
    if NEW in s:
        print("already matched; nothing to do")
    else:
        if OLD not in s:
            i = s.find("by ~7.6 orders")
            print("ANCHOR NOT FOUND. Actual text:")
            print(repr(s[max(0, i - 160):i + 80]) if i >= 0 else "  (phrase absent)")
            raise SystemExit(1)
        if s.count(OLD) != 1:
            raise SystemExit(f"ambiguous ({s.count(OLD)} matches)")
        s = s.replace(OLD, NEW, 1)
        io.open(PATH, "w", encoding="utf-8", newline="").write(s)
        print("wording matched to the registry's declared exempt phrasing")

    # Verify the family now sees nothing live in this file.
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "chn", "debug/qa/check_headline_numbers.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    e = [x for x in mod.REGISTRY
         if x["id"] == "p11-h2plus-0002pct-retired"][0]
    pat, req, exm = (re.compile(e["pattern"]), re.compile(e["require_nearby"]),
                     re.compile(e["exempt_if_nearby"]))
    win = getattr(mod, "WINDOW", 3)
    lines = io.open(PATH, encoding="utf-8").read().splitlines()
    live = []
    for k, ln in enumerate(lines):
        if not pat.search(ln):
            continue
        ctx = "\n".join(lines[max(0, k - win):k + win + 1])
        if req.search(ctx) and not exm.search(ctx):
            live.append(k + 1)
    print(f"  {PATH}: live occurrences now -> {live or 'none'}")


if __name__ == "__main__":
    main()
