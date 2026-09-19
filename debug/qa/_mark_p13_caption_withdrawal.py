r"""Append the id-carrying withdrawal marker to Paper 13's new table caption.

WHY.  The caption I wrote for `tab:hierarchy` explains the correction and
therefore QUOTES the retired number:

    ...the $0.0002\%$ this cell carried until then understated the method by
    roughly seven orders...

C17 family `p11-h2plus-0002pct-retired` correctly flags that as a live
occurrence -- it cannot tell "asserted" from "explained".  The registry's
answer to exactly this is the per-entry marker
`[retracted YYYY-MM-DD: <entry-id>]`, which exempts THIS entry and no other
(redesigned 2026-09-04 after a bare token silenced an unrelated live defect
in the next sentence).

WHY A FILE AND NOT A HEREDOC.  Two attempts at this edit through a bash
heredoc failed on `\%` escape mangling -- the second tripped its own anchor
assertion and wrote nothing, which is the guard working.  The standing rule
(memory: no-heredoc-backslashes) is that ANY backslash-bearing edit goes
through a Write-tool script with raw strings.  Third violation this session;
writing it down so the pattern is visible rather than repeated.

Run:  python debug/qa/_mark_p13_caption_withdrawal.py
"""
from __future__ import annotations

import io
import sys

PATH = "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex"

OLD = (r"(\textbf{[MEASURED 2026-09-19]}; the $0.0002\%$ this cell carried until"
       "\n"
       r"then understated the method by roughly seven orders).")

NEW = (r"(\textbf{[MEASURED 2026-09-19]}"
       "\n"
       r"[retracted 2026-09-19: p11-h2plus-0002pct-retired]; the $0.0002\%$ this"
       "\n"
       r"cell carried until then understated the method by roughly seven orders).")


def main() -> None:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

    s = io.open(PATH, encoding="utf-8").read()
    if "[retracted 2026-09-19: p11-h2plus-0002pct-retired]" in s:
        print("marker already present; nothing to do")
        return
    if OLD not in s:
        # Report the real text so the next attempt is anchored, not guessed.
        i = s.find("this cell carried until")
        print("ANCHOR NOT FOUND. Actual surrounding text:")
        print(repr(s[max(0, i - 260):i + 160]) if i >= 0
              else "  (phrase 'this cell carried until' absent entirely)")
        raise SystemExit(1)
    if s.count(OLD) != 1:
        raise SystemExit(f"anchor is ambiguous ({s.count(OLD)} matches)")

    io.open(PATH, "w", encoding="utf-8", newline="").write(s.replace(OLD, NEW, 1))
    print("marker added to Paper 13 caption")

    # Confirm the C17 family now treats this locus as withdrawn, not live.
    import importlib.util
    import re
    spec = importlib.util.spec_from_file_location(
        "chn", "debug/qa/check_headline_numbers.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    e = [x for x in mod.REGISTRY
         if x["id"] == "p11-h2plus-0002pct-retired"][0]
    pat, req, exm = (re.compile(e["pattern"]),
                     re.compile(e["require_nearby"]),
                     re.compile(e["exempt_if_nearby"]))
    win = getattr(mod, "WINDOW", 3)
    lines = io.open(PATH, encoding="utf-8").read().splitlines()
    still_live = []
    for k, ln in enumerate(lines):
        if not pat.search(ln):
            continue
        ctx = "\n".join(lines[max(0, k - win):k + win + 1])
        if req.search(ctx) and not exm.search(ctx):
            still_live.append(k + 1)
    print(f"  paper_13 live occurrences now: {still_live or 'none'}")


if __name__ == "__main__":
    main()
