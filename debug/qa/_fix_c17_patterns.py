r"""Rewrite the two C17 families I just registered with CORRECT regexes.

WHY THIS FILE EXISTS.  I wrote both patterns inside a bash heredoc, so every
backslash was halved/doubled on the way in and the stored patterns were
`0\\.0002` -- two literal backslashes, matching nothing.  The gate then
reported `[ok] ... clean` on 27 live loci.

This is the THIRD instance of that exact bug in this corpus: the
`brown-surjection-attribution` pattern was killed by doubled escapes on
2026-08-22, and a backspace-corrupted `pauli-advantage-floor` regex sat dead
for an unknown period.  Both were found only because someone tested them.
The standing rule (memory: no-heredoc-backslashes) says any backslash edit
goes through a Write-tool script file -- which is what this is.

Run:  python debug/qa/_fix_c17_patterns.py
Then: fire-test BOTH directions before trusting either family.
"""
from __future__ import annotations

import io
import re

PATH = "debug/qa/check_headline_numbers.py"

# The corpus writes the H2+ claim as `$0.0002\%$`, `0.0002\%`, and bare
# `0.0002` in table cells.  So: the literal digits, then optional space,
# then an optional LaTeX-escaped or bare percent sign.
GOOD_P11_PATTERN = r'r"0\.0002\s*\\?%?"'
GOOD_P11_REQUIRE = (
    r'r"H\$?_2\^?\{?\+|H2\+|prolate|spectral Laguerre|Laguerre|'
    r'paper11|loutey_paper11|Bates|0\.6026"'
)
GOOD_P11_EXEMPT = (
    r'r"retired|superseded|RETIRED|withdrawn|corrected 2026-09-19|'
    r'machine precision|reference-limited|understates|was wrong by"'
)

GOOD_P13_PATTERN = r'r"0\.019\s*\\?%"'
GOOD_P13_REQUIRE = (
    r'r"He|helium|2D variational|cusp|hyperspherical|paper13|'
    r'loutey_paper13|Track DI"'
)
GOOD_P13_EXEMPT = (
    r'r"retired|superseded|withdrawn|corrected 2026-09-19|'
    r'matches no He result|0\.022|0\.004"'
)


def _swap(src: str, entry_id: str, field: str, good: str) -> str:
    """Replace one field of one registry entry, located by its id."""
    i = src.find(f'"id": "{entry_id}"')
    if i < 0:
        raise SystemExit(f"entry {entry_id} not found")
    # end of this dict entry
    j = src.find('\n    },\n', i)
    block = src[i:j]
    m = re.search(rf'"{field}":\s*(r?"(?:[^"\\]|\\.)*"(?:\s*\n\s*r?"(?:[^"\\]|\\.)*")*)',
                  block)
    if not m:
        raise SystemExit(f"field {field} not found in {entry_id}")
    new_block = block[:m.start(1)] + good + block[m.end(1):]
    return src[:i] + new_block + src[j:]


def main() -> None:
    s = io.open(PATH, encoding="utf-8").read()
    for eid, field, good in (
        ("p11-h2plus-0002pct-retired", "pattern", GOOD_P11_PATTERN),
        ("p11-h2plus-0002pct-retired", "require_nearby", GOOD_P11_REQUIRE),
        ("p11-h2plus-0002pct-retired", "exempt_if_nearby", GOOD_P11_EXEMPT),
        ("p13-he-0019pct-nonexistent", "pattern", GOOD_P13_PATTERN),
        ("p13-he-0019pct-nonexistent", "require_nearby", GOOD_P13_REQUIRE),
        ("p13-he-0019pct-nonexistent", "exempt_if_nearby", GOOD_P13_EXEMPT),
    ):
        s = _swap(s, eid, field, good)
    io.open(PATH, "w", encoding="utf-8", newline="").write(s)

    # prove the stored regexes are what we meant, by importing and compiling
    import importlib.util
    spec = importlib.util.spec_from_file_location("chn", PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    for e in mod.REGISTRY:
        if e["id"] in ("p11-h2plus-0002pct-retired", "p13-he-0019pct-nonexistent"):
            pat = e["pattern"]
            print(f'  {e["id"]}: pattern={pat!r}')
            assert "\\\\." not in pat, "DOUBLED ESCAPE STILL PRESENT"
            rx = re.compile(pat)
            # the two directions, checked here so a silent family cannot ship
            live = {"p11-h2plus-0002pct-retired": r"achieving $0.0002\%$ accuracy with $N_b = 20$.",
                    "p13-he-0019pct-nonexistent": r"$0.019\%$ error (Track~DI v2.6.0)"}[e["id"]]
            fixed = {"p11-h2plus-0002pct-retired": "reproduced to machine precision at $N_b = 20$.",
                     "p13-he-0019pct-nonexistent": r"$0.022\%$ raw error (Track~DI)"}[e["id"]]
            assert rx.search(live), f"{e['id']} does NOT fire on the retired wording"
            assert not rx.search(fixed), f"{e['id']} fires on the CORRECTED wording"
            print("    FIRES on retired wording: yes;  silent on corrected: yes")
    print("both families rewritten and discrimination-checked")


if __name__ == "__main__":
    main()
