r"""Narrow both new C17 families' `exempt_if_nearby` to per-entry markers.

THE DEFECT THIS FIXES, measured rather than suspected.  My first version of
both families used a generic exemption alternation
(`retired|superseded|withdrawn|...`).  A probe over the declared files found it
excusing two LIVE loci on words that withdraw nothing:

    docs/paper_notes_archive.md:108    excused on 'withdrawn'
    docs/project_closeout_plan.md:26   excused on 'superseded'

Both are incidental vocabulary in NEIGHBOURING rows of archive tables -- within
WINDOW = 3 lines of the hit, about entirely different claims.

This is the 2026-09-04 DELTA #5 incident replayed exactly.  That run found a
bare `[retracted YYYY-MM-DD]` token silencing every entry in its window,
including a live volume-quotient defect in the next sentence, and the registry
was redesigned so the marker CARRIES THE ENTRY ID it withdraws.  The same
record states the general rule: authoring exemption vocabulary is the failure
mode, because correct text is exactly what surrounds a defect.  Three entries
written on 2026-09-03 exempted on words drawn from their own corrected prose
("misnomer", "corrected 2026-09-03", "bound|deficit|saturat") and each
reported clean on a live locus.

I wrote generic vocabulary anyway, one turn after reading that rule.

FIX.  Exempt on (a) the standardized per-entry marker, and (b) the specific
corrected phrasings this sweep will introduce -- nothing generic, nothing that
could plausibly appear near an unrelated claim.

Run:  python debug/qa/_narrow_c17_exemptions.py
"""
from __future__ import annotations

import io
import re
import sys

PATH = "debug/qa/check_headline_numbers.py"

# (a) the per-entry marker, id-carrying, per withdrawal_marker() in C16;
# (b) the exact corrected phrasings, which name THIS quantity and no other.
NARROW = {
    "p11-h2plus-0002pct-retired":
        r'r"\[retracted \d{4}-\d{2}-\d{2}:\s*p11-h2plus-0002pct-retired\]|'
        r'reproduced to machine precision|reference-limited residual|'
        r'wrong by ~?7\.6 orders|retired 0\.0002"',
    "p13-he-0019pct-nonexistent":
        r'r"\[retracted \d{4}-\d{2}-\d{2}:\s*p13-he-0019pct-nonexistent\]|'
        r'matches no He result|retired 0\.019"',
}


def _swap(src: str, entry_id: str, good: str) -> str:
    i = src.find(f'"id": "{entry_id}"')
    if i < 0:
        raise SystemExit(f"entry {entry_id} not found")
    j = src.find("\n    },\n", i)
    block = src[i:j]
    m = re.search(
        r'"exempt_if_nearby":\s*(r?"(?:[^"\\]|\\.)*"(?:\s*\n\s*r?"(?:[^"\\]|\\.)*")*)',
        block)
    if not m:
        raise SystemExit(f"exempt_if_nearby not found in {entry_id}")
    return src[:i] + block[:m.start(1)] + good + block[m.end(1):] + src[j:]


def main() -> None:
    # stdout may be cp1252 here; the corpus contains H₂⁺ and similar.
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

    s = io.open(PATH, encoding="utf-8").read()
    for eid, good in NARROW.items():
        s = _swap(s, eid, good)
    io.open(PATH, "w", encoding="utf-8", newline="").write(s)

    import importlib.util
    spec = importlib.util.spec_from_file_location("chn", PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    win = getattr(mod, "WINDOW", 3)

    for eid in NARROW:
        e = [x for x in mod.REGISTRY if x["id"] == eid][0]
        assert "\\\\." not in e["exempt_if_nearby"], f"DOUBLED ESCAPE in {eid}"
        pat = re.compile(e["pattern"])
        req = re.compile(e["require_nearby"])
        exm = re.compile(e["exempt_if_nearby"])
        live, exempt = [], []
        for f in e["files"]:
            try:
                lines = io.open(f, encoding="utf-8",
                                errors="replace").read().splitlines()
            except FileNotFoundError:
                continue
            for k, ln in enumerate(lines):
                if not pat.search(ln):
                    continue
                ctx = "\n".join(lines[max(0, k - win):k + win + 1])
                if not req.search(ctx):
                    continue
                (exempt if exm.search(ctx) else live).append(f"{f}:{k+1}")
        print(f"{eid}:  live={len(live)}  exempt={len(exempt)}")
        for locus in exempt:
            print(f"    still exempt: {locus}   <-- must be a REAL withdrawal")
        # the marker must still work, or the sweep cannot quiet the gate
        probe = f"foo [retracted 2026-09-19: {eid}] bar"
        assert exm.search(probe), f"{eid}: per-entry marker does not exempt!"
    print("exemptions narrowed; per-entry markers verified functional")


if __name__ == "__main__":
    main()
