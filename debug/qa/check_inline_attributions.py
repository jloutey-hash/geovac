"""C20 -- inline attribution resolvability (ratchet gate).

WHY THIS EXISTS
---------------
The 2026-08-28/29 QA arc found that every external-citation defect in Paper 34
lived in the same place: inline "Author Year" attributions carrying Layer-2
input values with NO corresponding \\bibitem.  A reader cannot resolve them
from the paper alone, and a reviewer cannot check them without guessing which
work is meant.  That layer produced a phantom date ("Eides 2024" -- no such
edition), a wrong section anchor recurring eight times, numbers credited to a
paper that does not contain them, and an ambiguity between two different works
by the same author pair.

The class is mechanical, so it gets a mechanical guard.

WHAT IT DOES
------------
For each gated paper: extract inline attributions of the form `Surname~YYYY`
or `Surname--Surname~YYYY`, and mark each RESOLVED if some \\bibitem in the
same file mentions that surname and that year (or a key encoding both), else
UNRESOLVED.

RATCHET SEMANTICS
-----------------
Fixing ~30 attributions at once is not realistic, and a gate that always fails
is a gate people learn to ignore.  So the check stores a per-paper baseline
SET of unresolved tags in `debug/qa/inline_attribution_baseline.json` and
FAILS only when a paper acquires an attribution that is NOT in its
baseline.  Progress is locked in (re-baseline when you improve); a new
unresolved attribution is caught immediately and named explicitly.

Usage:
  python debug/qa/check_inline_attributions.py --gate group6
  python debug/qa/check_inline_attributions.py --gate group6 --update-baseline
  python debug/qa/check_inline_attributions.py --selftest
"""
from __future__ import annotations

import json
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
PAPERS = ROOT / "papers"
BASELINE = ROOT / "debug" / "qa" / "inline_attribution_baseline.json"

# `Surname~1999`, `Surname--Surname~1999`, `Surname et al.~1999`
ATTR = re.compile(
    r"\b([A-Z][a-zA-Z]{2,}(?:--[A-Z][a-zA-Z]{2,})*)"      # surname(s)
    r"(?:~| )(?:\\emph\{et~al\.\}[~ ])?"
    r"(1[89]\d{2}|20\d{2})\b"
)

# things that look like attributions but are not external citations
MONTHS = {"January", "February", "March", "April", "May", "June", "July",
          "August", "September", "October", "November", "December"}
NOT_AUTHORS = MONTHS | {
    "Sprint", "Track", "Paper", "Section", "Table", "Figure", "Appendix",
    "CLAUDE", "GeoVac", "Corrected", "Withdrawn", "Retired", "Added",
    "Verified", "Note", "Status", "Updated", "Since", "During", "Phase",
    "CODATA", "NIST", "Zenodo", "Observation", "Prediction", "Remark",
    # correction-note verbs that precede a date in this corpus
    "Noted", "Noted", "RETRACTED", "Recounted", "Retracted", "Withdrawn",
    "Diagnosed", "Corrected", "Sourced", "Resolved", "Flagged", "Scoped",
    "Promoted", "Demoted", "Reviewed", "Certified", "Applied", "Fixed",
}


def bibitem_index(text: str):
    """Return (keys, bodies) for every \\bibitem in the file."""
    keys, bodies = [], []
    for m in re.finditer(r"\\bibitem\{([^}]*)\}", text):
        keys.append(m.group(1))
        nxt = text.find("\\bibitem{", m.end())
        bodies.append(text[m.end(): nxt if nxt > 0 else len(text)])
    return keys, bodies


def resolved(surnames: str, year: str, keys, bodies) -> bool:
    parts = [p for p in surnames.split("--") if p]
    for k, b in zip(keys, bodies):
        kl = k.lower()
        if year in b and any(p in b for p in parts):
            return True
        if year in kl and any(p.lower() in kl for p in parts):
            return True
    return False


def _strip_comments(text: str) -> str:
    """Remove LaTeX comments; an unescaped % starts one.

    Comments do not render, so they cannot carry an attribution.  Without
    this the gate read editorial notes as citations -- it FAILed on
    `% Reconstructed 2026-08-29 ...`, a provenance note on line 1 of a
    paper.  An escaped \\% is preserved, since that one does render.
    """
    out = []
    for line in text.split("\n"):
        i, n = 0, len(line)
        cut = None
        while i < n:
            if line[i] == "\\":
                i += 2
                continue
            if line[i] == "%":
                cut = i
                break
            i += 1
        out.append(line if cut is None else line[:cut])
    return "\n".join(out)


def scan(path: pathlib.Path):
    text = _strip_comments(path.read_text(encoding="utf-8", errors="replace"))
    body = text.split("\\begin{thebibliography}")[0]     # exclude the bib itself
    keys, bodies = bibitem_index(text)
    # classify each DISTINCT attribution once (an earlier version
    # incremented the unresolved count on every repeat occurrence without
    # re-checking, so any multi-use attribution looked unresolved)
    seen_ok, seen_bad = set(), set()
    for m in ATTR.finditer(body):
        sur, yr = m.group(1), m.group(2)
        if sur.split("--")[0] in NOT_AUTHORS:
            continue
        tag = f"{sur}~{yr}"
        if tag in seen_ok or tag in seen_bad:
            continue
        (seen_ok if resolved(sur, yr, keys, bodies) else seen_bad).add(tag)
    # count uses of the unresolved ones only
    unresolved = {}
    for tag in seen_bad:
        sur, yr = tag.split("~")
        unresolved[tag] = len(re.findall(
            re.escape(sur) + r"(?:~| )(?:\\emph\{et~al\.\}[~ ])?" + yr, body))
    return unresolved


def _gate_substr(argv):
    if "--gate" in argv:
        i = argv.index("--gate")
        if i + 1 < len(argv):
            return argv[i + 1]
    return None


def selftest() -> int:
    doc = (r"Per Smith~1999 and Jones--Brown~2004 the value is 3." "\n"
           r"Smith~1999 again, and Smith~1999 once more." "\n"
           r"\begin{thebibliography}{9}" "\n"
           r"\bibitem{smith1999} A. Smith, Phys. Rev. A 1, 1 (1999)." "\n"
           r"\end{thebibliography}")
    tmp = ROOT / "debug" / "qa" / "_c20_selftest.tex"
    tmp.write_text(doc, encoding="utf-8")
    try:
        u = scan(tmp)
        ok = ("Jones--Brown~2004" in u) and ("Smith~1999" not in u)
        print("selftest:", "PASS" if ok else "FAIL", "->", sorted(u))
        return 0 if ok else 1
    finally:
        tmp.unlink(missing_ok=True)


def main() -> int:
    if "--selftest" in sys.argv:
        return selftest()
    gate = _gate_substr(sys.argv)
    update = "--update-baseline" in sys.argv
    base = json.loads(BASELINE.read_text(encoding="utf-8")) if BASELINE.exists() else {}

    files = sorted(p for p in PAPERS.rglob("*.tex") if "archive" not in p.parts)
    if gate:
        files = [p for p in files if gate in str(p).replace("\\", "/")]

    print(f"\ninline attribution check over {len(files)} paper(s)"
          f"   [scope: {gate or 'ALL'}]")
    failed, new_base = [], dict(base)
    for p in files:
        rel = str(p.relative_to(ROOT)).replace("\\", "/")
        u = scan(p)
        tags = sorted(u)
        prev = base.get(rel)
        prev_set = set(prev) if isinstance(prev, list) else None
        new_base[rel] = tags
        flag, added = "", []
        if prev_set is not None:
            added = sorted(set(tags) - prev_set)
            gone = len(prev_set - set(tags))
            if added:
                flag = f"  *** {len(added)} NEW ***"
                failed.append((rel, added))
            elif gone:
                flag = f"  (improved: {gone} resolved since baseline)"
        if tags:
            print(f"  {rel}: {len(tags)} unresolved{flag}")
            for tag, cnt in sorted(u.items(), key=lambda kv: -kv[1])[:6]:
                print(f"      {tag}  ({cnt} use{'s' if cnt > 1 else ''})")
            if len(tags) > 6:
                print(f"      ... and {len(tags) - 6} more")

    if update:
        BASELINE.write_text(json.dumps(new_base, indent=2), encoding="utf-8")
        print(f"\nbaseline written -> {BASELINE.relative_to(ROOT)}")
        return 0

    if failed:
        print("\n*** NEW UNRESOLVED ATTRIBUTIONS (not in baseline) ***")
        for rel, added in failed:
            for tag in added:
                print(f"  {rel}:  {tag}  has no resolvable \\bibitem")
        print("\nRESULT: FAIL (a new inline attribution has no resolvable \\bibitem)")
        return 1
    print("\nRESULT: PASS (no paper carries an unresolved attribution outside its baseline)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
