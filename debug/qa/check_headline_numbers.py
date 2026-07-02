"""C17 -- headline-number registry gate (added 2026-07-02, PI direction).

The 7th group4 cert's meta-lesson: the MATERIAL classes that kept surviving
judgment review are MECHANICAL --
  (a) second-locus propagation: a decided headline value corrected at one
      locus and stale at another (Z=1-36 vs 56 across five loci; the
      memo-listed-but-never-applied KH fix), and
  (b) number-vs-source drift: a stated value contradicting its own cited
      source (the "190x" floor vs the cited table's 51x; the stale 33.3
      1-norm vs the live 32.6).
Both are registry-checkable. Each entry holds a headline FAMILY: either a
set of known-wrong variant patterns (C16 style) or a capture pattern plus
the CANONICAL value (any capture that disagrees is a live hit). Exempt
markers cover legitimately historical/withdrawn mentions.

MAINTENANCE RULE (mirrors C16): when a cert run corrects or demotes a
headline number, ADD/UPDATE its family here so the wrong value can never
silently re-surface at any locus.

Usage: python debug/qa/check_headline_numbers.py [--gate <branch>] [--all]
Exit 0 = PASS. Mirror test: tests/test_headline_numbers_check.py.
"""
from __future__ import annotations

import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
WINDOW = 3  # +/- lines for the exemption window

GROUP4_FILES = [
    "papers/group4_quantum_computing/*.tex",
    "papers/synthesis/group4_quantum_computing_synthesis.tex",
]

REGISTRY = [
    {
        "id": "library-z-span",
        "scope": "group4",
        "severity": "fail",
        "canonical_note": "library span = Z=1--56 (H through Ba; SrH Z=38, BaH "
                          "Z=56 registry-probed). Decided v4.58.0 M-C; the 6th "
                          "cert found 5 stale Z=1--36 loci (second-locus class).",
        "capture": r"Z\s*=?\s*1\s*\$?\s*--\s*(\d{2})",
        "canonical": "56",
        # only the LIBRARY-span statements are in this family; bare periodic-row
        # prose ("First-row (Z=1--10) atoms...") is a different, legitimate quantity
        "require_nearby": r"librar|spanning|systems|H\s+through\s+Ba",
        "exempt_if_nearby": r"historical|previously|was\s+corrected|stale",
        "files": GROUP4_FILES,
    },
    {
        "id": "pauli-advantage-floor",
        "scope": "group4",
        "severity": "fail",
        "canonical_note": "raw-JW Pauli advantage floor = 51x (equal-qubit "
                          "table: 51/746/1712; P20: '51--1,712x'). The 7th cert "
                          "found a drifted '190x--1,712x' floor (M2).",
        "capture": r"(\d{2,4})\s*\$?\\times\$?\s*--\s*1\{?,\}?712",
        "canonical": "51",
        "exempt_if_nearby": r"historical|previously|stale",
        "files": GROUP4_FILES,
    },
    {
        "id": "library-count",
        "scope": "group4",
        "severity": "fail",
        "canonical_note": "library = 37 systems (35 composed + He + H2), "
                          "decided PI 2026-06-28; retired wrong counts 28/30/38/40.",
        "pattern": r"\b(?:28|30|38|40)\s+systems\b",
        "exempt_if_nearby": r"was|stale|historical|corrected|retired|previously",
        "files": GROUP4_FILES,
    },
    {
        "id": "balanced-lih-binds-at-3015",
        "scope": "group4",
        "severity": "fail",
        "canonical_note": "balanced LiH binds at the COMPUTED R_eq=3.227 bohr "
                          "(7.0% above the experimental 3.015); 'binds at 3.015' "
                          "was the v4.56.0 M2 finding (recurred v4.57.0 in the "
                          "synthesis).",
        "pattern": r"binds[^.\n]{0,60}3\.015",
        "exempt_if_nearby": r"experimental|7\.0\s*\\?%|above|computed",
        "files": GROUP4_FILES,
    },
    {
        "id": "lih-onenorm-stale",
        "scope": "group4",
        "severity": "fail",
        "canonical_note": "composed LiH 1-norm = 32.6 Ha live (0.95x vs STO-3G "
                          "34.3); the stale 33.3 / 0.97x pair retired v4.60.0 "
                          "(PI-directed corpus-wide 2026-07-01). 0.97 is scoped "
                          "to 1-norm proximity (the l-parity Pauli-ratio 0.97 "
                          "cells are a different, legitimate quantity).",
        "pattern": r"\b33\.3\b\s*~?Ha|\\lambda\s*=\s*33\.3"
                   r"|1-norm[^.\n]{0,40}\b0\.97\b|\b0\.97\b\$?\\times\$?[^.\n]{0,25}1-norm",
        "exempt_if_nearby": r"historical|stale|rested\s+on|retired",
        "files": GROUP4_FILES,
    },
]


def _gate_substr(argv: "list[str]") -> "str | None":
    for i, a in enumerate(argv):
        if a.startswith("--gate="):
            return a.split("=", 1)[1]
        if a == "--gate" and i + 1 < len(argv):
            return argv[i + 1]
    return None


def _resolve(globs: "list[str]") -> "list[pathlib.Path]":
    out: "list[pathlib.Path]" = []
    for g in globs:
        out.extend(sorted(ROOT.glob(g)))
    seen, uniq = set(), []
    for p in out:
        if p not in seen and p.is_file():
            seen.add(p)
            uniq.append(p)
    return uniq


def scan_entry(entry: dict, text_override: "str | None" = None):
    """Return (live_hits, exempt_hits); item = (relpath, line_no, snippet).

    text_override: scan the given text as a single pseudo-file (self-test hook).
    """
    exempt = re.compile(entry["exempt_if_nearby"], re.IGNORECASE)
    require = (re.compile(entry["require_nearby"], re.IGNORECASE)
               if "require_nearby" in entry else None)
    if "pattern" in entry:
        pat = re.compile(entry["pattern"], re.IGNORECASE)
        is_wrong = lambda m: True  # any match of a wrong-variant pattern
    else:
        pat = re.compile(entry["capture"], re.IGNORECASE)
        canonical = entry["canonical"]
        is_wrong = lambda m: m.group(1) != canonical

    def scan_lines(lines, rel):
        live, ok = [], []
        for i, line in enumerate(lines):
            m = pat.search(line)
            if not m or not is_wrong(m):
                continue
            lo, hi = max(0, i - WINDOW), min(len(lines), i + WINDOW + 1)
            window_txt = "\n".join(lines[lo:hi])
            if require is not None and not require.search(window_txt):
                continue  # outside this family's context (different quantity)
            snip = re.sub(r"\s+", " ", line.strip())[:160]
            (ok if exempt.search(window_txt) else live).append((rel, i + 1, snip))
        return live, ok

    if text_override is not None:
        return scan_lines(text_override.splitlines(), "<override>")

    live_all, ok_all = [], []
    for path in _resolve(entry["files"]):
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        live, ok = scan_lines(lines, path.relative_to(ROOT))
        live_all.extend(live)
        ok_all.extend(ok)
    return live_all, ok_all


def main() -> int:
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    gate = _gate_substr(sys.argv)
    scope = f"scope '{gate}'" if gate else "ALL entries"

    def selected(e: dict) -> bool:
        return gate is None or e["scope"] == "all" or gate in e["scope"]

    n_live, n_exempt = 0, 0
    print(f"headline-number registry gate (C17)   [{scope}]\n")
    for e in REGISTRY:
        if not selected(e):
            continue
        live, ok = scan_entry(e)
        n_live += len(live)
        n_exempt += len(ok)
        status = "clean" if not live else f"{len(live)} LIVE"
        print(f"  [{'FAIL' if live else 'ok'}] {e['id']}: {status}"
              + (f"  (exempt: {len(ok)})" if ok else ""))
        for rel, ln, snip in live:
            print(f"      {rel}:{ln}  {snip}")

    if n_live:
        print(f"\nRESULT: FAIL ({n_live} live wrong-headline occurrence(s) in {scope})")
        return 1
    print(f"\nRESULT: PASS (no live wrong headline value in {scope}"
          + (f"; {n_exempt} exempt/historical mention(s))" if n_exempt else ")"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
