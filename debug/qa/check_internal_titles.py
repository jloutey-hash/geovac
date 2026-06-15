#!/usr/bin/env python3
"""
Deterministic internal-title consistency check  --  the /qa "internal-title" criterion.

Verifies that every internal GeoVac citation names the cited paper by its CURRENT
real ``\title{}``. This REPLACES the unreliable LLM-reviewer check for the
internal-title class (the /qa run-#1 blind spot, 2026-06-14): title drift is a
string-comparison problem, so a deterministic script catches it more reliably than
a prose reviewer, and the LLM reviewers are thereby narrowed to genuine judgment
calls.

Method: build a title map from each active paper's ``\title{}`` (full title and
main-title-without-subtitle); scan every active .tex for internal GeoVac
bibitem citations (a quoted title immediately before "GeoVac Paper(s) N (YEAR");
normalise (expand GeoVac geometry macros, strip LaTeX, lower, alphanumerics) and
match the cited title against the full OR main title.

Exit 0 = clean. Exit 1 = at least one MISMATCH outside the known exceptions.
Exclusions (by design): house-style capsule files (deliberate short descriptors);
archived papers (out of scope). The descope-pending propinquity cluster is
FLAGGED, not failed.

Usage:  python debug/qa/check_internal_titles.py
"""
from __future__ import annotations

import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
PAPERS = ROOT / "papers"

HOUSE_STYLE = {
    "geovac_field_guide.tex",
    "paper_55_periods_of_geovac.tex",
    "paper_57_forced_free_seam.tex",
    "group1_operator_algebras_synthesis.tex",
}
PROPINQUITY = {39, 40, 45, 46, 47, 48, 49}   # descope-pending titles -> flag, don't fail
ARCHIVED = {3, 4, 5, 6, 10, 21}              # out of scope

# GeoVac title macros that carry semantic content lost by a blind \word strip.
MACROS = [
    (r"\sthree", "S3"), (r"\sfive", "S5"), (r"\Uone", "U1"),
    (r"\SU", "SU"), (r"\SL", "SL"), (r"\Ga", "Ga"), (r"\Aut", "Aut"),
    (r"\R", "R"), (r"\N", "N"),
]


def norm(s: str) -> str:
    s = re.sub(r"\\\\", " ", s)                         # \\ line break -> space
    for mac, rep in MACROS:
        s = re.sub(re.escape(mac) + r"(?![a-zA-Z])", rep, s)
    s = re.sub(r"\\[a-zA-Z]+", " ", s)                  # remaining control sequences
    s = s.replace("^", "").replace("$", "")             # S^3 -> S3 (no space)
    s = re.sub(r"[{}~\\]", " ", s)                      # braces, ties, stray backslash
    s = re.sub(r"[^a-z0-9]+", " ", s.lower())
    return " ".join(s.split())


def main_part(title: str) -> str:
    """Title with a trailing {\\large ...} / {\\Large ...} subtitle removed."""
    return re.split(r"\{\s*\\[lL]arge\b", title, maxsplit=1)[0]


def active_tex():
    return sorted([*PAPERS.glob("group*/*.tex"), *PAPERS.glob("synthesis/*.tex")])


def build_title_map() -> dict:
    """{paper_number_or_FCI: (full_norm, main_norm, raw_title, filename)}."""
    tmap = {}
    for f in active_tex():
        txt = f.read_text(encoding="utf-8", errors="replace")
        m = re.search(r"\\title\{((?:[^{}]|\{[^{}]*\})*)\}", txt)
        if not m:
            continue
        title = m.group(1)
        entry = (norm(title), norm(main_part(title)), " ".join(title.split()), f.name)
        num = re.search(r"[Pp]aper[_ ](\d+)[_ ]", f.name)
        if num:
            tmap[int(num.group(1))] = entry
        elif "fci_atoms" in f.name:
            tmap["FCI-A"] = entry
        elif "fci_molecules" in f.name:
            tmap["FCI-M"] = entry
    return tmap


# Title quote that cannot span across another '' (prevents an in-text ``quote''
# from being glued to a later bibitem). Must sit immediately (whitespace only,
# optional \newblock) before "GeoVac Paper(s) N (YEAR" -- the year-paren excludes
# in-text "GeoVac Paper N~\cite{}" references.
QUOTE = r"``(?P<title>(?:(?!'')[\s\S])*?),?''"
NORMAL = re.compile(
    QUOTE + r"\s*(?:\\newblock\s*)?GeoVac\s+Papers?~?\s*(?P<n>\d+)(?:--\d+)?\s*\(",
    re.DOTALL,
)
ZENODO = re.compile(r"``\s*GeoVac\s+Paper~?\s*(?P<n>\d+):\s*(?P<title>(?:(?!'')[\s\S])*?),?''", re.DOTALL)


def check_one(fname, n, cited_title, tmap, mismatches, flagged):
    if n in ARCHIVED or n not in tmap:
        return
    full_norm, main_norm, raw, _ = tmap[n]
    cn = norm(cited_title)
    if cn == full_norm or cn == main_norm:
        return
    rec = {"file": fname, "paper": n, "cited": " ".join(cited_title.split())[:80], "real": raw[:80]}
    (flagged if n in PROPINQUITY else mismatches).append(rec)


def main() -> int:
    tmap = build_title_map()
    mismatches, flagged = [], []
    for f in active_tex():
        if f.name in HOUSE_STYLE:
            continue
        txt = f.read_text(encoding="utf-8", errors="replace")
        for m in ZENODO.finditer(txt):
            check_one(f.name, int(m.group("n")), m.group("title"), tmap, mismatches, flagged)
        for m in NORMAL.finditer(txt):
            t = m.group("title")
            if t.strip().startswith("GeoVac Paper"):
                continue
            check_one(f.name, int(m.group("n")), t, tmap, mismatches, flagged)

    print(f"title map: {len([k for k in tmap if isinstance(k,int)])} numbered papers + "
          f"{[k for k in tmap if not isinstance(k,int)]}")
    if flagged:
        print(f"\nFLAGGED (propinquity cluster, descope-pending -- not a failure): {len(flagged)}")
        for r in flagged:
            print(f"  [P{r['paper']}] {r['file']}: cited \"{r['cited']}\"")
    if mismatches:
        print(f"\n*** MISMATCH ({len(mismatches)}) -- internal cite does not match the paper's \\title: ***")
        for r in mismatches:
            print(f"  [P{r['paper']}] {r['file']}")
            print(f"      cited: \"{r['cited']}\"")
            print(f"      real : \"{r['real']}\"")
        print("\nRESULT: FAIL (internal-title criterion)")
        return 1
    print("\nRESULT: PASS -- every internal GeoVac citation matches the cited paper's \\title.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
