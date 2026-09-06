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
archived papers (out of scope). The descope-pending Lorentzian propinquity cluster
(45-49) is FLAGGED, not failed; Papers 38/39/40 are now descoped + cite-converged
and are ENFORCED (a title/cite mismatch on them FAILS).

Usage:  python debug/qa/check_internal_titles.py
"""
from __future__ import annotations

import pathlib
import re
import sys
import os

# Shared --gate scope resolution (see debug/qa/qa_scopes.py): named
# scopes resolve to an explicit file list and every RESULT line carries
# the file count, so a gate can never report PASS on an empty scope.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qa_scopes  # noqa: E402


ROOT = pathlib.Path(__file__).resolve().parents[2]
PAPERS = ROOT / "papers"

HOUSE_STYLE = {
    "geovac_field_guide.tex",
    "paper_55_periods_of_geovac.tex",
    "paper_57_forced_free_seam.tex",
    "group1_operator_algebras_synthesis.tex",
}
# The propinquity cluster (38/39/40 Euclidean + 45-49 Lorentzian) is now fully
# descoped, cite-consistent, and ENFORCED (empty exception): any title/cite drift on
# these papers FAILS rather than flags. Papers 38/39/40 (clean "Latremoliere
# propinquity -> van Suijlekom state-space GH" descopes) were enforced across
# v4.14.x / v4.20.0-2; the Lorentzian leg 45-49 (degeneracy / partial descopes)
# became enforceable once the U(1)/\Uone notation-normalization unification above let
# their literal cites match the macro-form titles (v4.20.4). Re-populate this set only
# for a NEW paper whose title is mid-descope.
PROPINQUITY: set = set()
ARCHIVED = {3, 4, 5, 6, 10, 21}              # out of scope

# GeoVac title macros that carry semantic content lost by a blind \word strip.
MACROS = [
    (r"\sthree", "S3"), (r"\sfive", "S5"), (r"\Uone", "U 1"),
    (r"\SU", "SU"), (r"\SL", "SL"), (r"\Ga", "Ga"), (r"\Aut", "Aut"),
    (r"\R", "R"), (r"\N", "N"),
]
# \Uone -> "U 1" (SEPARATED) so the macro normalises identically to the literal
# parenthetical form "U(1)" -> "u 1"; this unifies the three notations that the
# same paper's title and its citers use ($\SU(2)\times\Uone_T$ macro,
# "SU(2) U(1)_T" literal, "\mathrm{SU}(2)\times\mathrm{U}(1)_T"). Before this fix
# the three forms normalised to "su 2 1 t" / "su 2 u 1 t" / "su 2 u 1 t" and the
# Lorentzian-cluster cites (Papers 45-49) could never match their macro titles.


def norm(s: str) -> str:
    s = re.sub(r"\\\\", " ", s)                         # \\ line break -> space
    # Strip math-font wrappers, keeping their (space-padded) contents, so
    # \mathrm{SU} normalises the same as the bare \SU macro and the literal "SU".
    s = re.sub(
        r"\\(?:mathrm|mathbf|mathit|mathsf|mathcal|text|emph|operatorname)\s*"
        r"\{([^{}]*)\}", r" \1 ", s)
    # Canonical macros. Replacements are SPACE-PADDED so that an adjacent control
    # sequence (e.g. "\times\Uone" -> "\times U 1") is NOT swallowed by the generic
    # control-sequence strip below (the bug that dropped the "U" of \Uone, giving
    # "su 2 1 t" instead of "su 2 u 1 t").
    for mac, rep in MACROS:
        s = re.sub(re.escape(mac) + r"(?![a-zA-Z])", f" {rep} ", s)
    s = re.sub(r"\\[a-zA-Z]+", " ", s)                  # remaining control sequences
    s = s.replace("^", "").replace("$", "")             # S^3 -> S3 (no space)
    s = re.sub(r"[{}~\\]", " ", s)                      # braces, ties, stray backslash
    s = re.sub(r"[^a-z0-9]+", " ", s.lower())
    return " ".join(s.split())


def main_part(title: str) -> str:
    """Main title with a trailing subtitle removed.

    Handles both subtitle styles in the corpus: a ``{\\large ...}`` block and a
    ``Main Title:\\ Subtitle`` / ``Main Title: Subtitle`` colon separator. A cite
    by the main title alone (subtitle dropped) is legitimate, so main_norm must
    not carry the subtitle. Stripping only makes main_norm shorter, never adds a
    false match (the full title is still matched separately)."""
    t = re.split(r"\{\s*\\[lL]arge\b", title, maxsplit=1)[0]
    t = re.split(r":(?=[\s\\])", t, maxsplit=1)[0]  # ": " / ":\\ " / ":\ " subtitle separators
    return t


# Two in-corpus bibitem formats, and reading only one is how this class
# stayed invisible: "GeoVac Paper~7 (2026)." (house style) and the older
# "Paper~7 (2026)." without the "GeoVac" (Paper 32's early entries).
# Restricted to the bibliography so an in-text "Paper 7 (2026)" cannot
# false-positive.
BIBYEAR = re.compile(r"(?:GeoVac )?Papers?~?\s?(\d{1,2})\s*\((20\d\d)\)")

# FULL run #7 (2026-09-06) scope fix.  BIBYEAR only sees "Paper~N (YEAR)".
# A bibitem keyed {paperN}/{loutey_paperN} (or labelled "internal Paper~N")
# whose year lives in "GeoVac internal preprint (YEAR)" -- with no "Paper~N"
# adjacent to the year-paren -- escaped it entirely: Paper 38 cited Paper 7 as
# a 2024 preprint against its 2026 date and C11 passed clean (the completeness
# critic caught it, not the gate).  This block scanner reads each bibitem, takes
# the paper number from the KEY or the label, and any "(YEAR)" from the body.
BIBITEM_BLOCK = re.compile(
    r"\\bibitem(?:\[[^\]]*\])?\{([^}]+)\}(.*?)"
    r"(?=\\bibitem|\\end\{thebibliography\})", re.S)
_INTERNAL_KEY = re.compile(r"(?:loutey_)?paper[_]?(\d{1,2})\b", re.I)
_INTERNAL_LABEL = re.compile(r"internal Paper~?\s?(\d{1,2})", re.I)
_ANYYEAR = re.compile(r"\((20\d\d)\)")


def build_year_map() -> dict:
    """{paper_number: year} from each ACTIVE paper's own \\date{}."""
    ymap = {}
    for f in active_tex():
        num = re.search(r"[Pp]aper[_ ](\d+)[_ ]", f.name)
        if not num:
            continue
        txt = f.read_text(encoding="utf-8", errors="replace")
        d = re.search(r"\\date\{([^}]*)\}", txt)
        if not d:
            continue
        y = re.search(r"(20\d\d)", d.group(1))
        if y:
            ymap.setdefault(int(num.group(1)), y.group(1))
    return ymap


YEAR_BASELINE = pathlib.Path(__file__).with_name("internal_year_baseline.json")


def load_year_baseline() -> set:
    """Pre-existing year debt that PASSES; anything new FAILS.

    Recorded 2026-09-01 when this criterion was added: 30 bibitems citing a
    GeoVac paper with a year its own date contradicts, across five groups.
    Papers are corrected in place, so the citing text kept the original
    year.  Mass-editing them out of an unrelated sprint would be the
    out-of-scope edit this apparatus exists to catch;  each group's
    re-certification fixes and deletes its own entries.
    """
    if not YEAR_BASELINE.exists():
        return set()
    import json
    data = json.loads(YEAR_BASELINE.read_text(encoding="utf-8"))
    return {(e["path"], e["paper"], e["cited"])
            for e in data.get("entries", [])}


def scan_years(ymap: dict) -> list:
    """Internal bibitem years that contradict the cited paper's own date."""
    out = []
    for f in active_tex():
        txt = f.read_text(encoding="utf-8", errors="replace")
        i = txt.rfind(r"\begin{thebibliography}")
        if i < 0:
            continue
        self_num = re.search(r"[Pp]aper[_ ](\d+)[_ ]", f.name)
        self_num = int(self_num.group(1)) if self_num else None
        path = str(f.relative_to(ROOT)).replace(chr(92), '/')
        for m in BIBYEAR.finditer(txt, i):
            n, y = int(m.group(1)), m.group(2)
            if n == self_num or n not in ymap or ymap[n] == y:
                continue
            out.append({"file": f.name, "path": path,
                        "paper": n, "cited": y, "real": ymap[n]})
        # Second pass: the "internal preprint (YEAR)" dodge (see BIBITEM_BLOCK).
        for bm in BIBITEM_BLOCK.finditer(txt, i):
            key, body = bm.group(1), bm.group(2)
            km = _INTERNAL_KEY.search(key) or _INTERNAL_LABEL.search(body)
            if not km:
                continue
            n = int(km.group(1))
            if n == self_num or n not in ymap:
                continue
            years = _ANYYEAR.findall(body)
            if not years or ymap[n] in years:
                continue
            if any(o["path"] == path and o["paper"] == n for o in out):
                continue  # already reported by BIBYEAR
            out.append({"file": f.name, "path": path,
                        "paper": n, "cited": years[0], "real": ymap[n]})
    return out


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
# KEYED: bibitems whose paper number lives in the \bibitem KEY rather than in a
# trailing "GeoVac Paper N (YEAR)" marker -- e.g. P22's
#   \bibitem{GeoVac_Paper7} J.~Loutey, ``Title,'' GeoVac Technical Report (2025).
# Resolve N from the key (GeoVac_Paper7 / paper7 / Paper_7), then take the first
# quoted title before the next \bibitem. Closes the run-#6 C11 blind spot where
# the "Technical Report" format left an entire paper's internal cites uncertified.
KEYED = re.compile(
    r"\\bibitem\{(?:geovac[_-]?)?paper[_-]?(?P<n>\d+)[a-z_]*\}"
    r"(?P<body>(?:(?!\\bibitem)[\s\S])*?)" + QUOTE,
    re.DOTALL | re.IGNORECASE,
)


def check_one(fname, relpath, n, cited_title, tmap, mismatches, flagged):
    if n in ARCHIVED or n not in tmap:
        return
    full_norm, main_norm, raw, _ = tmap[n]
    cn = norm(cited_title)
    if cn == full_norm or cn == main_norm:
        return
    rec = {"file": fname, "path": relpath, "paper": n,
           "cited": " ".join(cited_title.split())[:80], "real": raw[:80]}
    (flagged if n in PROPINQUITY else mismatches).append(rec)


eol_hdr_yaudit = ("\n--- AUDIT (advisory, out-of-gate-scope) (%d) -- "
                  "bibitem YEAR vs the cited paper's own date: ---")
eol_hdr_ygated = ("\n*** YEAR MISMATCH (%d) -- an internal bibitem year "
                  "contradicts the cited paper's own date: ***")
yline = "  [P%s] %s: cited (%s) != that paper's date (%s)"
ybase_msg = ("bibitem-year ratchet: %d pre-existing mismatch(es) forgiven from a %d-entry baseline (debug/qa/internal_year_baseline.json) -- shrink it, never grow it")
ymsg = "\nRESULT: FAIL (internal-bibitem-year criterion%s)"


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    gate = None
    if "--gate" in argv:
        i = argv.index("--gate")
        gate = argv[i + 1] if i + 1 < len(argv) else None

    tmap = build_title_map()
    _year_all = scan_years(build_year_map())
    _year_base = load_year_baseline()
    year_known = [r for r in _year_all
                  if (r["path"], r["paper"], r["cited"]) in _year_base]
    year_bad = [r for r in _year_all
                if (r["path"], r["paper"], r["cited"]) not in _year_base]
    mismatches, flagged = [], []
    for f in active_tex():
        if f.name in HOUSE_STYLE:
            continue
        # Repo-relative: qa_scopes.make_predicate matches endswith()
        # against repo-relative paths, so a papers/-relative key here
        # makes the gate predicate False for EVERY finding and the gate
        # can never fail in scope (proven by probe, 2026-09-01).
        relpath = str(f.relative_to(ROOT)).replace("\\", "/")
        txt = f.read_text(encoding="utf-8", errors="replace")
        seen = set()  # (n, normalized cited) -> dedupe a bibitem matched by >1 pattern

        def _check(n: int, title: str, _rel=relpath, _fn=f.name) -> None:
            if title.strip().startswith("GeoVac Paper"):
                return
            key = (n, norm(title))
            if key in seen:
                return
            seen.add(key)
            check_one(_fn, _rel, n, title, tmap, mismatches, flagged)

        for m in ZENODO.finditer(txt):
            _check(int(m.group("n")), m.group("title"))
        for m in NORMAL.finditer(txt):
            _check(int(m.group("n")), m.group("title"))
        for m in KEYED.finditer(txt):
            _check(int(m.group("n")), m.group("title"))

    # Partition by --gate: a mismatch FAILs only if its CITING file's path matches
    # the gate token (e.g. 'group3'); out-of-gate mismatches are advisory AUDIT
    # (mirrors check_k_label.py / check_paper_test_refs.py). No --gate => corpus-wide.
    if gate:
        _in_scope, _scope_files, _scope_warnings = qa_scopes.make_predicate(gate)
        qa_scopes.emit_warnings(_scope_warnings)
        gated = [r for r in mismatches if _in_scope(r["path"])]
        audit = [r for r in mismatches if not _in_scope(r["path"])]
        ygated = [r for r in year_bad if _in_scope(r["path"])]
        yaudit = [r for r in year_bad if not _in_scope(r["path"])]
        scope_desc = qa_scopes.describe(gate, _scope_files)
    else:
        gated, audit = mismatches, []
        ygated, yaudit = year_bad, []
        scope_desc = "ALL"

    print(f"title map: {len([k for k in tmap if isinstance(k,int)])} numbered papers + "
          f"{[k for k in tmap if not isinstance(k,int)]}   [gated scope: {scope_desc}]")
    if flagged:
        print(f"\nFLAGGED (propinquity cluster, descope-pending -- not a failure): {len(flagged)}")
        for r in flagged:
            print(f"  [P{r['paper']}] {r['file']}: cited \"{r['cited']}\"")
    if audit:
        print(f"\n--- AUDIT (advisory, out-of-gate-scope) ({len(audit)}) -- stale titles in other branches: ---")
        for r in audit:
            print(f"  [P{r['paper']}] {r['path']}: cited \"{r['cited']}\"  != \"{r['real']}\"")
    print(ybase_msg % (len(year_known), len(_year_base)))
    if yaudit:
        print(eol_hdr_yaudit % len(yaudit))
        for r in yaudit:
            print(yline % (r['paper'], r['path'], r['cited'], r['real']))
    if ygated:
        print(eol_hdr_ygated % len(ygated))
        for r in ygated:
            print(yline % (r['paper'], r['file'], r['cited'], r['real']))
    if gated:
        print(f"\n*** MISMATCH ({len(gated)}) -- internal cite does not match the paper's \\title: ***")
        for r in gated:
            print(f"  [P{r['paper']}] {r['file']}")
            print(f"      cited: \"{r['cited']}\"")
            print(f"      real : \"{r['real']}\"")
        print(f"\nRESULT: FAIL (internal-title criterion{f', gated scope {gate}' if gate else ''})")
        return 1
    if ygated:
        print(ymsg % ((', gated scope ' + gate) if gate else ''))
        return 1
    print(f"\nRESULT: PASS -- every internal GeoVac citation in {scope_desc} "
          f"matches the cited paper's \\title, and no bibitem year "
          f"outside the recorded baseline contradicts that paper's own \\date.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
