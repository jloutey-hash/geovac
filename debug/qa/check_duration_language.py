"""C18 -- Project-course duration-language check (deterministic).

Rule (PI direction, 2026-07-04): paper prose must never attach wall-clock
units (years / months / weeks / days) to the PROJECT'S OWN course or events
-- neither historical ("Three years ago the project ...") nor as forward
effort estimates ("multi-month program", "multi-year frontier").  LLM
drafting is unreliable about elapsed project time (the incident: the field
guide claimed "Three years ago the project produced transcendentals ..."
when the entire project spans well under a year).  Sequence language
("earlier", "subsequently", "an earlier draft"), real dates ("2026-06-10",
"May 2026"), version anchors ("v3.19.0"), and the project's own unit-free
work vocabulary ("sprint-scale", "beyond sprint scale") are all fine --
they are verifiable.  External-world history ("a decade-long experimental
anomaly" about the proton radius puzzle) is exempt: the rule targets the
project's course, and the FAIL patterns below are written narrowly enough
that ordinary physics usage (log "decades" of a ratio, "light-year",
per-year clock-drift specs) does not trip them.

Scope: resolved by debug/qa/qa_scopes.py from the pre-registered DoD
(papers/archive/ is historical and exempt).  Exit 0 = PASS, 1 = FAIL.

Convention home: docs/authoring_conventions.md project-wide rule 12;
severity: fix-on-sight NIT unless it distorts a result's provenance
(docs/qa/criteria.md).  Registered as C18 in qa.md step 1.

Usage:
    python debug/qa/check_duration_language.py [--gate <branch>] [--selftest]
"""

from __future__ import annotations

import re
import sys
from pathlib import Path
from typing import List, Tuple
import os

# Shared --gate scope resolution (see debug/qa/qa_scopes.py): named
# scopes resolve to an explicit file list and every RESULT line carries
# the file count, so a gate can never report PASS on an empty scope.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qa_scopes  # noqa: E402


ROOT = Path(__file__).resolve().parents[2]

NUM = r"(?:\d+|one|two|three|four|five|six|seven|eight|nine|ten|several|a few|many)"
UNIT = r"(?:year|month|week|day)"

# All duration classes are FAIL tier (2026-07-05 promotion: the corpus-wide
# sweep retired the ~190-instance forward effort-estimate debt, so the former
# ADVISORY tier is promoted per the standing plan; the C14 debug/-refs
# precedent). ADVISORY_PATTERNS kept as an (empty) hook for future staging.
FAIL_PATTERNS: List[Tuple[str, str]] = [
    ("same-day-project-events",
     r"\bsame-day\s+(?:diagnostic\s+)?(?:sprints?|arcs?|tracks?|passes?)\b"),
    ("duration-ago",
     rf"\b{NUM}\s+{UNIT}s?\s+ago\b"),
    # allow up to two intervening adjectives ("weeks of FOCUSED work",
    # "months of INTENSIVE effort") -- the bare "of work" form missed the
    # 2026-08-24 group4 completeness-critic find ("one to several weeks of
    # focused work" in Paper 23)
    ("units-of-project-work",
     rf"\b{UNIT}s?\s+of\s+(?:\w+\s+){{0,2}}(?:work|effort|development|iteration|iterations|"
     rf"research|sprints|refinement|investigation)\b"),
    ("over-the-past-unit",
     rf"\bover\s+the\s+(?:past|last)\s+(?:{NUM}\s+)?{UNIT}s?\b"),
    ("across-n-units",
     rf"\b(?:across|within|in)\s+{NUM}\s+{UNIT}s?\b"),
    ("compressed-into-units",
     rf"\binto\s+(?:days|weeks|months)\b"),
    ("project-took-units",
     rf"\btook\s+(?:{NUM}\s+)?{UNIT}s?\b"),
    # -- promoted from ADVISORY 2026-07-05 --
    ("multi-unit-estimate",
     rf"\bmulti[- ]{UNIT}\b"),
    ("unit-scale-estimate",
     rf"\b{UNIT}s?-scale\b"),
    ("unit-long",
     rf"\b{NUM}?[- ]?{UNIT}s?-long\b"),
    # -- added 2026-07-05: numeric/adjectival forms the sweep found by
    #    judgment ("1-week sprint", "$\sim 6$--$12$ months", "2--3 days") --
    ("numeric-duration",
     rf"\b\d+\$?(?:\s*--?\s*\$?\d+\$?)?[-\s~]+{UNIT}s?\b"),
    # 4b. symbolic math-mode counts attached to a unit ("$N$-month effort",
    #     "$k$-week"): a variable standing in for a duration is still a
    #     duration (cert-3 owed fix, 2026-08-29)
    ("symbolic-math-duration",
     rf"\$[A-Za-z]\$[-~\s]+{UNIT}s?\b"),
    # -- added 2026-07-05 (sweep-delta finding): word-number adjectivals
    #    ("a focused two-day session") escape the digit-based pattern --
    ("wordnum-duration",
     rf"\b(?:one|two|three|four|five|six|seven|eight|nine|ten)-{UNIT}\b"),
]

ADVISORY_PATTERNS: List[Tuple[str, str]] = []

# Lines matching this are EXEMPT from all duration patterns: wall-clock
# figures describing COMPUTATIONAL RUNTIME are legitimate technical
# measurements, not project-course chronology (rule 12 targets the latter).
# Added 2026-08-29 when the LaTeX-tie extension started catching
# "$\sim 3.3$~days per PES point" (Paper 17).
EXEMPT_CONTEXT = (r"per\s+(PES\s+)?point|per\s+iteration|runtime|wall[- ]?clock"
                  r"|wall\s+time|CPU|compute\s+time|solver|diagonaliz")

# Per-occurrence exemptions: (filename-substring, regex) pairs whose matches
# are allowed (add sparingly; each entry needs a comment saying why).
ALLOWLIST: List[Tuple[str, str]] = [
    # (none currently -- physics "decades" (log) and external-world history
    #  do not trip the narrow patterns above)
]

# Scope tables RETIRED 2026-08-31.  GATED_DIRS/BRANCH_DIRS were
# directory-based, which made `--gate trunk` scan all 27 papers of
# group1+group3 instead of trunk's 6, dropped each group's synthesis
# (a different directory), and silently widened any single-paper target
# to the whole corpus.  The scope declarations now live once in
# debug/qa/qa_scopes.py, keyed to docs/qa/<target>.done.md.


def strip_comments(text: str) -> str:
    """Drop LaTeX comment tails (unescaped %) so commented-out prose can't trip."""
    out = []
    for line in text.splitlines():
        m = re.search(r"(?<!\\)%", line)
        out.append(line[: m.start()] if m else line)
    return "\n".join(out)


def scan_file(path: Path, patterns: List[Tuple[str, str]]) -> List[Tuple[str, int, str, str]]:
    """Scan with newlines treated as spaces (LaTeX semantics) so phrases
    spanning a source line break are still caught (the 'several years\\nof
    investigation' miss), mapping each hit back to its source line."""
    hits = []
    text = strip_comments(path.read_text(encoding="utf-8", errors="replace"))
    lines = text.splitlines()
    starts: List[int] = []          # char offset of each line in the joined text
    pos = 0
    for line in lines:
        starts.append(pos)
        pos += len(line) + 1        # +1 for the space that replaces the newline
    joined = " ".join(lines)
    from bisect import bisect_right
    for name, pat in patterns:
        for m in re.finditer(pat, joined, re.IGNORECASE):
            lineno = bisect_right(starts, m.start())          # 1-indexed line
            context = lines[lineno - 1].strip()[:120]
            # exemption window spans the hit line and its neighbour (LaTeX
            # wraps mid-phrase: '...$~days per' / 'PES point...')
            window = lines[lineno - 1]
            if lineno < len(lines):
                window += " " + lines[lineno]
            if lineno >= 2:
                window = lines[lineno - 2] + " " + window
            if re.search(EXEMPT_CONTEXT, window, re.IGNORECASE):
                continue
            allowed = any(
                sub in str(path).replace("\\", "/") and re.search(apat, context, re.IGNORECASE)
                for sub, apat in ALLOWLIST
            )
            if not allowed:
                hits.append((name, lineno, m.group(0), context))
    return hits


def selftest() -> int:
    all_patterns = FAIL_PATTERNS + ADVISORY_PATTERNS
    positives = [
        "Three years ago the project produced transcendentals",
        "after months of iteration on the solver",
        "over the past year the arc closed",
        "is the named multi-year frontier",
        "an explicitly open, multi-month program",
        "remains open; multi-month scale.",
        "has compressed multi-year content into days",
        "the derivation took three months",
        "a week-long diagnostic arc",
        "19 sub-sprints across two days",
        # phrase spanning a LaTeX source line break (join semantics)
        "Over several years\nof investigation, it acquired",
        # numeric/adjectival forms (post-promotion additions)
        "completed in a focused two-day session",
        "a 1-week sprint at highest priority",
        "commitment ($\\sim 6$--$12$ months estimated)",
        "spanning $12$~months of effort",       # LaTeX tie escape (2026-08-29)
        "an $N$-month effort",                  # symbolic math count (2026-08-29)
        "NotImplementedError ($\\sim 2$-$3$ days, bundled)",
        "multi-week+ architectural lifts",
        "Three same-day diagnostic sprints tested the wall",
    ]
    negatives = [
        "Three diagnostic sprints in immediate succession tested the wall",
        "closed POSITIVE-THIN in May 2026",
        "corrected 2026-06-18 per the register",
        "an earlier draft reported a non-match",
        "may close at sprint-scale once the machinery is in place",
        "five decades of nonzero rho_M across the sweep",
        "a decade-long experimental anomaly (proton radius, external)",
        "subsequently descoped by the degeneracy theorem",
        "13 CPU-h at 30 dps never returned",
        "J. Math. Phys. 59 (2018), 062303",
        "Fields Inst. Commun. 12",
        "the long-range Stage-2 program at sprint scale",
    ]
    ok = True
    for s in positives:
        s_joined = s.replace("\n", " ")     # the scan_file join semantics
        if not any(re.search(p, s_joined, re.IGNORECASE) for _, p in all_patterns):
            print(f"  [selftest FAIL] should match but did not: {s!r}")
            ok = False
    for s in negatives:
        matched = [n for n, p in all_patterns if re.search(p, s, re.IGNORECASE)]
        # 'a decade-long ...' is external-world; 'decade' is deliberately NOT
        # in UNIT, so it must not match.
        if matched:
            print(f"  [selftest FAIL] should NOT match but did ({matched}): {s!r}")
            ok = False
    # exemption leg (2026-08-29): computational-runtime durations match the
    # raw patterns but must be EXEMPT at scan level via EXEMPT_CONTEXT.
    exempt_cases = [
        r"the resulting ${\sim}3.3$~days per PES point",
        "a wall-clock time of 2 days for the full sweep",
        "solver needs 3 hours per iteration",
    ]
    nonexempt_cases = [
        "the arc took $6$--$12$~months of work",
    ]
    for c in exempt_cases:
        if not re.search(EXEMPT_CONTEXT, c, re.IGNORECASE):
            ok = False
            print(f"  [selftest FAIL] runtime case NOT exempt: {c!r}")
    for c in nonexempt_cases:
        if re.search(EXEMPT_CONTEXT, c, re.IGNORECASE):
            ok = False
            print(f"  [selftest FAIL] project-course case wrongly exempt: {c!r}")

    print(f"selftest: {'PASS' if ok else 'FAIL'} "
          f"({len(positives)} positives, {len(negatives)} negatives)")
    return 0 if ok else 1


def main(argv: List[str]) -> int:
    if "--selftest" in argv:
        return selftest()

    # Scope resolution moved to the shared resolver (debug/qa/qa_scopes.py).
    # The old BRANCH_DIRS was DIRECTORY-based, so `--gate trunk` scanned all 27
    # papers of group1+group3 rather than trunk's 6, group runs silently
    # omitted their synthesis (a different directory), and any single-paper
    # target fell through to the whole corpus while still printing its name.
    gate = ""
    if "--gate" in argv:
        gate = argv[argv.index("--gate") + 1]
    tex_files, _scope_warnings = qa_scopes.resolve(gate)
    qa_scopes.emit_warnings(_scope_warnings)

    fail_total = 0
    advisory_by_file: dict = {}
    files = 0
    for tex_path in tex_files:
        tex = Path(tex_path)
        files += 1
        rel = tex.relative_to(ROOT)
        for name, lineno, frag, line in scan_file(tex, FAIL_PATTERNS):
            fail_total += 1
            print(f"  [FAIL:{name}] {rel}:{lineno}  '{frag}'  in: {line}")
        adv = scan_file(tex, ADVISORY_PATTERNS)
        if adv:
            advisory_by_file[str(rel)] = len(adv)

    adv_total = sum(advisory_by_file.values())
    if advisory_by_file:
        print("\n  advisory (forward effort-estimate vocabulary; standing debt "
              "until the corpus sweep):")
        for rel, n in sorted(advisory_by_file.items(), key=lambda kv: -kv[1]):
            print(f"    {n:4d}  {rel}")

    print(f"\nduration-language check over "
          f"{qa_scopes.describe(gate, tex_files)}")
    print(f"  fail-tier hits: {fail_total};  advisory hits: {adv_total}")
    if fail_total:
        print(f"RESULT: FAIL ({fail_total} historical wall-clock duration(s) "
              f"attached to the project's course; fix per "
              f"authoring_conventions.md rule 12)")
        return 1
    print("RESULT: PASS -- no historical wall-clock duration attached to the "
          "project's course in the gated papers (advisory debt reported above).")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
