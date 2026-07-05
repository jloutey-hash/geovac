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

Scope: all .tex in papers/group*/ and papers/synthesis/ (papers/archive/
is historical and exempt).  Exit 0 = PASS, 1 = FAIL.

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

ROOT = Path(__file__).resolve().parents[2]

NUM = r"(?:\d+|one|two|three|four|five|six|seven|eight|nine|ten|several|a few|many)"
UNIT = r"(?:year|month|week|day)"

# Two severity tiers:
#   FAIL     -- historical wall-clock durations about the project's course
#               (the fabricated-elapsed-time class; zero tolerance).
#   ADVISORY -- forward effort-estimate vocabulary (multi-month / multi-year /
#               week-scale ...). Also banned by rule 12, but ~190 pre-existing
#               instances live in certified group papers; counted as STANDING
#               DEBT until the corpus-wide sweep sprint retires them, then
#               promote to FAIL (the C14 debug/-refs precedent).
FAIL_PATTERNS: List[Tuple[str, str]] = [
    ("duration-ago",
     rf"\b{NUM}\s+{UNIT}s?\s+ago\b"),
    ("units-of-project-work",
     rf"\b{UNIT}s?\s+of\s+(?:work|effort|development|iteration|iterations|"
     rf"research|sprints|refinement|investigation)\b"),
    ("over-the-past-unit",
     rf"\bover\s+the\s+(?:past|last)\s+(?:{NUM}\s+)?{UNIT}s?\b"),
    ("across-n-units",
     rf"\b(?:across|within|in)\s+{NUM}\s+{UNIT}s?\b"),
    ("compressed-into-units",
     rf"\binto\s+(?:days|weeks|months)\b"),
    ("project-took-units",
     rf"\btook\s+(?:{NUM}\s+)?{UNIT}s?\b"),
]

ADVISORY_PATTERNS: List[Tuple[str, str]] = [
    ("multi-unit-estimate",
     rf"\bmulti[- ]{UNIT}\b"),
    ("unit-scale-estimate",
     rf"\b{UNIT}s?-scale\b"),
    ("unit-long",
     rf"\b{NUM}?[- ]?{UNIT}s?-long\b"),
]

# Per-occurrence exemptions: (filename-substring, regex) pairs whose matches
# are allowed (add sparingly; each entry needs a comment saying why).
ALLOWLIST: List[Tuple[str, str]] = [
    # (none currently -- physics "decades" (log) and external-world history
    #  do not trip the narrow patterns above)
]

GATED_DIRS = [
    "papers/synthesis",
    "papers/group1_operator_algebras",
    "papers/group2_quantum_chemistry",
    "papers/group3_foundations",
    "papers/group4_quantum_computing",
    "papers/group5_qed_gauge",
    "papers/group6_precision_observations",
]

BRANCH_DIRS = {
    "trunk": ["papers/group1_operator_algebras", "papers/group3_foundations"],
    "group1": ["papers/group1_operator_algebras"],
    "group2": ["papers/group2_quantum_chemistry"],
    "group3": ["papers/group3_foundations"],
    "group4": ["papers/group4_quantum_computing"],
    "group5": ["papers/group5_qed_gauge"],
    "group6": ["papers/group6_precision_observations"],
    "synthesis": ["papers/synthesis"],
}


def strip_comments(text: str) -> str:
    """Drop LaTeX comment tails (unescaped %) so commented-out prose can't trip."""
    out = []
    for line in text.splitlines():
        m = re.search(r"(?<!\\)%", line)
        out.append(line[: m.start()] if m else line)
    return "\n".join(out)


def scan_file(path: Path, patterns: List[Tuple[str, str]]) -> List[Tuple[str, int, str, str]]:
    hits = []
    text = strip_comments(path.read_text(encoding="utf-8", errors="replace"))
    for lineno, line in enumerate(text.splitlines(), 1):
        for name, pat in patterns:
            for m in re.finditer(pat, line, re.IGNORECASE):
                allowed = any(
                    sub in str(path).replace("\\", "/") and re.search(apat, line, re.IGNORECASE)
                    for sub, apat in ALLOWLIST
                )
                if not allowed:
                    hits.append((name, lineno, m.group(0), line.strip()[:120]))
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
    ]
    negatives = [
        "closed POSITIVE-THIN in May 2026",
        "corrected 2026-06-18 per the register",
        "an earlier draft reported a non-match",
        "may close at sprint-scale once the machinery is in place",
        "five decades of nonzero rho_M across the sweep",
        "a decade-long experimental anomaly (proton radius, external)",
        "subsequently descoped by the degeneracy theorem",
        "13 CPU-h at 30 dps never returned",
    ]
    ok = True
    for s in positives:
        if not any(re.search(p, s, re.IGNORECASE) for _, p in all_patterns):
            print(f"  [selftest FAIL] should match but did not: {s!r}")
            ok = False
    for s in negatives:
        matched = [n for n, p in all_patterns if re.search(p, s, re.IGNORECASE)]
        # 'a decade-long ...' is external-world; 'decade' is deliberately NOT
        # in UNIT, so it must not match.
        if matched:
            print(f"  [selftest FAIL] should NOT match but did ({matched}): {s!r}")
            ok = False
    print(f"selftest: {'PASS' if ok else 'FAIL'} "
          f"({len(positives)} positives, {len(negatives)} negatives)")
    return 0 if ok else 1


def main(argv: List[str]) -> int:
    if "--selftest" in argv:
        return selftest()

    dirs = GATED_DIRS
    if "--gate" in argv:
        branch = argv[argv.index("--gate") + 1]
        dirs = BRANCH_DIRS.get(branch, GATED_DIRS)

    fail_total = 0
    advisory_by_file: dict = {}
    files = 0
    for d in dirs:
        base = ROOT / d
        if not base.is_dir():
            continue
        for tex in sorted(base.glob("*.tex")):
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

    print(f"\nduration-language check over {files} papers  [dirs: {', '.join(dirs)}]")
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
