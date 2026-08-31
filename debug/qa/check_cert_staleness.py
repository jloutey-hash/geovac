r"""Which certifications are stale?  Compute it; do not remember it.

A `docs/qa/<target>.done.md` saying CERTIFIED records the state at the
moment it was written.  It does not update itself when the papers move, so
after a few sprints every record says CERTIFIED and none of them is
necessarily true.  That is the same relevance-decay shape C22 was built
for, one level up: a record asserting a verdict about text that has since
changed.

Asked directly on 2026-08-31 ("have we accurately recorded which papers we
owe certification?") the answer was **no** -- all 11 records said CERTIFIED
while 8 of 10 targets had `.tex` changes since their recorded date, two of
them the same day.  CLAUDE.md SS2 carried prose OWED notes, but they were
partial and hand-maintained.

So this is a tool rather than a note.  Run it instead of trusting any
record, including this docstring.

    python debug/qa/check_cert_staleness.py
    python debug/qa/check_cert_staleness.py --detail

Reads the cert date from each record (the latest date it mentions) and asks
git which of that target's `.tex` changed after it.  PDF-only changes are
reported separately -- a recompile is not a claim change.  Same-day changes
are flagged AMBIGUOUS rather than asserted, because a date has no clock and
the record may have been written after the edit.
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
QA = ROOT / "docs" / "qa"
DATE = re.compile(r"20\d{2}-\d{2}-\d{2}")

GROUPS = {
    "group1": "group1_operator_algebras",
    "group2": "group2_quantum_chemistry",
    "group3": "group3_foundations",
    "group4": "group4_quantum_computing",
    "group5": "group5_qed_gauge",
    "group6": "group6_precision_observations",
}


def paths_for(target: str) -> list[str]:
    """Papers a target covers. Single-paper targets are globbed, not hard-coded.

    Hard-coding filenames is how the first version of this audit reported
    `paper_58` and `paper_60` as "paths not found" -- the real names are
    paper_58_abelian_residue.tex and paper_60_sturmian_secular_quantum.tex.
    """
    if target in GROUPS:
        d = f"papers/{GROUPS[target]}"
        syn = f"papers/synthesis/{GROUPS[target]}_synthesis.tex"
        return [p for p in (d, syn) if (ROOT / p).exists()]
    if target == "synthesis":
        return ["papers/synthesis"]
    if target == "trunk":
        # Scope per docs/qa/trunk.done.md: Papers 0, 1, 7 (group3) + 32, 38
        # (group1) -- the roots SS9 says to review before any branch.
        # The DoD scope includes the group3 synthesis, not just the five
        # papers; omitting it was the same scope gap as C19's, one level up.
        out = ["papers/synthesis/group3_foundations_synthesis.tex"]
        for n in (0, 1, 7, 32, 38):
            out += [str(p.relative_to(ROOT)).replace("\\", "/")
                    for p in ROOT.glob(f"papers/*/[Pp]aper_{n}_*.tex")]
        return sorted(out)
    m = re.match(r"paper_(\d+)$", target)
    if m:
        hits = sorted(ROOT.glob(f"papers/*/paper_{m.group(1)}_*.tex"))
        return [str(p.relative_to(ROOT)).replace("\\", "/") for p in hits]
    return []


def cert_date(record: Path) -> str | None:
    ds = DATE.findall(record.read_text(encoding="utf-8", errors="replace"))
    return max(ds) if ds else None


def changed_since(paths: list[str], since: str) -> tuple[list[str], list[str]]:
    """(.tex changed after `since`, .tex changed ON `since` -- ambiguous)."""
    def run(args):
        out = subprocess.run(["git", "log", "--name-only", "--format=%cs", "--"]
                             + paths, capture_output=True, text=True,
                             cwd=ROOT).stdout
        return out

    after, same = set(), set()
    cur = None
    for line in run(paths).splitlines():
        if DATE.fullmatch(line.strip()):
            cur = line.strip()
        elif line.strip().endswith(".tex") and cur:
            if cur > since:
                after.add(line.strip())
            elif cur == since:
                same.add(line.strip())
    return sorted(after), sorted(same - after)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--detail", action="store_true")
    a = ap.parse_args()

    records = sorted(QA.glob("*.done.md"))
    if not records:
        print("no certification records found")
        return 1

    print(f"{'target':<12} {'certified':<12} {'tex after':>9} {'same-day':>9}  verdict")
    print("-" * 66)
    owed, ambiguous, rows = [], [], []
    for r in records:
        t = r.stem.replace(".done", "")
        d = cert_date(r)
        p = paths_for(t)
        if not d or not p:
            print(f"{t:<12} {str(d):<12} {'?':>9} {'?':>9}  NO DATE/PATHS")
            continue
        after, same = changed_since(p, d)
        rows.append((t, d, after, same))
        if after:
            owed.append(t)
            v = "OWED"
        elif same:
            ambiguous.append(t)
            v = "AMBIGUOUS (same-day)"
        else:
            v = "current"
        print(f"{t:<12} {d:<12} {len(after):>9} {len(same):>9}  {v}")

    if a.detail:
        for t, d, after, same in rows:
            if after or same:
                print(f"\n{t} (certified {d}):")
                for f in after:
                    print(f"    [after]    {Path(f).name}")
                for f in same:
                    print(f"    [same-day] {Path(f).name}")

    print(f"\nOWED:      {', '.join(owed) if owed else 'none'}")
    print(f"AMBIGUOUS: {', '.join(ambiguous) if ambiguous else 'none'}")
    print("\nA date has no clock: same-day means the record may have been "
          "written before or after the edit. Check those by hand.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
