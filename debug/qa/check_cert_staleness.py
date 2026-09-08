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


# A date is only a CERTIFICATION date if it sits on a line that asserts one.
# Before 2026-09-07 this was `max(every date in the file)`, so a record whose
# newest line said "NOT re-certified" was reported as freshly certified -- the
# audit against stale certifications, fooled by a date, which is the exact
# failure it exists to prevent.
# (A same-line "CERTIFIED" test was tried and removed -- see cert_date.)
# An EXPLICIT decline, written by a run that chose not to certify.  Kept narrow
# on purpose: the generated staleness banner says "historical" and
# "RE-CERTIFICATION OWED" about every stale record, and reading those as
# declines would relabel the whole table.
DISCLAIM = re.compile(
    r"NOT\s+(?:re-)?certified"
    r"|NEVER\s+(?:been\s+)?CERTIFIED", re.I)


def _dated_lines(record: Path):
    text = record.read_text(encoding="utf-8", errors="replace")
    for line in text.splitlines():
        for d in DATE.findall(line):
            yield d, line


ASSERTS_CERT = re.compile(r"CERTIFIED|certifying pass|cert\s*#?\d*\s*=?\s*PASS",
                          re.I)


def cert_date(record: Path) -> str | None:
    """Newest date on a line that ASSERTS certification; else newest date at all.

    Both simpler rules are wrong, in opposite directions, and both were tried
    on 2026-09-07:

      * max-date-anywhere reads ANY later date as the certification date.  It
        overstated group3 by five days -- that record's certifying line says
        2026-08-24 and its banner said 2026-08-29, picked up from
        "Post-certification touch" notes -- and it read a scope note added the
        same day as a fresh certification.
      * same-line-only lost trunk completely, whose status line is
        "FROZEN -- certified PASS" with the date on a different line.

    So: prefer an asserting line, fall back only when the record has none.
    `_used_fallback` records which records relied on the fallback, so that is
    visible rather than silent.
    """
    text = record.read_text(encoding="utf-8", errors="replace")
    good = []
    for line in text.splitlines():
        if ASSERTS_CERT.search(line) and not DISCLAIM.search(line):
            good.extend(DATE.findall(line))
    if good:
        return max(good)
    _used_fallback.add(record.stem.replace(".done", ""))
    ds = DATE.findall(text)
    return max(ds) if ds else None


_used_fallback: set[str] = set()


STATUS_LINE = re.compile(r"STATUS\s*:", re.I)


def status_declines(record: Path) -> "str | None":
    """The record's own STATUS line, if it declines certification.

    This is the record's current-state declaration and outranks any date
    arithmetic.  Keying on dates alone got group6 wrong: its STATUS says
    "NOT CERTIFIED -- superseded 2026-08-22" while an unrelated later line
    dated 2026-08-24 mentions a certifying pass, so a `declined >= cert`
    comparison read it as merely OWED.
    """
    for line in record.read_text(encoding="utf-8", errors="replace").splitlines():
        if STATUS_LINE.search(line):
            return line.strip() if DISCLAIM.search(line) else None
    return None


def declined_date(record: Path) -> "str | None":
    """Newest date on a line that explicitly DECLINES certification.

    Secondary trigger, for a record whose STATUS line has not been updated but
    which carries a dated decline newer than its certification.
    """
    bad = [d for d, line in _dated_lines(record) if DISCLAIM.search(line)]
    return max(bad) if bad else None


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


def _selftest() -> int:
    """Pin the four ways this gate has actually been wrong (2026-09-07)."""
    import tempfile

    cases = [
        # (name, body, expect_cert_date, expect_declines)
        ("max_date_overstates",
         "> **STATUS: CERTIFIED 2026-08-24**\n\n"
         "## Post-certification touch (2026-08-29)\nsome later note\n",
         "2026-08-24", False),
        ("no_asserting_line_falls_back",
         "> **STATUS: FROZEN -- certified PASS** (run #4)\n\n"
         "note dated 2026-09-06\n",
         "2026-09-06", False),
        ("status_declines_wins_over_later_cert_line",
         "> **STATUS: NOT CERTIFIED -- superseded 2026-08-22**\n\n"
         "**FULL certifying pass = PASS, 2026-08-24.**\n",
         "2026-08-24", True),
        ("never_certified_is_a_decline",
         "> **STATUS: NEVER CERTIFIED -- first DoD, 2026-09-07**\n",
         "2026-09-07", True),
    ]

    ok = True
    for name, body, want_date, want_decline in cases:
        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / f"{name}.done.md"
            f.write_text(body, encoding="utf-8")
            _used_fallback.clear()
            got_date = cert_date(f)
            got_decline = bool(status_declines(f))
            if got_date != want_date:
                ok = False
                print(f"[FAIL] {name}: cert_date {got_date!r}, want {want_date!r}")
            elif got_decline != want_decline:
                ok = False
                print(f"[FAIL] {name}: declines {got_decline}, want {want_decline}")
            else:
                print(f"[ok]   {name}")

    # The fallback must be REPORTED, not silent -- a record with no asserting
    # line is exactly where the date can be overstated.
    with tempfile.TemporaryDirectory() as d:
        f = Path(d) / "nofallbackflag.done.md"
        f.write_text("no status line here, just 2026-01-02\n", encoding="utf-8")
        _used_fallback.clear()
        cert_date(f)
        if "nofallbackflag" not in _used_fallback:
            ok = False
            print("[FAIL] fallback used but not recorded in _used_fallback")
        else:
            print("[ok]   fallback is recorded, not silent")

    print("\nSELFTEST:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--detail", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        return _selftest()

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
        declined = declined_date(r)
        if status_declines(r):
            owed.append(t)
            v = (f"NOT CERTIFIED (record says so{f'; {declined}' if declined else ''})")
        elif declined and declined >= d:
            owed.append(t)
            v = f"NOT CERTIFIED (re-run {declined})"
        elif after:
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

    if _used_fallback:
        print(f"\nNOTE: no CERTIFIED-asserting dated line in "
              f"{', '.join(sorted(_used_fallback))} -- date taken as the "
              f"newest in the record, which can overstate it.")
    print(f"\nOWED:      {', '.join(owed) if owed else 'none'}")
    print(f"AMBIGUOUS: {', '.join(ambiguous) if ambiguous else 'none'}")
    print("\nA date has no clock: same-day means the record may have been "
          "written before or after the edit. Check those by hand.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
