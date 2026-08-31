r"""C22 -- test -> claim backing integrity.  C13 run backwards.

C13 asserts: every test a paper cites exists.
C22 asserts: every claim a test backs still exists, and is not retracted.

The reverse direction is where relevance decays.  A test pinned to a
withdrawn claim keeps passing, keeps costing wall time, and nothing can
tell it is dead -- the 2026-08-30 pass hit exactly that, a
`NOTE (FLAG, do not "fix" here)` in tests/test_paper22_density.py pointing
at a test name that no longer existed and calling a dissolved
carry-forward open.

Four checks:

  A. claim_test_matrix.md rows cite tests that EXIST.  Nothing covered
     this -- C13 gates papers, not the matrix.  Currently clean: the two
     rows first reported here were an unbounded regex matching
     `test_lih.py` inside `retest_lih.py`, and a row NAMING a stale cite
     in order to disown it.
  B. test_paperNN_*.py names a paper that exists.  Both filename
     conventions are honoured (`paper_14_*.tex` and `Paper_7_*.tex`);
     checking only the lower-case one would fire falsely on Papers 0, 7
     and 8, which is the fires-on-the-wrong-thing failure.
  C. No matrix row backs a claim whose wording is in C16's retracted
     registry without a withdrawal flag.  Reuses C16 rather than
     duplicating its patterns, so a claim retired in one place cannot stay
     live in the other.
  D. Paper-backing tests do not depend on debug/, which SS9 defines as
     the transient clean-room directory that is pruned over time.  Three
     paper-backing tests currently do (P26, P27 x2) and are baselined --
     fixing them means relocating the imported modules into geovac/ or
     tests/, a real refactor, so the ratchet holds the line meanwhile.

A and D are RATCHETS against `test_claim_backing_baseline.json`, following
C20: they fail on anything outside the recorded baseline, so today's known
debt does not block the gate while new debt cannot enter.  The baseline
count is printed on every run -- a ratchet that hides its own size is how
debt becomes permanent.

    python debug/qa/check_test_claim_backing.py
    python debug/qa/check_test_claim_backing.py --selftest
    python debug/qa/check_test_claim_backing.py --update-baseline
"""
from __future__ import annotations

import argparse
import ast
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "debug" / "qa"))

import check_retracted_terms as C16          # noqa: E402

TESTS = ROOT / "tests"
MATRIX = ROOT / "docs" / "claim_test_matrix.md"
PAPERS = ROOT / "papers"
BASELINE = Path(__file__).with_name("test_claim_backing_baseline.json")

# Left boundary is load-bearing: without it this matches the substring
# `test_lih.py` inside `chemistry_solver_retest_lih.py` and reports a
# dangling row that does not exist.
TESTFILE = re.compile(r"(?<![A-Za-z0-9_])test_[a-z0-9_]+\.py")

# A row may NAME a test in order to disown it ("Replaces stale ... cite").
# Such a name is a historical note, not a backing claim, and must not be
# read as one.  The marker is always on the same line as the name, so no
# context window is involved.
DISOWNED = re.compile(r"stale|retired|replaced|superseded|archived|no longer|withdrawn", re.I)
PAPER_NAMED = re.compile(r"^test_paper(\d+)[_.]")
DEBUG_IMPORT = re.compile(r"(?:from|import)\s+(debug\.[A-Za-z0-9_.]+)")


def existing_tests() -> set[str]:
    return {p.name for p in TESTS.glob("test_*.py")}


def known_papers() -> set[str]:
    """Paper numbers with a .tex, under BOTH naming conventions.

    `paper_14_qubit_encoding.tex` and `Paper_7_Dimensionless_Vacuum.tex`
    both exist; matching only the first reports Papers 0, 7 and 8 as
    missing.  papers/archive/ counts -- an archived paper still exists.
    """
    nums: set[str] = set()
    for tex in PAPERS.rglob("*.tex"):
        m = re.match(r"[Pp]aper_(\d+)_", tex.name)
        if m:
            nums.add(str(int(m.group(1))))
    return nums


def matrix_rows() -> list[tuple[int, str, set[str]]]:
    """(line number, row text, test files cited) for every matrix row."""
    if not MATRIX.exists():
        return []
    out = []
    for i, line in enumerate(
            MATRIX.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
        names = set(TESTFILE.findall(line))
        if names and DISOWNED.search(line):
            # The row disowns the name it mentions; treat it as prose.
            names = {x for x in names if (ROOT / "tests" / x).exists()}
        if names:
            out.append((i, line, names))
    return out


def debug_imports(path: Path) -> set[str]:
    try:
        src = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return set()
    return set(DEBUG_IMPORT.findall(src))


def load_baseline() -> dict:
    if BASELINE.exists():
        return json.loads(BASELINE.read_text(encoding="utf-8"))
    return {"dangling_matrix_rows": [], "paper_backing_debug_imports": []}


# --------------------------------------------------------------- the checks

def check_a(base, verbose=True):
    """Matrix rows cite tests that exist (ratcheted)."""
    live = existing_tests()
    allowed = set(base.get("dangling_matrix_rows", []))
    new, known = [], []
    for ln, _text, names in matrix_rows():
        for n in sorted(names - live):
            (known if n in allowed else new).append((ln, n))
    if verbose:
        print("\nA. claim-matrix rows cite live tests")
        for ln, n in new:
            print(f"   [NEW] {MATRIX.name}:{ln}  cites {n} -- does not exist")
        if known:
            print(f"   [baseline] {len(known)} known dangling row(s): "
                  f"{', '.join(sorted({n for _, n in known}))}")
        if not new:
            print("   [ok] no new dangling row")
    return len(new)


def check_b(verbose=True):
    """test_paperNN_* names a paper that exists."""
    papers = known_papers()
    bad = []
    for p in sorted(TESTS.glob("test_paper*.py")):
        m = PAPER_NAMED.match(p.name)
        if m and m.group(1).lstrip("0") not in papers and m.group(1) not in papers:
            bad.append((p.name, m.group(1)))
    if verbose:
        print("\nB. paper-backing test filenames name a real paper")
        for name, num in bad:
            print(f"   [FAIL] {name} claims paper {num}, which has no .tex")
        if not bad:
            print(f"   [ok] every test_paperNN_* names one of "
                  f"{len(papers)} known papers")
    return len(bad)


def check_c(verbose=True):
    """No matrix row backs a retracted claim (patterns reused from C16)."""
    entries = [e for e in C16.REGISTRY if e.get("severity") == "fail"]
    hits = []
    for ln, text, names in matrix_rows():
        for e in entries:
            if not re.search(e["pattern"], text):
                continue
            ex = e.get("exempt_if_nearby")
            if ex and re.search(ex, text, re.I):
                continue
            hits.append((ln, e["id"], sorted(names)))
    if verbose:
        print("\nC. no matrix row backs a retracted claim")
        for ln, eid, names in hits:
            print(f"   [FAIL] {MATRIX.name}:{ln} matches retired '{eid}'")
            print(f"          backing test(s): {', '.join(names)}")
        if not hits:
            print(f"   [ok] {len(entries)} retracted patterns checked, none live")
    return len(hits)


def check_d(base, verbose=True):
    """Paper-backing tests do not depend on the prunable debug/ tree."""
    allowed = set(base.get("paper_backing_debug_imports", []))
    new, known = [], []
    for p in sorted(TESTS.glob("test_paper*.py")):
        mods = debug_imports(p)
        if mods:
            (known if p.name in allowed else new).append((p.name, sorted(mods)))
    if verbose:
        print("\nD. paper-backing tests avoid the prunable debug/ tree")
        for name, mods in new:
            print(f"   [NEW] {name} imports {', '.join(mods)}")
            print("         debug/ is pruned over time (CLAUDE.md SS9); this "
                  "breaks paper backing silently")
        if known:
            print(f"   [baseline] {len(known)} known: "
                  f"{', '.join(n for n, _ in known)}")
        if not new:
            print("   [ok] no new debug/ dependency in a paper-backing test")
    return len(new)


# ------------------------------------------------------------------ selftest

def selftest() -> int:
    """Two-way discrimination: each check must FIRE and must go SILENT.

    Mandatory.  A gate proven only to stay quiet is indistinguishable from
    one that cannot fire at all, which is the failure this whole QA arc
    keeps finding.
    """
    ok = True

    # A fires on a SYNTHETIC dangling row, and is silent once baselined.
    # Deliberately synthetic: a selftest that proves a check can fire by
    # pointing at real debt stops working the moment the debt is fixed,
    # which is precisely backwards.
    probe_name = "test_c22_selftest_absent_probe.py"
    probe_row = f"| 14 | synthetic claim | `{probe_name}` | mod | BACKED |"
    seen = set(TESTFILE.findall(probe_row))
    fires = int(bool(seen - existing_tests()))
    silent = int(bool(seen - existing_tests() - {probe_name}))
    # and the real corpus must be clean under the shipped baseline
    live_a = check_a(load_baseline(), verbose=False)
    print(f"  A  fires={fires} (want >0)   baselined={silent} (want 0)   "
          f"live corpus={live_a} (want 0)   "
          f"{'PASS' if fires and not silent and not live_a else 'FAIL'}")
    ok &= bool(fires) and not silent and not live_a

    # B fires on a fabricated paper number.
    probe = TESTS / "test_paper9999_selftest_probe.py"
    probe.write_text("# C22 selftest probe\n", encoding="utf-8")
    try:
        fires_b = check_b(verbose=False)
    finally:
        probe.unlink()
    silent_b = check_b(verbose=False)
    print(f"  B  fires={fires_b} (want >0)   after removal={silent_b} (want 0)"
          f"   {'PASS' if fires_b and not silent_b else 'FAIL'}")
    ok &= bool(fires_b) and not silent_b

    # C fires on a synthetic row carrying a retired phrase.
    entries = [e for e in C16.REGISTRY if e.get("severity") == "fail"]
    probe_row = "| 26 | angular ERI density is $2.76\\% | test_probe.py |"
    hit = any(re.search(e["pattern"], probe_row) for e in entries)
    print(f"  C  fires on a synthetic retired row={hit} (want True)   "
          f"live rows={check_c(verbose=False)} (want 0)   "
          f"{'PASS' if hit and not check_c(verbose=False) else 'FAIL'}")
    ok &= hit and not check_c(verbose=False)

    # D fires on the real debug/-importing paper tests when un-baselined.
    fires_d = check_d({"paper_backing_debug_imports": []}, verbose=False)
    names = [p.name for p in TESTS.glob("test_paper*.py") if debug_imports(p)]
    silent_d = check_d({"paper_backing_debug_imports": names}, verbose=False)
    print(f"  D  fires={fires_d} (want >0)   baselined={silent_d} (want 0)   "
          f"{'PASS' if fires_d and not silent_d else 'FAIL'}")
    ok &= bool(fires_d) and not silent_d

    print(f"\nSELFTEST: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


def update_baseline() -> int:
    dangling = sorted({n for _, _, names in matrix_rows()
                       for n in names - existing_tests()})
    debug_users = sorted(p.name for p in TESTS.glob("test_paper*.py")
                         if debug_imports(p))
    BASELINE.write_text(json.dumps(
        {"_note": "Known debt at the time C22 was introduced. The gate FAILS "
                  "on anything outside these lists. Shrink them; do not grow "
                  "them.",
         "dangling_matrix_rows": dangling,
         "paper_backing_debug_imports": debug_users}, indent=2) + "\n",
        encoding="utf-8")
    print(f"wrote {BASELINE.name}: {len(dangling)} dangling row(s), "
          f"{len(debug_users)} debug-importing paper test(s)")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--update-baseline", action="store_true")
    a = ap.parse_args()

    if a.selftest:
        return selftest()
    if a.update_baseline:
        return update_baseline()

    base = load_baseline()
    n = check_a(base) + check_b() + check_c() + check_d(base)
    print()
    if n:
        print(f"RESULT: FAIL ({n} new test->claim backing defect(s))")
        return 1
    print("RESULT: PASS (every test-backed claim resolves and is not "
          "retracted; no new debug/ dependency)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
