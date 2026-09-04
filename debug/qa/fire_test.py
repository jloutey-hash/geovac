"""fire_test.py -- prove a guard can fail, at the moment you write it.

WHY THIS EXISTS
---------------
/qa trunk DELTA #4 and FULL run #4 (2026-09-03) found EIGHT guards in this
corpus that could not fail.  Two were written the day before they were caught.
The pattern behind all eight is the same: a test written right after
establishing a claim reproduces the claim's own reasoning, so it restates the
claim instead of testing it.  Examples, all real:

    assert max(a / 2) == max(a) / 2          # an identity, any input
    assert wrong.shape != direct.shape or ...  # short-circuits, never compares
    assert gap[1] < 1e-9                     # true whenever the line above holds
    h_vee = 2.0; assert 4.0 / h_vee == 2.0   # the literal is the answer

None of these was found by inspection.  Every one was found by BREAKING THE
SUBJECT and seeing the guard stay green -- which is what the reviewers do and
what the author should have done.  This makes that a single command.

USAGE
-----
    python debug/qa/fire_test.py tests/test_foo.py \\
        --plant "expected = 2.0=>expected = 3.0"

    python debug/qa/fire_test.py tests/test_foo.py \\
        --plant-in geovac/bar.py --plant "return 4=>return 5" \\
        -k test_the_specific_guard

Each --plant is "OLD=>NEW".  Multiple --plant flags apply in order.  Exit code
is 0 when every plant made the selected tests FAIL (the guard discriminates)
and 1 when any plant left them green (the guard is asleep).

The repo is never modified: plants are applied to a scratch copy, the original
is restored from a byte copy in a finally block, and the scratch file is
removed.  A plant that does not match its anchor is an error, not a silent
pass -- a no-op plant would otherwise read as "guard did not fire".
"""
from __future__ import annotations

import argparse
import io
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _apply(text: str, plants: list[tuple[str, str]], where: str) -> str:
    for old, new in plants:
        n = text.count(old)
        if n == 0:
            raise SystemExit(
                f"ERROR: plant anchor not found in {where}:\n  {old!r}\n"
                f"A plant that matches nothing looks exactly like a guard that "
                f"did not fire.  Fix the anchor.")
        text = text.replace(old, new, 1)
    return text


def run(test_path: str, plants: list[tuple[str, str]],
        plant_in: str | None = None, selector: str | None = None,
        quiet: bool = False) -> bool:
    """Return True when the planted defect makes the selected tests fail."""
    target = os.path.join(ROOT, plant_in or test_path)
    backup = target + ".firetest.bak"
    cmd = [sys.executable, "-m", "pytest", test_path, "-q", "-x"]
    if selector:
        cmd += ["-k", selector]

    shutil.copy(target, backup)
    try:
        src = io.open(target, encoding="utf-8").read()
        io.open(target, "w", encoding="utf-8", newline="\n").write(
            _apply(src, plants, plant_in or test_path))
        proc = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    finally:
        shutil.move(backup, target)

    fired = proc.returncode != 0
    if not quiet:
        tail = [l for l in proc.stdout.splitlines() if l.strip()][-1:]
        print(f"  plant -> {'FIRED' if fired else 'DID NOT FIRE'}"
              + (f"   ({tail[0].strip()})" if tail else ""))
    return fired


def _selftest() -> int:
    """Prove the helper itself discriminates, on a throwaway test file.

    A tool that reports FIRED for everything is worse than no tool, so it is
    checked in both directions: a real guard must fire, and a tautological one
    must not.
    """
    fd, path = tempfile.mkstemp(suffix="_firetest_selftest.py", dir=os.path.join(ROOT, "tests"))
    os.close(fd)
    rel = os.path.join("tests", os.path.basename(path))
    io.open(path, "w", encoding="utf-8", newline="\n").write(
        "VALUE = 2.0\n\n\n"
        "def test_real_guard():\n"
        "    assert VALUE == 2.0\n\n\n"
        "def test_tautological_guard():\n"
        "    assert VALUE / 2 == VALUE * 0.5\n")
    try:
        print("fire_test self-test")
        real = run(rel, [("VALUE = 2.0", "VALUE = 3.0")],
                   selector="test_real_guard")
        taut = run(rel, [("VALUE = 2.0", "VALUE = 3.0")],
                   selector="test_tautological_guard")
        ok = real and not taut
        print(f"  real guard fired: {real}   (expected True)")
        print(f"  tautological guard fired: {taut}   (expected False)")
        print(f"\nRESULT: {'PASS' if ok else 'FAIL'} "
              f"(the helper separates a guard that tests something from one "
              f"that restates itself)")
        return 0 if ok else 1
    finally:
        os.remove(path)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("test_path", nargs="?", help="test file to run")
    ap.add_argument("--plant", action="append", default=[],
                    metavar="OLD=>NEW", help='mutation, e.g. "x = 2=>x = 3"')
    ap.add_argument("--plant-in", default=None,
                    help="file to mutate (default: the test file itself)")
    ap.add_argument("-k", dest="selector", default=None,
                    help="pytest -k selector")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    if args.selftest:
        return _selftest()
    if not args.test_path or not args.plant:
        ap.error("give a test_path and at least one --plant (or --selftest)")

    plants = []
    for spec in args.plant:
        if "=>" not in spec:
            ap.error(f"--plant must be OLD=>NEW, got {spec!r}")
        old, new = spec.split("=>", 1)
        plants.append((old, new))

    print(f"fire-testing {args.test_path}"
          + (f" (mutating {args.plant_in})" if args.plant_in else "")
          + (f" [-k {args.selector}]" if args.selector else ""))
    fired = run(args.test_path, plants, args.plant_in, args.selector)
    print(f"\nRESULT: {'PASS -- the guard fires' if fired else
                       'FAIL -- the guard did NOT fire; it does not test its subject'}")
    return 0 if fired else 1


if __name__ == "__main__":
    raise SystemExit(main())
