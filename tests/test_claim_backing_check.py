"""Mirror test for C22 (`debug/qa/check_test_claim_backing.py`).

C22 is C13 run backwards: C13 asserts every test a paper cites exists;
C22 asserts every claim a test backs still exists and is not retracted.

The gate lives in `debug/`, which is not collected by pytest, so without
this file nothing in the ordinary test run would notice if it broke. That
is the same shape as the failure C22 itself guards against -- a check
whose failure has no consequence -- so the mirror is not optional.

Provenance: INTERNAL THEOREM (the checks are exact set operations over
the matrix, the papers tree and the test tree; no tolerances, no fitting).
"""
from __future__ import annotations

import json
import pathlib
import re
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "debug" / "qa"))

import check_test_claim_backing as C22          # noqa: E402


def test_gate_passes_on_the_current_corpus():
    """The shipped baseline plus the current tree is clean."""
    base = C22.load_baseline()
    assert C22.check_a(base, verbose=False) == 0
    assert C22.check_b(verbose=False) == 0
    assert C22.check_c(verbose=False) == 0
    assert C22.check_d(base, verbose=False) == 0


def test_testfile_regex_has_a_left_boundary():
    """`test_lih.py` must not match inside `retest_lih.py`.

    This exact substring bug reported two phantom dangling matrix rows,
    which were relayed as findings twice before being caught. The boundary
    is the fix; this pins it.
    """
    assert not C22.TESTFILE.findall("debug/x/chemistry_solver_retest_lih.py")
    assert C22.TESTFILE.findall("backed by `test_lih.py` here") == ["test_lih.py"]


def test_disowned_names_are_not_read_as_backing():
    """A row may NAME a test in order to retire it; that is prose.

    The row that triggered this reads "Replaces stale
    `debug/test_harmonic_phase_lock.py` cite" while its real backing,
    test_paper8_sigma_bond_selection.py, exists and passes.
    """
    assert C22.DISOWNED.search("Replaces stale `debug/test_foo.py` cite")
    assert C22.DISOWNED.search("archived `debug/.../test_bar.py`")
    assert not C22.DISOWNED.search("| 14 | live claim | `test_baz.py` | ok |")


def test_both_paper_filename_conventions_are_known():
    """`paper_14_*.tex` AND `Paper_7_*.tex` both exist in the corpus.

    Matching only the lower-case form makes check B report Papers 0, 7 and
    8 as nonexistent -- a gate firing confidently on the wrong thing.
    """
    papers = C22.known_papers()
    for n in ("0", "7", "8", "14", "60"):
        assert n in papers, f"paper {n} not recognised: {sorted(papers)[:12]}"


@pytest.mark.parametrize("check", ["a", "b", "c", "d"])
def test_each_check_can_fire(check):
    """Two-way discrimination, per check.

    A gate proven only to stay quiet is indistinguishable from one that
    cannot fire at all. Every check must be shown to fire on a defect it
    is meant to catch.
    """
    if check == "a":
        row = "| 14 | c | `test_c22_absent_probe.py` | m | BACKED |"
        assert set(C22.TESTFILE.findall(row)) - C22.existing_tests()
    elif check == "b":
        probe = ROOT / "tests" / "test_paper9999_probe.py"
        probe.write_text("# probe\n", encoding="utf-8")
        try:
            assert C22.check_b(verbose=False) > 0
        finally:
            probe.unlink()
        assert C22.check_b(verbose=False) == 0
    elif check == "c":
        entries = [e for e in C22.C16.REGISTRY if e.get("severity") == "fail"]
        row = "| 26 | angular ERI density is $2.76\\% | test_probe.py |"
        assert any(re.search(e["pattern"], row) for e in entries)
    else:
        assert C22.check_d({"paper_backing_debug_imports": []},
                           verbose=False) > 0


def test_baseline_is_declared_and_shrinking_is_the_intent():
    """A ratchet that hides its own size is how debt becomes permanent.

    Corrected 2026-09-12.  This asserted 3 while the baseline had held 4 since
    2026-09-07, so it had been RED and unnoticed -- the guard meant to watch the
    ratchet was itself stale.

    The count grew for a reason its own failure message misdiagnosed.  Growth
    has TWO causes and they are opposite in meaning:

      (a) NEW DEBT -- a freshly written paper-backing test reaches into the
          prunable tree.  That is what the ratchet exists to stop.
      (b) BETTER DETECTION -- the checker learned to see a dependency that was
          always there.  `test_paper58_headline_numbers.py` puts `debug/` on
          `sys.path` and then imports a BARE module name, carrying no `debug.`
          token for the original regex to find; adding DEBUG_SYSPATH surfaced
          it.  The debt did not grow, the blind spot shrank.

    So the assertion pins the MEMBERS, not just the count: a changed roster is
    reported by name, which distinguishes (a) from (b) on sight.
    """
    base = C22.load_baseline()
    assert "_note" in base, "the baseline must say what it is for"
    assert isinstance(base["paper_backing_debug_imports"], list)
    known = {
        "test_paper26_entanglement.py",
        "test_paper27_entropy.py",
        "test_paper27_entropy_locus.py",
        "test_paper58_headline_numbers.py",   # sys.path route, found 2026-09-07
    }
    got = set(base["paper_backing_debug_imports"])
    assert got == known, (
        "baseline roster changed. ADDED: "
        f"{sorted(got - known)} -- new debt, unless the checker's reach grew; "
        f"REMOVED: {sorted(known - got)} -- debt genuinely paid down, update "
        "this roster."
    )
