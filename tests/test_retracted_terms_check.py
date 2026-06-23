"""Regression backing for the C16 retracted-claims / zombie-drift screen.

Guards `debug/qa/check_retracted_terms.py` (the deterministic gate that FAILs when
a withdrawn claim re-surfaces without an inline withdrawal flag). The two things
worth locking in:
  1. the registry is well-formed (every entry parses; patterns compile);
  2. the screen is NON-VACUOUS -- it genuinely fires on a re-introduced zombie and
     genuinely exempts a properly-withdrawn one. (If someone later broadens an
     exemption to the point of vacuity, this test fails.)
"""
from __future__ import annotations

import importlib.util
import pathlib
import re

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
CHECK = ROOT / "debug" / "qa" / "check_retracted_terms.py"


def _load():
    spec = importlib.util.spec_from_file_location("check_retracted_terms", CHECK)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_registry_well_formed():
    m = _load()
    assert m.REGISTRY, "registry is empty"
    for e in m.REGISTRY:
        for key in ("id", "scope", "severity", "pattern", "exempt_if_nearby", "files"):
            assert key in e, f"{e.get('id', '?')} missing {key}"
        assert e["severity"] in ("fail", "advisory")
        re.compile(e["pattern"])             # compiles
        re.compile(e["exempt_if_nearby"])    # compiles
        assert isinstance(e["files"], list) and e["files"]


@pytest.mark.parametrize("entry_id,zombie,withdrawn", [
    ("withdrawn-pythagorean-mechanism",
     "the C_3 Pythagorean refinement gives C_3 < 1 at every cutoff",          # live zombie
     "the Pythagorean refinement was WITHDRAWN 2026-06-18 as false"),         # flagged
    ("s7-structural-negative",
     "Free scalar on $S^7$: PSLQ fails on simple ring; UNKNOWN",
     "an earlier draft reported $S^7$ as a structural non-match; that was a "
     "false negative (Erratum), the ladder generates in-ring closed forms"),
    ("batch2-false-closure",
     "the wedge Dirac $D_W$ does not close the period closure",
     "the wedge Dirac $D_W$ does not close the period closure (corrected: operator-level -I)"),
])
def test_non_vacuous(entry_id, zombie, withdrawn):
    """Each fail-severity entry must FIRE on its zombie and EXEMPT its withdrawn form."""
    m = _load()
    entry = next(e for e in m.REGISTRY if e["id"] == entry_id)
    pat = re.compile(entry["pattern"], re.IGNORECASE)
    ex = re.compile(entry["exempt_if_nearby"], re.IGNORECASE)
    # zombie: pattern matches AND no exemption in the (same-line) window -> live
    assert pat.search(zombie), f"{entry_id}: pattern failed to match its own zombie"
    assert not ex.search(zombie), f"{entry_id}: exemption vacuously clears the zombie"
    # withdrawn: pattern matches AND exemption present -> cleared
    assert pat.search(withdrawn), f"{entry_id}: pattern failed to match the withdrawn form"
    assert ex.search(withdrawn), f"{entry_id}: withdrawal flag not recognised as exemption"


def test_current_corpus_passes_group1():
    """The gate is clean on the current corpus (all known zombies are flagged/fixed)."""
    m = _load()
    import sys
    argv = sys.argv
    try:
        sys.argv = ["check_retracted_terms.py", "--gate", "group1"]
        rc = m.main()
    finally:
        sys.argv = argv
    assert rc == 0, "C16 retracted-terms gate is dirty on the current group1 corpus"
