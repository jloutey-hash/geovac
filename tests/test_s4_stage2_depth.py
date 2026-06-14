"""Falsifier for the S^(4) stage-2 depth verdict (v4.10.0).

Freezes the theorem-grade result: realized depth <= 3 for S^(4).
- Structural fact: depth-4 atoms enter S^(4) ONLY through the C4 core
  (so the verdict reduces to C4's depth-4 reduction).
- w >= 9: stuffle relations span the FULL depth-4 word space (rank =
  #depth-4 words), so C4's depth-4 content collapses to depth <= 3 —
  verified by exact mod-p rank at two independent primes.
- w = 5, 7: stuffle is rank-deficient (documented); settled separately by
  the low-weight depth filtration (not a code assertion).

Run:  pytest tests/test_s4_stage2_depth.py            (fast: w=5,7,9 + structure)
      pytest tests/test_s4_stage2_depth.py --slow      (adds w=11,13 full-rank)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DEBUG = ROOT / "debug"
DATA = DEBUG / "data"
sys.path.insert(0, str(DEBUG))


def test_depth4_only_in_C4():
    """The depth-4 atoms appear ONLY in the C4 core table — the structural
    fact that reduces the S^(4) depth verdict to C4's depth-4 reduction."""
    d = json.loads((DATA / "s4_decomp_tables.json").read_text())
    for name, tab in d["tables"].items():
        d4 = [ks for ks in tab
              if (s := eval(ks))[0] == "t" and len(s[1]) == 4]  # noqa: S307
        if name == "C4":
            assert d4, "C4 must carry the depth-4 atoms"
        else:
            assert not d4, "depth-4 atom leaked into core %s" % name


def _verdict(w, prime):
    import s4_stage2_assembly as A
    A.P = prime
    return A.analyze(w, A.load_C4())


@pytest.mark.parametrize("w,full_rank", [(9, 35)])
def test_w9_full_rank_collapse_two_primes(w, full_rank):
    """w=9: stuffle spans the full depth-4 space (rank = #words) and C4's
    depth-4 part reduces — at BOTH primes."""
    import s4_stage2_assembly as A
    for pr in A.PRIMES:
        res = _verdict(w, pr)
        assert res["n_depth4_words"] == full_rank
        assert res["stuffle_rank_on_D4"] == full_rank, (w, pr)
        assert res["C4_depth4_reduces_by_stuffle"] is True, (w, pr)


def test_w5_w7_stuffle_rank_deficient():
    """w=5,7: stuffle is rank-deficient on the depth-4 space (documented;
    settled by the low-weight depth filtration, not by code)."""
    import s4_stage2_assembly as A
    r5 = _verdict(5, A.PRIMES[0])
    r7 = _verdict(7, A.PRIMES[0])
    assert r5["stuffle_rank_on_D4"] == 0 and r5["n_depth4_words"] == 1
    assert r7["stuffle_rank_on_D4"] == 7 and r7["n_depth4_words"] == 10
    assert r5["C4_depth4_reduces_by_stuffle"] is False
    assert r7["C4_depth4_reduces_by_stuffle"] is False


@pytest.mark.slow
@pytest.mark.parametrize("w,full_rank", [(11, 84), (13, 165)])
def test_high_weight_full_rank_collapse_two_primes(w, full_rank):
    """w=11,13: full-rank depth-4 collapse, both primes (the bulk of the
    proof; slower rank computations)."""
    import s4_stage2_assembly as A
    for pr in A.PRIMES:
        res = _verdict(w, pr)
        assert res["n_depth4_words"] == full_rank
        assert res["stuffle_rank_on_D4"] == full_rank, (w, pr)
        assert res["C4_depth4_reduces_by_stuffle"] is True, (w, pr)
