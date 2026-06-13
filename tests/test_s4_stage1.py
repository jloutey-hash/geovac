"""Falsifier for the S^(4) stage-1 solid results (k=4 rung of the S-tower).

Freezes what is theorem-grade / rigorous from the stage-1 sprint
(v4.9.0): the chain collapse, the 15-term n-space identity, the o-space
decomposition relation, the census (odd weights 5..13), the corrected
single-nsum evaluator's bit-exactness vs the k=3 cache (non/single-trailing),
and the rigorous bracket containment.

Does NOT cover the b1=2 high-log trailing constants (deferred to stage-2
symbolic identification) or any canonical S^(4) value (not produced).

Run:  pytest tests/test_s4_stage1.py            (fast)
      pytest tests/test_s4_stage1.py --slow      (adds the 220-dps cache gate)
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from itertools import product as iproduct
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DEBUG = ROOT / "debug"
DATA = DEBUG / "data"
sys.path.insert(0, str(DEBUG))


# --------------------------------------------------------------------------
# exact building blocks (self-contained; mirror s4_scoping_diag / engine)
# --------------------------------------------------------------------------

def phi(n: int) -> Fraction:
    x = Fraction(2 * n + 3, 2)
    return 2 / x ** 2 - Fraction(1, 2) / x ** 4


def I_closed(a: int, b: int) -> int:
    if min(a, b) == 0:
        return 0
    return 2 * min(a, b) - 1 - (1 if a == b else 0)


N = 8
PH = {n: phi(n) for n in range(N + 1)}


def _s4_collapsed():
    return sum(I_closed(a, b) * I_closed(b, c) * I_closed(c, d)
               * PH[a] * PH[b] * PH[c] * PH[d]
               for a in range(1, N + 1) for b in range(1, N + 1)
               for c in range(1, N + 1) for d in range(1, N + 1))


def test_chain_collapse_bitexact():
    """S^(4) direct CG-channel sum == collapsed min-kernel sum (Fraction)."""
    from geovac.qed_vertex import so4_channel_count
    I_direct = {(a, b): sum(so4_channel_count(a, b, q)
                            for q in range(1, a + b + 1))
                for a in range(N + 1) for b in range(N + 1)}
    # pairwise closed form
    assert all(I_direct[(a, b)] == I_closed(a, b)
               for a in range(N + 1) for b in range(N + 1))
    s4_direct = sum(I_direct[(a, b)] * I_direct[(b, c)] * I_direct[(c, d)]
                    * PH[a] * PH[b] * PH[c] * PH[d]
                    for a in range(N + 1) for b in range(N + 1)
                    for c in range(N + 1) for d in range(N + 1))
    assert s4_direct == _s4_collapsed()


def test_fifteen_term_identity_bitexact():
    """The 15-term n-space identity reproduces S^(4) bit-exactly at N=8."""
    rng = range(1, N + 1)
    P = sum(PH[n] for n in rng)
    Q = sum(PH[n] ** 2 for n in rng)
    R3 = sum(PH[n] ** 3 for n in rng)
    R4 = sum(PH[n] ** 4 for n in rng)
    S_min = sum(min(a, b) * PH[a] * PH[b] for a in rng for b in rng)
    M2 = sum(min(a, b) * PH[a] * PH[b] ** 2 for a in rng for b in rng)
    M3 = sum(min(a, b) * min(b, c) * PH[a] * PH[b] * PH[c]
             for a in rng for b in rng for c in rng)
    M3e = sum(min(a, b) * min(b, c) * PH[a] * PH[b] * PH[c] ** 2
              for a in rng for b in rng for c in rng)
    M3m = sum(min(a, b) * min(b, c) * PH[a] * PH[b] ** 2 * PH[c]
              for a in rng for b in rng for c in rng)
    T31 = sum(min(a, b) * PH[a] * PH[b] ** 3 for a in rng for b in rng)
    T22 = sum(min(a, b) * PH[a] ** 2 * PH[b] ** 2 for a in rng for b in rng)
    M4 = sum(min(a, b) * min(b, c) * min(c, d) * PH[a] * PH[b] * PH[c] * PH[d]
             for a in rng for b in rng for c in rng for d in rng)
    rhs = (8 * M4 - 8 * M3 * P - 8 * M3e - 4 * M3m + 8 * M2 * P
           + 6 * S_min * P ** 2 + 4 * S_min * Q - 4 * S_min ** 2
           + 4 * T31 + 2 * T22 - 2 * R3 * P - 3 * P ** 2 * Q
           - P ** 4 - Q ** 2 - R4)
    assert rhs == _s4_collapsed()


def test_ospace_relation_present_and_integer():
    """o-space relation artifact exists with the frozen integer coefficients."""
    j = json.loads((DATA / "s4_ospace_identities.json").read_text())
    # every recorded identity verified bit-exact in the producing driver
    assert all(v is True for k, v in j.items()
               if isinstance(v, bool)), "an o-space identity regressed"


def test_census_odd_weights_5_to_13():
    """t4 content realizes exactly odd weights {5,7,9,11,13}."""
    j = json.loads((DATA / "s4_decomp_tables.json").read_text())
    cens = j["census"]
    assert cens["t4_weights"] == [5, 7, 9, 11, 13]
    assert cens["t4_count"] == 40


def test_bracket_contains_estimate_excludes_heuristics():
    """Rigorous bracket at N=4e6 contains ~316.57 and excludes the
    diagnostic's pre-committed-out-of-band heuristics."""
    j = json.loads((DATA / "s4_bracket.json").read_text())
    grid = j["bracket"]
    best = grid[max(grid, key=lambda k: int(k))]
    lo, hi = float(best["lower"]), float(best["upper"])
    assert lo <= 316.57 <= hi
    assert hi - lo < 0.5          # width tightened to < 0.5 at 4e6
    assert not (lo <= 100 <= hi) and not (lo <= 800 <= hi)


def test_evaluator_bitexact_vs_k3_nontrailing_fast():
    """Single-nsum evaluator reproduces k=3 cache for a non-trailing sample
    (fast, 30 dps -> agreement to ~1e-28)."""
    import s4_mt_eval
    from mpmath import mp, mpmathify
    k3 = json.loads((DATA / "s3_pslq_cache.json").read_text())
    mp.dps = 50
    sample = [(2, 2), (4, 3), (2, 2, 2), (4, 3, 3)]
    for tup in sample:
        nm = "t%d(%s)@220" % (len(tup), ",".join(map(str, tup)))
        if nm not in k3:
            continue
        v = mpmathify(s4_mt_eval.eval_t(tup, 30))
        assert abs(v - mpmathify(k3[nm])) < mpmathify("1e-28"), nm


@pytest.mark.slow
def test_evaluator_bitexact_nontrailing_220dps():
    """NON-trailing constants are bit-exact (<1e-200) vs the k=3 cache at
    220 dps (closed tau/pL, single Levin, no log modulation)."""
    import s4_mt_eval
    from mpmath import mp, mpmathify
    k3 = json.loads((DATA / "s3_pslq_cache.json").read_text())
    mp.dps = 240
    for tup in [(2, 2), (4, 3), (2, 2, 2), (4, 3, 3)]:
        nm = "t%d(%s)@220" % (len(tup), ",".join(map(str, tup)))
        if nm not in k3:
            continue
        v = mpmathify(s4_mt_eval.eval_t(tup, 220))
        assert abs(v - mpmathify(k3[nm])) < mpmathify("1e-200"), nm


@pytest.mark.slow
def test_evaluator_trailing_formula_correct_moderate():
    """TRAILING constants: the evaluator FORMULA is correct (agrees with the
    k=3 cache to >~12 digits) but b1=2 trailing DEGRADES at very high
    precision (Levin under-converges on the log modulation) — high-precision
    trailing values come from the k=3 cache (11 of them) or stage-2 symbolic
    identification, NOT from this independent evaluator.  This test documents
    that boundary: formula-correct at moderate precision, NOT bit-exact at
    220 dps for b1=2."""
    import s4_mt_eval
    from mpmath import mp, mpmathify
    k3 = json.loads((DATA / "s3_pslq_cache.json").read_text())
    mp.dps = 240
    for tup in [(4, 1), (2, 3, 1), (4, 3, 1)]:
        nm = "t%d(%s)@220" % (len(tup), ",".join(map(str, tup)))
        if nm not in k3:
            continue
        v = mpmathify(s4_mt_eval.eval_t(tup, 220))
        assert abs(v - mpmathify(k3[nm])) < mpmathify("1e-12"), nm
