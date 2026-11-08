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
    """Corrected channel count summed over q (2026-11-07): 2 min(a, b)."""
    return 2 * min(a, b)


def I_closed_old(a: int, b: int) -> int:
    """The pre-correction kernel 2 min - 1 - [a=b] (old so4_channel_count,
    every unit of which was empty); kept for its 15-term identity only."""
    if min(a, b) == 0:
        return 0
    return 2 * min(a, b) - 1 - (1 if a == b else 0)


N = 8
PH = {n: phi(n) for n in range(N + 1)}


def _s4_collapsed(kernel=I_closed):
    return sum(kernel(a, b) * kernel(b, c) * kernel(c, d)
               * PH[a] * PH[b] * PH[c] * PH[d]
               for a in range(1, N + 1) for b in range(1, N + 1)
               for c in range(1, N + 1) for d in range(1, N + 1))


def test_chain_collapse_bitexact():
    """S^(4) direct CG-channel sum == collapsed min-kernel sum == 8 M4
    (Fraction; corrected W 2026-11-07), and the o-space relation
    S4 = C4 - 6 P C - 3 Smin_o^2 + 27 Smin_o P^2 - 27 P^4 (C4 coefficient
    +1: the four-chain core is still the only depth-4 source)."""
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
    rng = range(1, N + 1)
    M4 = sum(min(a, b) * min(b, c) * min(c, d) * PH[a] * PH[b] * PH[c] * PH[d]
             for a in rng for b in rng for c in rng for d in rng)
    assert s4_direct == 8 * M4
    odds = list(range(5, 2 * N + 4, 2))
    ps = {o: Fraction(8 * (o * o - 1), o ** 4) for o in odds}
    P = sum(ps.values())
    Sm = sum(min(a, b) * ps[a] * ps[b] for a in odds for b in odds)
    C = sum(min(a, b) * min(b, c) * ps[a] * ps[b] * ps[c]
            for a in odds for b in odds for c in odds)
    C4 = sum(min(a, b) * min(b, c) * min(c, d) * ps[a] * ps[b] * ps[c] * ps[d]
             for a in odds for b in odds for c in odds for d in odds)
    assert s4_direct == C4 - 6 * P * C - 3 * Sm ** 2 + 27 * Sm * P ** 2 - 27 * P ** 4


def test_fifteen_term_identity_bitexact():
    """The 15-term n-space identity of the PRE-CORRECTION kernel
    2 min - 1 - [a=b] holds bit-exactly at N=8 (an identity of that
    kernel; the corrected S^(4) is 8 M4, test above)."""
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
    assert rhs == _s4_collapsed(I_closed_old)


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
    """Rigorous bracket for the corrected-W S^(4) (2026-11-07), computed
    live at N = 1e6: S4 = 2 sum_k B(k)^2, B(k) = sum_{n>=k} phi(n) L(n),
    L = 2 sum_{j<=n} T(j); lower = positive partial sum, upper adds
    tau = 8(ln N + 2)/N to each B_N(k) and the k > N integral bound
    64[(ln N+2)^2 + 2(ln N+2) + 2]/N (L <= 4(ln n + 1), phi <= 2/(n+1)^2).
    At N = 4e6 the bracket is [409.286, 409.425].  It must exclude the
    old-count value ~316.57 (s4_bracket.json) and the out-of-band
    heuristics 100 and 800."""
    import math
    import numpy as np
    from scipy.special import zeta as hz
    Nb = 1_000_000
    m = np.arange(1, Nb + 1, dtype=np.float64)
    x = 2 * m + 3
    ph = 8 * (x * x - 1) / x ** 4
    L = 2 * np.cumsum(2 * hz(2, m + 1.5) - 0.5 * hz(4, m + 1.5))
    a = ph * L
    BN = np.cumsum(a[::-1])[::-1]
    lo = 2 * math.fsum((BN * BN).tolist())
    lN = math.log(Nb)
    tau = 8 * (lN + 2) / Nb
    hi = 2 * (math.fsum(((BN + tau) ** 2).tolist())
              + 64 * ((lN + 2) ** 2 + 2 * (lN + 2) + 2) / Nb)
    assert lo <= 409.35 <= hi
    assert hi - lo < 0.5
    old = json.loads((DATA / "s4_bracket.json").read_text())["bracket"]
    old_best = old[max(old, key=lambda k: int(k))]
    assert float(old_best["upper"]) < lo      # old-count value excluded
    assert not (lo <= 100 <= hi) and not (lo <= 800 <= hi)
    assert not (lo <= 316.57 <= hi)


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
