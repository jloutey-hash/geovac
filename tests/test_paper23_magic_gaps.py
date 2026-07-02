"""Paper 23 magic-number GAP STRUCTURE pins (7th group4 cert, 2026-07-02).

The scan helper `verify_magic_ordering` uses min_gap=1e-10 -- a NECESSARY
condition only (each magic number falls on a positive-gap level boundary),
nine orders of magnitude below the real ~0.1-0.9 hw gap scale, so it could
never catch a drift collapsing a magic gap (the C2 finding of the 7th cert).
These tests pin the actual gap structure at the paper's optimal parameters
(eq:optimal, v_ls/hw = 0.1709, d_ll/hw = 0.0211):

  - the six lower magic numbers 2, 8, 20, 28, 50, 82 are exactly the SIX
    LARGEST finite gaps in the single-particle spectrum (dominance), and
  - 126 appears as a clear but NON-dominant gap (~0.107 hw), tied with the
    known sub-shell boundary at 40 and below several intermediate sub-shell
    gaps (14, 56, 104, 118) -- the honest structure disclosed in the paper.

Gap values are pinned in units of hw (scale-free). Production code is
untouched: `find_magic_numbers` / `find_optimal_vls` keep their scan
semantics; these pins guard the paper's headline against silent drift.
"""
from __future__ import annotations

import math

from geovac.nuclear.spin_orbit import find_magic_numbers

HW = 10.0
V_LS = 1.709   # v_ls/hw = 0.1709 (paper eq:optimal)
D_LL = 0.211   # d_ll/hw = 0.0211

MAGIC = [2, 8, 20, 28, 50, 82, 126]

# pinned gap/hw values at the optimal parameters (2026-07-02 probe)
EXPECTED_GAP_OVER_HW = {
    2: 0.8723,
    8: 0.5738,
    20: 0.3607,
    28: 0.3819,
    50: 0.4663,
    82: 0.3607,
    126: 0.1075,
}


def _gaps_over_hw() -> dict:
    """Finite gaps in units of hw, restricted to the physically meaningful
    region cum <= 126: boundaries above 126 (e.g. 136, 154) sit at the
    n_max=7 truncation edge, where the N=7 intruder levels that would fill
    those gaps in the real Mayer-Jensen scheme are absent -- edge artifacts,
    outside the paper's claims (which end at 126)."""
    res = find_magic_numbers(n_max=7, hw=HW, v_ls=V_LS, d_ll=D_LL)
    return {c: g / HW for c, g in res['gaps']
            if not math.isinf(g) and c <= 126}


def test_all_seven_magic_gaps_pinned() -> None:
    """Every magic number sits on a boundary whose gap magnitude is pinned
    (a collapse toward the old min_gap=1e-10 floor now FAILS)."""
    gaps = _gaps_over_hw()
    for m in MAGIC:
        assert m in gaps, f"magic number {m} not on a level boundary"
        assert abs(gaps[m] - EXPECTED_GAP_OVER_HW[m]) < 0.005, (m, gaps[m])


def test_six_lower_magics_are_the_six_dominant_gaps() -> None:
    """Dominance: the six largest finite gaps are exactly {2,8,20,28,50,82}."""
    gaps = _gaps_over_hw()
    top6 = sorted(gaps, key=gaps.get, reverse=True)[:6]
    assert sorted(top6) == [2, 8, 20, 28, 50, 82], sorted(top6)


def test_126_present_but_non_dominant() -> None:
    """126's gap is real (~0.107 hw, far above noise) but non-dominant --
    a collapse OR a silent structural promotion both break these pins."""
    gaps = _gaps_over_hw()
    g126 = gaps[126]
    assert 0.09 < g126 < 0.13, g126
    larger_nonmagic = [c for c, g in gaps.items()
                       if g > g126 and c not in MAGIC]
    # measured 10 (2026-07-02): 6,14,16,32,40,56,78,90,104,118
    assert 8 <= len(larger_nonmagic) <= 12, sorted(larger_nonmagic)
