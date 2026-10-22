"""Paper 59 sec:density_ring, the correlated triangle Y at gamma = 3/2.

Record docs/predictions/2026-10-22_correlated_triangle_g32.md (GT1-GT3, committed 4b1be301 before
computing).  Written as a separate step from the result; no debug/ import.  The exact pieces found by
search (the order-nine annihilator L9, the chain ra < R2 < V4 < V5 < K6 < L9, the order-three piece F
and the map A from F_3's operator) are loaded from tests/data/paper59_ladder_2026_10_22.json and
re-verified here.  The chain and F are theta-form coefficient lists in x = t^2 with theta = t d/dt.
Not re-run here: that the four order-one pieces have algebraic solutions and that the order-two piece
is irreducible (local-exponent bounds, in the record).
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest
import sympy as sp

import test_paper59_density_ring as DR
import test_paper59_generic_slice as GS

t, u = GS.t, GS.u
x = sp.Symbol("x")
DATA = json.loads((Path(__file__).parent / "data" / "paper59_ladder_2026_10_22.json").read_text())["g32"]
CHAIN = ["chain_ra", "chain_R2", "chain_V4", "chain_V5", "chain_K6"]
L3 = GS._coeff_list([4 * u**3, -(u + 1) ** 3])     # F_3's operator 4 theta^3 - t^2 (theta + 1)^3


def _op(key):
    return [sp.expand(sp.sympify(s.replace("^", "**"), locals={"x": x}).subs(x, t**2)) for s in DATA[key]]


def _L9(sym=u):
    return [sp.sympify(q, locals={"u": sym}) for q in DATA["L9_Q"]]


def _map():
    N = [sp.sympify(s.replace("^", "**"), locals={"x": x}).subs(x, t**2) for s in DATA["A_numerators_theta_0_1_2"]]
    return [sp.cancel(t * n / (t**2 - 4) ** 4) for n in N]


def _zero(rem):
    return all(sp.cancel(r) == 0 for r in rem)


def _theta_to_Q(cs):
    """theta-form coefficients c_j(t) (even in t) -> Q-form [Q_i(u)], sum_i t^(2i) Q_i(theta), in DR's u."""
    polys = [sp.Poly(sp.expand(c), t) for c in cs]
    deg = max(p.degree() for p in polys)
    return [sp.expand(sum(p.coeff_monomial(t ** (2 * i)) * DR.u**j for j, p in enumerate(polys)))
            for i in range(deg // 2 + 1)]


@pytest.mark.slow
def test_g32_composition_series():
    """Rejects 'Y's annihilator has a piece of order above three' and a wrong chain: each of ra, R2, V4,
    V5, K6 right-divides the next exactly, with orders 1, 2, 3, 4, 6, and L9 = (F up to a left factor)
    o K6 exactly.  The pieces have orders 1, 1, 1, 1, 2, 3."""
    ops = [_op(c) for c in CHAIN]
    assert [len(o) - 1 for o in ops] == [1, 2, 3, 4, 6]
    for a, b in zip(ops, ops[1:]):
        assert _zero(GS._rdiv(b, a)[1])
    L9 = GS._coeff_list(_L9())
    assert len(L9) - 1 == 9
    q, rem = GS._rdiv(L9, ops[-1])
    F = _op("factor_F")
    r = sp.cancel(q[-1] / F[-1])
    assert _zero(rem) and len(q) == len(F) and all(sp.cancel(a - r * b) == 0 for a, b in zip(q, F))


@pytest.mark.slow
def test_g32_k3_piece_is_exact():
    """Rejects 'the order-three piece is not F_3's' and a wrong stored map: F o A = 0 modulo F_3's
    operator exactly, A is nonzero of order two (so injective on the irreducible F_3 solutions, an
    equivalence), and a 1/1000 change in one coefficient of A breaks the identity."""
    F, A = _op("factor_F"), _map()
    assert _zero(GS._rdiv(GS._compose(F, A), L3)[1])
    assert len(A) == 3 and any(sp.cancel(a) != 0 for a in A)
    A2 = list(A)
    A2[1] = A2[1] * sp.Rational(1001, 1000)
    assert not _zero(GS._rdiv(GS._compose(F, A2), L3)[1])


def test_g32_y_needs_the_k3_piece():
    """Rejects 'Y sits below the K3 piece' (K6 Y = 0): on Y itself, from the kite evaluation of the
    gamma = 1 guard, L9 vanishes to 1e-9 of its largest t-power group at t = 0.8, 1.3, 1.7, a 1e-3
    change in one coefficient of L9 does not, and K6, the right factor below F, leaves more than half
    the largest group.  Normalisation: Y(0.8) against the record's 32-digit grid."""
    ts = (0.8, 1.3, 1.7)
    thY = DR._theta_Y(1.5, ts, J=9)
    assert abs(thY[0.8][0] / 0.001088348769854764451 - 1) < 1e-10
    L9 = _L9(DR.u)
    bad = list(L9)
    bad[3] = L9[3] * sp.Rational(1001, 1000)
    K6 = _theta_to_Q(_op("chain_K6"))
    for tv in ts:
        assert DR._apply_Q_num(L9, thY[tv], tv) < 1e-9
        assert DR._apply_Q_num(bad, thY[tv], tv) > 1e-7
        assert DR._apply_Q_num(K6, thY[tv], tv) > 0.5
