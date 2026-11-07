"""Paper 63 prop:tree_bound -- the rung bound for trees of s-type densities on the one-scale slice.

Claim: an even section x^(2k) K0(x)^a (x K1(x))^b, a + b = n - 1, integrated against K0(tx), equals R(t, theta_t)
applied to F_n(t) = int K0(tx) K0(x)^(n-1) dx, with R built by the recursion eq:tree_bound_recursion and the
integration-by-parts rules theta_x -> -(1 + theta_t), x^2 -> t^-2 theta_t^2. So the s-type tree period is a
differential operator applied to F_n and sits at or below rung n. Record docs/predictions/2026-11-06_rung_bound_trees.md.

Guards, each naming the wrong answer it rejects:
(1) exact: for n = 3 the product K0 (x K1) maps to R = (1 + theta)/2. Rejects the untransposed -theta/2, the
    denominator 1 for 2, and the sign-flipped -(1 + theta)/2.
(2) numeric: g_T = R F_n at two values of t for n = 3 products with b = 1, 2 and an x^2 weight, and for n = 4 with
    b = 3. Rejects a wrong recursion denominator (misses by 25-80%), the x^2 rule with t^-1 for t^-2, and the
    transpose without its -1.
(3) the scaled leaf: int (tx) K1(tx) K0^2 dx = -theta_t F_3.
"""
from __future__ import annotations

from fractions import Fraction
from typing import Dict, Tuple

import mpmath as mp
import pytest

Op = Dict[Tuple[int, int], Fraction]  # {(p, m): c} = c t^p theta_t^m, t-power on the left


def _add(a: Op, b: Op, s: Fraction = Fraction(1)) -> Op:
    r = dict(a)
    for k, v in b.items():
        r[k] = r.get(k, Fraction(0)) + s * v
        if r[k] == 0:
            del r[k]
    return r


def _theta_left(a: Op) -> Op:
    r: Op = {}
    for (p, m), c in a.items():
        r = _add(r, {(p, m): c * p, (p, m + 1): c})
    return r


def _x2(a: Op) -> Op:
    a2 = _theta_left(_theta_left(a))
    return {(p - 2, m): c for (p, m), c in a2.items()}


def _theta_x(a: Op) -> Op:
    return _add(_add({}, a, Fraction(-1)), _theta_left(a), Fraction(-1))


def _section(n: int, k: int, a: int, b: int) -> Op:
    if b == 0:
        r: Op = {(0, 0): Fraction(1)}
        for _ in range(k):
            r = _x2(r)
        return r
    s1 = _section(n, k, a + 1, b - 1)
    r = _add({}, s1, Fraction(2 * k))
    r = _add(r, _theta_x(s1), Fraction(-1))
    if b >= 2:
        r = _add(r, _section(n, k + 1, a + 2, b - 2), Fraction(-(b - 1)))
    return {key: c / (a + 1) for key, c in r.items() if c != 0}


def _theta_pow_K0(m: int) -> Dict[Tuple[int, int], Fraction]:
    f = {(0, 0): Fraction(1)}
    for _ in range(m):
        g: Dict[Tuple[int, int], Fraction] = {}
        for (i, kind), c in f.items():
            terms = [((i, 0), c * i), ((i + 1, 1), -c)] if kind == 0 else [((i, 1), c * (i - 1)), ((i + 1, 0), -c)]
            for key, v in terms:
                g[key] = g.get(key, Fraction(0)) + v
        f = {k2: v for k2, v in g.items() if v != 0}
    return f


_BREAKS = [0, mp.mpf(1) / 8, mp.mpf(1) / 2, 1, 2, 4, 8, 16, 32, 64, mp.inf]


def _theta_F(n: int, m: int, t: mp.mpf) -> mp.mpf:
    terms = _theta_pow_K0(m)

    def f(x):
        y = t * x
        v = sum(mp.mpf(c.numerator) / c.denominator * y**i * mp.besselk(kind, y) for (i, kind), c in terms.items())
        return v * mp.besselk(0, x) ** (n - 1)

    return mp.quad(f, _BREAKS)


def _apply(n: int, R: Op, t: mp.mpf) -> mp.mpf:
    return sum(mp.mpf(c.numerator) / c.denominator * t**p * _theta_F(n, m, t) for (p, m), c in R.items())


def _direct(n: int, k: int, a: int, b: int, t: mp.mpf) -> mp.mpf:
    return mp.quad(lambda x: mp.besselk(0, t * x) * x ** (2 * k) * mp.besselk(0, x) ** a
                   * (x * mp.besselk(1, x)) ** b, _BREAKS)


def test_k0_xk1_maps_to_one_plus_theta_over_two() -> None:
    assert _section(3, 0, 1, 1) == {(0, 0): Fraction(1, 2), (0, 1): Fraction(1, 2)}


def test_x2_weight_maps_to_t_minus2_theta2() -> None:
    """x^2 K0^2 -> t^-2 theta^2 exactly; rejects t^-1 theta^2 and theta^2 without the t-power."""
    assert _section(3, 1, 2, 0) == {(-2, 2): Fraction(1)}
    # (xK1)^2 = -theta(K0 xK1) - x^2 K0^2 against K0(tx): a fixed exact operator with a t^-2 part
    assert _section(3, 0, 0, 2) == {(0, 0): Fraction(1, 2), (0, 1): Fraction(1), (0, 2): Fraction(1, 2),
                                    (-2, 2): Fraction(-1)}


@pytest.mark.parametrize("n,k,a,b", [(3, 0, 1, 1), pytest.param(3, 0, 0, 2, marks=pytest.mark.slow),
                                     pytest.param(3, 1, 1, 1, marks=pytest.mark.slow),
                                     pytest.param(4, 0, 0, 3, marks=pytest.mark.slow)])
def test_section_is_operator_on_F_n(n: int, k: int, a: int, b: int) -> None:
    mp.mp.dps = 20
    R = _section(n, k, a, b)
    for ts in ("0.5", "2.5"):
        t = mp.mpf(ts)
        lhs = _direct(n, k, a, b, t)
        rhs = _apply(n, R, t)
        assert abs(lhs - rhs) / abs(lhs) < mp.mpf("1e-15"), (n, k, a, b, ts)


def test_scaled_leaf_is_theta_image() -> None:
    mp.mp.dps = 20
    t = mp.mpf("1.5")
    lhs = mp.quad(lambda x: (t * x) * mp.besselk(1, t * x) * mp.besselk(0, x) ** 2, _BREAKS)
    rhs = -_theta_F(3, 1, t)
    assert abs(lhs - rhs) / abs(lhs) < mp.mpf("1e-15")
    assert abs(lhs - _theta_F(3, 1, t)) / abs(lhs) > mp.mpf("0.5")  # the sign matters
