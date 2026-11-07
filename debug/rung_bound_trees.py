"""Rung bound for s-type trees on the one-scale slice (record docs/predictions/2026-11-06_rung_bound_trees.md).

Claim: an even section T(x) = x^(2k) K0(x)^a (x K1(x))^b, a + b = n - 1, is Q(x^2, theta) K0^(n-1) with Q built by
a triangular recursion; integrating by parts against K0(tx) turns theta_x into -(1 + theta_t) and x^2 into
t^-2 theta_t^2, so g_T(t) = int K0(tx) T(x) dx = R(t, theta_t) F_n(t), F_n = int K0(tx) K0(x)^(n-1) dx.

Operators on t-functions are represented as dicts {(p, m): c}: the term c * t^p * theta_t^m (t-power on the LEFT).
Applying R to F_n needs theta_t^m F_n = int (theta^m K0)(tx) K0^(n-1) dx, with theta K0 = -y K1, theta (y K1) = -y^2 K0.
"""
from __future__ import annotations

import sys
from fractions import Fraction
from functools import lru_cache
from typing import Dict, Tuple

import mpmath as mp

Op = Dict[Tuple[int, int], Fraction]


def op_add(a: Op, b: Op, s: Fraction = Fraction(1)) -> Op:
    r = dict(a)
    for k, v in b.items():
        r[k] = r.get(k, Fraction(0)) + s * v
        if r[k] == 0:
            del r[k]
    return r


def theta_left(a: Op) -> Op:
    """theta_t o (t^p theta^m) = t^p (p + theta) theta^m."""
    r: Op = {}
    for (p, m), c in a.items():
        r = op_add(r, {(p, m): c * p})
        r = op_add(r, {(p, m + 1): c})
    return r


def x2_image(a: Op) -> Op:
    """g_{x^2 T} = t^-2 theta_t^2 g_T: left-multiply by t^-2 theta^2."""
    a2 = theta_left(theta_left(a))
    return {(p - 2, m): c for (p, m), c in a2.items()}


def theta_image(a: Op) -> Op:
    """g_{theta_x T} = -(1 + theta_t) g_T."""
    r = op_add({}, a, Fraction(-1))
    return op_add(r, theta_left(a), Fraction(-1))


def section_op(n: int, k: int, a: int, b: int, plant: str = "") -> Op:
    """Operator R with int K0(tx) x^(2k) K0^a (xK1)^b dx = R F_n."""
    assert a + b == n - 1 and a >= 0 and b >= 0 and k >= 0
    if b == 0:
        r: Op = {(0, 0): Fraction(1)}
        for _ in range(k):
            r = x2_image(r)
        return r
    # theta(S1) = 2k S1 - (a+1) T - (b-1) x^(2k+2) K0^(a+2) (xK1)^(b-2), S1 = x^(2k) K0^(a+1) (xK1)^(b-1)
    s1 = section_op(n, k, a + 1, b - 1, plant)
    r = op_add({}, s1, Fraction(2 * k))
    r = op_add(r, theta_image(s1), Fraction(-1))
    if b >= 2:
        r = op_add(r, section_op(n, k + 1, a + 2, b - 2, plant), Fraction(-(b - 1)))
    denom = a + 1 if plant != "denom" else a + 2
    return {key: c / denom for key, c in r.items()}


def theta_pow_K0(m: int):
    """theta^m K0(y) = sum c_i y^i K_{i mod 2}(y): returned as {(i, kind): c}, kind 0 = K0, 1 = K1."""
    f = {(0, 0): Fraction(1)}
    for _ in range(m):
        g: Dict[Tuple[int, int], Fraction] = {}
        for (i, kind), c in f.items():
            # theta (y^i K0) = i y^i K0 - y^(i+1) K1 ; theta (y^i K1) = (i-1) y^i K1 - y^(i+1) K0
            if kind == 0:
                terms = [((i, 0), c * i), ((i + 1, 1), -c)]
            else:
                terms = [((i, 1), c * (i - 1)), ((i + 1, 0), -c)]
            for key, v in terms:
                g[key] = g.get(key, Fraction(0)) + v
        f = {k2: v for k2, v in g.items() if v != 0}
    return f


def quad(f) -> mp.mpf:
    return mp.quad(f, [0, mp.mpf(1) / 8, mp.mpf(1) / 2, 1, 2, 4, 8, 16, 32, 64, mp.inf])


@lru_cache(maxsize=None)
def theta_F(n: int, m: int, t_str: str) -> mp.mpf:
    t = mp.mpf(t_str)
    terms = theta_pow_K0(m)

    def integrand(x):
        y = t * x
        val = mp.mpf(0)
        for (i, kind), c in terms.items():
            val += mp.mpf(c.numerator) / c.denominator * y**i * (mp.besselk(0, y) if kind == 0 else mp.besselk(1, y))
        return val * mp.besselk(0, x) ** (n - 1)

    return quad(integrand)


def apply_op(n: int, R: Op, t_str: str) -> mp.mpf:
    t = mp.mpf(t_str)
    return sum(mp.mpf(c.numerator) / c.denominator * t**p * theta_F(n, m, t_str) for (p, m), c in R.items())


def direct(n: int, k: int, a: int, b: int, t_str: str) -> mp.mpf:
    t = mp.mpf(t_str)
    return quad(lambda x: mp.besselk(0, t * x) * x ** (2 * k) * mp.besselk(0, x) ** a * (x * mp.besselk(1, x)) ** b)


CASES = {3: [(0, 2, 0), (0, 1, 1), (0, 0, 2), (1, 2, 0), (1, 1, 1)],
         4: [(0, 2, 1), (0, 1, 2), (0, 0, 3), (1, 2, 1)]}
TS = {3: ["0.5", "1.5", "2.5"], 4: ["0.5", "2", "2.5"]}


def run(dps: int, plant: str = "") -> float:
    mp.mp.dps = dps
    theta_F.cache_clear()
    worst = mp.mpf(0)
    for n, cases in CASES.items():
        for (k, a, b) in cases:
            R = section_op(n, k, a, b, plant)
            for ts in TS[n]:
                lhs = direct(n, k, a, b, ts)
                rhs = apply_op(n, R, ts)
                rel = abs(lhs - rhs) / abs(lhs)
                worst = max(worst, rel)
                print(f"n={n} T=x^{2*k} K0^{a} (xK1)^{b} t={ts}: rel {mp.nstr(rel, 3)}  ops {len(R)}")
    print(f"WORST {mp.nstr(worst, 3)}")
    return float(worst)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"
    if mode == "controls":
        mp.mp.dps = 40
        for n, ts_list in TS.items():
            for ts in ts_list:
                f40 = theta_F(n, 0, ts)
                mp.mp.dps = 60
                theta_F.cache_clear()
                f60 = theta_F(n, 0, ts)
                mp.mp.dps = 40
                theta_F.cache_clear()
                print(f"control quad n={n} t={ts}: F={mp.nstr(f60, 20)}  |F40-F60|/F={mp.nstr(abs(f40-f60)/abs(f60), 3)}")
        print("planted control (denominator a+2 for a+1):")
        run(30, plant="denom")
    else:
        run(int(sys.argv[2]) if len(sys.argv) > 2 else 40)
