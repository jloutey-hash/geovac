"""FD1-FD3 step O: evaluate operators (function-left, delta = w d/dw = theta_t/2, w = t^2) on S from
the cached theta_t^j S of step D.  Residual = |V S| / max_k |c_k(t^2) delta^k S|."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import mpmath as mp
import sympy as sp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_star_ops as O  # noqa: E402

DFILE = HERE / "data" / "four_density_star_D.json"


def thetas(tv: str, dps: int):
    d = json.loads(DFILE.read_text())
    return [mp.mpf(x) for x in d[f"{tv}@{dps}"]]


def evaluate(V, tv: str, dps: int):
    mp.mp.dps = dps
    th = thetas(tv, dps)
    if len(V) > len(th):
        raise ValueError("order exceeds cached theta powers")
    wv = mp.mpf(tv) ** 2
    terms = []
    for k, c in enumerate(V):
        if c == 0:
            continue
        num = c.numer.as_expr()
        den = c.denom.as_expr()
        cv = mp.mpf(sp.Rational(sp.Poly(num, O.w_sym).eval(sp.Rational(tv) ** 2))) / mp.mpf(sp.Rational(sp.Poly(den, O.w_sym).eval(sp.Rational(tv) ** 2)))
        terms.append(cv * th[k] / mp.mpf(2) ** k)
    tot = sum(terms)
    big = max(abs(x) for x in terms)
    return tot, big, abs(tot) / big
