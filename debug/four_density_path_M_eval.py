"""FP step M: residuals of operators on P from the cached theta_t^j P (step E).  Residual = |V P| / max_k |c_k delta^k P|
(V function-left in delta = theta/2, w = t^2)."""
from __future__ import annotations
import json, sys
from pathlib import Path
import mpmath as mp
import sympy as sp
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_star_ops as O  # noqa: E402
u = O.u
DFILE = HERE / "data" / "four_density_path_E.json"


def thetas(tv, dps):
    d = json.loads(DFILE.read_text())
    return [mp.mpf(x) for x in d[f"{tv}@{dps}"]]


def avail(dps):
    d = json.loads(DFILE.read_text())
    return [k.split("@")[0] for k in d if k.endswith(f"@{dps}")]


def evaluate(V, tv, dps):
    mp.mp.dps = dps
    th = thetas(tv, dps)
    assert len(V) <= len(th)
    wv = sp.Rational(tv) ** 2
    terms = []
    for k, c in enumerate(V):
        if c == 0:
            continue
        cv = sp.Rational(c.as_expr().subs(O.w_sym, wv))
        terms.append(mp.mpf(cv.p) / cv.q * th[k] / mp.mpf(2) ** k)
    tot = sum(terms)
    big = max(abs(x) for x in terms)
    return tot, big, abs(tot) / big


def loadQ(path, key, sub=None):
    d = json.loads((HERE / "data" / path).read_text())
    d = d[key] if sub is None else d[key][sub]
    return O.from_Q([sp.sympify(q).subs(sp.Symbol("u"), u) for q in d])


if __name__ == "__main__":
    dps = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    ops = {"W12": loadQ("four_density_path_D.json", "W12", "Q"),
           "theta-2": O.from_Q([u - 2]),
           "L": loadQ("four_density_path_A.json", "L", "Q")}
    for tv in avail(dps):
        print(tv, {k: mp.nstr(evaluate(v, tv, dps)[2], 4) for k, v in ops.items()}, flush=True)
