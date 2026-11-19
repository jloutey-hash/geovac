"""Review P3: phi(a) T(H) = T(aH) for the stars of three and four in mpmath, phi applied under the integral sign.

Own fixed-node Gauss-Legendre panels; chi = int_0^x s K1 by cumulative panel quadrature (no Struve functions).
a and g are rebuilt with the memo's analyse(); P2 of tree_attain_review_main.py checks a H = g exactly, so what
this adds is the integration by parts on [0, oo) (the transfer, boundary terms at 0 and infinity) at n = 4.
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Dict

import mpmath as mp
import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "debug"))
T0 = time.time()
mp.mp.dps = 40
NG, NS = 36, 10


def log(*a) -> None:
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


def gl(n):
    # explicit Gauss-Legendre nodes by Newton on P_n
    nodes, weights = [], []
    for i in range(1, n + 1):
        x = mp.cos(mp.pi * (i - mp.mpf(1) / 4) / (n + mp.mpf(1) / 2))
        for _ in range(100):
            p0, p1 = mp.mpf(1), x
            for k in range(2, n + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            dp = n * (x * p1 - p0) / (x * x - 1)
            dx = p1 / dp
            x -= dx
            if abs(dx) < mp.mpf(10) ** (-(mp.mp.dps - 3)):
                break
        p0, p1 = mp.mpf(1), x
        for k in range(2, n + 1):
            p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
        dp = n * (x * p1 - p0) / (x * x - 1)
        nodes.append(x)
        weights.append(2 / ((1 - x * x) * dp * dp))
    return nodes, weights


def build_grid(xmax: int):
    gx, gw = gl(NG)
    sx, sw = gl(NS)
    edges = [mp.mpf(2) ** (-j) for j in range(24, 0, -1)] + [mp.mpf(1)]
    v = 1
    while v < xmax:
        v += 1 if v < 8 else 2
        edges.append(mp.mpf(v))

    def sk1(s):
        return s * mp.besselk(1, s)

    def sub(a, b):
        h, c = (b - a) / 2, (a + b) / 2
        return h * mp.fsum(w * sk1(c + h * x) for x, w in zip(sx, sw))

    # chi at the first edge: s K1 = 1 + O(s^2 log s) -> chi(e0) by quadrature on geometric panels below it
    chi_a = mp.mpf(0)
    lo = edges[0]
    for j in range(30):
        chi_a += sub(lo / 2, lo)
        lo /= 2
    chi_a += lo                                   # int_0^lo s K1 ~ lo, error lo^3 log lo ~ 1e-80
    X, W, CHI = [], [], []
    for a, b in zip(edges[:-1], edges[1:]):
        h, c = (b - a) / 2, (a + b) / 2
        prev, run = a, chi_a                      # nodes ascending: integrate node to node (NS-point rule each)
        for x, w in sorted(zip(gx, gw)):
            xi = c + h * x
            run += sub(prev, xi)
            prev = xi
            X.append(xi)
            W.append(h * w)
            CHI.append(run)
        chi_a = run + sub(prev, b)
    return X, W, CHI, chi_a


def theta_K_tables(mmax: int):
    y = sp.symbols("y")
    A, B = [sp.Integer(1)], [sp.Integer(0)]
    for _ in range(mmax):
        al, be = A[-1], B[-1]
        A.append(sp.expand(y * sp.diff(al, y) - y * be))
        B.append(sp.expand(y * sp.diff(be, y) - y * al - be))
    return y, A, B


def to_mp(c):
    c = sp.Rational(c)
    return mp.mpf(int(c.p)) / mp.mpf(int(c.q))


def identity(n: int, a: Dict, g: Dict, tq: Fraction, X, W, CHI) -> Dict:
    th = sp.symbols("theta")
    tS = sp.Rational(tq.numerator, tq.denominator)
    pol = sp.Poly(sp.expand(sum(sp.Integer(int(c)) * tS ** (-2 * k) * sp.prod([(th - 2 * i) ** 2 for i in range(k)]) * (-1 - th) ** j
                                for (k, j), c in a.items())), th)
    coeffs = pol.all_coeffs()[::-1]
    y, A, B = theta_K_tables(len(coeffs) - 1)
    alc = [to_mp(c) for c in sp.Poly(sp.expand(sum(c * A[m] for m, c in enumerate(coeffs))), y).all_coeffs()]
    bec = [to_mp(c) for c in sp.Poly(sp.expand(sum(c * B[m] for m, c in enumerate(coeffs))), y).all_coeffs()]
    t = mp.mpf(tq.numerator) / tq.denominator
    zz = sp.symbols("z")
    gp = {i: [mp.mpf(int(c)) for c in sp.Poly(p, zz).all_coeffs()] for i, p in g.items()}
    lhs_terms, rhs_terms = [], []
    for x, w, chi in zip(X, W, CHI):
        yv = t * x
        k0t, k1t = mp.besselk(0, yv), mp.besselk(1, yv)
        H = chi ** (n - 1) * x ** (3 - n)
        lhs_terms.append(w * (mp.polyval(alc, yv) * k0t + mp.polyval(bec, yv) * k1t) * H)
        k0v, D = mp.besselk(0, x), -x * mp.besselk(1, x)
        s = mp.fsum(mp.polyval(cs, x * x) * k0v ** (n - 1 - i) * D ** i for i, cs in gp.items())
        rhs_terms.append(w * k0t * s)
    lhs, rhs = mp.fsum(lhs_terms), mp.fsum(rhs_terms)
    res = {"n": n, "t": str(tq), "dps": mp.mp.dps, "lhs_phi_a_T_H": mp.nstr(lhs, 30), "rhs_T_aH": mp.nstr(rhs, 30),
           "rel_diff": mp.nstr(abs(lhs - rhs) / abs(rhs), 5),
           "largest_lhs_term": mp.nstr(max(abs(v) for v in lhs_terms), 5),
           "tail_last_panel": mp.nstr(abs(lhs_terms[-NG // 2]) / abs(lhs), 5)}
    log(f"n={n} t={tq}: lhs {res['lhs_phi_a_T_H']}  rhs {res['rhs_T_aH']}  rel {res['rel_diff']}  "
        f"largest term {res['largest_lhs_term']}  last-panel/result {res['tail_last_panel']}")
    return res


if __name__ == "__main__":
    X, W, CHI, chi_end = build_grid(200)
    log(f"grid: {len(X)} nodes; chi(200) - pi/2 = {mp.nstr(chi_end - mp.pi / 2, 5)}")
    # spot check of chi against the Struve closed form at one small-x node
    xi = X[len(X) // 5]
    ref = -xi * mp.besselk(0, xi) + (mp.pi * xi / 2) * (mp.besselk(0, xi) * mp.struvel(-1, xi) + mp.besselk(1, xi) * mp.struvel(0, xi))
    log(f"chi check at x = {mp.nstr(xi, 8)}: rel diff vs Struve form {mp.nstr(abs(CHI[len(X) // 5] / ref - 1), 5)}")
    import tree_attain_algebra as TA  # noqa: E402
    out = {"chi_end_minus_half_pi": mp.nstr(chi_end - mp.pi / 2, 5), "runs": []}
    for name in ("star3_root_centre", "star4"):
        n, edges = TA.TREES[name]
        res, a, g = TA.analyse(TA.Tree(n, edges), name)
        for tq in (Fraction(3, 2), Fraction(5, 2)):
            r = identity(n, a, g, tq, X, W, CHI)
            r["tree"] = name
            out["runs"].append(r)
            (ROOT / "debug" / "data" / "tree_attain_review_p3.json").write_text(json.dumps(out, indent=1))
    # wrong-answer control: drop the lowest-order term of a (an operator that does NOT send H to g)
    n, edges = TA.TREES["star4"]
    res, a, g = TA.analyse(TA.Tree(n, edges), "star4")
    kmin = min(a, key=lambda kj: (kj[1], kj[0]))
    a_bad = {k: v for k, v in a.items() if k != kmin}
    r = identity(n, a_bad, g, Fraction(3, 2), X, W, CHI)
    r["tree"] = "star4 CONTROL: a with one term dropped"
    out["runs"].append(r)
    (ROOT / "debug" / "data" / "tree_attain_review_p3.json").write_text(json.dumps(out, indent=1))
    log("done")
