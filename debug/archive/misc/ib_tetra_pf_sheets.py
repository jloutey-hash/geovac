"""PF2 Y clause, second arm, on chosen sheets: x0 is genuine for Y's minimal annihilator M if some analytic
continuation of Y changes on a small loop around x0 (every continuation of Y solves M).  Exact order-19 operator L
(data/ib_tetra_pf_op_Q.json, form_D on S = Y/t^3 in x = t^2), continued by Taylor stepping at high precision from
Y's values at t = t_s on a grid, through the upper half-plane (where Y, a convergent integral for Re t > 0, is
analytic) to the base point x_b = 0.12.
A sheet is a word of loops applied at the base point, e.g. "0,49/4": once around 0, then once around 49/4; each
loop around x0 != 0 goes along the corridor Im x = H to x0 + iH, down to x0 + i r, round the circle |x - x0| = r
counterclockwise, and back the same way.  The empty word "-" is Y itself (the control).
For each target the script reports the change of the scaled vector over the loop, relative to its size.
usage: python ib_tetra_pf_sheets.py <grid json> <t_s> <dps> <sheet word> <targets, comma-separated> <out tag>"""
from __future__ import annotations
import json, sys, time
from fractions import Fraction
from pathlib import Path
import mpmath as mp
from flint import fmpz_poly
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_8_onY as E
from sympy.functions.combinatorial.numbers import stirling as _st

DATA = Path(__file__).parent / "data"
grid, ts, DPS, word, targ, tag = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5], sys.argv[6]
mp.mp.dps = DPS
d = json.loads((DATA / "ib_tetra_pf_op_Q.json").read_text())
polys = [[int(v) for v in p] for p in d["polys"]]
R = len(polys) - 1
roots = [mp.mpc(complex(r)) for r, _ in fmpz_poly(polys[-1]).complex_roots()]
XB = mp.mpf("0.12")
H = mp.mpf("0.12")
_FF = {}


def ff(n, k):
    v = _FF.get((n, k))
    if v is None:
        v = _FF[(n, k)] = mp.ff(n, k)
    return v


def shift(p, z):
    c = [mp.mpc(v) for v in p]
    n = len(c)
    for k in range(n):
        for j in range(n - 2, k - 1, -1):
            c[j] += z * c[j + 1]
    return c


def step(a, z, h, N):
    P = [shift(p, z) for p in polys]
    lead0 = P[R][0]
    c = list(a) + [mp.mpc(0)] * (N - R + 1)
    for m in range(0, N - R + 1):
        terms = []
        for k in range(R + 1):
            pk = P[k]
            for i in range(len(pk)):
                n = m - i + k
                if n < k or n >= m + R:
                    continue
                terms.append(pk[i] * ff(n, k) * c[n])
        c[m + R] = -mp.fsum(terms) / (lead0 * ff(m + R, R))
    return [mp.fsum(c[n] * mp.binomial(n, k) * h ** (n - k) for n in range(k, N + 1)) for k in range(R)]


def dist(z):
    return min([abs(z - r) for r in roots] + [abs(z)])


def path(a, pts, frac=mp.mpf("0.3")):
    z = pts[0]
    for w in pts[1:]:
        while abs(w - z) > 0:
            h = w - z
            r = dist(z)
            if abs(h) > frac * r:
                h = h / abs(h) * frac * r
            N = R + 20 + int(mp.log(mp.mpf(10) ** (-DPS)) / mp.log(abs(h) / r))
            a = step(a, z, h, N)
            z = z + h
    return a


def circle(c, r, a0, n=24):
    return [c + r * mp.expj(a0 + 2 * mp.pi * k / n) for k in range(n + 1)]


def radius(x0):
    return mp.mpf("0.08") if x0 == mp.mpf(1) / 4 else H


def loop_at_base(a, x0):
    """apply one loop around x0, based at x_b, to the vector a at x_b."""
    if x0 == 0:
        return path(a, circle(mp.mpc(0), XB, 0))
    r = radius(x0)
    out = [mp.mpc(XB), mp.mpc(XB, H), mp.mpc(x0, H), mp.mpc(x0, r)]
    a = path(a, out)
    a = path(a, circle(mp.mpc(x0), r, mp.pi / 2))
    return path(a, out[::-1])


def init(tstr):
    rec = json.loads(Path(grid).read_text())
    t = mp.mpf(tstr); x = t ** 2
    th = E.thetaY(rec, tstr, R - 1)
    thx = []
    for j in range(R):
        cf = [mp.binomial(j, i) * mp.mpf(-3) ** (j - i) / mp.mpf(2) ** j for i in range(j + 1)]
        thx.append(x ** mp.mpf(-1.5) * mp.fsum(cf[i] * th[i] for i in range(j + 1)))
    a = [mp.mpc(mp.fsum(mp.mpf(int(_st(k, j, kind=1, signed=True))) * thx[j] for j in range(k + 1))
                / x ** k / mp.factorial(k)) for k in range(R)]
    return a, x


def q(s):
    f = Fraction(s)
    return mp.mpf(f.numerator) / f.denominator


def main():
    t0 = time.time()
    a, xs = init(ts)
    a = path(a, [mp.mpc(xs), mp.mpc(xs, H), mp.mpc(XB, H), mp.mpc(XB)])
    print(f"Y at x_b ({time.time()-t0:.0f}s)", flush=True)
    loops = [] if word == "-" else [q(s) for s in word.split(",")]
    for x0 in loops:
        a = loop_at_base(a, x0)
        print(f"applied loop around {mp.nstr(x0, 6)} ({time.time()-t0:.0f}s)", flush=True)
    res = {}
    z = mp.mpc(XB)
    cur = path(a, [z, mp.mpc(XB, H)])
    z = mp.mpc(XB, H)
    for s in sorted(targ.split(","), key=lambda s: Fraction(s)):
        x0 = q(s)
        r = radius(x0)
        cur = path(cur, [z, mp.mpc(x0, H)])
        z = mp.mpc(x0, H)
        top = path(cur, [z, mp.mpc(x0, r)])
        lp = path(top, circle(mp.mpc(x0), r, mp.pi / 2))
        sc = [r ** k for k in range(R)]
        rel = max(abs(lp[k] - top[k]) * sc[k] for k in range(R)) / max(abs(top[k]) * sc[k] for k in range(R))
        res[s] = mp.nstr(rel, 4)
        print(f"sheet [{word}] x0 = {s:>5s}: loop changes the scaled vector by {mp.nstr(rel, 4)} "
              f"({time.time()-t0:.0f}s)", flush=True)
    out = {"grid": grid, "t_s": ts, "dps": DPS, "sheet": word, "results": res}
    (DATA / f"ib_tetra_pf_sheets_{tag}.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
