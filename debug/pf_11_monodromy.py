"""PF1, PF3, PF4: monodromy of the exact order-19 operator L of Y/t^3 (data/ib_tetra_pf_op_Q.json, form_D in x = t^2),
independent of the probes of pf_7_chain.py.  Each of the 19 unit initial vectors (Taylor data at the base point
B = x_b + iH, x_b = H = 0.12) is carried by the Taylor stepping of ib_tetra_pf_sheets.py (copied) once around x = 0
(from B down to x_b, the circle |x| = x_b, back up) and, along the corridor Im x = H, to each target x0: down to x0 + ir,
round the circle |x - x0| = r counterclockwise.  Saved: the vectors before and after each circle; the analysis is
pf_11b_monodromy_analyse.py.  Targets: the eight cuts and the ordinary points 25/4 and 10 (controls).
usage: python pf_11_monodromy.py <dps> <out tag>"""
from __future__ import annotations
import json, sys, time
from fractions import Fraction
from pathlib import Path
import mpmath as mp
from flint import fmpz_poly

DATA = Path(__file__).parent / "data"
_d = json.loads((DATA / "ib_tetra_pf_op_Q.json").read_text())
POLYS = [[int(v) for v in p] for p in _d["polys"]]
R = len(POLYS) - 1
TARGETS = ["1/4", "1", "9/4", "4", "25/4", "9", "10", "49/4", "16", "81/4"]
ST = {}


def setup(dps):
    mp.mp.dps = dps
    ST["roots"] = [mp.mpc(complex(r)) for r, _ in fmpz_poly(POLYS[-1]).complex_roots()]
    ST["ff"] = {}
    ST["XB"] = mp.mpf("0.12")
    ST["H"] = mp.mpf("0.12")


def ff(n, k):
    v = ST["ff"].get((n, k))
    if v is None:
        v = ST["ff"][(n, k)] = mp.ff(n, k)
    return v


def shift(p, z):
    c = [mp.mpc(v) for v in p]
    n = len(c)
    for k in range(n):
        for j in range(n - 2, k - 1, -1):
            c[j] += z * c[j + 1]
    return c


def step(a, z, h, N):
    P = [shift(p, z) for p in POLYS]
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
    return min([abs(z - r) for r in ST["roots"]] + [abs(z)])


def path(a, pts, frac=mp.mpf("0.3")):
    z = pts[0]
    for w in pts[1:]:
        while abs(w - z) > 0:
            h = w - z
            r = dist(z)
            if abs(h) > frac * r:
                h = h / abs(h) * frac * r
            N = R + 20 + int(mp.log(mp.mpf(10) ** (-mp.mp.dps)) / mp.log(abs(h) / r))
            a = step(a, z, h, N)
            z = z + h
    return a


def circle(c, r, a0, n=24):
    return [c + r * mp.expj(a0 + 2 * mp.pi * k / n) for k in range(n + 1)]


def q(s):
    f = Fraction(s)
    return mp.mpf(f.numerator) / f.denominator


def radius(x0):
    return mp.mpf("0.08") if x0 == mp.mpf(1) / 4 else ST["H"]


def enc(v):
    return [[str(z.real), str(z.imag)] for z in v]


def column(args):
    i, dps = args
    setup(dps)
    XB, H = ST["XB"], ST["H"]
    t0 = time.time()
    B = mp.mpc(XB, H)
    e = [mp.mpc(0)] * R
    e[i] = mp.mpc(1)
    out = {}
    a = path(e, [B, mp.mpc(XB)])
    a = path(a, circle(mp.mpc(0), XB, 0))
    a = path(a, [mp.mpc(XB), B])
    out["0"] = enc(a)
    cur, z = e, B
    for s in TARGETS:
        x0 = q(s)
        r = radius(x0)
        cur = path(cur, [z, mp.mpc(x0, H)])
        z = mp.mpc(x0, H)
        top = path(cur, [z, mp.mpc(x0, r)])
        lp = path(top, circle(mp.mpc(x0), r, mp.pi / 2))
        out[s] = {"top": enc(top), "after": enc(lp)}
    print(f"column {i} done ({time.time()-t0:.0f}s)", flush=True)
    return i, out


def main():
    dps, tag = int(sys.argv[1]), sys.argv[2]
    from multiprocessing import Pool
    t0 = time.time()
    with Pool(R) as pool:
        res = pool.map(column, [(i, dps) for i in range(R)], chunksize=1)
    cols = {str(i): o for i, o in res}
    rec = {"dps": dps, "base": "0.12+0.12i", "targets": TARGETS, "order": R, "columns": cols,
           "seconds": time.time() - t0}
    (DATA / f"pf_11_monodromy_{tag}.json").write_text(json.dumps(rec))
    print(f"saved ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
