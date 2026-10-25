"""PF2 Y clause, second arm: is each of the eight points genuine for Y's MINIMAL annihilator M?
Every analytic continuation of Y lies in M's solution space, so x0 is genuine for M as soon as one continuation of Y
is not invariant under a small loop around x0 (log or half-integer component).  Numerical continuation of the exact
order-19 operator (data/ib_tetra_pf_op_Q.json, form_D on S = Y/t^3 in x = t^2) by Taylor stepping at high precision,
from initial data S, S', ..., S^(18) at x_b = t_b^2 built from Y's 32-digit grid.
Sheets: 'principal' = Y itself (the convergent integral, analytic for Re t > 0, so the control: invariant at all
eight); 'around0' = Y continued once counterclockwise around x = 0 (t -> -t).  Paths stay in the upper half-plane.
Controls on the around0 sheet: loops around ordinary points 25/4, 25 and 3.
usage: python ib_tetra_pf_monodromy.py <grid json> <t_b> <dps>"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import mpmath as mp
from flint import fmpz_poly
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_8_onY as E

DATA = Path(__file__).parent / "data"
rec = json.loads(Path(sys.argv[1]).read_text())
tb = sys.argv[2]
mp.mp.dps = int(sys.argv[3])
d = json.loads((DATA / "ib_tetra_pf_op_Q.json").read_text())
polys = [[int(v) for v in p] for p in d["polys"]]
R = len(polys) - 1
lead = fmpz_poly(polys[-1])
roots = [mp.mpc(complex(r)) for r, _ in lead.complex_roots()]
print("order", R, "leading roots", len(roots), flush=True)


def shift(p, z):
    """coefficients of p(z + s) in s."""
    c = [mp.mpc(v) for v in p]
    n = len(c)
    for k in range(n):
        for j in range(n - 2, k - 1, -1):
            c[j] += z * c[j + 1]
    return c


def step(a, z, h, N):
    """a = [S^(k)(z)/k!], k < R; return the same at z + h from the Taylor series at z (N terms)."""
    P = [shift(p, z) for p in polys]
    c = list(a) + [mp.mpc(0)] * (N - R + 1)
    for m in range(0, N - R + 1):
        acc = mp.mpc(0)
        for k in range(R + 1):
            pk = P[k]
            for i in range(len(pk)):
                n = m - i + k
                if n < 0 or n >= m + R or (k == R and i == 0):
                    continue
                if n < k:
                    continue
                acc += pk[i] * c[n] * mp.ff(n, k)
        c[m + R] = -acc / (P[R][0] * mp.ff(m + R, R))
    out = []
    for k in range(R):
        out.append(mp.fsum(c[n] * mp.binomial(n, k) * h ** (n - k) for n in range(k, N + 1)))
    return out


def dist(z):
    return min([abs(z - r) for r in roots] + [abs(z)])


def go(a, z, w, frac=mp.mpf("0.35")):
    """continue along the straight segment z -> w, subdividing so each step is < frac * distance to a root."""
    while abs(w - z) > 0:
        h = w - z
        r = dist(z)
        if abs(h) > frac * r:
            h = h / abs(h) * frac * r
        N = R + 20 + int(mp.log(mp.mpf(10) ** (-mp.mp.dps)) / mp.log(abs(h) / r))
        a = step(a, z, h, N)
        z = z + h
    return a, z


def path(a, pts):
    z = pts[0]
    for w in pts[1:]:
        a, z = go(a, z, w)
    return a


def circle(c, r, start_angle, n=24):
    return [c + r * mp.expj(start_angle + 2 * mp.pi * k / n) for k in range(n + 1)]


from sympy.functions.combinatorial.numbers import stirling as _st


def init(tstr):
    """theta_t^k Y at t from the grid -> [S^(k)(x)/k!], x = t^2."""
    t = mp.mpf(tstr); x = t ** 2
    th = E.thetaY(rec, tstr, R - 1)
    # theta_x^j S = x^(-3/2) ((theta_t - 3)/2)^j Y
    thx = []
    for j in range(R):
        cf = [mp.binomial(j, i) * mp.mpf(-3) ** (j - i) / mp.mpf(2) ** j for i in range(j + 1)]
        thx.append(x ** mp.mpf(-1.5) * mp.fsum(cf[i] * th[i] for i in range(j + 1)))
    # x^k S^(k) = sum_j s(k, j) theta_x^j S (Stirling numbers of the first kind, signed)
    a = []
    for k in range(R):
        v = mp.fsum(mp.mpf(int(_st(k, j, kind=1, signed=True))) * thx[j] for j in range(k + 1))
        a.append(mp.mpc(v / x ** k / mp.factorial(k)))
    return a, x


a0, xs = init(tb)
print("start x =", mp.nstr(xs, 8), "S =", mp.nstr(a0[0], 20), flush=True)
H = mp.mpf("0.12")  # height of the upper-half-plane corridor
t0 = time.time()
# check the stepping and the conversion: principal Y from t_b to t = 1.7 through the upper half-plane
a_chk, xc = init("1.7")
a_cont = path(a0, [mp.mpc(xs), mp.mpc(xs, H), mp.mpc(xc, H), mp.mpc(xc)])
err = max(abs(a_cont[k] - a_chk[k]) * xc ** k for k in range(R)) / max(abs(a_chk[k]) * xc ** k for k in range(R))
print("check: continued vs grid at t = 1.7, scaled relative difference", mp.nstr(err, 4),
      f"({time.time()-t0:.0f}s)", flush=True)
# both sheets start at x_b = 0.12 on the real axis, reached through the upper half-plane (Y analytic there)
xb = mp.mpf("0.12")
a_p = path(a0, [mp.mpc(xs), mp.mpc(xs, H), mp.mpc(xb, H), mp.mpc(xb)])
# sheet 'around0': one counterclockwise loop around 0 at radius x_b (encloses no other singular point)
a_r = path(a_p, circle(mp.mpc(0), xb, 0))
print("around 0 done", f"{time.time()-t0:.0f}s; S changes by", mp.nstr(abs(a_r[0] - a_p[0]) / abs(a_p[0]), 5),
      flush=True)
sheets = {"principal": a_p, "around0": a_r}
targets = [(mp.mpf(1) / 4, mp.mpf("0.08")), (mp.mpf(1), mp.mpf("0.12")), (mp.mpf(9) / 4, mp.mpf("0.12")),
           (mp.mpf(4), mp.mpf("0.12")), (mp.mpf(9), mp.mpf("0.12")), (mp.mpf(49) / 4, mp.mpf("0.12")),
           (mp.mpf(16), mp.mpf("0.12")), (mp.mpf(81) / 4, mp.mpf("0.12")),
           (mp.mpf(3), mp.mpf("0.12")), (mp.mpf(25) / 4, mp.mpf("0.12")), (mp.mpf(25), mp.mpf("0.12"))]
res = {"check_t1.7": mp.nstr(err, 4)}
for name, a in sheets.items():
    # walk the corridor once, from x_b up to x_b + iH, then right; at each target branch off and loop
    z = mp.mpc(xb)
    a_cur, z = go(a, z, mp.mpc(xb, H))
    for x0, r in sorted(targets, key=lambda q: q[0]):
        a_cur, z = go(a_cur, z, mp.mpc(x0, H))
        # down to x0 + i r, loop around x0
        a_top, _ = go(a_cur, z, mp.mpc(x0, r))
        a_loop = path(a_top, circle(mp.mpc(x0), r, mp.pi / 2))
        scale = [r ** k for k in range(R)]
        num = max(abs(a_loop[k] - a_top[k]) * scale[k] for k in range(R))
        den = max(abs(a_top[k]) * scale[k] for k in range(R))
        rel = num / den
        res[f"{name} x0={mp.nstr(x0, 6)}"] = mp.nstr(rel, 4)
        print(f"{name:9s} x0 = {mp.nstr(x0, 6):>7s}: loop changes the scaled vector by {mp.nstr(rel, 4)}"
              f"  ({time.time()-t0:.0f}s)", flush=True)
out = {"grid": sys.argv[1], "t_start": tb, "dps": mp.mp.dps, "results": res}
(DATA / f"ib_tetra_pf_monodromy_tb{tb}_d{mp.mp.dps}.json").write_text(json.dumps(out, indent=1))
