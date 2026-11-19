"""exact9_review_numeric2.py -- the two parts of exact9_review_numeric.py that did not finish inside the cap
(N3 at an irrational ratio, N4 on K), with f0 by cumulative quadrature (one tail integral + short intervals).
"""
import json
import sys
import time
from fractions import Fraction

import mpmath as mp

sys.path.insert(0, "tests")
import _pachucki_extension as PE  # noqa: E402

T0 = time.time()
mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 50
ONLY_N3 = len(sys.argv) > 2
out = {}


def log(*a):
    print(*a, "[%.1fs]" % (time.time() - T0), flush=True)


def gfun(a, b):
    c = 1 / (2 * (a * a - b * b))
    LC = mp.log((a + b) / (b - a))
    g = {
        2: lambda R: mp.ei(-(a + b) * R),
        3: lambda R: mp.exp(-(a + b) * R),
        4: lambda R: mp.exp(-2 * a * R),
        5: lambda R: mp.exp(-2 * b * R),
        6: lambda R: mp.exp(2 * a * R) * mp.ei(-2 * (a + b) * R),
        7: lambda R: mp.exp(2 * b * R) * mp.ei(-2 * (a + b) * R),
        8: lambda R: mp.exp(-2 * a * R) * (mp.ei(2 * (a - b) * R) - 2 * mp.ei((a - b) * R) - LC),
        9: lambda R: mp.exp(-2 * b * R) * (mp.ei(2 * (b - a) * R) - 2 * mp.ei((b - a) * R) - LC),
    }
    fm1 = lambda R: c * (g[7](R) - g[6](R) - g[9](R) + g[8](R))  # noqa: E731
    return c, g, fm1


def f0_cumulative(fm1, pts):
    """f0 at the ascending points pts: tail from the last point, then short intervals downward."""
    Rm = pts[-1]
    acc = mp.quad(lambda s: fm1(s) / s, [Rm, Rm + 1, Rm + 2, Rm + 4, Rm + 10, Rm + 30, mp.inf])
    vals = [None] * len(pts)
    vals[-1] = acc
    for k in range(len(pts) - 2, -1, -1):
        acc += mp.quad(lambda s: fm1(s) / s, [pts[k], pts[k + 1]])
        vals[k] = acc
    return vals


def svals(rows_):
    A = mp.matrix(rows_)
    for j in range(A.cols):
        nrm = mp.sqrt(sum(A[i, j] ** 2 for i in range(A.rows)))
        for i in range(A.rows):
            A[i, j] /= nrm
    S = mp.svd_r(A, compute_uv=False)
    return sorted([S[i] for i in range(len(S))])


NP, DEG = 90, 3
pts = sorted(mp.mpf("0.15") + (mp.mpf(8) - mp.mpf("0.15")) * (1 - mp.cos(mp.pi * (k + mp.mpf(1) / 2) / NP)) / 2
             for k in range(NP))

# ---- N3 at b/a = sqrt 2
a, b = mp.mpf(1), mp.sqrt(2)
c, g, fm1 = gfun(a, b)
f0 = f0_cumulative(fm1, pts)
log("f0 sqrt2 done")
vals = [[f0[k]] + [g[i](R) for i in range(2, 10)] for k, R in enumerate(pts)]
rows_ = [[v * R ** d for v in vals[k] for d in range(DEG + 1)] for k, R in enumerate(pts)]
s9 = svals(rows_)
ctrl = [R ** 2 * vals[k][1] + 3 * vals[k][5] - R * vals[k][0] for k, R in enumerate(pts)]
s10 = svals([rows_[k] + [ctrl[k] * R ** d for d in range(DEG + 1)] for k, R in enumerate(pts)])
out["N3_sqrt2"] = {"smallest_3": [mp.nstr(x, 5) for x in s9[:3]], "largest": mp.nstr(s9[-1], 5),
                   "control_smallest_4": [mp.nstr(x, 5) for x in s10[:4]]}
log("N3 sqrt2", out["N3_sqrt2"])

if ONLY_N3:
    json.dump(out, open("debug/data/exact9_review_numeric2_dps%d.json" % mp.mp.dps, "w"), indent=1)
    sys.exit(0)
# ---- N4: K, ..., K^(8) at (1, 5/2), and the same with K^(9) added at degree 16 is NOT attempted (170 unknowns)
tq = Fraction(5, 2)
a, b = mp.mpf(1), mp.mpf(5) / 2
c, g, fm1 = gfun(a, b)
f0 = f0_cumulative(fm1, pts)
chk = abs(f0[10] - PE.f0_quadrature(Fraction(1), tq, pts[10]))
cq = PE.Ctx("Q", 1, tq)
calq, fS, Fs, P2, nfree = PE.solve_all(cq, PE.case_downset("1s"))
Kq = PE.assemble_K(cq, "1s", fS)
ders = calq.derivatives(Kq, 8)
rowsK = []
for k, R in enumerate(pts):
    dv = [PE.evaluate(Fraction(1), tq, d, R, f0=f0[k]) for d in ders]
    rowsK.append([v * R ** d for v in dv for d in range(DEG + 1)])
sK = svals(rowsK)
# derivative sanity: K' from the formal calculus vs a numerical derivative of the closed form
Rt = pts[40]
with mp.workdps(50):
    def Kval(x):
        return PE.evaluate(Fraction(1), tq, Kq, x, f0=PE.f0_quadrature(Fraction(1), tq, x))
    num = mp.diff(Kval, Rt, h=mp.mpf("1e-8"))
    form = PE.evaluate(Fraction(1), tq, ders[1], Rt, f0=f0[40])
Kn, _ = PE.neumann_integral("1s", 1, tq, Fraction(3, 2), dps=20)
Kc = PE.evaluate(Fraction(1), tq, Kq, mp.mpf(3) / 2, f0=PE.f0_quadrature(Fraction(1), tq, mp.mpf(3) / 2))
out["N4"] = {"smallest_3": [mp.nstr(x, 5) for x in sK[:3]], "largest": mp.nstr(sK[-1], 5),
             "f0_cumulative_vs_helper": mp.nstr(chk, 5),
             "Kprime_formal_minus_numeric": mp.nstr(form - num, 5),
             "K_closed_minus_neumann_R=3/2_t=5/2": mp.nstr(Kc - Kn, 5), "K": mp.nstr(Kc, 15)}
log("N4", out["N4"])
out["time_s"] = round(time.time() - T0, 1)
json.dump(out, open("debug/data/exact9_review_numeric2.json", "w"), indent=1)
log("done")
