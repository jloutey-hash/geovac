"""exact9_review_numeric.py -- adversarial numerics for debug/exact9_independence_memo.md (attacks 1, 2, 7).

N1  log coefficient of f0 at R = 0, reviewer's own derivation: f0 + Lam ln R must be smooth at 0, with
    Lam = 2c[Cinh(2bR) - Cinh(2aR)].  Second difference D(h) = [B(2h) - 2B(h) + B(0)]/h^2 must converge as h -> 0
    (a leftover eps R^2 ln R would make it drift by 2 eps ln 10 per decade).  Wrong candidates are scored too.
N2  the reviewer's own FM1 = c(g7 - g6 - g9 + g8) (from the memo's g's, mp.ei) equals the helper's label-60 Func.
N3  collocation: sum_i r_i(R) g_i(R) = 0 with deg r_i <= 3 at 90 points -> singular values (null space must be
    trivial); positive control with a tenth, dependent function (must show 2 null vectors).
    (a, b) = (1, 5/2) with the helper's f0 quadrature, and (1, sqrt 2) with the reviewer's own.
N4  collocation on K, ..., K^(8) (1s, (1, 5/2)), polynomial coefficients of degree <= 3: null space trivial.
"""
import json
import sys
import time
from fractions import Fraction

import mpmath as mp

sys.path.insert(0, "tests")
import _pachucki_extension as PE  # noqa: E402

T0 = time.time()
mp.mp.dps = 60
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


def cinh(x):
    s, term, k = mp.mpf(0), mp.mpf(1), 0
    while True:
        k += 1
        term = term * x * x / ((2 * k - 1) * (2 * k))
        add = term / (2 * k)
        s += add
        if abs(add) < mp.mpf(10) ** (-mp.mp.dps - 8):
            return s


# ---------------------------------------------------------------- N1, N2
n1 = {}
for (a, b, lab) in ((mp.mpf(1), mp.mpf(5) / 2, "1,5/2"), (mp.mpf(1), mp.sqrt(2), "1,sqrt2"), (mp.mpf(1), mp.mpf(3), "1,3")):
    c, g, fm1 = gfun(a, b)
    cands = {
        "memo: -2c[Cinh(2bR)-Cinh(2aR)]": lambda R: -2 * c * (cinh(2 * b * R) - cinh(2 * a * R)),
        "sign flipped": lambda R: 2 * c * (cinh(2 * b * R) - cinh(2 * a * R)),
        "zero": lambda R: mp.mpf(0),
        "half": lambda R: -c * (cinh(2 * b * R) - cinh(2 * a * R)),
    }
    res = {}
    for nm, A in cands.items():
        Ds = []
        for h in (mp.mpf("1e-2"), mp.mpf("1e-3"), mp.mpf("1e-4"), mp.mpf("1e-5")):
            # f0(2h) - 2 f0(h) + f0(0) = -int_h^{2h} FM1/s + int_0^h FM1/s   (f0 = int_R^oo FM1/s)
            i1 = mp.quad(lambda s: fm1(s) / s, [h, 2 * h])
            i0 = mp.quad(lambda s: fm1(s) / s, [0, h / 1000, h / 10, h])
            d2 = -i1 + i0 - (A(2 * h) * mp.log(2 * h) - 2 * A(h) * mp.log(h))      # B = f0 - A ln R, A(0) = 0
            Ds.append(d2 / h ** 2)
        res[nm] = {"D(h), h=1e-2..1e-5": [mp.nstr(x, 12) for x in Ds],
                   "drift_last_decade": mp.nstr(Ds[-1] - Ds[-2], 5)}
    n1[lab] = res
    log("N1", lab, {k: v["drift_last_decade"] for k, v in res.items()})
out["N1_f0_log_coefficient"] = n1

n2 = {}
for tq in (Fraction(5, 2), Fraction(3), Fraction(7, 5)):
    a, b = mp.mpf(1), mp.mpf(tq.numerator) / tq.denominator
    c, g, fm1 = gfun(a, b)
    hf = PE.Calc(PE.Ctx("Q", 1, tq)).FM1
    n2[str(tq)] = mp.nstr(max(abs(fm1(R) - PE.evaluate(Fraction(1), tq, hf, R)) for R in (mp.mpf("0.3"), mp.mpf("1.7"), mp.mpf(4))), 5)
out["N2_FM1_reviewer_vs_helper_maxdiff"] = n2
log("N2", n2)


# ---------------------------------------------------------------- N3
def svals(rows_):
    A = mp.matrix(rows_)
    for j in range(A.cols):                       # column scaling
        nrm = mp.sqrt(sum(A[i, j] ** 2 for i in range(A.rows)))
        for i in range(A.rows):
            A[i, j] /= nrm
    S = mp.svd_r(A, compute_uv=False)
    return sorted([S[i] for i in range(len(S))])


NP, DEG = 90, 3
pts = [mp.mpf("0.15") + (mp.mpf(8) - mp.mpf("0.15")) * (1 - mp.cos(mp.pi * (k + mp.mpf(1) / 2) / NP)) / 2 for k in range(NP)]
n3 = {}
for lab in ("1,5/2", "1,sqrt2"):
    if lab == "1,5/2":
        a, b = mp.mpf(1), mp.mpf(5) / 2
        c, g, fm1 = gfun(a, b)
        f0 = [PE.f0_quadrature(Fraction(1), Fraction(5, 2), R) for R in pts]
    else:
        a, b = mp.mpf(1), mp.sqrt(2)
        c, g, fm1 = gfun(a, b)
        f0 = [mp.quad(lambda s: fm1(s) / s, [R, R + 1, R + 2, R + 4, R + 10, R + 30, mp.inf]) for R in pts]
    log("N3 f0 values done", lab)
    vals = [[f0[k]] + [g[i](R) for i in range(2, 10)] for k, R in enumerate(pts)]
    rows_ = [[v * R ** d for v in vals[k] for d in range(DEG + 1)] for k, R in enumerate(pts)]
    s9 = svals(rows_)
    # positive control: tenth function in the span with coefficients of degree <= 3
    ctrl = [R ** 2 * vals[k][1] + 3 * vals[k][5] - R * vals[k][0] for k, R in enumerate(pts)]
    rows_c = [rows_[k] + [ctrl[k] * R ** d for d in range(DEG + 1)] for k, R in enumerate(pts)]
    s10 = svals(rows_c)
    n3[lab] = {"points": NP, "coeff_degree": DEG, "unknowns": 9 * (DEG + 1),
               "smallest_3_singular_values": [mp.nstr(x, 5) for x in s9[:3]], "largest": mp.nstr(s9[-1], 5),
               "control_smallest_4": [mp.nstr(x, 5) for x in s10[:4]]}
    log("N3", lab, n3[lab])
out["N3_collocation_nine_functions"] = n3

# ---------------------------------------------------------------- N4
tq = Fraction(5, 2)
cq = PE.Ctx("Q", 1, tq)
calq, fS, Fs, P2, nfree = PE.solve_all(cq, PE.case_downset("1s"))
Kq = PE.assemble_K(cq, "1s", fS)
ders = calq.derivatives(Kq, 8)
ptsK = pts[::1]
f0K = [PE.f0_quadrature(Fraction(1), tq, R) for R in ptsK]
rowsK = []
for k, R in enumerate(ptsK):
    dv = [PE.evaluate(Fraction(1), tq, d, R, f0=f0K[k]) for d in ders]
    rowsK.append([v * R ** d for v in dv for d in range(DEG + 1)])
sK = svals(rowsK)
# sanity: K itself against an independent Neumann evaluation at one point
Kn, _ = PE.neumann_integral("1s", 1, tq, Fraction(3, 2), dps=20)
Kc = PE.evaluate(Fraction(1), tq, Kq, mp.mpf(3) / 2, f0=PE.f0_quadrature(Fraction(1), tq, mp.mpf(3) / 2))
out["N4_collocation_K_order8_deg3"] = {"smallest_3": [mp.nstr(x, 5) for x in sK[:3]], "largest": mp.nstr(sK[-1], 5),
                                        "K_closed_minus_neumann_at_R=3/2": mp.nstr(Kc - Kn, 5), "K": mp.nstr(Kc, 15)}
log("N4", out["N4_collocation_K_order8_deg3"])
out["time_s"] = round(time.time() - T0, 1)
json.dump(out, open("debug/data/exact9_review_numeric.json", "w"), indent=1)
log("done")
