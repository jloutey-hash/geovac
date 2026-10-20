"""Pre-registered EFT checks X4-1 and X3-1 (docs/predictions/2026-10-19_eft_checks.md).

X4-1: n-dependence of the s-wave quantum defect at x^4 (uniform sphere, Z = 1).
X3-1: the profile-independent constant of the hyperfine x^3 ln x coefficient.
Reuses the bound-state machinery of debug/zemach_door.py.
"""
import json
import sys
from pathlib import Path

import mpmath as mp

sys.path.insert(0, str(Path(__file__).parent))
import zemach_door as zd  # noqa: E402

DATA = Path(__file__).parent / "data"
DPS = 60


def logspace(a, b, k):
    la, lb = mp.log(a), mp.log(b)
    return [mp.e**(la + (lb - la) * i / (k - 1)) for i in range(k)]


def f(n):
    n = mp.mpf(n)
    return mp.digamma(n) - mp.log(n) + 1 / (2 * n)


def x4_pred(n, m):
    d0sq = mp.mpf(4) / 25
    r4 = mp.mpf(3) / 7
    return -d0sq * (f(n) - f(m)) - r4 / 15 * (1 / mp.mpf(n)**2 - 1 / mp.mpf(m)**2)


def stage_x4(out):
    tol = mp.mpf(10)**(-65)
    xs = logspace(mp.mpf("1e-4"), mp.mpf("3e-3"), 14)
    mus = {}
    for n in (1, 2, 3):
        mus[n] = []
        for x in xs:
            nu = zd.solve_nu(n, x, tol)
            mus[n].append(n - nu)
            print("x4", n, mp.nstr(x, 5), mp.nstr(n - nu, 20), flush=True)
    res = {"xs": [str(x) for x in xs], "mu": {n: [str(v) for v in mus[n]] for n in mus}}
    # control: mu(1) against the threshold series through x^3, at the x^4 level
    x0 = mp.mpf("1e-3")
    mu1 = 1 - zd.solve_nu(1, x0, tol)
    res["control_mu1_minus_series_x1e-3"] = str(mu1 - (-mp.mpf(2) / 5 * x0**2 + mp.mpf(32) / 63 * x0**3))
    cols5 = [lambda x: 1, lambda x: x * mp.log(x), lambda x: x, lambda x: x**2 * mp.log(x), lambda x: x**2]
    for (n, m) in ((1, 2), (1, 3)):
        ys = [(mus[n][i] - mus[m][i]) / xs[i]**4 for i in range(len(xs))]
        fits = {}
        for name, cols in (("5col", cols5), ("4col", cols5[:4]),
                           ("5col_plus_lnx", [lambda x: mp.log(x)] + cols5)):
            sol, r = zd.lsq(xs, ys, cols)
            fits[name] = [str(s) for s in sol]
        pred = x4_pred(n, m)
        c5 = mp.mpf(fits["5col"][0])
        c4 = mp.mpf(fits["4col"][0])
        lnx = mp.mpf(fits["5col_plus_lnx"][0])
        const_plus = mp.mpf(fits["5col_plus_lnx"][1])
        res["D%d%d" % (n, m)] = {
            "pred": str(pred), "fits": fits,
            "rel_err_5col": str(abs(c5 / pred - 1)), "rel_err_4col": str(abs(c4 / pred - 1)),
            "lnx_over_const": str(abs(lnx / const_plus)),
        }
        print("D%d%d" % (n, m), "pred", mp.nstr(pred, 12), "fit5", mp.nstr(c5, 12), "fit4", mp.nstr(c4, 12),
              "lnx/const", mp.nstr(abs(lnx / const_plus), 3), flush=True)
    out["X4"] = res


def stage_x3(out):
    tol = mp.mpf(10)**(-65)
    xs = logspace(mp.mpf("1e-5"), mp.mpf("2e-3"), 16)
    C = zd.COLS
    cols8 = [C["x"], C["x2lnx"], C["x2"], C["x3lnx"], C["x3"], C["x4lnx"], C["x4"],
             lambda x: x**5 * mp.log(x)]
    res = {"xs": [str(x) for x in xs]}
    target = -mp.mpf(64) / 63
    for n in (1, 2):
        H = {"point": [], "ball": [], "shell": []}
        for x in xs:
            bs = zd.bound_state(n, x, tol)
            h = zd.hyperfine(n, x, bs)
            for k in H:
                H[k].append(h[k])
            print("x3", n, mp.nstr(x, 5), {k: mp.nstr(h[k], 15) for k in h}, flush=True)
        for prof in H:
            entry = {"H": [str(v) for v in H[prof]]}
            for name, cols in (("8col", cols8), ("7col", cols8[:7])):
                sol, r = zd.lsq(xs, H[prof], cols)
                c = sol[3] - mp.mpf(4) / 5 * sol[0]
                entry[name] = {"coefs": [str(s) for s in sol], "c": str(c),
                               "rel_err": str(abs(c / target - 1))}
            res["n%d_%s" % (n, prof)] = entry
            print(n, prof, "x coef", mp.nstr(mp.mpf(entry["8col"]["coefs"][0]), 12),
                  "x2lnx", mp.nstr(mp.mpf(entry["8col"]["coefs"][1]), 10),
                  "c8", entry["8col"]["c"][:16], "c7", entry["7col"]["c"][:16],
                  "rel8", mp.nstr(mp.mpf(entry["8col"]["rel_err"]), 3),
                  "rel7", mp.nstr(mp.mpf(entry["7col"]["rel_err"]), 3), flush=True)
    out["X3"] = res


if __name__ == "__main__":
    mp.mp.dps = DPS
    stages = sys.argv[1:] or ["x4", "x3"]
    for st in stages:
        out = {"dps": DPS}
        {"x4": stage_x4, "x3": stage_x3}[st](out)
        (DATA / ("eft_checks_%s.json" % st)).write_text(json.dumps(out, indent=1, default=str))
