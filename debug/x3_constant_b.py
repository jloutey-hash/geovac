"""Pre-registered X3-1b (docs/predictions/2026-10-21_x3_constant.md).

Stages (run as separate processes):
  python x3_constant_b.py raw N     compute H(x) for n = N, all three profiles, dps 120, tol 1e-125
  python x3_constant_b.py fit       fit, controls, scoring; writes debug/data/x3_constant_b.json
Reuses debug/zemach_door.py unmodified.
"""
import json
import sys
import time
from pathlib import Path

import mpmath as mp

sys.path.insert(0, str(Path(__file__).parent))
import zemach_door as zd  # noqa: E402

DATA = Path(__file__).parent / "data"
DPS = 120
TOL_EXP = 125
NPTS = 24
TARGET = -mp.mpf(64) / 63


def setup():
    mp.mp.dps = DPS


def xs_grid():
    a, b = mp.mpf("1e-6"), mp.mpf("1e-3")
    la, lb = mp.log(a), mp.log(b)
    return [mp.e ** (la + (lb - la) * i / (NPTS - 1)) for i in range(NPTS)]


def col(k, withlog, log2=False):
    if log2:
        return lambda x: x**k * mp.log(x)**2
    if withlog:
        return lambda x: x**k * mp.log(x)
    return lambda x: x**k


# Set A: x, x^2 ln x, x^2, x^3 ln x, x^3, x^4 ln x, x^4, x^5 ln x, x^5, x^6 ln x
COLS_A = [col(1, 0), col(2, 1), col(2, 0), col(3, 1), col(3, 0), col(4, 1), col(4, 0),
          col(5, 1), col(5, 0), col(6, 1)]
COLS_B = COLS_A[:-1]
COLS_C = COLS_A + [col(4, 0, True), col(5, 0, True)]
SETS = {"A": COLS_A, "B": COLS_B, "C": COLS_C}
# index of x coefficient = 0, x^3 ln x coefficient = 3 in all sets
ZEMACH = {"point": -mp.mpf(3) / 2, "ball": -mp.mpf(72) / 35, "shell": -mp.mpf(12) / 5}


def c_of(sol):
    return sol[3] - mp.mpf(4) / 5 * sol[0]


def stage_raw(n):
    setup()
    tol = mp.mpf(10) ** (-TOL_EXP)
    xs = xs_grid()
    out = {"n": n, "dps": DPS, "xs": [], "H": {"point": [], "ball": [], "shell": []}, "nu": []}
    t0 = time.time()
    for x in xs:
        bs = zd.bound_state(n, x, tol)
        h = zd.hyperfine(n, x, bs)
        out["xs"].append(mp.nstr(x, DPS))
        out["nu"].append(mp.nstr(bs["nu"], DPS))
        for k in h:
            out["H"][k].append(mp.nstr(h[k], DPS))
        print("raw n", n, "x", mp.nstr(x, 6), {k: mp.nstr(h[k], 15) for k in h},
              "t=%.0fs" % (time.time() - t0), flush=True)
    (DATA / ("x3_constant_b_raw_n%d.json" % n)).write_text(json.dumps(out, indent=1))


def stage_fit():
    setup()
    res = {"dps": DPS, "tol": "1e-%d" % TOL_EXP, "window": "1e-6..1e-3, %d log-spaced" % NPTS}
    cvals = {}  # (n,prof,set) -> c
    xcoef = {}
    sols = {}
    for n in (1, 2):
        raw = json.loads((DATA / ("x3_constant_b_raw_n%d.json" % n)).read_text())
        xs = [mp.mpf(s) for s in raw["xs"]]
        for prof in ("point", "ball", "shell"):
            ys = [mp.mpf(s) for s in raw["H"][prof]]
            for name, cols in SETS.items():
                sol, r = zd.lsq(xs, ys, cols)
                cvals[(n, prof, name)] = c_of(sol)
                xcoef[(n, prof, name)] = sol[0]
                sols[(n, prof, name)] = sol
                print("fit", n, prof, name, "c", mp.nstr(c_of(sol), 20), "xcoef", mp.nstr(sol[0], 15),
                      "resid", mp.nstr(r, 3), flush=True)
    res["c"] = {"n%d_%s_%s" % k: mp.nstr(v, 30) for k, v in cvals.items()}
    res["xcoef"] = {"n%d_%s_%s" % k: mp.nstr(v, 30) for k, v in xcoef.items()}
    # spread: largest difference of c among A,B,C over the six cases, relative to 64/63
    mag = mp.mpf(64) / 63
    spread = max(abs(cvals[(n, p, a)] - cvals[(n, p, b)])
                 for n in (1, 2) for p in ("point", "ball", "shell")
                 for a in "ABC" for b in "ABC") / mag
    # six-case agreement within each set
    case_spread = {nm: max(abs(cvals[k1] - cvals[k2]) for k1 in cvals for k2 in cvals
                           if k1[2] == nm and k2[2] == nm) for nm in "ABC"}
    gaps = {k: abs(cvals[k] + mag) / mag for k in cvals}
    maxgap = max(gaps.values())
    mingap = min(gaps.values())
    res["spread_rel"] = mp.nstr(spread, 10)
    res["gap_rel"] = {"n%d_%s_%s" % k: mp.nstr(v, 10) for k, v in gaps.items()}
    res["gap_rel_max"] = mp.nstr(maxgap, 10)
    res["gap_rel_min"] = mp.nstr(mingap, 10)
    res["six_case_agreement_abs"] = {nm: mp.nstr(v, 10) for nm, v in case_spread.items()}

    # control 1: x coefficients vs Zemach to 1e-10 (relative)
    zerr = max(abs(xcoef[k] / ZEMACH[k[1]] - 1) for k in xcoef)
    res["ctrl_zemach_max_rel_err"] = mp.nstr(zerr, 10)
    res["ctrl_zemach_pass"] = bool(zerr < mp.mpf("1e-10"))
    zerr_A = max(abs(xcoef[k] / ZEMACH[k[1]] - 1) for k in xcoef if k[2] == "A")
    res["ctrl_zemach_max_rel_err_setA"] = mp.nstr(zerr_A, 10)

    # control 3: six-case agreement to 1e-8 (absolute difference of c, per column set)
    res["ctrl_sixcase_pass"] = bool(max(case_spread.values()) < mp.mpf("1e-8"))

    # control 2: synthetic planted shift, case ball n=1, column set A
    raw = json.loads((DATA / "x3_constant_b_raw_n1.json").read_text())
    xs = [mp.mpf(s) for s in raw["xs"]]
    sol = list(sols[(1, "ball", "A")])
    sol[3] = sol[3] + mp.mpf("1e-4")
    c_planted = c_of(sol)
    import random
    rng = random.Random(12345)
    syn = []
    for x in xs:
        v = mp.fsum(sol[j] * COLS_A[j](x) for j in range(len(COLS_A)))
        noise = mp.mpf(rng.uniform(-1, 1)) * mp.mpf("1e-100")
        syn.append(v + noise)
    sol2, r2 = zd.lsq(xs, syn, COLS_A)
    c_rec = c_of(sol2)
    res["ctrl_synthetic"] = {"c_planted": mp.nstr(c_planted, 30), "c_recovered": mp.nstr(c_rec, 30),
                             "abs_err": mp.nstr(abs(c_rec - c_planted), 10),
                             "pass": bool(abs(c_rec - c_planted) < mp.mpf("1e-8"))}

    # scoring, literal
    pass_ = (spread < mp.mpf("1e-7")) and (maxgap < mp.mpf("1e-6"))
    fals = (spread < mp.mpf("1e-7")) and (maxgap > mp.mpf("1e-5"))
    # falsifier: "|c + 64/63| > 1e-5 x 64/63 in any case" -> any case/set exceeding
    res["verdict"] = "HELD" if pass_ else ("FAILED" if fals else "UNCLEAR")
    res["verdict_inputs"] = {"spread_lt_1e-7": bool(spread < mp.mpf("1e-7")),
                             "max_gap_lt_1e-6": bool(maxgap < mp.mpf("1e-6")),
                             "max_gap_gt_1e-5": bool(maxgap > mp.mpf("1e-5"))}
    print("spread", mp.nstr(spread, 6), "maxgap", mp.nstr(maxgap, 6), "mingap", mp.nstr(mingap, 6))
    print("ctrl zemach", res["ctrl_zemach_max_rel_err"], "synthetic", res["ctrl_synthetic"],
          "six-case", res["six_case_agreement_abs"])
    print("VERDICT", res["verdict"], res["verdict_inputs"])
    (DATA / "x3_constant_b.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    if sys.argv[1] == "raw":
        stage_raw(int(sys.argv[2]))
    else:
        stage_fit()
