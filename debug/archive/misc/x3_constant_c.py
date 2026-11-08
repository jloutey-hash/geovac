"""Pre-registered X3-1c (docs/predictions/2026-10-22_x3_constant_c.md).

  python debug/x3_constant_c.py raw N   H(x) for n = N, three profiles, dps 120, tol 1e-125
  python debug/x3_constant_c.py fit     fit, controls, scoring; writes debug/data/x3_constant_c.json
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
DPS, TOL_EXP, NPTS = 120, 125, 30
NS = (1, 2, 3)
PROFS = ("point", "ball", "shell")
D = [(1, 0), (2, 1), (2, 0), (3, 1), (3, 0), (4, 2), (4, 1), (4, 0), (5, 2), (5, 1), (5, 0),
     (6, 2), (6, 1), (6, 0)]
SETS = {"D": D, "E": D + [(7, 1), (7, 0)], "F": D + [(3, 2), (6, 3)]}
ZEMACH = {"point": -mp.mpf(3) / 2, "ball": -mp.mpf(72) / 35, "shell": -mp.mpf(12) / 5}


def xs_grid():
    la, lb = mp.log(mp.mpf("3e-7")), mp.log(mp.mpf("3e-4"))
    return [mp.e ** (la + (lb - la) * i / (NPTS - 1)) for i in range(NPTS)]


def cols(spec):
    return [(lambda x, k=k, j=j: x**k * mp.log(x)**j) for k, j in spec]


def c_of(sol, spec):
    return sol[spec.index((3, 1))] - mp.mpf(4) / 5 * sol[spec.index((1, 0))]


def stage_raw(n):
    mp.mp.dps = DPS
    tol = mp.mpf(10) ** (-TOL_EXP)
    out = {"n": n, "dps": DPS, "xs": [], "H": {p: [] for p in PROFS}}
    t0 = time.time()
    for x in xs_grid():
        bs = zd.bound_state(n, x, tol)
        h = zd.hyperfine(n, x, bs)
        out["xs"].append(mp.nstr(x, DPS))
        for k in PROFS:
            out["H"][k].append(mp.nstr(h[k], DPS))
        print("raw n", n, "x", mp.nstr(x, 6), "t=%.0fs" % (time.time() - t0), flush=True)
    (DATA / ("x3_constant_c_raw_n%d.json" % n)).write_text(json.dumps(out, indent=1))


def stage_fit():
    mp.mp.dps = DPS
    mag = mp.mpf(64) / 63
    cv, xc, sols = {}, {}, {}
    for n in NS:
        raw = json.loads((DATA / ("x3_constant_c_raw_n%d.json" % n)).read_text())
        xs = [mp.mpf(s) for s in raw["xs"]]
        for p in PROFS:
            ys = [mp.mpf(s) for s in raw["H"][p]]
            for nm, spec in SETS.items():
                sol, r = zd.lsq(xs, ys, cols(spec))
                cv[(n, p, nm)] = c_of(sol, spec)
                xc[(n, p, nm)] = sol[0]
                sols[(n, p, nm)] = (xs, sol)
                print("fit", n, p, nm, "c+64/63", mp.nstr(cv[(n, p, nm)] + mag, 5), "resid", mp.nstr(r, 3),
                      flush=True)
    spread = max(abs(cv[(n, p, a)] - cv[(n, p, b)]) for n in NS for p in PROFS
                 for a in SETS for b in SETS) / mag
    gaps = {k: abs(v + mag) / mag for k, v in cv.items()}
    case_agree = {nm: max(abs(cv[a] - cv[b]) for a in cv for b in cv if a[2] == nm and b[2] == nm)
                  for nm in SETS}
    zem = max(abs(xc[k] - ZEMACH[k[1]]) / abs(ZEMACH[k[1]]) for k in xc)
    # synthetic control: case (ball, n=1, D), c shifted by 1e-6, noise 1e-100
    xs, sol = sols[(1, "ball", "D")]
    sol2 = list(sol)
    sol2[D.index((3, 1))] += mp.mpf("1e-6")
    f = cols(D)
    ys2 = [sum(s * g(x) for s, g in zip(sol2, f)) + mp.mpf("1e-100") * mp.sin(17 * i) for i, x in enumerate(xs)]
    s3, _ = zd.lsq(xs, ys2, f)
    syn_err = abs(c_of(s3, D) - (c_of(sol, D) + mp.mpf("1e-6")))
    passed = spread < mp.mpf("1e-9") and max(gaps.values()) < mp.mpf("1e-9")
    fals = spread < mp.mpf("1e-9") and max(gaps.values()) > mp.mpf("1e-8")
    verdict = "HELD" if passed else ("FALSIFIED" if fals else "UNCLEAR")
    res = {"dps": DPS, "window": "3e-7..3e-4, 30 log-spaced", "spread_rel": mp.nstr(spread, 6),
           "gap_rel": {"n%d_%s_%s" % k: mp.nstr(v, 6) for k, v in gaps.items()},
           "gap_rel_max": mp.nstr(max(gaps.values()), 6),
           "c": {"n%d_%s_%s" % k: mp.nstr(v, 30) for k, v in cv.items()},
           "case_agreement_abs": {k: mp.nstr(v, 6) for k, v in case_agree.items()},
           "ctrl_zemach_max_rel_err": mp.nstr(zem, 6), "ctrl_zemach_pass": bool(zem < mp.mpf("1e-12")),
           "ctrl_synthetic_err": mp.nstr(syn_err, 6), "ctrl_synthetic_pass": bool(syn_err < mp.mpf("1e-12")),
           "ctrl_cases_pass": bool(max(case_agree.values()) < mp.mpf("1e-10")), "verdict": verdict}
    (DATA / "x3_constant_c.json").write_text(json.dumps(res, indent=1))
    print(json.dumps({k: v for k, v in res.items() if k not in ("c", "gap_rel")}, indent=1))


if __name__ == "__main__":
    if sys.argv[1] == "raw":
        stage_raw(int(sys.argv[2]))
    else:
        stage_fit()
