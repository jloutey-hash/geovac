"""Cross-check of a second-pipeline run against the UNPATCHED geovac evaluator (this process never imports
hybrid_ratio_engine, so geovac.two_center_eri runs with its own normalised radial helpers).

The parametric atom form, specialised at u = u0, is compared with `hybrid_closed_form` at the corresponding
exponents at three R; the ratio must be one constant (the dropped normalisation).  Same method as
hybrid_ratio_check.py; reads debug/data/hybrid_ratio2_runs/, writes debug/data/hybrid_ratio2_checks/<run>.json.

Usage: python debug/hybrid_ratio2_check.py <case> <region> <u0> [<u0> ...]
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from geovac.two_center_eri import hybrid_closed_form, R_s  # noqa: E402

import os
mp.mp.dps = int(os.environ.get("HR2_DPS", "40"))
u = sp.Symbol("u", positive=True)
ORB = {"1s": (1, 0, 0), "2s": (2, 0, 0), "3s": (3, 0, 0), "2p0": (2, 1, 0), "2p1": (2, 1, 1),
       "2p-1": (2, 1, -1), "3p0": (3, 1, 0), "3p1": (3, 1, 1), "3d0": (3, 2, 0), "3d1": (3, 2, 1),
       "3d-1": (3, 2, -1), "3d2": (3, 2, 2), "3d-2": (3, 2, -2)}
RVALS = (sp.Rational(9, 10), sp.Rational(17, 10), sp.Rational(13, 5))


def atom_val(c, kind, kap, Rv):
    e = mp.e ** (c * Rv)
    if kind == "1":
        return e
    if kind == "log":
        return e * mp.log(Rv)
    return e * mp.ei(kap * Rv)


def const_val(consts, u0):
    v = sp.Integer(1)
    for n, e in consts:
        e = sp.Rational(e)
        if n == "PI":
            v *= sp.pi ** e
        elif n == "GAM":
            v *= sp.EulerGamma ** e
        elif n.startswith("P:"):
            v *= sp.Integer(n[2:]) ** e
        elif n.startswith("LOG:"):
            v *= sp.log(sp.sympify(n[4:], locals={"u": u}).subs(u, u0)) ** e
        else:
            raise ValueError(n)
    return v


def main():
    case, region = sys.argv[1], sys.argv[2]
    stem = case.replace(",", "_") + "__" + region
    run = json.loads((ROOT / "debug/data/hybrid_ratio2_runs" / (stem + ".json")).read_text())
    forms, res = run["form"], run["result"]
    info = {k: sp.sympify(v, locals={"u": u}) for k, v in res["info"].items()}
    Rp = sp.Symbol("R", positive=True)
    parsed = []
    for ks, cs in forms.items():
        c, kind, kap, consts = eval(ks)  # noqa: S307 -- our own json
        parsed.append((sp.sympify(c, locals={"u": u}), kind, sp.sympify(kap, locals={"u": u}), consts,
                       sp.sympify(cs, locals={"u": u, "R": Rp})))
    out = {}
    a, b, c_, d = (ORB[x] for x in case.split(","))
    for u0s in sys.argv[3:]:
        t0 = time.time()
        u0 = sp.Rational(u0s)
        ZA, ZB = info["ZA"].subs(u, u0), info["ZB"].subs(u, u0)
        ref = hybrid_closed_form(Fraction(int(ZA.p), int(ZA.q)), a, b, c_, Fraction(int(ZB.p), int(ZB.q)), d)
        tref = time.time() - t0
        ratios = []
        for Rv in RVALS:
            par = mp.mpf(0)
            for c, kind, kap, consts, cf in parsed:
                cv = mp.mpf(sp.N(c.subs(u, u0), mp.mp.dps + 10))
                kv = mp.mpf(sp.N(kap.subs(u, u0), mp.mp.dps + 10))
                val = sp.N(cf.subs({u: u0, Rp: Rv}) * const_val(consts, u0), mp.mp.dps + 10)
                par += mp.mpf(val) * atom_val(cv, kind, kv, mp.mpf(sp.N(Rv, mp.mp.dps + 10)))
            rv = mp.mpf(sp.N(ref.subs(R_s, Rv), mp.mp.dps + 10))
            ratios.append(rv / par)
        spread = max(abs(r / ratios[0] - 1) for r in ratios)
        out[u0s] = {"ZA": str(ZA), "ZB": str(ZB), "ratio": mp.nstr(ratios[0], 20), "rel_spread": mp.nstr(spread, 5),
                    "ref_build_s": round(tref, 1)}
        print(u0s, out[u0s], flush=True)
    od = ROOT / "debug/data/hybrid_ratio2_checks"
    od.mkdir(parents=True, exist_ok=True)
    (od / (stem + ".json")).write_text(json.dumps({"key": f"{case}|{region}", "points": out}, indent=1))


if __name__ == "__main__":
    main()
