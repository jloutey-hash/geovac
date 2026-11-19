"""Independent-of-the-patch check: the parametric atom form, specialised at u = u0, against the UNPATCHED
geovac evaluator at the corresponding exponents (normalised orbitals), at several R.  The ratio must be one
constant (the dropped normalisation, times pi powers).

Usage: python debug/hybrid_ratio_check.py <case> <region> <u0> [<u0> ...]
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from geovac.two_center_eri import hybrid_closed_form, R_s  # noqa: E402  (unpatched: engine not imported)

mp.mp.dps = 40
u = sp.Symbol("u", positive=True)
R = R_s
ORB = {"1s": (1, 0, 0), "2s": (2, 0, 0), "3s": (3, 0, 0), "2p0": (2, 1, 0), "2p1": (2, 1, 1),
       "2p-1": (2, 1, -1), "3p0": (3, 1, 0), "3p1": (3, 1, 1), "3d0": (3, 2, 0), "3d1": (3, 2, 1),
       "3d-1": (3, 2, -1), "3d2": (3, 2, 2), "3d-2": (3, 2, -2)}


def atom_val(lab, Rv):
    c, kind, kap = lab
    e = mp.e ** (c * Rv)
    if kind == "1":
        return e
    if kind == "log":
        return e * mp.log(Rv)
    return e * mp.ei(kap * Rv)


def main():
    case, region = sys.argv[1], sys.argv[2]
    key = f"{case}|{region}"
    run = json.loads((ROOT / "debug/data/hybrid_ratio_runs" /
                      (key.replace(",", "_").replace("|", "__") + ".json")).read_text())
    forms, res = run["form"], run["result"]
    info = {k: sp.sympify(v, locals={"u": u}) for k, v in res["info"].items()}
    Rp = sp.Symbol("R", positive=True)
    parsed = []
    for ks, cs in forms.items():
        c, kind, kap, consts = eval(ks)  # noqa: S307 -- our own json
        parsed.append((sp.sympify(c, locals={"u": u}), kind, sp.sympify(kap, locals={"u": u}), consts,
                       sp.sympify(cs, locals={"u": u, "R": Rp})))

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
    out = {}
    for u0s in sys.argv[3:]:
        u0 = sp.Rational(u0s)
        a, b, c, d = (ORB[x] for x in case.split(","))
        ZA, ZB = info["ZA"].subs(u, u0), info["ZB"].subs(u, u0)
        ref = hybrid_closed_form(Fraction(int(ZA.p), int(ZA.q)), a, b, c, Fraction(int(ZB.p), int(ZB.q)), d)
        ratios = []
        for Rv in (sp.Rational(9, 10), sp.Rational(17, 10), sp.Rational(13, 5)):
            par = mp.mpf(0)
            for c, kind, kap, consts, cf in parsed:
                lab = (mp.mpf(sp.N(c.subs(u, u0), 50)), kind, mp.mpf(sp.N(kap.subs(u, u0), 50)))
                val = sp.N(cf.subs({u: u0, Rp: Rv}) * const_val(consts, u0), 50)
                par += mp.mpf(val) * atom_val(lab, mp.mpf(sp.N(Rv, 50)))
            rv = mp.mpf(sp.N(ref.subs(R, Rv), 50))
            ratios.append(rv / par)
        spread = max(abs(r / ratios[0] - 1) for r in ratios)
        out[u0s] = {"ZA": str(ZA), "ZB": str(ZB), "ratio": mp.nstr(ratios[0], 20), "rel_spread": mp.nstr(spread, 5)}
        print(u0s, out[u0s], flush=True)
    p = ROOT / "debug/data/hybrid_ratio_checks.json"
    allc = json.loads(p.read_text()) if p.exists() else {}
    allc[key] = out
    p.write_text(json.dumps(allc, indent=1))


if __name__ == "__main__":
    main()
