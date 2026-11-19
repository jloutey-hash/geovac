"""Review, attacks 2 and 7: the stored parametric atom forms at FRESH u (not among the memo's 27 points), against the
current (unpatched-by-the-engine) geovac `hybrid_closed_form`, at 4 R values incl. small and large, 100-digit atoms.
Also checks that the ratio ref/param equals the product of the four radial normalisations exactly (the dropped
factor), not just "a constant".

Usage: python debug/hybrid_ratio_review_fresh.py <case> <region> <u0> [...]
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
from geovac.two_center_eri import hybrid_closed_form, radial_norm, R_s  # noqa: E402

mp.mp.dps = 100
u = sp.Symbol("u", positive=True)
Rp = sp.Symbol("R", positive=True)
ORB = {"1s": (1, 0, 0), "2s": (2, 0, 0), "3s": (3, 0, 0), "2p0": (2, 1, 0), "2p1": (2, 1, 1),
       "3d0": (3, 2, 0), "3d1": (3, 2, 1)}
RVALS = [sp.Rational(1, 7), sp.Rational(6, 5), sp.Rational(7, 2), sp.Rational(11, 1)]


def main():
    case, region = sys.argv[1], sys.argv[2]
    key = f"{case}|{region}"
    run = json.loads((ROOT / "debug/data/hybrid_ratio_runs" /
                      (key.replace(",", "_").replace("|", "__") + ".json")).read_text())
    info = {k: sp.sympify(v, locals={"u": u}) for k, v in run["result"]["info"].items()}
    parsed = []
    for ks, cs in run["form"].items():
        c, kind, kap, consts = eval(ks)  # noqa: S307
        parsed.append((sp.sympify(c, locals={"u": u}), kind, sp.sympify(kap, locals={"u": u}), consts,
                       sp.sympify(cs, locals={"u": u, "R": Rp})))
    orbs = [ORB[x] for x in case.split(",")]
    out = {}
    for u0s in sys.argv[3:]:
        u0 = sp.Rational(u0s)
        ZA, ZB = sp.Rational(info["ZA"].subs(u, u0)), sp.Rational(info["ZB"].subs(u, u0))
        t0 = time.time()
        ref = hybrid_closed_form(Fraction(int(ZA.p), int(ZA.q)), *orbs[:3], Fraction(int(ZB.p), int(ZB.q)), orbs[3])
        tb = time.time() - t0
        norm = 1
        for o, Z in zip(orbs, (ZA, ZA, ZA, ZB)):
            norm *= radial_norm(Fraction(int(Z.p), int(Z.q)), o[0], o[1])
        normv = mp.mpf(sp.N(norm, 110))
        rels = []
        for Rv in RVALS:
            par = mp.mpf(0)
            for c, kind, kap, consts, cf in parsed:
                cv = sp.Rational(c.subs(u, u0))
                kv = sp.Rational(kap.subs(u, u0))
                coef = cf.subs({u: u0, Rp: Rv})
                assert coef.is_Rational, coef
                val = mp.mpf(int(coef.p)) / int(coef.q)
                for n, e in consts:
                    e = sp.Rational(e)
                    if n == "PI":
                        val *= mp.pi ** (mp.mpf(e.p) / e.q)
                    elif n == "GAM":
                        val *= mp.euler ** int(e)
                    elif n.startswith("P:"):
                        val *= mp.mpf(int(n[2:])) ** (mp.mpf(e.p) / e.q)
                    elif n.startswith("LOG:"):
                        arg = sp.Rational(sp.sympify(n[4:], locals={"u": u}).subs(u, u0))
                        assert arg > 0, (n, arg)
                        val *= mp.log(mp.mpf(arg.p) / arg.q) ** int(e)
                    else:
                        raise ValueError(n)
                x = mp.mpf(Rv.p) / Rv.q
                at = mp.e ** (mp.mpf(cv.p) / cv.q * x)
                if kind == "log":
                    at *= mp.log(x)
                elif kind == "Ei":
                    at *= mp.ei(mp.mpf(kv.p) / kv.q * x)
                par += val * at
            rv = mp.mpf(sp.N(ref.subs(R_s, Rv), 90))
            rels.append(rv / (par * normv) - 1)
        worst = max(abs(r) for r in rels)
        out[u0s] = {"ZA": str(ZA), "ZB": str(ZB), "build_s": round(tb, 1),
                    "max_rel_dev_from_norm_product": mp.nstr(worst, 5)}
        print(key, u0s, out[u0s], flush=True)
    p = ROOT / "debug/data/hybrid_ratio_review_fresh.json"
    allc = json.loads(p.read_text()) if p.exists() else {}
    allc[key] = out
    p.write_text(json.dumps(allc, indent=1))


if __name__ == "__main__":
    main()
