"""Review, attacks 4 and 5: from the STORED atom form of one run, independently of debug/hybrid_ratio_analyse.py,
  * extract the E-coefficients per shift, test proportionality, and compute gcd over Q[u] of q's R-coefficients;
  * build L (Wronskian operator of {q_- e^{-a_d R}, q_+ e^{+a_d R}}) by my own cofactor expansion, apply to f,
    and for each surviving rate compute the gcd over Q[u] of its numerator's coefficients (w.r.t. R and constants);
  * also list the positive real roots in u of every such gcd and of every u-factor of a denominator.
Arithmetic: sympy Poly / cancel on Q(u, R, consts) (no expand-cancel on the full form: per coefficient only).

Usage: python debug/hybrid_ratio_review_zero.py <case> <region>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
u = sp.Symbol("u", positive=True)
R = sp.Symbol("R", positive=True)


def main():
    case, region = sys.argv[1], sys.argv[2]
    key = f"{case}|{region}"
    run = json.loads((ROOT / "debug/data/hybrid_ratio_runs" /
                      (key.replace(",", "_").replace("|", "__") + ".json")).read_text())
    info = {k: sp.sympify(v, locals={"u": u}) for k, v in run["result"]["info"].items()}
    ad = info["a_d"]
    cs: dict = {}
    f: dict = {}
    for ks, val in run["form"].items():
        c, kind, kap, consts = eval(ks)  # noqa: S307
        m = sp.Integer(1)
        for n, e in consts:
            e = sp.Rational(e)
            if e.q != 1:          # common fractional factor (sqrt pi etc.) -- drop, but record
                cs.setdefault("_frac", set()).add((n, str(e)))
                continue
            if n not in cs:
                cs[n] = sp.Symbol(f"K{len([k for k in cs if k != '_frac'])}")
            m *= cs[n] ** int(e)
        lab = (sp.sympify(c, locals={"u": u}), kind, sp.sympify(kap, locals={"u": u}))
        f[lab] = f.get(lab, 0) + sp.sympify(val, locals={"u": u, "R": R}) * m
    f = {k: sp.cancel(v) for k, v in f.items()}
    f = {k: v for k, v in f.items() if v != 0}
    consts = [v for k, v in cs.items() if k != "_frac"]

    def D(d):
        o: dict = {}
        for (c, kind, kap), p in d.items():
            o[(c, kind, kap)] = o.get((c, kind, kap), 0) + sp.diff(p, R) + c * p
            if kind == "log":
                k2 = (c, "1", sp.Integer(0))
                o[k2] = o.get(k2, 0) + p / R
            elif kind == "Ei":
                k2 = (sp.expand(c + kap), "1", sp.Integer(0))
                o[k2] = o.get(k2, 0) + p / R
        o = {k: sp.cancel(v) for k, v in o.items()}
        return {k: v for k, v in o.items() if v != 0}

    E = {k: v for k, v in f.items() if k[1] != "1"}
    shifts = sorted({k[0] for k in E}, key=str)
    rep = {"shifts": {}}
    qs = {}
    for c in shifts:
        vals = [v for k, v in E.items() if k[0] == c]
        q0 = vals[0]
        ratios = [sp.cancel(v / q0) for v in vals]
        prop = all(not r.has(R) for r in ratios)
        num, den = sp.fraction(sp.cancel(q0))
        P = sp.Poly(num, R, *consts)
        g = 0
        for cf in P.coeffs():
            g = sp.gcd(g, cf) if g != 0 else cf
        qs[c] = q0
        rep["shifts"][str(c)] = {"n_atoms": len(vals), "ratios": [str(r) for r in ratios], "proportional": prop,
                                 "q_has_consts": bool(q0.free_symbols & set(consts)),
                                 "gcd_of_q_R_coeffs": str(sp.factor(g)),
                                 "den_u_factors": [str(fa) for fa, _ in sp.factor_list(den)[1] if fa.has(u)]}
    (cm, qm), (cp, qp) = sorted(qs.items(), key=lambda t: str(t[0]))[:2]
    rep["q_plus_over_q_minus_reflected"] = str(sp.cancel(qp / qm.subs(R, -R)))
    rep["shift_sum"] = str(sp.expand(cm + cp))

    def tw(q, c, j):
        for _ in range(j):
            q = sp.diff(q, R) + c * q
        return q
    g1 = [tw(qm, cm, j) for j in range(3)]
    g2 = [tw(qp, cp, j) for j in range(3)]
    # L y = det [[g1_0, g1_1, g1_2], [g2_0, g2_1, g2_2], [y, y', y'']] (exp(cm + cp) R) = 1 dropped)
    A0 = sp.cancel(g1[1] * g2[2] - g1[2] * g2[1])
    A1 = sp.cancel(-(g1[0] * g2[2] - g1[2] * g2[0]))
    A2 = sp.cancel(g1[0] * g2[1] - g1[1] * g2[0])
    f1 = D(f)
    f2 = D(f1)
    Lf: dict = {}
    for A, d in ((A0, f), (A1, f1), (A2, f2)):
        for k, v in d.items():
            Lf[k] = Lf.get(k, 0) + A * v
    Lf = {k: sp.cancel(v) for k, v in Lf.items()}
    Lf = {k: v for k, v in Lf.items() if v != 0}
    rep["Lf_E_leftover"] = [str(k) for k in Lf if k[1] != "1"]
    rep["Lf_rates"] = {}
    for k, v in Lf.items():
        num, den = sp.fraction(v)
        P = sp.Poly(num, R, *consts)
        g = 0
        for cf in P.coeffs():
            g = sp.gcd(g, cf) if g != 0 else cf
        gf = sp.factor(g)
        pr = [str(r) for r in sp.Poly(g, u).real_roots() if r > 0] if sp.sympify(g).has(u) else []
        rep["Lf_rates"][str(k[0])] = {"gcd": str(gf), "pos_roots": pr,
                                      "consts": sorted(str(s) for s in num.free_symbols & set(consts)),
                                      "den_u_factors": [str(fa) for fa, _ in sp.factor_list(den)[1] if fa.has(u)]}
    rep["k"] = len(rep["Lf_rates"])
    rep["order"] = 2 + rep["k"]
    print(json.dumps(rep, indent=1))
    p = ROOT / "debug/data/hybrid_ratio_review_zero.json"
    allc = json.loads(p.read_text()) if p.exists() else {}
    allc[key] = rep
    p.write_text(json.dumps(allc, indent=1))


if __name__ == "__main__":
    main()
