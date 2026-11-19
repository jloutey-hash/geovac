"""Collect debug/data/hybrid_ratio_runs/*.json into debug/data/hybrid_ratio_results.json and add the label checks
the order argument needs at every u > 0 of a region:

  * distinct exponential rates (f's, the closure's c + kappa, and L f's) never collide for u > 0;
  * E-labels with the same shift have distinct kappa for u > 0, and no kappa vanishes for u > 0 (a vanishing
    kappa is a hidden ln R coincidence);
  * the u-factors of every coefficient denominator, of every zero_factor (the u at which a rate of L f or an
    E-coefficient could vanish identically in R), have no root u > 0.

A root in (0, oo) of any of these would be an exceptional exponent ratio inside the open region; listed.
"""
from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
u = sp.Symbol("u", positive=True)


def pos_roots(expr):
    e = sp.sympify(expr, locals={"u": u})
    if not e.has(u):
        return [] if e != 0 else ["identically 0"]
    num = sp.numer(sp.together(e))
    rr = [r for r in sp.Poly(num, u).real_roots() if r > 0]
    return [str(r) for r in rr]


def check(res, form):
    labs = [eval(k) for k in form]  # noqa: S307 -- our own json
    exp_rates = {lab[0] for lab in labs if lab[1] == "1"}
    E = [(lab[0], lab[1], lab[2]) for lab in labs if lab[1] != "1"]
    closure = {str(sp.expand(sp.sympify(c, locals={"u": u}) + sp.sympify(k, locals={"u": u})))
               for c, kind, k in E if kind == "Ei"} | {c for c, kind, _k in E if kind == "log"}
    lf = set(res.get("Lf_rates", {}))
    rates = sorted(exp_rates | closure | lf)
    issues = []
    for i in range(len(rates)):
        for j in range(i + 1, len(rates)):
            rt = pos_roots(f"({rates[i]}) - ({rates[j]})")
            if rt:
                issues.append(f"rates {rates[i]} = {rates[j]} at u = {rt}")
    shifts = sorted({c for c, _k, _x in E})
    for i in range(len(shifts)):
        for j in range(i + 1, len(shifts)):
            rt = pos_roots(f"({shifts[i]}) - ({shifts[j]})")
            if rt:
                issues.append(f"shifts collide at u = {rt}")
    for c in shifts:
        ks = sorted({k for cc, kind, k in E if cc == c and kind == "Ei"})
        for k in ks:
            if pos_roots(k):
                issues.append(f"kappa {k} vanishes at u = {pos_roots(k)}")
        for i in range(len(ks)):
            for j in range(i + 1, len(ks)):
                rt = pos_roots(f"({ks[i]}) - ({ks[j]})")
                if rt:
                    issues.append(f"kappa {ks[i]} = {ks[j]} at u = {rt}")
    for fac in res.get("coeff_u_denominator_factors", []):
        if pos_roots(fac):
            issues.append(f"coefficient denominator {fac} vanishes at u = {pos_roots(fac)}")
    for c, g in res.get("shifts", {}).items():
        if pos_roots(g["q_first_zero_factor"]):
            issues.append(f"E-coefficient at shift {c} vanishes identically at u = {pos_roots(g['q_first_zero_factor'])}")
    for c, g in res.get("Lf_rates", {}).items():
        if pos_roots(g["zero_factor"]):
            issues.append(f"L f rate {c} coefficient vanishes identically at u = {pos_roots(g['zero_factor'])}")
        if g["consts_in_coeff"]:
            issues.append(f"L f rate {c} coefficient carries constants {g['consts_in_coeff']}")
    return {"rates_all": rates, "n_rates_f": len({r for r in exp_rates} | {c for c, kind, _k in E if kind == 'Ei'}),
            "issues": issues}


def main():
    out = {}
    for fn in sorted((ROOT / "debug/data/hybrid_ratio_runs").glob("*.json")):
        run = json.loads(fn.read_text())
        res, form = run["result"], run["form"]
        chk = check(res, form)
        shifts_c = list(res.get("shifts", {}))
        nrates = len(set(res.get("exp_rates_of_f", [])) | set(shifts_c))
        out[run["key"]] = {"order": res.get("order"), "d_E": res.get("d_E"), "k": res.get("k"),
                           "rates_of_f": sorted(set(res.get("exp_rates_of_f", [])) | set(shifts_c)),
                           "n_rates": nrates, "order_minus_rates": (res.get("order") or 0) - nrates,
                           "Lf_rates": list(res.get("Lf_rates", {})), "constants": res.get("constants"),
                           "issues": chk["issues"], "seconds": res.get("seconds"), "info": res.get("info")}
        print(run["key"], out[run["key"]]["order"], "rates", nrates, "dE", res.get("d_E"), "k", res.get("k"),
              "issues", chk["issues"])
    (ROOT / "debug/data/hybrid_ratio_results.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
