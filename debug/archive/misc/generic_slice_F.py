"""GS4 + GS5 (record docs/predictions/2026-10-21_generic_slice.md).
Local exponents (indicial roots) of R_V and W7 = R o E2 o (theta+2) at t = 0, 1, 2, 3 and infinity, the
pole/degree bound they set for a splitting gauge (leading-term argument: ord_t0 b_j >= min exp(W7) -
max exp(R_V), numerator degree from infinity), then the mod-p split search: does Sol(R_V) map into
Sol(W7) itself?  Positive control: E3 o R_V (split by construction) at the same box."""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import generic_slice_C as GC  # noqa: E402
import density_loop_D4_factors as F4  # noqa: E402
from generic_slice_B import compose, coeff_list, cth, t, u  # noqa: E402

HERE = Path(__file__).parent
rho = sp.Symbol("rho")


def indicial(Lc, t0):
    """Indicial polynomial of sum Lc[k] theta^k at t = t0 (t0 finite) or t0 = 'oo'."""
    if t0 == "oo":
        s = sp.Symbol("s")
        f = t ** (-rho)
    else:
        f = (t - t0) ** rho
    acc, g = 0, f
    for c in Lc:
        acc += sp.cancel(c) * g
        g = sp.expand(t * sp.diff(g, t))
    if t0 == "oo":
        e = sp.cancel(acc / t ** (-rho))
        e = sp.cancel(e.subs(t, 1 / sp.Symbol("w")))
        w = sp.Symbol("w")
        num, den = sp.fraction(sp.together(e))
        nl = sp.Poly(num, w).monoms()[-1][0] if True else 0
        dl = sp.Poly(den, w).monoms()[-1][0]
        pn = sp.Poly(num, w)
        lowest = min(m[0] for m in pn.monoms())
        return sp.factor(pn.coeff_monomial(w**lowest) if lowest else pn.as_expr().subs(w, 0))
    e = sp.cancel(acc / (t - t0) ** rho)
    x = sp.Symbol("x")
    e = sp.cancel(e.subs(t, t0 + x))
    num, den = sp.fraction(sp.together(e))
    pn = sp.Poly(sp.expand(num), x)
    lowest = min(m[0] for m in pn.monoms())
    return sp.factor(sum(c for m, c in zip(pn.monoms(), pn.coeffs()) if m[0] == lowest))


def main():
    rec = {}
    Dd = json.loads((HERE / "data" / "generic_slice_D2.json").read_text())
    Cd = json.loads((HERE / "data" / "generic_slice_C.json").read_text())
    R = [sp.sympify(c, locals={"t": t}) for c in Dd["R"]]
    E2 = [sp.sympify(c, locals={"t": t}) for c in Cd["E2"]]
    E3 = compose(E2, [sp.Integer(2), sp.Integer(1)])
    W7 = compose(R, E3)
    RV = coeff_list(GC.QV)
    for name, Lc in [("R_V", RV), ("R", R), ("W7", W7)]:
        for t0 in [0, 1, 2, 3, "oo"]:
            ip = indicial(Lc, t0)
            roots = sp.roots(sp.Poly(ip, rho)) if ip.has(rho) else {}
            rec[f"{name}|{t0}"] = {str(k): v for k, v in roots.items()}
            print(f"{name:4s} t={t0}: {sp.factor(ip)}  roots {dict(roots)}", flush=True)
    (HERE / "data" / "generic_slice_F.json").write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
