"""GY: the Newton polygon at x = infinity of an operator L = sum_k c_k(x) theta^k (theta = theta_t = 2 x d/dx, x = t^2;
library convention of correlated_triangle_g32_lib), from the degrees of its coefficients.

Points (k, d_k), d_k = deg_x c_k.  On e^{c x^sigma} the term k contributes x^{d_k + k sigma} (theta x^s ~ s x^s and
theta e^{c x^sigma} ~ 2 c sigma x^sigma e^{..}), so formal solutions of exponential type x^sigma, sigma > 0, come from the edges
of the UPPER convex hull of the points with slope -sigma; the operator is regular singular (Fuchsian) at infinity iff
max_k d_k is attained at k = order (every upper-hull slope >= 0; height := max d_k - d_order = 0).  For an edge of slope -1
(e^{c x}) the characteristic polynomial in C = 2c is  sum_{k on the edge} lc(c_k) C^k = 0  (theta e^{cx} = 2 c x e^{cx}),
whose roots give the exponents c = C/2 of e^{c x} = e^{c t^2}.  At x = 0 the same construction on the valuations (lower hull)
reads the regularity there.
Controls (run with 'controls'): the Slater L9 at gamma = 2 (debug/data/cf_g2_5_op_g2_L9.json) must be Fuchsian at infinity
(height 0); F_3's operator 4 th^3 - x (th+1)^3 Fuchsian; GA2's avatar operator L_G = x^2 (th+1)^2 - x th^2 (4 th + 1)
+ 4 th^2 (th-2)^2 (in theta_t, x = t^2) must give slope -1 with c = 1/4 (double) [GA2: J ~ e^{t^2/4}].
usage: python gy_newton.py controls | python gy_newton.py <op json> <key|Q|theta_form> <out tag>"""
from __future__ import annotations
import json, sys
from fractions import Fraction
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G

DATA = Path(__file__).parent / "data"
u = G.u_sym


def upper_hull(pts):
    """upper convex hull of points (k, d) sorted by k (Andrew's monotone chain, upper part)"""
    pts = sorted(pts)
    hull = []
    for p_ in pts:
        while len(hull) >= 2:
            (x1, y1), (x2, y2) = hull[-2], hull[-1]
            # remove hull[-1] if it is below or on the segment hull[-2] -> p_
            if (x2 - x1) * (p_[1] - y1) - (y2 - y1) * (p_[0] - x1) >= 0:
                hull.pop()
            else:
                break
        hull.append(p_)
    return hull


def newton_infinity(L):
    r = len(L) - 1
    pts = [(k, L[k].degree()) for k in range(r + 1) if not L[k].is_zero()]
    hull = upper_hull(pts)
    dmax = max(d for _, d in pts)
    height = dmax - L[r].degree()
    edges = []
    for (k1, d1), (k2, d2) in zip(hull, hull[1:]):
        slope = Fraction(d2 - d1, k2 - k1)
        edge = {"from": (k1, d1), "to": (k2, d2), "slope": str(slope), "sigma": str(-slope)}
        if slope < 0:
            sigma = -slope
            # characteristic polynomial: points on this edge: d_k = d1 + slope (k - k1)
            on = [k for k, d in pts if d == d1 + slope * (k - k1)]
            C = sp.Symbol("C")
            poly = sum(int(L[k].coeffs()[-1]) * C ** k for k in on)
            poly = sp.Poly(poly, C)
            # roots: exponent c of e^{c x^sigma} with theta e^{c x^sigma} ~ 2 c sigma x^sigma
            roots = sp.roots(poly.as_expr(), C)
            edge["points_on_edge"] = on
            edge["char_poly_in_C"] = str(sp.factor(poly.as_expr()))
            edge["exponents_c_of_exp_c_x_sigma"] = {str(sp.nsimplify(rt / (2 * sigma))): m for rt, m in roots.items()}
            edge["note"] = "a root C = 0 belongs to the lower-slope part, not to an exponential; the roots are exact if the polynomial factors over Q, else numerical"
            if sum(roots.values()) < poly.degree():
                edge["unresolved_factor_degree"] = poly.degree() - sum(roots.values())
                edge["numeric_roots_C"] = [str(sp.N(rt, 15)) for rt in sp.Poly(poly, C).nroots()]
            if sigma == 1:
                edge["reading"] = "formal solutions exp(c x) = exp(c t^2) with c = C/2"
        edges.append(edge)
    return {"order": r, "points": pts, "upper_hull": hull, "height": height, "fuchsian_at_infinity": height == 0,
            "edges": edges, "growth_exponents_if_fuchsian": [str(e) for e in G.exps_of(G.indicial_inf(L))[0]] if height == 0 else None}


def newton_zero(L):
    r = len(L) - 1
    pts = [(k, G.valuation(L[k])) for k in range(r + 1) if not L[k].is_zero()]
    vmin = min(v for _, v in pts)
    return {"points": pts, "regular_at_0": L[r].is_zero() is False and G.valuation(L[r]) == vmin}


def report(L, name):
    inf = newton_infinity(L)
    z0 = newton_zero(L)
    print(f"== {name}: order {inf['order']}, points (k, deg c_k) = {inf['points']}", flush=True)
    print(f"   upper hull {inf['upper_hull']}, height {inf['height']}, Fuchsian at infinity: {inf['fuchsian_at_infinity']}", flush=True)
    for e in inf["edges"]:
        print("   edge", e, flush=True)
    print(f"   at x = 0: regular {z0['regular_at_0']}", flush=True)
    return {"infinity": inf, "zero": z0}


def load(path, key):
    d = json.loads(Path(path).read_text())
    if key == "Q":
        return G.from_Q([sp.sympify(q, locals={"u": u}) for q in d["Q"]])
    if key == "theta_form":
        return G.parse_op(d["theta_form"])
    return G.parse_op(d[key])


if __name__ == "__main__":
    if sys.argv[1] == "controls":
        out = {}
        L9 = load(DATA / "cf_g2_5_op_g2_L9.json", "Q")
        out["L9_slater_g2"] = report(L9, "Slater L9 (gamma = 2, record CF)")
        out["F3"] = report(G.L3_op(), "F_3 operator")
        LG = G.from_Q([4 * u**2 * (u - 2) ** 2, -u**2 * (4 * u + 1), (u + 1) ** 2])
        out["GA2_L_G"] = report(LG, "GA2 avatar L_G")
        ok = (out["L9_slater_g2"]["infinity"]["fuchsian_at_infinity"] and out["F3"]["infinity"]["fuchsian_at_infinity"]
              and not out["GA2_L_G"]["infinity"]["fuchsian_at_infinity"]
              and any({k: v for k, v in e.get("exponents_c_of_exp_c_x_sigma", {}).items() if k != "0"} == {"1/4": 2}
                      for e in out["GA2_L_G"]["infinity"]["edges"]))
        out["controls_pass"] = bool(ok)
        print("controls pass (L9 Fuchsian, F_3 Fuchsian, L_G slope -1 with c = 1/4 double):", ok, flush=True)
        (DATA / "gy_newton_controls.json").write_text(json.dumps(out, indent=1, default=str))
    else:
        L = load(sys.argv[1], sys.argv[2])
        out = report(L, sys.argv[3])
        (DATA / f"gy_newton_{sys.argv[3]}.json").write_text(json.dumps(out, indent=1, default=str))
