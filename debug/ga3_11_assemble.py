"""GA3: assemble debug/data/ga3_pieces.json: the exact operator L9lin at gamma = 2, its composition factors and chain
(strings), the order-one solutions, the Legendre map (order-two factor) and the F_3 map (order-three factor), the local
data, the controls and the on-Y residuals, as recorded by the ga3_* stages (copy of ga1_11_assemble.py with ga3 names;
a 1/1000 perturbation re-verification of the F_3 map and an exact comparison of the order-two factor against the
exponential shape's C2b at gamma = 2 are done here)."""
import json, sys
from pathlib import Path
import flint
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib2 as H

DATA = Path(__file__).parent / "data"
J = lambda n: json.loads((DATA / n).read_text())
x = sp.Symbol("x")
P = lambda s: sp.expand(sp.sympify(s.replace("^", "**"), locals={"x": x}))
pieces = J("ga3_pieces_lin_g2.json")
out = {k: v for k, v in pieces.items()}
out["operator_L9lin"] = J("ga3_5_op_lin_g2_L9lin.json")
out["local_data"] = J("ga3_8_pieces_lin_g2.json")
out["order_one_solutions"] = J("ga3_8b_order1_lin_g2.json")
out["survey"] = {k: v for k, v in J("ga3_6_survey_lin_g2.json").items() if k in ("L", "apparent", "hyperexp_L", "hyperexp_Ladj")}
orders = {k: len(v) - 1 for k, v in pieces.items() if k.startswith("factor_")}
k2 = [k for k, o in orders.items() if o == 2]
k3 = [k for k, o in orders.items() if o == 3]
assert len(k2) == 1 and len(k3) == 1, orders
f2, f3 = k2[0], k3[0]
km = J(f"ga3_10_k3map_lin_g2_{f3}.json")
out["map_F3_to_order_three_factor"] = km
if km.get("identity_exact"):
    F = G.parse_op(pieces[f3])
    A = G.parse_op(km["A_numerators_theta_0_1_2"])
    s, b = km["t_power"], -km["x_minus_4_power"]
    nu, a = (s, 0) if s > 0 else (s % 2, (s % 2 - s) // 2)
    den = [int(c) for c in (flint.fmpz_poly([-4, 1]) ** b).coeffs()]
    out["map_F3_to_order_three_factor"]["reverified_exact"] = H.verify_map(F, nu, a, den, A)
lm = J(f"ga3_map_lin_g2_{f2}.json")
out["map_Legendre_to_order_two_factor"] = {
    "form": f"A = (N0 + N1 d/dx)/E from L_Leg(gamma = 2) (degree -1/2 in z = (t^2 - 3)/(2t), twisted by t^(-1/2), in x = t^2; "
            f"'twist' x^(1/2) means a further twist by t^(1/2)) into {f2}; {f2}(A y) = 0 for y in Sol(source), exact over Q",
    "source_coefficients_d_dx_form": lm.get("source"), "map": lm.get("map"), "twists_searched": lm.get("twists_searched"),
    "boxes_searched": [{k: r[k] for k in ("twist", "k", "m", "a", "dN", "nullity")} for r in lm.get("runs", [])]}
# the order-two factor against the exponential shape's C2b at gamma = 2
C2b = J("cf_g2_pieces.json")["factor_C2b"]
rs = {sp.cancel(P(u) / P(v)) for u, v in zip(pieces[f2], C2b)}
out["order_two_factor_vs_exponential_C2b_g2"] = {
    "same_up_to_scale": bool(len(rs) == 1 and list(rs)[0].is_number),
    "x_degree_lin": max(sp.degree(P(v), x) for v in pieces[f2]), "x_degree_C2b": max(sp.degree(P(v), x) for v in C2b),
    "lead_factors_lin": str(sp.factor(P(pieces[f2][-1]))), "lead_factors_C2b": str(sp.factor(P(C2b[-1])))}
out["on_Y"] = J("ga3_9_onY_L9lin.json")
out["on_Y_control_dominant_monomial"] = J("ga3_9b_onY_control_L9lin.json")
out["controls"] = {
    "loop_and_kite_checks": J("ga3_1_checks.json"),
    "radius_check": (DATA / "ga3_series_radius_check.log").read_text().strip(),
    "modp_series_validation": {"exp_vs_CF5_float_series": J("ga3_5d_phi_modp_exp_N400_p2305843009213693951.json").get("validation_mismatches"),
                               "lin_vs_ga3_float_series": J("ga3_5d_phi_modp_lin_N60_p2305843009213693951.json").get("validation_mismatches")},
    "modp_guesser_reproduces_L9g": J("ga3_5f_exp_control.json"),
    "peeler_on_L9g": {"log_tail": (DATA / "ga3_7_chain_ctrl_L9g.log").read_text().strip().splitlines()[-9:],
                      "factors": J("ga3_pieces_ctrl_L9g.json")["factors"],
                      "order_two_factor_Legendre_map": J("ga3_map_ctrl_L9g_factor_5.json").get("map")},
    "F3_map_control_on_L9g_factor_F": J("ga3_10_k3map_ctrl_L9g.json"),
    "Legendre_map_control_CF5": {k: v for k, v in J("ga3_map_ctrl_C2b_g2.json").items() if k in ("map", "equals_CF5_map_up_to_scale")},
    "minimality": {"float_search_tol330": J("ga3_5_search_lin_g2_N800_tol330.json"),
                   "float_search_tol350": J("ga3_5_search_lin_g2_N800_tol350.json"),
                   "modp_nullity_N520": [l for l in (DATA / "ga3_5e_nullity_lin.log").read_text().strip().splitlines() if l.startswith("r =")],
                   "modp_nullity_N570": ([l for l in (DATA / "ga3_5e_nullity_lin_N570.log").read_text().strip().splitlines() if l.startswith("r =")]
                                         if (DATA / "ga3_5e_nullity_lin_N570.log").exists() else "pending")},
}
out["note"] = ("GA3 (docs/predictions/2026-10-23_geminal_shapes_g2.md): Y_lin, f = r exp(-gamma r) on both edges, gamma = 2, scales (1,1,t), apex, "
               "fingerprint fibre. Operators sum_k c_k(x) theta^k, x = t^2, theta = t d/dt, coefficient lists as flint polynomial strings. "
               "factor_i: composition factors, i = 1 rightmost; chain_Ri = (1/g) factor_i o ... o factor_1 exact right divisors of L9lin.")
(DATA / "ga3_pieces.json").write_text(json.dumps(out, indent=1, default=str))
print("written debug/data/ga3_pieces.json")
