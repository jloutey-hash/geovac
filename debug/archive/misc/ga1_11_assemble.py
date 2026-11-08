"""GA1: assemble debug/data/ga1_pieces.json: the exact operator L9lin, its composition factors and chain (strings),
the order-one solutions, the Legendre map (factor_5) and the F_3 map (factor_6), the local data, the controls and
the on-Y residuals, as recorded by the ga1_* stages.  Also a 1/1000 perturbation check of the F_3 map here."""
import json, sys
from pathlib import Path
import flint
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib2 as H

DATA = Path(__file__).parent / "data"
J = lambda n: json.loads((DATA / n).read_text())
pieces = J("ga1_pieces_lin_g1.5.json")
out = {k: v for k, v in pieces.items()}
out["operator_L9lin"] = J("ga1_5_op_lin_g1.5_L9lin.json")
out["local_data"] = J("ga1_8_pieces_lin_g1.5.json")
out["order_one_solutions"] = J("ga1_8b_order1_lin_g1.5.json")
out["survey"] = {k: v for k, v in J("ga1_6_survey_lin_g1.5.json").items() if k in ("L", "apparent", "hyperexp_L", "hyperexp_Ladj")}
k3 = J("ga1_10_k3map_lin_g1.5.json")
out["map_F3_to_factor_6"] = k3
if k3.get("identity_exact"):
    F = G.parse_op(pieces["factor_6"])
    A = G.parse_op(k3["A_numerators_theta_0_1_2"])
    s, b = k3["t_power"], -k3["x_minus_4_power"]
    nu, a = (s, 0) if s > 0 else (s % 2, (s % 2 - s) // 2)
    den = [int(c) for c in (flint.fmpz_poly([-4, 1]) ** b).coeffs()]
    out["map_F3_to_factor_6"]["reverified_exact"] = H.verify_map(F, nu, a, den, A)
lm = J("ga1_map_lin_g1.5_factor_5.json")
out["map_Legendre_to_factor_5"] = {"form": "A = (N0 + N1 d/dx)/E from L_Leg(gamma = 3/2) (degree -1/2 in z = (t^2 + 1 - gamma^2)/(2t), twisted by t^(-1/2), in x = t^2; 'twist' x^(1/2) means a further twist by t^(1/2)) into factor_5; factor_5(A y) = 0 for y in Sol(source), exact over Q",
                                   "source_coefficients_d_dx_form": lm.get("source"), "map": lm.get("map"), "twists_searched": lm.get("twists_searched"),
                                   "boxes_searched": [{k: r[k] for k in ("twist", "k", "m", "a", "dN", "nullity")} for r in lm.get("runs", [])]}
out["on_Y"] = J("ga1_9_onY_L9lin.json")
out["controls"] = {
    "loop_and_kite_checks": J("ga1_1_checks.json"),
    "radius_check": (DATA / "ga1_series_radius_check.log").read_text().strip(),
    "L9_from_stored_exp_series": (DATA / "ga1_5b_control_compare.log").read_text().strip(),
    "modp_series_validation": {"exp": J("ga1_5d_phi_modp_exp_N400_p2305843009213693951.json").get("validation_mismatches"),
                               "lin": J("ga1_5d_phi_modp_lin_N420_p2305843009213693951.json").get("validation_mismatches")},
    "modp_guesser_reproduces_L9": (DATA / "ga1_5e_exp_control.log").read_text().strip(),
    "peeler_on_L9": (DATA / "ga1_7_chain_ctrl_L9.log").read_text().strip().splitlines()[-4:],
    "F3_map_control_on_L9_factor": J("ga1_10_k3map_ctrl_L9.json"),
    "Legendre_map_control_on_L9_factor": J("ga1_map_ctrl_L9_factor_5.json").get("map"),
    "Legendre_map_control_CF2": (DATA / "ga1_map_control.log").read_text().strip(),
    "minimality": {"float_search_tol250": J("ga1_5_search_lin_g1.5_N800_tol250.json"), "float_search_tol330": J("ga1_5_search_lin_g1.5_N800_tol330.json"),
                   "modp_nullity": (DATA / "ga1_5e_nullity_lin.log").read_text().strip().splitlines()},
}
out["note"] = ("GA1 (docs/predictions/2026-10-22_geminal_shapes.md): Y_lin, f = r exp(-gamma r) on both edges, gamma = 3/2, scales (1,1,t), apex, "
               "fingerprint fibre. Operators sum_k c_k(x) theta^k, x = t^2, theta = t d/dt, coefficient lists as flint polynomial strings. "
               "factor_i: composition factors, i = 1 rightmost; chain_Ri = (1/g) factor_i o ... o factor_1 exact right divisors of L9lin.")
(DATA / "ga1_pieces.json").write_text(json.dumps(out, indent=1, default=str))
print("written debug/data/ga1_pieces.json")
