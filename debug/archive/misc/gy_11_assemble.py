"""GY: assemble debug/data/gy_pieces.json in the format of debug/data/ga3_pieces.json: the exact operator L9G (Y_G,log's
minimal annihilator over Q(x), theta form and Q form), its composition factors and chain (strings), the kite's exact
operator and its factors, the Newton polygons, the local data, the Jordan structure at x = 0, the F_3 map searches
(controls and targets), the on-Y residuals with both perturbation normalisations, the series and loop controls, the
minimality bounds.  Everything exact as strings.  Run after the gy_* stages."""
import json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

DATA = Path(__file__).parent / "data"


def J(n, default=None):
    p = DATA / n
    if not p.exists():
        return default
    try:
        return json.loads(p.read_text())
    except Exception:
        return {"raw": p.read_text()[-4000:]}


def L(n):
    p = DATA / n
    return p.read_text().strip() if p.exists() else None


pieces = J("gy_pieces_L9G.json") or {}
out = {k: v for k, v in pieces.items()}
op = J("gy_5_op_ylog.json")
out["operator_L9G"] = op
kite_op = J("gy_5_op_kite.json")
out["kite_operator"] = kite_op
out["kite_pieces"] = J("gy_pieces_kite.json")
out["kite_local_data"] = J("gy_8_pieces_kite.json")
out["local_data"] = J("gy_8_pieces_L9G.json")
out["newton_polygon"] = {"controls": J("gy_newton_controls.json"), "L9G": J("gy_newton_L9G.json"), "kite": J("gy_newton_kite.json"),
                         "mellin_route": {k: v for k, v in (J("gy_6_mellin_main.json") or {}).items() if k not in ("LM_theta_form", "LM_Q")}}
out["jordan_at_x0"] = {"controls": J("gy_logfrob_controls.json"), "L9G": J("gy_logfrob_L9G.json"), "factor_1": J("gy_logfrob_F1.json"),
                       "factor_2": J("gy_logfrob_F2.json"), "logfree_counts": {"controls": J("gy_logfree_controls.json"), "L9G": J("gy_logfree_L9G.json"),
                                                                              "F1": J("gy_logfree_F1.json"), "F2": J("gy_logfree_F2.json")}}
out["order_one_twisted_searches"] = {"control_G4": J("gy_8c_hyperexp_control.json"), "factor_1": J("gy_8c_hyperexp_F1.json"),
                                    "factor_2": J("gy_8c_hyperexp_F2.json")}
out["split_parts"] = {"scans_each_part_alone_modp_N600": L("gy_5e_single_parts.log"), "E_part_annihilator_is_factor_1": L("gy_part_a_vs_F1.log"),
                      "phi_alone_float_operator_r5_d62_on_Y": J("gy_5f_float_op_r5d62.json"), "phi_alone_float_operator_r9_d13_on_Y_control": J("gy_5f_float_op_r9d13_ctrl.json")}
out["map_F3_note"] = ("No map from F_3's operator (twisted by an odd power of t) was found into factor_1, factor_2, L9G or their adjoints "
                      "within the box (t-power in {1,-1,...,-11,3,5}, (x-4)^-b with b <= 6, apparent factor^-c with c <= 2, numerator degree <= 128); "
                      "so there are no N0, N1, E to re-check. The same code recovers the stored Slater gamma = 2 map (control_cf_g2_factor_F).")
out["deeper_probes"] = {"factor_1": J("gy_8b_probe_F1_d220.json"), "factor_2": J("gy_8b_probe_F2_d220.json")}
out["map_F3"] = {"control_cf_g2_factor_F": J("gy_10_k3map_ctrl_cf_g2_F.json"),
                 "targets": {t: J(f"gy_10b_k3map_{t}.json") for t in ("F2", "F2adj", "F1", "F1adj", "L9G", "L9Gadj")}}
out["on_Y"] = {"grid_h_1_128": J("gy_9_onY_L9G_h128.json"), "grid_h_1_16_UNCONVERGED": J("gy_9_onY_L9G.json")}
out["controls"] = {
    "kernel_derivation": J("gy_0_kernel.json"),
    "routes_and_forms": J("gy_1_checks.json"),
    "radius_check": L("gy_radius_check.log"),
    "exp_pair_series_control": L("gy_exp_control.log"),
    "modp_vs_float_validation": (J("gy_5d_phi_modp_N80_p2305843009213693951.json") or {}).get("validation"),
    "float_recheck_of_L9G": (op or {}).get("float_recheck"),
    "float_recheck_of_kite_operator": (kite_op or {}).get("float_recheck"),
    "F3_map_control": (J("gy_10_k3map_ctrl_cf_g2_F.json") or {}).get("equals_stored_CF6_map_up_to_scale"),
    "newton_controls_pass": (J("gy_newton_controls.json") or {}).get("controls_pass"),
    "minimality": {"modp_scan_N600": {"ylog": L("gy_5e_scan_ylog_N600.log"), "kite": L("gy_5e_scan_kite_N600.log")},
                   "modp_scan_N1200": {"ylog": L("gy_5e_scan_ylog_N1200.log"), "kite": L("gy_5e_scan_kite_N1200.log")},
                   "float_search_phi_alone_tol300": J("gy_5_search_phi_N800_tol300.json"),
                   "float_search_phi_alone_tol340": J("gy_5_search_phi_N800_tol340.json"),
                   "float_search_stacked_ab_tol300": J("gy_5_search_ab_N800_tol300.json")},
}
out["note"] = ("GY1-GY3 (docs/predictions/2026-10-23_gaussian_geminal_triangle.md): Y_G, g = exp(-r^2) on both edges, gamma = 1, scales (1,1,t), "
               "apex, fingerprint fibre. Operators sum_k c_k(x) theta^k, x = t^2, theta = t d/dt, coefficient lists as flint polynomial strings; "
               "Q-form sum_i x^i Q_i(u), u = theta. factor_i: composition factors of L9G, i = 1 rightmost; chain_Ri exact right divisors. "
               "The kite operator is in z with theta_t = 2 z d/dz in 'Q' and theta_z in 'Q_theta_z'.")
(DATA / "gy_pieces.json").write_text(json.dumps(out, indent=1, default=str))
print("written debug/data/gy_pieces.json")
