"""PF1-PF4: assemble debug/data/pf_pieces.json from the pf_* stages (format of debug/data/ga3_pieces.json where an object
exists; here no annihilator of Y or of the kite was recovered, so the operator/factor/map keys are present and empty, with
the reason, and the controls and the principal-sheet discriminator are recorded as computed).

Object: the exponential-geminal correlated triangle Y(t) at the fingerprint fibre with vertex 1 = density rho_i of
scale a = 1, vertex 2 = density rho_j of scale b = 2, vertex 3 = apex density rho_k of scale t; geminal
exp(-gamma_1 r_13) with gamma_1 = 3/2, geminal exp(-gamma_2 r_23) with gamma_2 = 3, Coulomb 1/r_12 (QG conventions)."""
import json
from pathlib import Path
import mpmath as mp

DATA = Path(__file__).parent / "data"
J = lambda n: json.loads((DATA / n).read_text())
mp.mp.dps = 60


def betas(path, R, ns):
    c = [mp.mpf(s) for s in J(path)["phi_n"]]
    out = {}
    for n in ns:
        a = mp.log(abs(c[n] * mp.mpf(R) ** n))
        b = mp.log(abs(c[n + 1] * mp.mpf(R) ** (n + 1)))
        out[str(n)] = {"ratio_c_n_over_c_n+1": mp.nstr(c[n] / c[n + 1], 12),
                       "beta_n": mp.nstr(-(b - a) / (mp.log(n + 1) - mp.log(n)) - 1, 6)}
    return out


out = {
    "record": "docs/predictions/2026-10-23_price_formula_slots.md (PF1-PF4, committed b84e1231 before computing)",
    "object": ("Y(t) = <rho_i rho_j rho_k exp(-gamma_1 r_13) exp(-gamma_2 r_23) / r_12>, scales (a, b, t) = (1, 2, t), t on the apex k, "
               "gamma_1 = 3/2 on 1-3, gamma_2 = 3 on 2-3; QG conventions: rho~ = (c k^2+1)^(-3/2), scale = inverse length, "
               "f~ = 8 pi gamma/(k^2+gamma^2)^2, Y = (1/(4 pi^2)) t^3 int z^(1/2) (z+t^2)^(-3/2) phi(z) dz, x = t^2, theta = t d/dt"),
    "status": "STOP for feasibility (decision gate): no annihilator of the kite phi(z) or of Y(t) recovered within the run cap; PF1-PF4 UNSCORED",
    "note": ("Operators would be sum_k c_k(x) theta^k, x = t^2, theta = t d/dt, coefficient lists as flint polynomial strings, as in ga3_pieces.json. "
             "None was recovered, so the operator, factor, chain and map keys are empty."),
    "L": None, "factors": [], "chain": [], "divisions_exact": {}, "maps": {}, "order_one_solutions": {}, "on_Y": None,
    "kite_operator": None,
    "kite_formula": ("Phi~ = 2 g1 g2 a^3 b^3 int int du dv sqrt(u(1-u)) v(1-v) K/(M^3 N^5), M^2 = u a^2 + (1-u) b^2 + u(1-u) z, "
                     "N^2 = (1-v) g1^2 + v g2^2 + v(1-v) z, k2 = (1-u-v)^2 z, s = M + N, D = s^2 + k2, "
                     "K = 3/D + 6 s^2/D^2 + 2(N^2+3MN)(3s^2-k2)/D^3 + 24 M N^2 s (s^2-k2)/D^4 (u weights the a line, v the gamma_2 line)"),
    "controls": {
        "pf_1_checks": J("pf_1_checks.json"),
        "radius_and_slice_series": J("pf_series_radius_check.json"),
        "direct_route_vs_cauchy": J("pf_series_direct_test_N200_d120.json"),
        "legendre_map_search_control_mg_g3": {k: v for k, v in J("pf_map_ctrl_mg_g3.json").items() if k in ("hits", "equals_stored_mg_g3_map_up_to_scale", "general_source_equals_c2b_feyn_A_at_B1")},
        "map_recheck_controls": {"sources_agree": J("pf_map_recheck_source_control.json"), "stored_mg_g3_map": J("pf_map_recheck_ctrl_mg_g3.json")},
        "R_V_map_controls": J("pf_10_rvmap_control.json"),
        "kite_first_route_control_g3": {k: v for k, v in J("pf_5_op_ctrl_g3.json").items() if k != "Q"},
        "kite_first_route_control_g3_equals_stored_L9": True,
        "discriminator_script_control_g3": J("pf_discriminator_ctrl_g3.json"),
    },
    "series": {
        "route": "direct power series in z at every Feynman node (arb_series), substituted quadrature mu = M0, sigma = N0 (pf_4b_series_direct.py)",
        "main": "debug/data/pf_series_direct_explore_N300_d280_460x460_sub.json",
        "twin": "debug/data/pf_series_direct_explore_N300_d280_400x400_sub.json",
        "twin_agreement_by_n": {"0": "4.48e-280", "50": "5.98e-280", "100": "3.01e-280", "150": "3.07e-280", "200": "1.88e-254", "250": "1.61e-232", "299": "1.68e-216"},
        "phi0": J("pf_1_checks.json")["e_phi0"],
    },
    "guessing": {"search_tol150_N300": J("pf_5_search_kite_N300.json"), "scan_N300": J("pf_5s_scan_N300.json"),
                 "verdict": "residual falls smoothly with the number of unknowns at orders 6, 9, 12 (1e-29 at U = 63 to 1e-196 at U = 259); no drop to the data's precision: no genuine kite annihilator with (r+1)(d+1) <= 260"},
    "discriminator": {
        "question": "x = 25/4 (t = 5/2) and x = 25 (t = 5) for Y; z = -25/4 and z = -25 for the kite (the direct-pairing thresholds)",
        "operator_level": "UNDETERMINED: no annihilator of Y or of the kite in hand, so no local exponents at any point",
        "principal_sheet": {
            "true_kite": {"taylor_radius": "9 = (a+b)^2 (two radii 4.5 and 6 agree; ratio test)",
                          "asymptotics_R9": betas("pf_series_direct_explore_N300_d280_460x460_sub.json", 9, (50, 100, 150, 200, 250, 298)),
                          "z=-25/4": "analytic (|z| = 6.25 < 9)", "z=-25": "beyond the radius: not tested",
                          "z in {-1/4,-1,-9/4,-4}": "analytic (inside the radius)", "z=-9": "singular, leading exponent -> -1 (beta_n = -1.008 at n = 298)",
                          "z in {-49/4,-16,-81/4}": "beyond the radius: not tested"},
            "geminal_swapped_kite (the direct-pairing object)": {
                "series": "debug/data/pf_series_swapped_g_M256_d60_rho3.json (pf_4_series.py, config a=1, b=2, gamma_1=3, gamma_2=3/2)",
                "taylor_radius": "25/4 = (a + 3/2)^2",
                "asymptotics_R25/4": betas("pf_series_swapped_g_M256_d60_rho3.json", mp.mpf(25) / 4, (25, 50, 75, 100, 125)),
                "z=-25/4": "singular on the principal sheet, leading exponent -1/2"},
            "route_A_selects": "the derived kite (k2 = (1-u-v)^2 z): route A agrees to 0, 0, 1.7e-21 at w = 0.5, 0.9, 2.1; the swapped object differs by 1.2e-2, 3.4e-2, 1.2e-1",
        },
    },
    "predictions": {"PF1": "UNSCORED (STOP)", "PF2": "UNSCORED (STOP)", "PF3": "UNSCORED (STOP)", "PF4": "UNSCORED (STOP)"},
    "grid_Y_32_digits": {"file": "debug/data/pf_grid_d32_h0.0625.json", "nodes": 105,
                         "Y": {"0.8": "6.7944071820188723310970226184055e-4", "1.3": "1.6007486613677833864873338827962e-3",
                               "1.7": "2.3993024977063884575086159022422e-3", "2.5": "3.904454202484652037339603261408e-3"}},
}
(DATA / "pf_pieces.json").write_text(json.dumps(out, indent=1, default=str))
print("written debug/data/pf_pieces.json")
