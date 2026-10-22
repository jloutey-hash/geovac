"""FP: collect the exact pieces into one JSON (strings).  Q-form = list of Q_i(u), operator sum_i t^(2i) Q_i(theta_t)."""
import json
from pathlib import Path
HERE = Path(__file__).parent / "data"
ld = lambda n: json.loads((HERE / f"four_density_path_{n}.json").read_text())
A, D, K, S, S2, X, V, Z, O2, U, AA, G2 = (ld(n) for n in ["A", "D", "K", "S", "S2", "X", "V", "Z", "O2", "U", "AA", "G2"])
out = {
    "convention": "Q-form: [Q_0(u), Q_1(u), ...] = sum_i t^(2i) Q_i(theta_t); delta form: coefficients in w = t^2 of delta^k, delta = theta/2",
    "L_order74_Q": A["L"]["Q"],
    "W12_Q": D["W12"]["Q"],
    "factorization": "L = M62 o W12 (exact, flint); W12 = Y4 o Z2 o Z2' o L1 o R2 o (theta - 2) (exact)",
    "theta-2": ["u - 2"],
    "R2_Q": S["R2_Q"], "L1_delta_w": S["L1_delta_w"], "R3_Q": K["R3_Q"],
    "Z2p_Q": X["Z2p_Q"], "Z2_Q": X["Z2_Q"], "X4_Q": V["X4_Q"], "R8_Q": V["R8_Q"], "Y4_Q": V["Y4_Q"],
    "maps": {
        "Legendre theta^2 - t^2 (theta+1)^2 -> R2 (theta-form coefficients a_0, a_1)": S2["A"],
        "Legendre -> Z2'": Z["Z2p"]["A"], "Legendre -> Z2": Z["Z2"]["A"],
        "F_4 Q4 = [9u^2(u-2)^2, -2u^2(5u^2+1), (u+1)^4] -> Y4": AA["A_theta_form"],
    },
    "F4_route": {"A4_Q (Mellin operator of P_h4, h4 = x^6 K0^3)": O2["A4_Q"],
                 "B_Q (B P + q P_h4 = 0)": O2["B_Q"], "q": O2["q"],
                 "C_Q (R8 P = C P_h4)": U["C_Q"],
                 "G (Q4 -> A4)": "t^-3 (theta-4)^2 (theta-2)^2 theta^2"},
    "exact_checks": {"L_div_W12_zero": G2["L_div_W12_zero_flint"], "W12_x_identity": D["W12"]["exact_x_identity"],
                     "B_x_identity": O2["exact_x_identity"], "R8_x_identity": U["exact_x_identity"],
                     **{k: v for k, v in V.items() if isinstance(v, bool)},
                     "R2oA_mod_Leg": S2["R2oA_mod_Leg_zero"], "Z2p": Z["Z2p"]["ZoA_mod_Leg_zero"], "Z2": Z["Z2"]["ZoA_mod_Leg_zero"],
                     "Y4oA_mod_Q4": AA["Y4oA_mod_Q4_zero"], "image_ann_Y4": AA["image_ann_equals_Y4"]},
}
(HERE / "four_density_path_pieces.json").write_text(json.dumps(out, indent=1))
print("written", sorted(out["exact_checks"].items()))
