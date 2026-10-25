"""GA3 positive control for the exact mod-p route at gamma = 2: on the EXPONENTIAL shape's mod-p kite series
(ga3_5d_modp_series.py exp, validated against the CF5 float series), the (r, d) nullspace of the Y_log ansatz
must be empty at (8, 39) and (9, 24) and one-dimensional at (9, 25), and that null vector must equal L9g
(debug/data/cf_g2_5_op_g2_L9.json, the CF5-CF7 annihilator) reduced mod p.
usage: python ga3_5f_exp_control.py <exp prime json>"""
from __future__ import annotations
import json, sys
from pathlib import Path
sys.set_int_max_str_digits(10**7)
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import ga3_5e_modp_guess as Q5

rec = json.loads(Path(sys.argv[1]).read_text())
p, phi_p = rec["p"], rec["phi_over_pi_mod_p"]
assert rec["shape"] == "exp" and rec["gamma"] == "2"
out = {"p": p, "validation_mismatches_vs_CF5_float_series": rec.get("validation_mismatches")}
for r, d in ((8, 39), (9, 24), (9, 25)):
    vecs, labels = Q5.nullspace(phi_p, p, r, d)
    print(f"exp control: r = {r}, d = {d}: nullity {len(vecs)}", flush=True)
    out[f"nullity_r{r}_d{d}"] = len(vecs)
L9 = json.loads((Path(__file__).parent / "data" / "cf_g2_5_op_g2_L9.json").read_text())
op = G.from_Q([sp.sympify(q, locals={"u": G.u_sym}) for q in L9["Q"]])
vec = []
for (i, k) in labels:
    cs = op[k].coeffs()
    vec.append(int(cs[i]) % p if i < len(cs) else 0)
v = vecs[0]
idx = max(i for i, a in enumerate(vec) if a)
same = G.normalise(v, idx, p) == G.normalise(vec, idx, p) if v[idx] else False
print("mod-p null vector equals L9g mod p (after normalisation):", same, flush=True)
out["null_vector_equals_L9g_mod_p"] = bool(same)
(Path(__file__).parent / "data" / "ga3_5f_exp_control.json").write_text(json.dumps(out, indent=1))
