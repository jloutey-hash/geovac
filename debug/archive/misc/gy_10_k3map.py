"""GY stage 10 (copy of ga3_10_k3map.py): exact map A = t^s (x-4)^-b sum_{j<3} p_j(x) theta^j from Sol(L3) (F_3's operator
4 th^3 - t^2 (th+1)^3, twisted by the odd power t^s) into Sol(F), F a named order-three factor of a pieces json.
Mod p search over (s, b, M), CRT + rational reconstruction, exact identity F o A = 0 mod L3, 1/1000 perturbation.
Control mode: the same code on the exponential shape's order-three factor at gamma = 2 (debug/data/cf_g2_pieces.json,
key factor_F) must recover the stored map map_CF6_F3_to_F (same t-power, (x-4)-power, numerator degree, numerators up
to scale).
usage: python gy_10_k3map.py control | python gy_10_k3map.py <pieces json> <key> <out tag>"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import flint
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib2 as H

S_RANGE = (1, -1, -3, -5, -7, -9, -11, 3, 5)


def run(F, key, tag):
    assert len(F) - 1 == 3
    print("F: order", len(F) - 1, "x-degree", G.op_degree(F), flush=True)
    t0 = time.time()
    found = None
    box = []
    for M in (8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 96, 112, 128):
        for b in range(0, 9):
            for s in S_RANGE:
                nu, a = (s, 0) if s > 0 else (s % 2, (s % 2 - s) // 2)
                den_rest = [int(c) for c in (flint.fmpz_poly([-4, 1]) ** b).coeffs()]
                vecs, labels = H.map_search_mod(F, nu, a, den_rest, M)
                if vecs:
                    found = (s, nu, a, b, M, len(vecs))
                    break
            if found:
                break
        box.append(M)
        print(f"  M = {M}: {'found ' + str(found) if found else 'none'} ({time.time()-t0:.0f}s)", flush=True)
        if found:
            break
    rec = {"target": key, "searched_M": box, "s_range": list(S_RANGE), "b_range": [0, 8]}
    if found:
        s, nu, a, b, M, nul = found
        den_rest = [int(c) for c in (flint.fmpz_poly([-4, 1]) ** b).coeffs()]
        while M > 0 and H.map_search_mod(F, nu, a, den_rest, M - 1)[0]:
            M -= 1
        while b > 0:
            dr = [int(c) for c in (flint.fmpz_poly([-4, 1]) ** (b - 1)).coeffs()]
            if H.map_search_mod(F, nu, a, dr, M)[0]:
                b -= 1
                den_rest = dr
            else:
                break
        print(f"minimal box: t^{s} (x-4)^-{b}, numerator degree {M} (nullity {nul})", flush=True)
        res = H.map_exact(F, nu, a, den_rest, M)
        if res:
            Aop, ok = res
            print("exact identity F o A = 0 mod L3:", ok, f"({time.time()-t0:.0f}s)", flush=True)
            cs = [int(v) for v in Aop[0].coeffs()]
            j = max(range(len(cs)), key=lambda i: abs(cs[i]))
            cs2 = list(cs); cs2[j] = cs2[j] * 1001
            A2 = [flint.fmpz_poly(cs2)] + [c * 1000 for c in Aop[1:]]
            bad = H.verify_map(F, nu, a, den_rest, A2)
            print("with one coefficient times 1.001 the identity holds:", bad, flush=True)
            rec.update({"t_power": s, "x_minus_4_power": -b, "numerator_degree": M, "identity_exact": bool(ok),
                        "perturbed_1e-3_identity": bool(bad),
                        "A_numerators_theta_0_1_2": [str(c) for c in Aop],
                        "A": "A = t^%d (x-4)^(-%d) * sum_j N_j(x) theta^j, x = t^2, theta = t d/dt, from Sol(F_3 op) into Sol(%s)" % (s, b, key)})
    (G.DATA / f"gy_10_k3map_{tag}.json").write_text(json.dumps(rec, indent=1))
    return rec


if __name__ == "__main__":
    if sys.argv[1] == "control":
        d = json.loads((G.DATA / "cf_g2_pieces.json").read_text())
        rec = run(G.parse_op(d["factor_F"]), "factor_F (Slater gamma = 2)", "ctrl_cf_g2_F")
        ref = d["map_CF6_F3_to_F"]
        import ast
        refA = ast.literal_eval(ref["A_numerators_theta_0_1_2"]) if isinstance(ref["A_numerators_theta_0_1_2"], str) else ref["A_numerators_theta_0_1_2"]
        same = False
        if rec.get("identity_exact"):
            xs = sp.Symbol("x")
            P = lambda s_: sp.expand(sp.sympify(s_.replace("^", "**"), locals={"x": xs}))
            ratios = {sp.cancel(P(a_) / P(b_)) for a_, b_ in zip(rec["A_numerators_theta_0_1_2"], refA) if P(b_) != 0}
            same = len(ratios) == 1 and list(ratios)[0].is_number and ref["form"].startswith(f"A = t^{rec['t_power']} (x-4)^({rec['x_minus_4_power']})")
        print("control: recovers the stored CF6 map (form and numerators up to scale):", same, flush=True)
        rec["equals_stored_CF6_map_up_to_scale"] = bool(same)
        (G.DATA / "gy_10_k3map_ctrl_cf_g2_F.json").write_text(json.dumps(rec, indent=1))
    else:
        d = json.loads(Path(sys.argv[1]).read_text())
        run(G.parse_op(d[sys.argv[2]]), sys.argv[2], sys.argv[3])
