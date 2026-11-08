"""CF6 (gamma = 2; copy of debug/correlated_triangle_g32_8_k3map.py, target = cf_g2_pieces factor_F). GT2: exact map A = t^s (x-4)^-b sum_{j<3} p_j(x) theta^j from Sol(L3) (F_3's operator
4 th^3 - t^2 (th+1)^3; the twist is the odd power t^s) into Sol(F), F the order-3 composition factor
(exact left quotient L9 = F o K6).  Denominators only at F_3's and F's true singular points (t = 0, +-2):
Sol(F) = K6(Sol(L9)) is analytic away from them, so A = (images) x W_f^-1 has no other poles.
Mod p search over (s, b, M), then CRT + rational reconstruction and the exact identity F o A = 0 mod L3.
Same pipeline as the gamma = 1 control (correlated_triangle_g32_2_control_g1.py)."""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import flint
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib2 as H

d = json.loads((G.DATA / "cf_g2_pieces.json").read_text())
F = G.parse_op(d["factor_F"])
print("F: order", len(F) - 1, "x-degree", G.op_degree(F), flush=True)
t0 = time.time()
found = None
box = []
for M in (8, 16, 24, 32, 40, 48, 56, 64, 72, 80):
    for b in range(0, 7):
        for s in (1, -1, -3, -5, -7, -9, -11):
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
rec = {"searched_M": box, "s_range": [1, -1, -3, -5, -7, -9, -11], "b_range": [0, 6]}
if found:
    s, nu, a, b, M, nul = found
    den_rest = [int(c) for c in (flint.fmpz_poly([-4, 1]) ** b).coeffs()]
    # minimal M and b at this s
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
        A, ok = res
        print("exact identity F o A = 0 mod L3:", ok, f"({time.time()-t0:.0f}s)", flush=True)
        rec.update({"t_power": s, "x_minus_4_power": -b, "numerator_degree": M, "identity_exact": bool(ok),
                    "A_numerators_theta_0_1_2": [str(c) for c in A],
                    "A": "A = t^%d (x-4)^(-%d) * sum_j N_j(x) theta^j, x = t^2, theta = t d/dt" % (s, b)})
(G.DATA / "cf_g2_10_k3map.json").write_text(json.dumps(rec, indent=1))
