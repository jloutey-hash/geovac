"""GA1 (copy of debug/cf_g2_10_k3map.py, target = a named order-three factor of a ga1 pieces json): exact map
A = t^s (x-4)^-b sum_{j<3} p_j(x) theta^j from Sol(L3) (F_3's operator 4 th^3 - t^2 (th+1)^3; the twist is the
odd power t^s) into Sol(F).  Denominators only at F_3's and F's true singular points (t = 0, +-2): Sol(F) is
analytic away from them when F is a left quotient of L (its solutions are images of Sol(L)).
Mod p search over (s, b, M), then CRT + rational reconstruction and the exact identity F o A = 0 mod L3.
usage: python ga1_10_k3map.py <pieces json> <key> <out tag>"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import flint
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib2 as H

d = json.loads(Path(sys.argv[1]).read_text())
key, tag = sys.argv[2], sys.argv[3]
F = G.parse_op(d[key])
assert len(F) - 1 == 3
print("F: order", len(F) - 1, "x-degree", G.op_degree(F), flush=True)
t0 = time.time()
found = None
box = []
S_RANGE = (1, -1, -3, -5, -7, -9, -11, 3, 5)
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
        A, ok = res
        print("exact identity F o A = 0 mod L3:", ok, f"({time.time()-t0:.0f}s)", flush=True)
        cs = [int(v) for v in A[0].coeffs()]
        j = max(range(len(cs)), key=lambda i: abs(cs[i]))
        cs2 = list(cs); cs2[j] = cs2[j] * 1001
        A2 = [flint.fmpz_poly(cs2)] + [c * 1000 for c in A[1:]]
        bad = H.verify_map(F, nu, a, den_rest, A2)
        print("with one coefficient times 1.001 the identity holds:", bad, flush=True)
        rec.update({"t_power": s, "x_minus_4_power": -b, "numerator_degree": M, "identity_exact": bool(ok),
                    "perturbed_1e-3_identity": bool(bad),
                    "A_numerators_theta_0_1_2": [str(c) for c in A],
                    "A": "A = t^%d (x-4)^(-%d) * sum_j N_j(x) theta^j, x = t^2, theta = t d/dt, from Sol(F_3 op) into Sol(%s)" % (s, b, key)})
(G.DATA / f"ga1_10_k3map_{tag}.json").write_text(json.dumps(rec, indent=1))
