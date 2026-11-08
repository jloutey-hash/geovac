"""GY stage 10b: the F_3 map search of gy_10_k3map.py against a target of ANY order (F2, F1, L9G and their adjoints):
a nonzero A = t^s (x-4)^-b app(x)^-c sum_{j<3} p_j(x) theta^j with F o A = 0 mod L3 (F_3's operator 4 th^3 - t^2 (th+1)^3)
maps Sol(L3) into Sol(F); since L3 is irreducible the image is a right factor of F equivalent to F_3's operator twisted by t.
On the adjoint F* the same finds a LEFT factor (quotient) of F of that type (F_3's operator is self-adjoint up to sign).
app(x) = the apparent factors of F's leading coefficient (optional extra denominator for adjoints, whose solutions may have
poles there).  Box: s in S_RANGE, b <= bmax, c <= cmax, numerator degree M <= Mmax.
usage: python gy_10b_k3map_any.py <pieces json> <key> <adjoint 0|1> <Mmax> <bmax> <cmax> <out tag>"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import flint
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib2 as H
import cf_g2_lib as C

S_RANGE = (1, -1, -3, -5, -7, -9, -11, 3, 5)
d = json.loads(Path(sys.argv[1]).read_text())
key, adj, Mmax, bmax, cmax, tag = sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), sys.argv[7]
F = G.parse_op(d[key]) if key != "Q" else C.load_Q(sys.argv[1])
if adj:
    F = G.content_free(G.adjoint(F))
app = flint.fmpz_poly([1])
for f, e in C.lead_factors(F):
    if f.degree() >= 2 or (f.degree() == 1 and C.rational_roots([f])[0][0] != 4):
        app = app * f
print(f"target {key}{'*' if adj else ''}: order {len(F)-1}, x-degree {G.op_degree(F)}; apparent factor degree {app.degree()}", flush=True)
t0 = time.time()
found = None
Ms = [m for m in (8, 16, 24, 32, 40, 48, 56, 64, 80, 96, 112, 128, 160, 200) if m <= Mmax]
for M in Ms:
    for c in range(cmax + 1):
        for b in range(bmax + 1):
            den = flint.fmpz_poly([-4, 1]) ** b * app ** c
            den_rest = [int(v) for v in den.coeffs()]
            for s in S_RANGE:
                nu, a = (s, 0) if s > 0 else (s % 2, (s % 2 - s) // 2)
                vecs, labels = H.map_search_mod(F, nu, a, den_rest, M)
                if vecs:
                    found = (s, nu, a, b, c, M, len(vecs))
                    break
            if found:
                break
        if found:
            break
    print(f"  M = {M}: {'found ' + str(found) if found else 'none'} ({time.time()-t0:.0f}s)", flush=True)
    if found:
        break
rec = {"target": key, "adjoint": bool(adj), "order": len(F) - 1, "searched_M": Ms, "s_range": list(S_RANGE), "b_range": [0, bmax], "app_power_range": [0, cmax],
       "apparent_factor": str(app), "seconds": time.time() - t0}
if found:
    s, nu, a, b, c, M, nul = found
    den = flint.fmpz_poly([-4, 1]) ** b * app ** c
    den_rest = [int(v) for v in den.coeffs()]
    res = H.map_exact(F, nu, a, den_rest, M)
    if res:
        Aop, ok = res
        cs = [int(v) for v in Aop[0].coeffs()]
        j = max(range(len(cs)), key=lambda i: abs(cs[i]))
        cs2 = list(cs); cs2[j] = cs2[j] * 1001
        A2 = [flint.fmpz_poly(cs2)] + [cc * 1000 for cc in Aop[1:]]
        bad = H.verify_map(F, nu, a, den_rest, A2)
        rec.update({"t_power": s, "x_minus_4_power": -b, "app_power": -c, "numerator_degree": M, "identity_exact": bool(ok),
                    "perturbed_1e-3_identity": bool(bad), "A_numerators_theta_0_1_2": [str(cc) for cc in Aop],
                    "A": f"A = t^{s} (x-4)^(-{b}) app^(-{c}) * sum_j N_j(x) theta^j, from Sol(F_3 op) into Sol({key}{'*' if adj else ''})"})
        print("exact identity:", ok, " perturbed holds:", bad, flush=True)
(G.DATA / f"gy_10b_k3map_{tag}.json").write_text(json.dumps(rec, indent=1))
print("saved", flush=True)
