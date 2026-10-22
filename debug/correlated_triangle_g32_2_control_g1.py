"""GT positive control (gamma = 1): the same pipeline recovers G4 = R3 o R1 and the exact F_3 map
stored in debug/data/ring_ext_QG_7_k3map_g1.json."""
from __future__ import annotations
import json, sys, time
from fractions import Fraction
from pathlib import Path
import sympy as sp
import flint
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib2 as H

t, xs = sp.symbols("t x")
out = {}
G4 = G.load_G4()
print("G4 exponents: t0", G.exps_of(G.indicial_t0(G4)), " x=4", G.exps_of(G.indicial_at(G4, 4)), " inf", G.exps_of(G.indicial_inf(G4)))
# (1) right factor through the exponent -9/2 Frobenius solution at t = 2
fn = H.frob_series_fn(G4, Fraction(4), Fraction(-9, 2))
r, d, nul = H.find_order_degree(fn, Fraction(4), 3, 12)
R1 = H.exact_factor(fn, Fraction(4), r, d, L_check=G4)
print("R1:", [str(c) for c in R1])
m, Qt, rem = G.pseudo_rdiv(G4, R1)
F = G.content_free(Qt)
print("G4 = F o R1 exact:", G.is_zero_op(rem), " F order", len(F) - 1, "x-degree", G.op_degree(F))
# compare with stored R1, R3
k3 = json.loads((G.DATA / "ring_ext_QG_7_k3map_g1.json").read_text())
r_st = sp.sympify(k3["R1_r"], locals={"t": t})
r_me = sp.cancel(-sp.sympify(str(R1[0]).replace("^", "**"), locals={"x": xs}).subs(xs, t**2) /
                 sp.sympify(str(R1[1]).replace("^", "**"), locals={"x": xs}).subs(xs, t**2))
print("R1 matches stored R1 (theta - r):", sp.simplify(r_me - r_st) == 0)
R3s = [sp.sympify(c, locals={"t": t}) for c in k3["R3_monic"]]
Fc1 = G.compose(F, [R1[1]])
Fm = [sp.cancel(sp.sympify(str(c).replace("^", "**"), locals={"x": xs}).subs(xs, t**2) /
                sp.sympify(str(Fc1[-1]).replace("^", "**"), locals={"x": xs}).subs(xs, t**2)) for c in Fc1]
# stored R3 is the quotient by c*(theta - r) for QG's own normalisation c; compare up to the gauge g
import ring_ext_QG_ops as O
lg = sp.cancel((Fm[2] - R3s[2]) / 3)            # theta g / g
gsol = sp.simplify(sp.exp(sp.integrate(sp.apart(lg / t, t), t)))
R3g = O.monic(O.comp(R3s, [gsol]))
same = all(sp.simplify(a - b) == 0 for a, b in zip(R3g, Fm))
print("left quotient equals stored R3 up to the gauge g =", gsol, ":", same)
out["R3_match_up_to_gauge"] = {"g": str(gsol), "match": bool(same)}
out["R1"] = [str(c) for c in R1]
out["F"] = [str(c) for c in F]
print("F exponents: t0", G.exps_of(G.indicial_t0(F)), " x=4", G.exps_of(G.indicial_at(F, 4)), " inf", G.exps_of(G.indicial_inf(F)))
print("F leading factors:", [(str(f), e) for f, e in F[-1].factor()[1]])
# (2) the stored map: into Sol(G4*), den t^3 (9t^4+312t^2-12560)^2, numerator degree 5 in x
G4a = G.adjoint(G4)
ap = [-12560, 312, 9]
den = flint.fmpz_poly(ap) ** 2
res = H.map_exact(G4a, -5, 0, [int(c) for c in den.coeffs()], 6)
A, ok = res
Ast = [sp.sympify(c, locals={"t": t}) for c in k3["A"]]
Ame = [sp.sympify(str(c).replace("^", "**"), locals={"x": xs}).subs(xs, t**2) / (t**5 * (9 * t**4 + 312 * t**2 - 12560) ** 2) for c in A]
ratio = [sp.cancel(a / b) for a, b in zip(Ame, Ast) if b != 0]
print("stored map recovered (exact identity", ok, "); ratio to stored A constant and equal:", len(set(ratio)) == 1 and all(q.is_number for q in ratio))
out["stored_map_recovered"] = bool(ok and len(set(ratio)) == 1)
# (3) GT2-style map into the order-3 composition factor F, denominators at true singular points only
found = None
for nu in (-1, 1, -3, 3, -5, 5):
    for a in range(0, 6):
        for b in range(0, 5):
            den_rest = [int(c) for c in (flint.fmpz_poly([-4, 1]) ** b).coeffs()]
            for M in (4, 8, 12, 16):
                vecs, labels = H.map_search_mod(F, nu, a, den_rest, M)
                if vecs:
                    found = (nu, a, b, M, len(vecs))
                    break
            if found:
                break
        if found:
            break
    if found:
        break
print("map F_3 -> Sol(F) (mod p): t^nu x^-a (x-4)^-b, numerator degree M:", found)
if found:
    nu, a, b, M, _ = found
    # minimal M at this denominator
    while M > 0 and H.map_search_mod(F, nu, a, [int(c) for c in (flint.fmpz_poly([-4, 1]) ** b).coeffs()], M - 1)[0]:
        M -= 1
    res = H.map_exact(F, nu, a, [int(c) for c in (flint.fmpz_poly([-4, 1]) ** b).coeffs()], M)
    print("exact map into F:", res[1] if res else None, " (nu, a, b, M) =", (nu, a, b, M))
    out["map_into_F"] = {"nu": nu, "a": a, "b": b, "M": M, "identity_exact": bool(res and res[1]),
                         "A_numerators": [str(c) for c in res[0]] if res else None}
(G.DATA / "correlated_triangle_g32_control_g1.json").write_text(json.dumps(out, indent=1))
