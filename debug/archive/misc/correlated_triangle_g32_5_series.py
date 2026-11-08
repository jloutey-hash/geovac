"""GT stage 5 (gamma = 3/2): the right-factor chain below F.
ra (order 1, R2's rational solution) | R2 (stored, exact) | B2 (minimal annihilator of L9's t = 3 solution
with exponent -3/2) | V4 = LCLM(R2, B2) (minimal annihilator of R2's basis and that solution) | K6.
Each inclusion is an exact right division; the composition factors are the exact left quotients
rb = R2/ra, B2q = V4/R2, C2 = K6/V4, F = L9/K6."""
from __future__ import annotations
import json, sys, time
from fractions import Fraction
from pathlib import Path
import flint
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib2 as H
import correlated_triangle_g32_lib3 as Lc

u = G.u_sym
L9 = G.load_L9()
d = json.loads((G.DATA / "correlated_triangle_g32_K6_F.json").read_text())
F, K6 = G.parse_op(d["F"]), G.parse_op(d["K6"])
rec = json.loads((G.DATA / "ring_ext_QG_6e_order2_g1.5.json").read_text())
R2 = G.content_free(G.from_Q([sp.sympify(q, locals={"u": u}) for q in rec["R2_Q"]]))
out = {}


def divides(L, R):
    m, Qt, rem = G.pseudo_rdiv(L, R)
    return G.is_zero_op(rem), G.content_free(Qt) if G.is_zero_op(rem) else None


# ra: R2's rational solution s = (88 x + 939)/x ;  theta s = -1878/x  ->  ra = (88x+939) theta + 1878
ra = [flint.fmpz_poly([1878]), flint.fmpz_poly([939, 88])]
ok, rb = divides(R2, ra)
print("R2 = rb o ra exact:", ok)
ok2, _ = divides(K6, R2)
print("K6 = . o R2 exact:", ok2)
X9 = Fraction(9)
for e in (Fraction(-3, 2), Fraction(-5, 2)):
    fn = H.frob_series_fn(L9, X9, e)
    res = H.find_order_degree(fn, X9, 5, 40)
    print(f"L9 solution at t = 3, exponent {e}: minimal (order, x-degree, nullity) = {res}")
    out[f"t3_exp_{e}"] = res
fn = H.frob_series_fn(L9, X9, Fraction(-3, 2))
r, dd, _ = H.find_order_degree(fn, X9, 5, 40, verbose=False)
B2 = H.exact_factor(fn, X9, r, dd, L_check=L9)
print("B2 exact; L9 = . o B2:", B2 is not None, " K6 = . o B2:", divides(K6, B2)[0])
# V4 = LCLM(R2, B2): R2's Taylor basis at x = 9 (regular for R2) plus the B2 solution


def fnV(p, N):
    R2p = G.op_mod(R2, p)
    sl = [(Fraction(j), G.frob_mod(R2p, X9, j, N, p)) for j in (0, 1)]
    sl.append((Fraction(-3, 2), G.frob_mod(G.op_mod(L9, p), X9, Fraction(-3, 2), N, p)))
    return sl


rV, dV, _ = H.find_order_degree(fnV, X9, 6, 60, rmin=3)
V4 = H.exact_factor(fnV, X9, rV, dV, L_check=K6)
print(f"V4 = LCLM(R2, B2): order {rV}, x-degree {dV}; exact; K6 = . o V4:", V4 is not None)
okR, B2q = divides(V4, R2)
okB, _ = divides(V4, B2)
okK, C2 = divides(K6, V4)
print("V4 = B2q o R2:", okR, " V4 = . o B2:", okB, " K6 = C2 o V4:", okK)
okF, Fq = divides(L9, K6)
print("L9 = F o K6:", okF, " (F as stored):", all((a - b).is_zero() for a, b in zip(Fq, F)))
pieces = {"ra": ra, "rb": rb, "R2": R2, "B2": B2, "B2q": B2q, "V4": V4, "C2": C2, "K6": K6, "F": F}
for k, v in pieces.items():
    print(f"  {k}: order {len(v)-1}, x-degree {G.op_degree(v)}")
G.save_op(G.DATA / "correlated_triangle_g32_pieces.json", pieces,
          {"note": "operators sum_k c_k(x) theta^k, x = t^2, theta = t d/dt. Exact chain: Sol(ra) < Sol(R2) < Sol(V4) < Sol(K6) < Sol(L9); "
                   "R2 = rb o ra, V4 = B2q o R2, K6 = C2 o V4, L9 = F o K6 (each verified by exact pseudo-division); "
                   "V4 = LCLM(R2, B2), B2 = minimal annihilator of L9's t = 3 exponent -3/2 solution.",
           "divisions_exact": {"R2/ra": ok, "K6/R2": ok2, "V4/R2": okR, "V4/B2": okB, "K6/V4": okK, "L9/K6": okF}})
for k in ("rb", "B2q", "C2", "B2"):
    v = pieces[k]
    print(f"== {k} local data: t0 {Lc.exps_t0(v)[0]}, mu {Lc.exps_inf_growth(v)[0]}, " +
          ", ".join(f"x={x0}: {Lc.exps_finite_t(v, x0)[0]}" for x0 in Lc.XPTS))
    print("   leading factors:", [(str(f)[:50], e) for f, e in v[-1].factor()[1]])
    for x0 in [Fraction(0)] + Lc.XPTS + ["inf"]:
        print(f"   {x0}:", [c["exponents"] for c in Lc.monodromy_type(v, x0)])
(G.DATA / "correlated_triangle_g32_5_series.json").write_text(json.dumps(out, indent=1, default=str))
