"""GT stage 7: refine the order-3 piece C2 = K6/V4.  The split search found a hyperexponential solution
s = t^3 (x-1/4)^-4 (x-25/4)^-4 (x-9)^-9/2 P(x) of C2.  Here: P exactly (mod p nullspace, CRT, rational
reconstruction), rho1 = the order-1 right factor, C2 = C2b o rho1 exactly, split search on C2b (order 2),
and, if C2b splits, its order-1 right factor exactly.  Then the full chain is re-verified by exact
pseudo-division."""
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
import correlated_triangle_g32_lib4 as S

d = json.loads((G.DATA / "correlated_triangle_g32_pieces.json").read_text())
P_ = {k: G.parse_op(d[k]) for k in ("ra", "rb", "R2", "B2", "B2q", "V4", "C2", "K6", "F")}
L9 = G.load_L9()
C2 = P_["C2"]


def hyperexp_exact(L, f0, combo, N, maxprimes=200):
    """exact P (in t) for s = t^f0 prod (1 - t/a)^{e_a} P(t), via mod p nullspaces + CRT"""
    got, last = [], None
    for p in H.PRIMES[:maxprimes]:
        NT = 2 * (N + 1) + 60 + 2 * G.op_degree(L)
        pre = S.prefactor(combo, NT, p)
        Lt = S.t_series_ops(L, p)
        ep = G.frac_mod(Fraction(f0), p)
        cols = []
        for j in range(N + 1):
            s = [0] * j + pre[: NT - j]
            cols.append(S.apply_t(Lt, ep, s, p, NT))
        valid = max(NT - 2 * G.op_degree(L) - 8, N + 31)
        rows = [[col[n] for col in cols] for n in range(min(valid, NT))]
        vecs = G.nullspace_mod(rows, N + 1, p)
        if len(vecs) != 1:
            continue
        v = vecs[0]
        idx = max(i for i, x in enumerate(v) if x)
        got.append((p, G.normalise(v, idx, p)))
        if len(got) % 4:
            continue
        m, res = G.crt_vectors(got)
        fr = G.reconstruct(m, res)
        if fr is not None and fr == last:
            return fr
        last = fr
    return None


def order1_from_s(f0, combo, Pt):
    """rho = Dn theta - Dn (theta s / s), as a polynomial operator in x (s even/odd in t handled: P must be even)."""
    t, x = sp.symbols("t x")
    Ps = sum(sp.Rational(c.numerator, c.denominator) * t**i for i, c in enumerate(Pt))
    assert sp.expand(Ps - Ps.subs(t, -t)) == 0, "P not even in t"
    s_log = sp.Integer(f0) * sp.log(t) + sum(sp.Rational(e.numerator, e.denominator) * sp.log(1 - t / sp.Rational(a.numerator, a.denominator))
                                             for a, e in combo.items() if e != 0) + sp.log(Ps)
    r = sp.together(sp.expand(t * sp.diff(s_log, t)))
    num, den = sp.fraction(sp.cancel(r))
    num, den = sp.expand(num.subs(t, sp.sqrt(x))), sp.expand(den.subs(t, sp.sqrt(x)))
    pn = sp.Poly(num, x).all_coeffs()[::-1]
    pd = sp.Poly(den, x).all_coeffs()[::-1]
    lc = sp.ilcm(*[sp.Rational(c).q for c in pn + pd])
    op = [flint.fmpz_poly([-int(c * lc) for c in pn]), flint.fmpz_poly([int(c * lc) for c in pd])]
    return G.content_free(op), sp.factor(Ps.subs(t, sp.sqrt(x)))


out = {}
cands = S.candidates(C2, "sol")
hits = [(f0, combo, N) for f0, combo, N in cands if S.test_candidate(C2, f0, combo, N, "sol")]
f0, combo, N = hits[0]
Pt = hyperexp_exact(C2, f0, combo, N)
rho1, Px = order1_from_s(f0, combo, Pt)
print("C2's hyperexponential solution: s = t^%s * prod (1-t/a)^e_a * P,  P(x) =" % f0, Px)
print("   exponents:", {str(a): str(e) for a, e in combo.items() if e})
ok, C2b = (lambda m: (G.is_zero_op(m[2]), G.content_free(m[1])))(G.pseudo_rdiv(C2, rho1))
print("C2 = C2b o rho1 exact:", ok, " C2b order", len(C2b) - 1, "x-degree", G.op_degree(C2b))
out["s_C2"] = {"t_power": str(f0), "factors": {str(a): str(e) for a, e in combo.items() if e}, "P_x": str(Px)}
print("== C2b local data: t0", Lc.exps_t0(C2b)[0], "mu", Lc.exps_inf_growth(C2b)[0],
      {str(x0): Lc.exps_finite_t(C2b, x0)[0] for x0 in Lc.XPTS})
for x0 in [Fraction(0)] + Lc.XPTS + ["inf"]:
    print(f"   {x0}:", [c["exponents"] for c in Lc.monodromy_type(C2b, x0)])
n, hb = S.split_search(C2b, "sol")
out["C2b_split_search"] = {"admissible": n, "hits": hb}
chain = {"ra": P_["ra"], "R2": P_["R2"], "V4": P_["V4"], "V5": G.content_free(G.compose(rho1, P_["V4"]))}
pieces = {"ra": P_["ra"], "rb": P_["rb"], "B2q": P_["B2q"], "rho1": rho1}
if hb:
    f0b, combob, Nb = None, None, None
    for f0_, cb, N_ in S.candidates(C2b, "sol"):
        if S.test_candidate(C2b, f0_, cb, N_, "sol"):
            f0b, combob, Nb = f0_, cb, N_
            break
    Ptb = hyperexp_exact(C2b, f0b, combob, Nb)
    sig1, Pxb = order1_from_s(f0b, combob, Ptb)
    okb, sig2 = (lambda m: (G.is_zero_op(m[2]), G.content_free(m[1])))(G.pseudo_rdiv(C2b, sig1))
    print("C2b's hyperexponential solution: t^%s" % f0b, {str(a): str(e) for a, e in combob.items() if e}, "P(x) =", Pxb)
    print("C2b = sig2 o sig1 exact:", okb)
    out["s_C2b"] = {"t_power": str(f0b), "factors": {str(a): str(e) for a, e in combob.items() if e}, "P_x": str(Pxb)}
    chain["V6"] = G.content_free(G.compose(sig1, chain["V5"]))
    pieces["sig1"] = sig1
    pieces["sig2"] = sig2
else:
    pieces["C2b"] = C2b
chain["K6"] = P_["K6"]
chain["L9"] = L9
pieces["F"] = P_["F"]
# verify the chain by exact pseudo-division, consecutive members
names = list(chain)
divs = {}
for a, b in zip(names[:-1], names[1:]):
    m, Qt, rem = G.pseudo_rdiv(chain[b], chain[a])
    divs[f"{b}/{a}"] = G.is_zero_op(rem)
    print(f"   {b} = (order {len(chain[b]) - len(chain[a])}) o {a}: exact {G.is_zero_op(rem)}")
out["chain_orders"] = {k: len(v) - 1 for k, v in chain.items()}
out["chain_divisions_exact"] = divs
out["factor_orders"] = {k: len(v) - 1 for k, v in pieces.items()}
G.save_op(G.DATA / "correlated_triangle_g32_series_final.json", {**{"chain_" + k: v for k, v in chain.items() if k != "L9"},
                                                                  **{"factor_" + k: v for k, v in pieces.items()}},
          {"note": "composition series of L9 over Q(t): right-factor chain chain_* (each an exact right divisor of the next and of L9); "
                   "factor_* are representatives of the composition factors (exact left quotients). x = t^2, theta = t d/dt.",
           "summary": out})
(G.DATA / "correlated_triangle_g32_7_refine.json").write_text(json.dumps(out, indent=1, default=str))
