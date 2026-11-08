"""EC1 (docs/predictions/2026-10-22_c2b_identification.md): local data of C2b at every singular point.
Rational points exactly; roots of irreducible factors of the leading coefficient mod p (several primes)."""
import json, sys
from fractions import Fraction
from pathlib import Path
import flint, sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib3 as G3

D = json.load(open(Path(__file__).parent / "data" / "correlated_triangle_g32_series_final.json"))
L = G.parse_op(D["factor_C2b"])
x = sp.Symbol("x")
c2 = sp.Poly(sum(int(v) * x**i for i, v in enumerate(L[2].coeffs())), x)
fl = sp.factor_list(c2.as_expr())
out = {"leading_theta2_factors": [(str(f), m) for f, m in fl[1]], "points": {}}
print("theta^2 coefficient factors:")
for f, m in fl[1]:
    print("  ", sp.Poly(f, x).degree(), m, str(f)[:90])

# rational points: 0, roots of linear factors, infinity
rat = [Fraction(0)]
for f, m in fl[1]:
    pf = sp.Poly(f, x)
    if pf.degree() == 1:
        a, b = pf.all_coeffs()
        rat.append(Fraction(int(-b), int(a)))
for x0 in rat + ["inf"]:
    desc = G3.monodromy_type(L, x0)
    out["points"][str(x0)] = desc
    print("x0 =", x0, desc)

# irreducible factors of degree >= 2: mod-p roots
def ev(row, v, p):
    s = 0
    for c in reversed(row):
        s = (s * v + c) % p
    return s

cands = sorted({Fraction(n, d) for n in range(-80, 81) for d in (1, 2, 3, 4, 6)})
alg = {}
for f, m in fl[1]:
    pf = sp.Poly(f, x)
    if pf.degree() < 2:
        continue
    res = []
    for p in G.primes_from(1 << 40, 400):
        fp = flint.nmod_poly([int(c) % p for c in pf.all_coeffs()[::-1]], p)
        roots = [int(r) for r, _ in fp.roots()]
        if not roots:
            continue
        Lp = G.op_mod(L, p)
        for r0 in roots[:2]:
            B, _ = G.hform(Lp, r0)
            row = [int(c) % p for c in B[0]]
            exps = [e for e in cands if ev(row, G.frac_mod(e, p), p) == 0]
            info = []
            for e in exps:
                try:
                    G.frob_mod(Lp, Fraction(r0), e, 40, p)
                    ok = True
                except ZeroDivisionError:
                    ok = False
                info.append((str(e), ok))
            res.append({"p": p, "root": r0, "exps_logfree": info})
        if len(res) >= 6:
            break
    alg[str(f)[:60]] = {"degree": pf.degree(), "mult": m, "samples": res}
    print("factor deg", pf.degree(), "mult", m, "->", [s["exps_logfree"] for s in res])
out["algebraic_factors"] = alg
json.dump(out, open(Path(__file__).parent / "data" / "c2b_ident_A_local.json", "w"), indent=1, default=str)
