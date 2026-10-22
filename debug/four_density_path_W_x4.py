"""FP step W: structure of X4 = R8 / (R3 o (theta-2)) (step V): leading coefficient, exponents at 0 and infinity,
log-free local solutions and their minimal annihilator orders, same for the adjoint X4*."""
import json, sys
from pathlib import Path
import sympy as sp
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_path_tools as T  # noqa: E402
u, t = T.u, sp.Symbol("t")
V = json.loads((HERE / "data" / "four_density_path_V.json").read_text())
X4 = [sp.sympify(q).subs(sp.Symbol("u"), u) for q in V["X4_Q"]]
r = max(sp.Poly(q, u).degree() for q in X4)
print("X4 order", r, "t-terms", len(X4))
print("lead:", sp.factor(sum(sp.Poly(q, u).coeff_monomial(u**r) * t**(2*i) for i, q in enumerate(X4))))
print("t^0:", sp.factor(X4[0]), " top:", sp.factor(X4[-1]))
X4s = [sp.expand(q.subs(u, -u - 2 * i)) for i, q in enumerate(X4)]
for nm, Q in [("X4", X4), ("X4*", X4s)]:
    qc = T.qcoeffs(Q)
    roots = sorted(set(sp.roots(sp.Poly(Q[0], u).as_expr(), u)))
    for lam in roots:
        if not lam.is_integer:
            print(nm, lam, "non-integer"); continue
        ser, free = T.frob(qc, int(lam), 900)
        if ser is None:
            print(nm, lam, "log forced", flush=True); continue
        print(nm, "lam", lam, "free", free, "min order:", T.min_annihilator_order([(ser, int(lam) + 1800)], 4, 120), flush=True)
