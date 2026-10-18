"""QG step 9b (gamma = 3/2): quadric count on the exact order-2 right factor R2 of L9 (6e), and a
reducibility check: hyperexponential solutions y = t^a (t^2-4)^b P(t^2) from its local exponents."""
import json, sys
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_9_quadrics as Qd
import ring_ext_QG_ops as O
from ring_ext_QG_ops import u, tt as t

d = json.loads((Path(__file__).parent / "data" / "ring_ext_QG_6e_order2_g1.5.json").read_text())
QR = [sp.sympify(q, locals={"u": u}) for q in d["R2_Q"]]
Lc = O.Q_to_coeffs(QR)
for t0 in (sp.Rational(1, 2), sp.Rational(7, 5)):
    print("R2 quadrics at", t0, Qd.quadrics(Lc, t0))
# hyperexponential solutions: exponents at 0 {0,-2}, at 2 {0,-3/2}, at inf {0,-3} (t^-rho with rho in {0,3})
found = []
for a in (0, -2):
    for b in (0, sp.Rational(-3, 2)):
        for dg in range(0, 4):
            cs = sp.symbols(f"c0:{dg+1}")
            y = t**a * (t**2 - 4) ** b * sum(c * t ** (2 * k) for k, c in enumerate(cs))
            e = sp.simplify(O.apply_fn(Lc, y) / (t**a * (t**2 - 4) ** b))
            num = sp.numer(sp.together(e))
            eqs = sp.Poly(sp.expand(num), t).coeffs()
            ns = sp.Matrix([[sp.diff(q, c) for c in cs] for q in eqs]).nullspace()
            if ns:
                found.append((a, b, sp.factor(sum(ns[0][k] * t ** (2 * k) for k in range(dg + 1)))))
                break
print("hyperexponential solutions of R2:", found)
