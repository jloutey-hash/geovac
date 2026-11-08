"""Loop test PT2, step D6: peel the chain control's operator.  For M6 = L_C / (theta - 2) and for the
adjoints, list the Frobenius solutions at t = 0 (log-free; resonances handled by allowing a free
coefficient when numerator and denominator both vanish) and their minimal annihilator orders."""
from __future__ import annotations
import json, sys
from fractions import Fraction
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import density_loop_C_factor as C
import density_loop_D_gauge as D
import density_loop_D4_factors as F4
u = sp.Symbol("u"); P = C.P

def frob_free(Q, lam, NT=400):
    Qp = [sp.Poly(q, u) for q in Q]
    def ev(i, x):
        val = Fraction(0)
        for (m,), c in zip(Qp[i].monoms(), Qp[i].coeffs()):
            val += Fraction(int(sp.Rational(c).p), int(sp.Rational(c).q)) * Fraction(x) ** m
        return D.fr(val)
    c = [1]; free = []
    for n in range(1, NT // 2):
        s_ = sum(ev(i, lam + 2 * (n - i)) * c[n - i] for i in range(1, len(Q)) if n - i >= 0) % P
        d = ev(0, lam + 2 * n)
        if d == 0:
            if s_ != 0:
                return None, None
            c.append(0); free.append(lam + 2 * n)
        else:
            c.append((-s_) * pow(d, P - 2, P) % P)
    return {lam + 2 * n: cn for n, cn in enumerate(c) if cn}, free

def run(name, Q):
    roots = sorted(set(sp.roots(sp.Poly(Q[0], u).as_expr(), u)))
    print(name, " t^0 coefficient:", sp.factor(Q[0]), flush=True)
    for lam in roots:
        ser, free = frob_free(Q, lam)
        if ser is None:
            print(f"   lam {lam}: log forced (no log-free solution)", flush=True); continue
        mo = D.minimal_order(None, ser, len([q for q in Q]) * 0 + 7, 10)
        print(f"   lam {lam}: log-free solution, free slots {free}, minimal order {mo}", flush=True)

B = json.loads((Path(__file__).parent / "data" / "density_loop_B.json").read_text())
QC = [sp.sympify(q) for q in B["QC"]]
M6 = [sp.cancel(q / (u - 2)) for q in QC]
run("M6 = L_C/(th-2)", M6)
run("L_C* (adjoint)", F4.adjoint(QC))
QG = [sp.sympify(q) for q in B["QG"]]
run("L_G (triangle, for comparison)", QG)
run("L_G* (triangle adjoint)", F4.adjoint(QG))
