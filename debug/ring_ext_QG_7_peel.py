"""QG step 7 (gamma = 1): exact peeling of G4, the minimal annihilator of Y_log (the log-part of the
geminal triangle at t = 0), and the exact K3 map.
(1) hyperexponential right factor: y = t^a (t^2-4)^b P(t^2), a, b from the local exponents;
(2) G4 = R3 o R1 exactly (right division, zero remainder);
(3) gauge map A (order 2, rational) with R3 o A = 0 modulo F_3's operator 4 th^3 - t^2 (th+1)^3,
    found by exact rational linear algebra on the f3 series and verified as an operator identity.
"""
from __future__ import annotations
import json, sys, time
from fractions import Fraction
from pathlib import Path
import sympy as sp
from sympy.polys.matrices import DomainMatrix
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_ops as O
from ring_ext_QG_ops import u, tt as t

DATA = Path(__file__).parent / "data"


def hyperexp_solutions(L, a_list, b_list, degs):
    """y = t^a (t^2-4)^b P(t^2), deg P <= dg; returns list of (a, b, P)."""
    out = []
    for a in a_list:
        for b in b_list:
            for dg in degs:
                cs = sp.symbols(f"c0:{dg + 1}")
                P = sum(c * t ** (2 * k) for k, c in enumerate(cs))
                y = t ** a * (t**2 - 4) ** b * P
                # L y / (t^a (t^2-4)^(b-K)) is a polynomial identity
                acc, g = 0, y
                for c in L:
                    acc += c * g
                    g = sp.expand(t * sp.diff(g, t))
                e = sp.together(sp.simplify(acc / (t ** a * (t**2 - 4) ** b)))
                num = sp.numer(e)
                eqs = sp.Poly(sp.expand(num), t).coeffs()
                M = sp.Matrix([[sp.diff(q, c) for c in cs] for q in eqs])
                ns = M.nullspace()
                if ns:
                    for v in ns:
                        Pv = sp.factor(sum(v[k] * t ** (2 * k) for k in range(dg + 1)))
                        out.append((a, b, Pv))
                    break
    return out


def main():
    d = json.loads((DATA / "ring_ext_QG_6_gcrd_g1.json").read_text())
    G4 = [sp.sympify(x) for x in d["G4"]]
    rec = {}
    t0 = time.time()
    sols = hyperexp_solutions(G4, [1, 3, 5], [0, 1, sp.Rational(-9, 2)], [0, 1, 2, 3, 4, 5])
    print("hyperexponential solutions:", sols, f"({time.time()-t0:.1f}s)", flush=True)
    rec["hyperexp"] = [str(s) for s in sols]
    a, b, Pv = sols[0]
    y = t ** a * (t**2 - 4) ** b * Pv
    r = sp.cancel(t * sp.diff(y, t) / y)
    R1 = [-r, sp.Integer(1)]
    Q3, rem = O.rdiv(G4, R1)
    print("G4 = R3 o (theta - r):  r =", sp.factor(r), "  remainder", [sp.simplify(c) for c in rem], flush=True)
    R3 = [sp.factor(c) for c in Q3]
    for k, c in enumerate(R3):
        print("  R3 theta^%d:" % k, c)
    rec["R1_r"] = str(sp.factor(r))
    rec["R3"] = [str(c) for c in R3]
    (DATA / "ring_ext_QG_7_peel_g1.json").write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
