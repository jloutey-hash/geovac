"""ND1 + ND2 of docs/predictions/2026-10-27_kite_ninth_point_origin.md.
ND1: on-shell system in explicit 2D components (w = (om, 0), k = (k1, k2), q = (q1, q2), q^2 = 0), saturated by k2
     (non-collinear), lex Groebner basis eliminating k1, k2, q1, q2; eliminant in om, then om^2 -> z. No Gram identity.
     At P2 and at two seeded random rational rate points R1, R2.
ND2: at P2, the exact configuration at the ND1 root, then the loop equations for alpha (4 scalar equations), exact
     nullspace, U(alpha*), F(alpha*) reported.
usage: python nq_nd12.py nd1 | nd2"""
import random
import sys
import time

import sympy as sp

k1, k2, q1, q2, om, sat, z = sp.symbols("k1 k2 q1 q2 om sat z")


def nd1(A, B, G1, G2, label):
    t0 = time.time()
    eqs = [k1**2 + k2**2 + A, (om - k1)**2 + k2**2 + B, (k1 - q1)**2 + (k2 - q2)**2 + G1,
           (om - k1 + q1)**2 + (q2 - k2)**2 + G2, q1**2 + q2**2, sat * k2 - 1]
    gb = sp.groebner(eqs, sat, k1, k2, q1, q2, om, order="lex")
    elim = [g for g in gb.exprs if g.free_symbols <= {om}]
    model = z * (A - G1) * (B - G2) - (A * G2 - B * G1) * (A - B - G1 + G2)
    cut12 = sp.expand((z + A + B)**2 - 4 * A * B)
    print(f"[{label}] rates (A,B,G1,G2) = {(A, B, G1, G2)}  time {time.time() - t0:.1f}s", flush=True)
    for e in elim:
        ez = sp.Poly(e, om)
        assert all(d % 2 == 0 for (d,) in ez.monoms()), "eliminant not even in om"
        pz = sp.Poly(sum(c * z**(d // 2) for (d,), c in ez.terms()), z)
        fl = sp.factor_list(pz.as_expr())
        print(f"[{label}] eliminant in z factors:", fl, flush=True)
        others = []
        for f, mult in fl[1]:
            if sp.rem(sp.expand(model), f, z) == 0 and sp.degree(f, z) == 1:
                print(f"[{label}]   model factor present, power {mult}", flush=True)
            elif sp.rem(cut12, f, z) == 0:
                print(f"[{label}]   {{1,2}} cut factor {f}, power {mult}", flush=True)
            else:
                others.append((f, mult))
        print(f"[{label}]   other factors: {others}", flush=True)
    return gb


def nd2():
    A, B, G1, G2 = sp.Integer(1), sp.Rational(9, 4), sp.Rational(49, 9), sp.Rational(64, 9)
    z9 = sp.Rational(-111, 1120)
    w0 = sp.sqrt(z9)                      # om = i sqrt(111/1120)
    K1 = sp.simplify((w0**2 - A + B) / (2 * w0))   # from p^2 = z - 2 om k1 + k^2 = -B
    K2 = sp.sqrt(-A - K1**2)
    kv = sp.Matrix([K1, K2])
    wv = sp.Matrix([w0, 0])
    pv = wv - kv
    # q from the two on-shell differences (linear, with q^2 = 0 assumed only via these): k.q = (G1 - A)/2, p.q = (B - G2)/2
    Q1, Q2 = sp.symbols("Q1 Q2")
    qv = sp.Matrix([Q1, Q2])
    sol = sp.solve([kv.dot(qv) - (G1 - A) / 2, pv.dot(qv) - (B - G2) / 2], [Q1, Q2], dict=True)[0]
    qv = sp.Matrix([sp.simplify(sol[Q1]), sp.simplify(sol[Q2])])
    checks = {"k^2+A": kv.dot(kv) + A, "p^2+B": pv.dot(pv) + B, "q^2": qv.dot(qv),
              "(k-q)^2+G1": (kv - qv).dot(kv - qv) + G1, "(p+q)^2+G2": (pv + qv).dot(pv + qv) + G2}
    for n, v in checks.items():
        print("[ND2] on-shell residual", n, "=", sp.nsimplify(sp.simplify(sp.expand(v))), flush=True)
    al = sp.symbols("al1:6")
    L1 = al[0] * kv - al[1] * pv + al[2] * (kv - qv) - al[3] * (pv + qv)
    L2 = -al[2] * (kv - qv) + al[3] * (pv + qv) + al[4] * qv
    M = sp.Matrix([[sp.diff(e, a) for a in al] for e in list(L1) + list(L2)])
    M = M.applyfunc(lambda e: sp.radsimp(sp.simplify(e)))
    ns = M.nullspace(simplify=True)
    print("[ND2] rank of loop-equation matrix:", M.rank(simplify=True), " nullspace dim:", len(ns), flush=True)
    for v in ns:
        v = v / v[4] if v[4] != 0 else v
        v = v.applyfunc(lambda e: sp.radsimp(sp.simplify(e)))
        print("[ND2] alpha* =", list(v), flush=True)
        print("[ND2] alpha* numeric:", [sp.N(c, 20) for c in v], flush=True)
        a1, a2, a3, a4, a5 = v
        U = (a1 + a2) * (a3 + a4) + a5 * (a1 + a2 + a3 + a4)
        print("[ND2] U(alpha*) =", sp.nsimplify(sp.simplify(sp.expand(U))), " numeric", sp.N(U, 30), flush=True)
        # F(alpha*) from the same construction as nq_engine.build_G (massless, at z9)
        import nq_engine as E
        from fractions import Fraction as Fr
        Gd = E.build_G()
        Fv = 0
        for m, cd in Gd.items():
            if sum(m) != 3:
                continue
            c = sum(sp.Rational(cf.numerator, cf.denominator) * z9**dz for (dz, de), cf in cd.items() if de == 0)
            Fv += c * sp.prod([v[j]**m[j] for j in range(5)])
        print("[ND2] F(alpha*) (reported, not scored) numeric:", sp.N(Fv, 30), flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "nd1":
        rng = random.Random(20261027)
        def rr():
            return sp.Rational(rng.randint(1, 40), rng.randint(1, 12))**2
        pts = [("P2", (sp.Integer(1), sp.Rational(9, 4), sp.Rational(49, 9), sp.Rational(64, 9)))]
        pts += [("R1", tuple(rr() for _ in range(4))), ("R2", tuple(rr() for _ in range(4)))]
        for lab, r in pts:
            nd1(*r, lab)
    else:
        nd2()
