"""QG step 7 (gamma given): exact K3 piece of G4 (the minimal annihilator of the geminal triangle Y).
(1) exact rational A = sum_{j<3} a_j theta^j, a_j = p_j(t)/den, with G4* (A f3) = 0 on the
    f3 = 3F2(1/2,1/2,1/2;1,1;t^2/4) series (exact rational linear algebra, size fixed by the mod-p scan);
(2) operator identity: G4* o A = 0 modulo L3 = 4 th^3 - t^2 (th+1)^3, in exact arithmetic (no series);
(3) R = annihilator of A Sol(L3) (order 3), G4* = R1' o R exactly; so G4 = R* o R1'*, R* equivalent to
    F_3's operator (the K3 rung) and R1'* of order 1;
(4) the order-1 factor's solution in closed form.
usage: python ring_ext_QG_7_k3map.py <gamma> <den-expr> <Dg> <a>"""
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
NT = 200


def f3_series(N):
    out, term = {}, Fraction(1)
    for n in range(N // 2):
        out[2 * n] = term
        term = term * Fraction(2 * n + 1, 2) ** 3 / Fraction(n + 1) ** 3 / 4
    return out


def mul(a, b, top):
    out = {}
    for e1, c1 in a.items():
        for e2, c2 in b.items():
            if e1 + e2 < top:
                out[e1 + e2] = out.get(e1 + e2, 0) + c1 * c2
    return {e: c for e, c in out.items() if c}


def inv_poly_series(poly, N):
    cs = [Fraction(int(sp.Rational(c).p), int(sp.Rational(c).q)) for c in sp.Poly(poly, t).all_coeffs()[::-1]]
    out = [1 / cs[0]]
    for n in range(1, N):
        out.append(-sum(cs[k] * out[n - k] for k in range(1, min(n, len(cs) - 1) + 1)) / cs[0])
    return {n: c for n, c in enumerate(out) if c}


def apply_Q(Q, s):
    Qp = [sp.Poly(q, u) for q in Q]
    out = {}
    for i, q in enumerate(Qp):
        for e, c in s.items():
            v = q.eval(e)
            if v:
                out[e + 2 * i] = out.get(e + 2 * i, 0) + Fraction(int(sp.Rational(v).p), int(sp.Rational(v).q)) * c
    return {e: c for e, c in out.items() if c}


def theta_s(s, j):
    return {e: c * e**j for e, c in s.items()}


def main():
    g = sys.argv[1]
    den = sp.sympify(sys.argv[2], locals={"t": t})
    Dg, a = int(sys.argv[3]), int(sys.argv[4])
    d = json.loads((DATA / f"ring_ext_QG_6_gcrd_g{g}.json").read_text())
    QG4 = [sp.sympify(q) for q in d["QG4"]]
    G4 = [sp.sympify(c) for c in d["G4"]]
    adjQ = [sp.expand(q.subs(u, -u - 2 * i)) for i, q in enumerate(QG4)]
    t0 = time.time()
    w = inv_poly_series(sp.expand(den), NT)
    f = f3_series(NT)
    top = NT - 2 * len(QG4) - Dg - 6
    cols, labels = [], []
    for j in range(3):
        gj = mul(w, theta_s(f, j), NT)
        for dd in range(-a, Dg + 1):
            cols.append(apply_Q(adjQ, {e + dd: v for e, v in gj.items()}))
            labels.append((j, dd))
    exps = sorted(e for e in set().union(*cols) if e < top)
    rows = [[sp.Rational(col.get(e, 0).numerator, col.get(e, 0).denominator) if col.get(e, 0) else sp.Integer(0)
             for col in cols] for e in exps]
    M = DomainMatrix.from_list_sympy(len(rows), len(cols), rows).convert_to(sp.QQ)
    ns = M.nullspace().to_Matrix()
    print("exact nullity", ns.shape[0], f"({time.time()-t0:.1f}s, {len(rows)} rows x {len(cols)} cols)", flush=True)
    vec = list(ns.row(0))
    A = [sp.Integer(0)] * 3
    for v, (j, dd) in zip(vec, labels):
        A[j] += v * t**dd
    A = [sp.factor(sp.together(x / den)) for x in A]
    print("A =", A, flush=True)
    # (2) operator identity modulo L3: theta^3 f = S0 f + S1 th f + S2 th^2 f
    L3 = [sp.Integer(0), sp.Integer(0), sp.Integer(0), sp.Integer(0)]
    L3 = O.Q_to_coeffs([4 * u**3, -(u + 1) ** 3])
    S = [sp.cancel(-c / L3[3]) for c in L3[:3]]

    def th_vec(p):
        c0, c1, c2 = p
        return [sp.cancel(O.th_fn(c0) + c2 * S[0]), sp.cancel(c0 + O.th_fn(c1) + c2 * S[1]),
                sp.cancel(c1 + O.th_fn(c2) + c2 * S[2])]

    def apply_mod(Lc, p):
        acc, cur = [0, 0, 0], p
        for m in Lc:
            acc = [sp.cancel(x + m * y) for x, y in zip(acc, cur)]
            cur = th_vec(cur)
        return acc
    G4adj = O.Q_to_coeffs(adjQ)
    rem = apply_mod(G4adj, A)
    print("G4* o A mod L3 =", [sp.simplify(x) for x in rem], flush=True)
    # (3) R: annihilator of the image
    vecs = [A]
    for _ in range(3):
        vecs.append(th_vec(vecs[-1]))
    Mx = sp.Matrix([[vecs[k][i] for k in range(4)] for i in range(3)])
    nsv = Mx.nullspace()
    R = [sp.cancel(x / nsv[0][3]) for x in nsv[0]]
    Qo, rem2 = O.rdiv(G4adj, R)
    print("G4* = R1' o R: remainder", [sp.simplify(x) for x in rem2], " R1' order", len(Qo) - 1, flush=True)

    def adjoint_coeffs(Lc):
        """(sum c_k theta^k)* = sum (-theta)^k c_k  -> function-left coefficients."""
        out = [sp.Integer(0)] * len(Lc)
        for k, c in enumerate(Lc):
            # (-theta)^k o c
            term = O.comp([sp.Integer(0)] * k + [sp.Integer((-1) ** k)], [c])
            for i, x in enumerate(term):
                out[i] += x
        return [sp.cancel(x) for x in out]
    R1 = adjoint_coeffs(Qo)
    R3 = adjoint_coeffs(R)
    chk = O.comp(R3, R1)
    G4m = O.monic(G4)
    chk = O.monic(chk)
    print("G4 == R3 o R1 (monic):", all(sp.simplify(x - y) == 0 for x, y in zip(chk, G4m)), flush=True)
    r1 = sp.factor(sp.cancel(-R1[0] / R1[1]))
    print("R1 = theta - r,  r =", r1, flush=True)
    y = sp.exp(sp.integrate(sp.apart(r1 / t, t), t))
    y = sp.simplify(y)
    print("solution of R1:", y, flush=True)
    R3m = O.monic(R3)
    rec = {"A": [str(x) for x in A], "remainder_mod_L3": [str(sp.simplify(x)) for x in rem],
           "R": [str(x) for x in R], "R1_r": str(r1), "R1_solution": str(y), "R3_monic": [str(sp.factor(x)) for x in R3m]}
    (DATA / f"ring_ext_QG_7_k3map_g{g}.json").write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
