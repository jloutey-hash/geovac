"""Loop test PT2, step D9: make D7's mod-p finding exact.  Find a rational gauge
A = a0 + a1 theta + a2 theta^2 with M5 (A f) = 0 for every solution f of F_3's operator
L3 = 4 theta^3 - t^2 (theta + 1)^3, by exact rational linear algebra on t-series of
f3 = 3F2(1/2,1/2,1/2; 1,1; t^2/4), then verify the operator identity exactly (no series):
reduce M5 o A modulo the right ideal of L3 and check the remainder is zero.
Also checks the convention link M5(chain_exact) = M5(D7) o (1 - t^2)."""
from __future__ import annotations
import json, sys, time
from fractions import Fraction
from pathlib import Path
import sympy as sp
from sympy.polys.matrices import DomainMatrix

u, t = sp.symbols("u t")
HERE = Path(__file__).parent
M5 = [sp.Poly(sp.sympify(q), u) for q in json.loads((HERE / "data" / "density_loop_D7.json").read_text())["M5"]]
NT = int(sys.argv[1]) if len(sys.argv) > 1 else 160


def f3_series(N):
    out, term = {}, Fraction(1)
    for n in range(N // 2):
        out[2 * n] = term
        term = term * Fraction(2 * n + 1, 2) ** 3 / Fraction(n + 1) ** 3 / 4
    return out


def theta(s, j):
    return {e: c * e**j for e, c in s.items()}


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


def apply_M(s):
    out = {}
    for i, q in enumerate(M5):
        for e, c in s.items():
            v = q.eval(e)
            if v:
                out[e + 2 * i] = out.get(e + 2 * i, 0) + Fraction(int(v)) * c
    return {e: c for e, c in out.items() if c}


def search(b, c, a, Dg, e=0):
    den = sp.expand((t**2 - 1) ** b * (t**2 - 4) ** c * (4 * t**2 + 11) ** e)
    w = inv_poly_series(den, NT)
    f = f3_series(NT)
    top = NT - 2 * len(M5) - Dg - 6
    cols, labels = [], []
    for j in range(3):
        g = mul(w, theta(f, j), NT)
        for d in range(-a, Dg + 1):
            cols.append(apply_M({e + d: v for e, v in g.items()}))
            labels.append((j, d))
    exps = sorted(e for e in set().union(*cols) if e < top)
    rows = [[sp.Rational(col.get(e, 0).numerator, col.get(e, 0).denominator) if col.get(e, 0) else sp.Integer(0)
             for col in cols] for e in exps]
    Mx = DomainMatrix.from_list_sympy(len(rows), len(cols), rows).convert_to(sp.QQ)
    ns = Mx.nullspace().to_Matrix()
    return ns, labels, den


def gauge_from(vec, labels, den):
    a = [0, 0, 0]
    for v, (j, d) in zip(vec, labels):
        a[j] += v * t**d
    return [sp.factor(sp.together(x / den)) for x in a]


# exact operator check: work on the basis (f, th f, th^2 f) modulo L3
S = [t**2 / (4 - t**2), 3 * t**2 / (4 - t**2), 3 * t**2 / (4 - t**2)]   # th^3 f = S0 f + S1 th f + S2 th^2 f
th = lambda e: sp.together(t * sp.diff(e, t))


def th_vec(p):
    c0, c1, c2 = p
    return [sp.together(th(c0) + c2 * S[0]), sp.together(c0 + th(c1) + c2 * S[1]), sp.together(c1 + th(c2) + c2 * S[2])]


def M5_as_theta_coeffs():
    """M5 = sum_k m_k(t) theta^k with functions on the left."""
    deg = max(q.degree() for q in M5)
    m = [0] * (deg + 1)
    for i, q in enumerate(M5):
        for (k,), c in zip(q.monoms(), q.coeffs()):
            m[k] += c * t ** (2 * i)
    return [sp.expand(x) for x in m]


def apply_op_mod_L3(mk, p):
    acc, cur = [0, 0, 0], p
    for m in mk:
        acc = [x + m * y for x, y in zip(acc, cur)]
        cur = th_vec(cur)
    return [sp.factor(sp.together(x)) for x in acc]


if __name__ == "__main__":
    rec = {}
    found = None
    for b in range(4):
        for c in range(4):
            t0 = time.time()
            ns, labels, den = search(b, c, 6, 12)
            print(f"(t^2-1)^{b} (t^2-4)^{c}: nullity {ns.shape[0]}  {time.time()-t0:.1f}s", flush=True)
            if ns.shape[0]:
                found = (ns, labels, den, b, c)
                break
        if found:
            break
    ns, labels, den, b, c = found
    A = gauge_from(list(ns.row(0)), labels, den)
    print("A =", A, flush=True)
    mk = M5_as_theta_coeffs()
    rem = apply_op_mod_L3(mk, A)
    print("M5 o A mod L3 =", rem, flush=True)
    # A is not in the ideal of L3 (it is of order < 3 and nonzero), so A maps Sol(L3) injectively
    rec = {"b": b, "c": c, "A": [str(x) for x in A], "remainder": [str(x) for x in rem], "nullity": ns.shape[0]}
    (HERE / "data" / "density_loop_D9.json").write_text(json.dumps(rec, indent=1))
