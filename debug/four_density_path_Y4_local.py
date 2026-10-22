"""FP step Y4: logarithms at t = 1 and t = 3 for Z2 and Z2' (order 2), exactly.
theta-form L = sum_k c_k(t) theta^k; in s = t - t0, theta s^e = t0 e s^(e-1) + e s^e.  For exponents r0 < r1 = r0 + k
(integer difference) the solution s^r0 sum a_n s^n is log-free iff the recursion is consistent at n = k."""
import json
import sys
from fractions import Fraction
from pathlib import Path

import sympy as sp

HERE = Path(__file__).parent
u, t, s = sp.symbols("u t s")


def theta_form(Q):
    Q = [sp.sympify(q).subs(sp.Symbol("u"), u) for q in Q]
    n = max(sp.Poly(q, u).degree() for q in Q) + 1
    L = [sp.Integer(0)] * n
    for i, q in enumerate(Q):
        for k, c in enumerate(sp.Poly(sp.expand(q), u).all_coeffs()[::-1]):
            L[k] += c * t ** (2 * i)
    return [sp.expand(c) for c in L]


def F(x):
    x = sp.Rational(x)
    return Fraction(int(x.p), int(x.q))


def apply(cs, y, rho, t0, top):
    """y: dict offset -> coeff (exponent rho + offset); cs: list of dict power -> coeff (c_k(t0 + s))."""
    out = {}
    g = dict(y)
    for k, c in enumerate(cs):
        for pc, vc in c.items():
            for n, v in g.items():
                if n + pc < top:
                    out[n + pc] = out.get(n + pc, 0) + vc * v
        ng = {}
        for n, v in g.items():
            e = rho + n
            ng[n - 1] = ng.get(n - 1, 0) + t0 * e * v
            ng[n] = ng.get(n, 0) + e * v
        g = {n: v for n, v in ng.items() if v}
    return {n: v for n, v in out.items() if v}


def local(Q, t0):
    Lc = theta_form(Q)
    cs = []
    for c in Lc:
        p = sp.Poly(sp.expand(c.subs(t, t0 + s)), s)
        cs.append({m[0]: F(v) for m, v in zip(p.monoms(), p.coeffs())})
    # indicial: lowest offset of L(s^e) as a polynomial in e
    e = sp.Symbol("e")
    vals = {}
    for ev in range(-12, 13):
        r = apply(cs, {0: Fraction(1)}, Fraction(ev), Fraction(t0), 10**6)
        vals[ev] = r
    low = min(min(r) for r in vals.values() if r)
    pts = [(ev, vals[ev].get(low, 0)) for ev in vals]
    ind = sp.interpolate([(a, sp.Rational(b.numerator, b.denominator)) for a, b in pts], e)
    roots = sorted(r for r in sp.roots(sp.Poly(ind, e)).keys())
    out = {"indicial": str(sp.factor(ind)), "roots": [str(r) for r in roots]}
    if len(roots) == 2 and (roots[1] - roots[0]).is_integer:
        r0, k = F(roots[0]), int(roots[1] - roots[0])
        a = {0: Fraction(1)}
        for n in range(1, k + 1):
            cur = apply(cs, a, r0, Fraction(t0), low + n + 1)
            rhs = cur.get(low + n, 0)
            # coefficient of a_n at offset low + n: indicial(r0 + n)
            dn = Fraction(ind.subs(e, sp.Rational(r0.numerator, r0.denominator) + n))
            if n < k:
                a[n] = -rhs / dn
            else:
                out["obstruction_at_resonance"] = str(rhs)
                out["log"] = rhs != 0
    return out


def main():
    X = json.loads((HERE / "data" / "four_density_path_X.json").read_text())
    rec = {}
    for nm in ["Z2p_Q", "Z2_Q"]:
        for t0 in (1, 3):
            r = local(X[nm], t0)
            print(nm, "t0 =", t0, r, flush=True)
            rec[f"{nm}@{t0}"] = r
    (HERE / "data" / "four_density_path_Y4.json").write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
