"""Minimality audit, step 2: does the triangle G really carry the elliptic piece?

(theta-2) o L_G annihilates G (L_G G = -3 t^2 on the actual function), L_G = M o (theta - r),
M equivalent to Legendre (irreducible).  The rung claim 'G sits on the elliptic rung' needs M to
be a composition factor of G's MINIMAL operator.  It is not iff G is Liouvillian.  With
h = (theta - r) G one has M h = -3 t^2, and h is Liouvillian iff O' = (theta - 2) o M has a
hyperexponential solution y with M y != 0 (its Liouvillian solutions form a Galois-invariant
subspace meeting ker M trivially, so of dimension <= 1).  For a Fuchsian O' a hyperexponential
solution is prod_p (t-p)^{e_p} P(t), e_p a local exponent, P a polynomial of degree
-e_oo - sum e_p.  Exhaustive search, no degree bound.

Then a numerical cross-check on the function: h(t) = (theta - r) G at several t, 25+ digits,
and the residual M h + 3 t^2.
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import mpmath as mp
import sympy as sp

sys.path.insert(0, str(Path(__file__).parent))
from minimality_audit_fuchs import indicial, roots_multiset  # noqa: E402

t, u, lam = sp.symbols("t u lambda")
OUT = Path(__file__).parent / "data" / "minimality_audit_ring.json"
QG = [-4 * u * (u - 3) ** 2, 5 * u**3 - 7 * u**2 - 5 * u + 3, -u * (u + 2) * (u + 3)]


def _th(c):
    return sp.together(t * sp.diff(c, t))


def coeff_list(Q):
    n = max(sp.Poly(q, u).degree() for q in Q) + 1
    out = [sp.Integer(0)] * n
    for i, q in enumerate(Q):
        for k, c in enumerate(sp.Poly(sp.expand(q), u).all_coeffs()[::-1]):
            out[k] += c * t ** (2 * i)
    return [sp.expand(c) for c in out]


def compose(A, B):
    out = [0] * (len(A) + len(B) - 1)
    for k, a in enumerate(A):
        for j, b in enumerate(B):
            bi = b
            for i in range(k + 1):
                out[k - i + j] += a * sp.binomial(k, i) * bi
                bi = _th(bi)
    return [sp.together(c) for c in out]


def apply(Lc, f):
    acc, g = 0, f
    for c in Lc:
        acc += c * g
        g = _th(g)
    return sp.simplify(acc)


def theta_list_to_D(L):
    """function-left theta coefficients -> polynomial D-coefficients (denominators cleared)."""
    n = len(L) - 1
    a = [sp.Integer(0)] * (n + 1)
    for k, c in enumerate(L):
        for j in range(k + 1):
            a[j] += c * sp.functions.combinatorial.numbers.stirling(k, j) * t**j
    den = sp.lcm([sp.denom(sp.together(x)) for x in a])
    return [sp.expand(sp.cancel(x * den)) for x in a]


def divide_LG():
    LG = coeff_list(QG)
    y = (1 - t**2 / 4) ** sp.Rational(-3, 2)
    r = sp.together(_th(y) / y)
    R = [-r, sp.Integer(1)]
    rem = list(LG)
    M = [0, 0, 0]
    for k in (3, 2, 1):
        c = sp.factor(sp.together(rem[k]))
        M[k - 1] = c
        sub = compose([0] * (k - 1) + [c], R)
        rem = [sp.together(rem[i] - (sub[i] if i < len(sub) else 0)) for i in range(4)]
    assert all(sp.simplify(c) == 0 for c in rem)
    return LG, r, [sp.factor(c) for c in M]


def hyperexp_search(Lc):
    a = theta_list_to_D(Lc)
    n = len(a) - 1
    lead = sp.Poly(a[-1], t)
    fin = sorted(sp.roots(lead).keys(), key=lambda z: float(z))
    exps = {p: roots_multiset(indicial(a, p)) for p in fin + [sp.oo]}
    sols = []
    tried = 0
    for combo in itertools.product(*[sorted(set(exps[p])) for p in fin]):
        for einf in sorted(set(exps[sp.oo])):
            N = -einf - sum(combo)
            if not (N.is_integer and N >= 0):
                continue
            tried += 1
            cs = sp.symbols(f"c0:{int(N) + 1}")
            P = sum(c * t**i for i, c in enumerate(cs))
            base = sp.Mul(*[(t - p) ** e for p, e in zip(fin, combo)])
            yv = base * P
            expr = sum(aj * sp.diff(yv, t, j) for j, aj in enumerate(a))
            expr = sp.expand(sp.simplify(expr / base * sp.Mul(*[(t - p) ** n for p in fin])))
            eqs = sp.Poly(sp.numer(sp.together(expr)), t).coeffs()
            solv = sp.linsolve(eqs, cs)
            for sol in solv:
                if any(v != 0 for v in sol):
                    Pp = sp.expand(P.subs(dict(zip(cs, sol))))
                    free = sorted(Pp.free_symbols - {t}, key=str)
                    for f in free:
                        Pk = Pp.subs({g: (1 if g == f else 0) for g in free})
                        sols.append(sp.factor(base * Pk))
                    if not free:
                        sols.append(sp.factor(base * Pp))
    return exps, tried, sols


def numeric_h(ts, dps=30):
    """h = (theta - r) G on the momentum form T = (1/pi) int dp dq D3(p) D3(q) [k(p+q) - k(|p-q|)],
    k(w) = t w / sqrt(w^2 + t^2) (scales (1,1,t)), theta acting analytically on k."""
    mp.mp.dps = dps
    w = sp.Symbol("w", positive=True)
    kexpr = t * w / sp.sqrt(w**2 + t**2)
    kG = sp.lambdify((t, w), kexpr, "mpmath")
    kth = sp.lambdify((t, w), sp.simplify(t * sp.diff(kexpr, t)), "mpmath")
    D3 = lambda k: (k * k + 1) ** mp.mpf(-1.5)
    res = {}
    for tv in ts:
        tv = mp.mpf(tv)
        def F(kf):
            kk = lambda x: kf(tv, x) if x > 0 else mp.mpf(0)
            inner = lambda p: mp.quad(lambda q: D3(q) * (kk(p + q) - kk(abs(p - q))), [0, p, 2 * p + 1, mp.inf])
            return mp.quad(lambda p: D3(p) * inner(p), [0, 1, 4, mp.inf]) / mp.pi
        G, thG = F(kG), F(kth)
        rr = 3 * tv**2 / (4 - tv**2)
        res[str(tv)] = {"G": G, "thetaG": thG, "h": thG - rr * G}
    return res


if __name__ == "__main__":
    LG, r, M = divide_LG()
    print("M (function-left theta coefficients):", M)
    Op = compose([sp.Integer(-2), sp.Integer(1)], M)
    exps, tried, sols = hyperexp_search(Op)
    print("O' = (theta-2) o M exponents:", {str(k): [str(x) for x in v] for k, v in exps.items()})
    print("candidate exponent combinations with integer degree:", tried)
    print("hyperexponential solutions of O':", sols)
    imgs = [sp.simplify(apply(M, y)) for y in sols]
    print("M applied to them:", imgs)
    # positive control: the same search on L_G itself must find the algebraic solution
    exL, trL, solL = hyperexp_search(LG)
    print("control, hyperexponential solutions of L_G:", solL)
    # second control: Legendre theta^2 - t^2 (theta+1)^2 has none
    exLe, trLe, solLe = hyperexp_search([-t**2, -2 * t**2, 1 - t**2])
    print("control, hyperexponential solutions of Legendre:", solLe)
    rec = {"M": [str(c) for c in M], "Oprime_exponents": {str(k): [str(x) for x in v] for k, v in exps.items()},
           "combos_tried": tried, "hyperexp_Oprime": [str(s) for s in sols], "M_images": [str(i) for i in imgs],
           "control_LG": [str(s) for s in solL], "control_legendre": [str(s) for s in solLe]}
    OUT.write_text(json.dumps(rec, indent=1))
