"""Paper 63 owed item (A): what are the odd-lambda pieces phi(L_m(theta - lam)) of prop:tree_rung_l?

Exact checks (sympy, Q(z)[S] recurrence algebra of the banana-dictionary sprint, debug/archive/misc):
  E1  phi(L_m(theta - 1)) = Q_m . V_m exactly, V_m = Vanhove's equal-mass 2-D banana operator with m propagators
      (arXiv:1401.6438 Table 1, s = t^2), m = 2..5; the quotient Q_m has order 2.
  E2  Q_m splits into two order-one factors: Q_m = (theta + a(t)) (theta + b(t)) with a, b rational (found by
      an exact rational right factor theta - r(t) of Q_m, r from the indicial data), m = 2..5.
  E3  the parity-transport identities used in the proof: P_{lam+2} U = U P_lam (U = t^-2 theta^2) and its
      adjoint U* P*_{lam+2} = P*_lam U*, exact, lam = -3..5, m = 2..5.
  E4  control: the even piece phi(L_m) is NOT right-divisible by V_m (remainder nonzero), m = 2..5.
Numerical (mpmath, dps 30):
  N1  phi(L_m(theta - lam)) annihilates b_m^(lam)(t) = int x^lam I0(t x) K0(x)^m dx, lam = 1, 3 (|t| < m), m = 2..4;
  N2  V_m b_m^(1) = const (the inhomogeneous banana equation): theta_s (V_m b) = 0 at two t, m = 2..4.
Output: debug/data/p63odd_factor.json
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import mpmath as mp
import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "tests"))
import test_paper63_tree_rung as TR  # noqa: E402  (operator builders _sym_op, _shift, _phi, _t_lists)

u, z, S, T = sp.symbols("u z s t")  # u = theta_t; z = Mellin variable; s = t^2

VAN = {2: ["s-2", "s*(s-4)"],
       3: ["s-3", "3*s**2-20*s+9", "s*(s-1)*(s-9)"],
       4: ["s-4", "7*s**2-68*s+64", "6*s**3-90*s**2+192*s", "s**2*(s-4)*(s-16)"],
       5: ["s-5", "(3*s-5)*(5*s-57)", "25*s**3-518*s**2+1839*s-450",
           "10*s**4-280*s**3+1554*s**2-900*s", "s**2*(s-25)*(s-1)*(s-9)"]}


# ---------------------------------------------------------------- operator forms
def P_tlist(m: int, lam: int) -> list:
    """phi(L_m(theta - lam)) as [Q_0(u), Q_1(u), ...], operator = sum_i t^{2i} Q_i(theta)."""
    Q = TR._t_lists(TR._phi(TR._shift(TR._sym_op(m), lam)))
    return [sp.expand(q.subs(TR.th, u)) for q in Q]


def van_tlist(m: int, strip: bool = True) -> list:
    qs = [sp.sympify(q, locals={"s": S}) for q in VAN[m]]
    K = len(qs) - 1
    tot = 0
    for k, q in enumerate(qs):
        ff = sp.prod([u - j for j in range(k)])  # s^k d^k/ds^k = theta_s falling factorial, theta_s = theta_t/2
        tot += sp.expand(q * S ** (K - k)) * ff
    P = sp.Poly(sp.expand(tot), S)
    Q = [sp.expand(P.coeff_monomial(S**i).subs(u, u / 2)) for i in range(P.degree() + 1)]
    while strip and Q and Q[0] == 0:
        Q.pop(0)
    return Q


def to_rec(Q: list) -> list:
    return [sp.cancel(sp.sympify(q).subs(u, -z - 2 * i)) for i, q in enumerate(Q)]


def from_rec(R: list) -> list:
    return [sp.cancel(sp.sympify(r).subs(z, -u - 2 * i)) for i, r in enumerate(R)]


def mul(a: list, b: list) -> list:
    out = [sp.Integer(0)] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] += ai * sp.sympify(bj).subs(z, z + 2 * i)
    return [sp.cancel(c) for c in out]


def rdiv(a: list, b: list):
    a = [sp.cancel(x) for x in a]
    q = [sp.Integer(0)] * max(1, len(a) - len(b) + 1)
    while a and sp.simplify(a[-1]) == 0:
        a.pop()
    while len(a) >= len(b):
        k = len(a) - len(b)
        c = sp.cancel(a[-1] / sp.sympify(b[-1]).subs(z, z + 2 * k))
        q[k] += c
        sub = mul([0] * k + [c], b)
        a = [sp.cancel(x - (sub[i] if i < len(sub) else 0)) for i, x in enumerate(a + [0] * (len(sub) - len(a)))]
        while a and sp.simplify(a[-1]) == 0:
            a.pop()
    return q, a


# ---------------------------------------------------------------- theta-form with rational coefficients in t
def tlist_to_theta(Q: list) -> dict:
    """sum_i t^{2i} Q_i(theta) -> {k: coefficient of theta^k} (coefficients polynomials in t)."""
    out: dict = {}
    for i, q in enumerate(Q):
        for (k,), c in sp.Poly(sp.expand(q), u).terms():
            out[k] = sp.expand(out.get(k, 0) + c * T ** (2 * i))
    return out


def th_mul(A: dict, B: dict) -> dict:
    """(sum a_k theta^k)(sum b_j theta^j), theta f = f theta + t f'."""
    out: dict = {}
    for k, a in A.items():
        for j, b in B.items():
            # theta^k b(t) = sum_i C(k,i) (theta_t^i b) theta^{k-i}
            bi = b
            for i in range(k + 1):
                c = sp.binomial(k, i) * bi
                out[k - i + j] = sp.cancel(out.get(k - i + j, 0) + a * c)
                bi = sp.cancel(T * sp.diff(bi, T))
    return {k: v for k, v in out.items() if sp.simplify(v) != 0}


def th_rdiv_order1(Q: dict) -> list:
    """rational r(t) with theta - r a right factor of the order-2 operator Q: Q(theta - r)-quotient remainder
    = Q applied to a 'hyperexponential' y with theta y = r y; i.e. solve the Riccati equation
    q2 (theta r + r^2) + q1 r + q0 = 0 for r = c (constant) or r = c + d t^2/(t^2 - e) shapes. Returns list of r."""
    q0, q1, q2 = (sp.sympify(Q.get(k, 0)) for k in range(3))
    sols = []
    # constant r
    c = sp.symbols("c")
    eq = sp.Poly(sp.numer(sp.together(q2 * c**2 + q1 * c + q0)), T)
    cs = sp.solve(eq.coeffs(), c, dict=True)
    for d in cs:
        sols.append(d[c])
    if sols:
        return sols
    # r = a + b t^2 / (t^2 - e): poles of r sit at singular points of Q
    a, b = sp.symbols("a b")
    lead = sp.factor(q2)
    roots = [rt for rt in sp.roots(sp.Poly(sp.numer(sp.together(lead)), T)).keys() if rt != 0]
    for e in set(rt**2 for rt in roots):
        r = a + b * T**2 / (T**2 - e)
        expr = sp.numer(sp.together(q2 * (T * sp.diff(r, T) + r**2) + q1 * r + q0))
        eqs = sp.Poly(sp.expand(expr), T).coeffs()
        for d in sp.solve(eqs, [a, b], dict=True):
            if d:
                sols.append(sp.simplify(r.subs(d)))
    return sols


# ---------------------------------------------------------------- numerics
def residual_tlist(Q: list, f) -> mp.mpf:
    """apply sum t^{2i} Q_i(theta) to f at the point (f given as callable of t), via mp.diff in log t."""
    raise NotImplementedError


def b_mp(m: int, lam: int, tv) -> mp.mpf:
    return mp.quad(lambda x: x**lam * mp.besseli(0, tv * x) * mp.besselk(0, x) ** m, [0, 1, 4, 12, 30, mp.inf])


def apply_theta_ops(Q: list, fun, tv) -> mp.mpf:
    """sum_i t^{2i} Q_i(theta_t) f at t = tv; theta^k f via mp.diff in w = log t (exact derivative orders)."""
    K = max(sp.Poly(q, u).degree() for q in Q if q != 0)
    g = lambda w: fun(mp.exp(w))  # noqa: E731
    w0 = mp.log(tv)
    ders = [mp.diff(g, w0, k) for k in range(K + 1)]
    tot = mp.mpf(0)
    scale = mp.mpf(0)
    for i, q in enumerate(Q):
        for (k,), c in sp.Poly(q, u).terms():
            term = mp.mpf(sp.Rational(c).p) / sp.Rational(c).q * tv ** (2 * i) * ders[k]
            tot += term
            scale = max(scale, abs(term))
    return tot / scale


def e3() -> dict:
    Ug = {(1, 2): sp.Integer(1)}
    res = {}
    for m in (2, 3, 4, 5):
        L = TR._sym_op(m)
        for lam in (-3, -1, 1, 3, 5):
            lhs = TR._a_mul(TR._phi(TR._shift(L, lam + 2)), Ug, -1)
            rhs = TR._a_mul(Ug, TR._phi(TR._shift(L, lam)), -1)
            res[f"{m},{lam}"] = set(lhs) == set(rhs) and all(sp.expand(lhs[k] - rhs[k]) == 0 for k in lhs)
    return res


_TH: dict = {}


def _theta_I0(k: int):
    """theta^k I0(y) = a_k(y) I0(y) + b_k(y) y I1(y), a_k, b_k polynomials (sympy -> mpmath callables)."""
    if k not in _TH:
        Y = sp.Symbol("Y")
        a, b = sp.Integer(1), sp.Integer(0)
        for _ in range(k):
            a, b = sp.expand(Y * sp.diff(a, Y) + Y**2 * b), sp.expand(a + Y * sp.diff(b, Y))
        _TH[k] = (sp.lambdify(Y, a, "mpmath"), sp.lambdify(Y, b, "mpmath"))
    return _TH[k]


def theta_b(m: int, lam: int, tv, k: int):
    a, b = _theta_I0(k)
    f = lambda x: x**lam * (a(tv * x) * mp.besseli(0, tv * x) + b(tv * x) * tv * x * mp.besseli(1, tv * x)) * mp.besselk(0, x) ** m  # noqa: E731
    return mp.quad(f, [0, 1, 4, 12, 30, mp.inf])


def apply_tlist(Q: list, m: int, lam: int, tv):
    """(|sum t^{2i} Q_i(theta) b| / largest term, raw value) for b = int x^lam I0(t x) K0^m dx."""
    K = max(sp.Poly(q, u).degree() for q in Q if q != 0)
    th = [theta_b(m, lam, tv, k) for k in range(K + 1)]
    tot, scale = mp.mpf(0), mp.mpf(0)
    for i, q in enumerate(Q):
        for (k,), c in sp.Poly(q, u).terms():
            c = sp.Rational(c)
            term = mp.mpf(c.p) / c.q * tv ** (2 * i) * th[k]
            tot += term
            scale = max(scale, abs(term))
    return abs(tot) / scale, tot


def main() -> None:
    out: dict = {"E1": {}, "E2": {}, "E3": {}, "E4": {}, "N1": {}, "N2": {}}
    t0 = time.time()
    for m in (2, 3, 4, 5):
        V = to_rec(van_tlist(m))
        P1 = to_rec(P_tlist(m, 1))
        q, r = rdiv(P1, V)
        ok = len(r) == 0
        Qth = tlist_to_theta(from_rec(q)) if ok else {}
        out["E1"][m] = {"remainder_zero": ok, "quotient_order": max(Qth) if Qth else None,
                        "quotient": {k: str(sp.factor(v)) for k, v in Qth.items()}}
        # E2: order-one right factor of the quotient, then the left cofactor is order one automatically
        rs = th_rdiv_order1(Qth) if ok else []
        out["E2"][m] = {"rational_right_factors_theta_minus_r": [str(x) for x in rs], "splits": bool(rs)}
        # E4 control
        P0 = to_rec(P_tlist(m, 0))
        _, r0 = rdiv(P0, V)
        out["E4"][m] = {"even_piece_divisible_by_V": len(r0) == 0}
        print(m, out["E1"][m]["remainder_zero"], out["E1"][m]["quotient_order"], out["E2"][m], out["E4"][m], flush=True)
    # E3: intertwiners, exact, in the u = t^-2 form of the tree_rung guard (a left t-power normalisation would
    # break the identity, so the t-list form is not used here)
    out["E3"] = e3()
    print("E3", out["E3"], flush=True)
    # numerics: theta_t^k under the integral, theta I0(y) = y I1(y), theta (y I1) = y^2 I0 (no numerical derivatives)
    mp.mp.dps = 30
    for m in (2, 3, 4):
        out["N1"][m] = {}
        for lam in (1, 3):
            vals = []
            for tv in (mp.mpf("0.6"), mp.mpf(m) - mp.mpf("0.8")):
                vals.append(float(apply_tlist(P_tlist(m, lam), m, lam, tv)[0]))
            out["N1"][m][lam] = vals
        out["N1"][m]["control_even_op_on_odd_fn"] = float(apply_tlist(P_tlist(m, 0), m, 1, mp.mpf("0.6"))[0])
        Kv = len(VAN[m]) - 1
        v1 = apply_tlist(van_tlist(m, strip=False), m, 1, mp.mpf("0.5"))[1] / mp.mpf("0.5") ** (2 * Kv)
        v2 = apply_tlist(van_tlist(m, strip=False), m, 1, mp.mpf("0.9"))[1] / mp.mpf("0.9") ** (2 * Kv)
        out["N2"][m] = {"V_b_at_t0.5": float(v1), "V_b_at_t0.9": float(v2), "rel_diff": float(abs(v1 - v2) / abs(v1))}
        print("N", m, out["N1"][m], out["N2"][m], flush=True)
    out["seconds"] = time.time() - t0
    (HERE / "data").mkdir(exist_ok=True)
    json.dump(out, open(HERE / "data" / "p63odd_factor.json", "w"), indent=1, default=str)


def _apply_unnormalized(Q: list, fun, tv) -> mp.mpf:
    """sum t^{2i} Q_i(theta) f, divided by the t-power t^{2 i_max} so a constant stays comparable: Vanhove's
    s^N-scaled operator; we return (operator f)/s^{deg}, i.e. sum_k q_k(s) d^k f/ds^k up to the s^N factor."""
    K = max(sp.Poly(q, u).degree() for q in Q if q != 0)
    g = lambda w: fun(mp.exp(w))  # noqa: E731
    w0 = mp.log(tv)
    ders = [mp.diff(g, w0, k) for k in range(K + 1)]
    tot = mp.mpf(0)
    for i, q in enumerate(Q):
        for (k,), c in sp.Poly(q, u).terms():
            tot += mp.mpf(sp.Rational(c).p) / sp.Rational(c).q * tv ** (2 * i) * ders[k]
    return tot


if __name__ == "__main__":
    main()
