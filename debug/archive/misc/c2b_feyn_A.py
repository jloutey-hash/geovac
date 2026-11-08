"""CF1 (docs/predictions/2026-10-22_c2b_feynman.md, committed d68e6fe5 before computing).

H(t) = int_0^oo K0(t x) K0(x) exp(-g x) dx.  L_Leg: Legendre of degree -1/2 in z = (t^2 + 1 - g^2)/(2t),
gauged by t^(-1/2) (H = t^(-1/2) w(z)).
 (0) symbolic check of the hand algebra j = 64 (z^2+3)^3/(1-z^2)^2 = J1 at g = 3/2.
 (1) L_Leg in d/dt form with polynomial coefficients; numerics: residual on H at four t (one below 1/2),
     with H, H', H'' by quadrature (derivatives under the integral).  Controls: 1/200 perturbation,
     no twist, g = 0 (F_2).
 (2) route 2: Mellin recurrence of g(x) = K0(x) exp(-g x), shift-2 elimination, operator in theta, t^2;
     exact right division by L_Leg.
Usage: python debug/c2b_feyn_A.py"""
import json
from pathlib import Path

import mpmath as mp
import sympy as sp

t, x, s, u = sp.symbols("t x s u")
DATA = Path(__file__).parent / "data"
OUT = {}


def check_j():
    g = sp.Rational(3, 2)
    z2 = (x + 1 - g**2) ** 2 / (4 * x)
    j = sp.cancel(64 * (z2 + 3) ** 3 / (1 - z2) ** 2)
    P = 16 * x**2 + 152 * x + 25
    J1 = P**3 / (x * (4 * x - 1) ** 2 * (4 * x - 25) ** 2)
    ok = sp.cancel(j - J1) == 0
    lam = (g**2 - (t - 1) ** 2) / (4 * t)
    kallen = (t**2 - (1 + g) ** 2) * (t**2 - (1 - g) ** 2)
    ok2 = sp.cancel(lam * (1 - lam) + kallen / (16 * t**2)) == 0
    print("j = J1:", ok, "  lambda(1-lambda) = -Kallen/(16 t^2):", ok2)
    return ok and ok2


def leg_op(g, twist=True, nu_shift=0):
    """coefficients [c0, c1, c2] (polynomials in t) of c2 H'' + c1 H' + c0 H for H = t^(-1/2) w(z(t))"""
    h0, h1, h2 = sp.symbols("h0 h1 h2")
    z = (t**2 + 1 - g**2) / (2 * t)
    rt = sp.sqrt(t)
    if twist:
        W, W1, W2 = rt * h0, h0 / (2 * rt) + rt * h1, -h0 / (4 * t * rt) + h1 / rt + rt * h2
    else:
        W, W1, W2 = h0, h1, h2
    zp, zpp = sp.diff(z, t), sp.diff(z, t, 2)
    w1 = W1 / zp
    w2 = (W2 - zpp * w1) / zp**2
    nn = -sp.Rational(1, 4) * (1 + nu_shift)
    e = sp.together(sp.expand(((1 - z**2) * w2 - 2 * z * w1 + nn * W) / (rt if twist else 1)))
    num = sp.expand(sp.numer(e))
    cs = [num.coeff(h) for h in (h0, h1, h2)]
    g0 = sp.gcd_list([c for c in cs if c != 0])
    return [sp.expand(sp.cancel(c / g0)) for c in cs]


def Hvals(tv, g):
    f0 = lambda y: mp.besselk(0, tv * y) * mp.besselk(0, y) * mp.exp(-g * y)
    f1 = lambda y: -y * mp.besselk(1, tv * y) * mp.besselk(0, y) * mp.exp(-g * y)
    f2 = lambda y: y**2 * (mp.besselk(0, tv * y) + mp.besselk(1, tv * y) / (tv * y)) * mp.besselk(0, y) * mp.exp(-g * y)
    return [mp.quad(f, [0, 1, 4, 16, mp.inf]) for f in (f0, f1, f2)]


def residual(cs, g, tv):
    vals = Hvals(tv, g)
    terms = [sp.lambdify(t, c, "mpmath")(tv) * v for c, v in zip(cs, vals)]
    return abs(sum(terms)) / max(abs(v) for v in terms)


def mellin_route(g):
    """operator sum_k t^(2k) p_k(theta) annihilating H, from G(sig)'s recurrence
    sig^2 G(sig) - g (2 sig + 1) G(sig+1) + (g^2 - 1) G(sig+2) = 0, G = M[K0(x) e^(-g x)]."""
    sig = -3 - s
    a = [(sig + i) ** 2 for i in range(3)]
    b = [-g * (2 * (sig + i) + 1) for i in range(3)]
    c = [g**2 - 1] * 3
    q0 = -a[0] * a[1] * b[2]                                  # G(sig)   = G(-3-s) -> h(s+4)/A(s+4)
    q2 = b[0] * b[1] * b[2] - a[1] * b[2] * c[0] - a[2] * b[0] * c[1]   # G(sig+2) -> h(s+2)/A(s+2)
    q4 = -b[0] * c[1] * c[2]                                  # G(sig+4) = G(1-s) -> h(s)/A(s)
    # times A(s+4) = (s+2)^2 s^2 A(s), A(s+2) = s^2 A(s): coefficients of h(s+4), h(s+2), h(s)
    Q = {2: sp.expand(q0), 1: sp.expand(q2 * (s + 2) ** 2), 0: sp.expand(q4 * (s + 2) ** 2 * s**2)}
    # q_k(s) h(s+2k) = M[t^(2k) p_k(theta) H] with p_k(u) = q_k(-u-2k)
    return {k: sp.factor(Q[k].subs(s, -u - 2 * k)) for k in Q}


def theta_to_d(ops):
    """sum_k t^(2k) p_k(theta) -> d/dt coefficients"""
    f = sp.Function("f")
    expr = 0
    for k, p in ops.items():
        pp = sp.Poly(p, u)
        th = f(t)
        acc = 0
        powers = [f(t)]
        for _ in range(pp.degree()):
            powers.append(sp.expand(t * sp.diff(powers[-1], t)))
        acc = sum(pp.coeff_monomial(u**j) * powers[j] for j in range(pp.degree() + 1))
        expr += t ** (2 * k) * acc
    expr = sp.expand(expr)
    cs = []
    for k in range(12):
        d = sp.Derivative(f(t), (t, k)) if k else f(t)
        cs.append(sp.expand(expr.coeff(d)))
    while cs and cs[-1] == 0:
        cs.pop()
    return cs


def compose_D(a, j):
    """D^j o (sum_i a_i D^i) as coefficient list"""
    cur = list(a)
    for _ in range(j):
        new = [0] * (len(cur) + 1)
        for i, ai in enumerate(cur):
            new[i] += sp.diff(ai, t)
            new[i + 1] += ai
        cur = [sp.cancel(v) for v in new]
    return cur


def right_divide(L, R):
    L = [sp.cancel(v) for v in L]
    m = len(R) - 1
    while len(L) - 1 >= m and any(v != 0 for v in L):
        k = len(L) - 1
        c = sp.cancel(L[-1] / R[-1])
        CR = compose_D(R, k - m)
        L = [sp.cancel(L[i] - c * (CR[i] if i < len(CR) else 0)) for i in range(len(L))]
        while L and L[-1] == 0:
            L.pop()
    return L   # remainder


if __name__ == "__main__":
    mp.mp.dps = 50
    OUT["j_algebra"] = check_j()
    g = sp.Rational(3, 2)
    L = leg_op(g)
    print("L_Leg (d/dt):", [sp.factor(c) for c in L])
    OUT["L_Leg"] = [str(c) for c in L]
    gm = mp.mpf(3) / 2
    pts = [mp.mpf("0.3"), mp.mpf("0.8"), mp.mpf("1.7"), mp.mpf("3.2")]
    res = {str(p): float(residual(L, gm, p)) for p in pts}
    print("CF1 residuals:", res)
    Lp = leg_op(g, nu_shift=sp.Rational(1, 200))
    resp = {str(p): float(residual(Lp, gm, p)) for p in pts}
    print("1/200 perturbation:", resp)
    Ln = leg_op(g, twist=False)
    resn = {str(p): float(residual(Ln, gm, p)) for p in pts}
    print("no twist:", resn)
    L0 = leg_op(sp.Integer(0))
    res0 = {str(p): float(residual(L0, mp.mpf(0), p)) for p in pts}
    print("gamma = 0 control (F_2):", res0)
    OUT.update({"res": res, "res_perturbed": resp, "res_no_twist": resn, "res_gamma0": res0})
    ops = mellin_route(g)
    print("route 2 theta-form:", ops)
    LM = theta_to_d(ops)
    print("route 2 order:", len(LM) - 1)
    rem = right_divide(LM, L)
    print("route 2 operator mod L_Leg: remainder zero:", len(rem) == 0)
    OUT["route2_order"] = len(LM) - 1
    OUT["route2_divisible"] = len(rem) == 0
    json.dump(OUT, open(DATA / "c2b_feyn_A.json", "w"), indent=1)
