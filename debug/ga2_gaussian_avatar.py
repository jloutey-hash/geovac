"""GA2 (docs/predictions/2026-10-22_geminal_shapes.md, committed 5e40da9a before computing).

H_G(t) = int_0^oo K0(tx) K0(x) exp(-x^2) dx.  Derived in the record: with g = K0 e^{-x^2},
theta^2 g + x^2 (4 theta + 3) g + 4 x^4 g = 0, so s^2 G(s) - (4s+5) G(s+2) + 4 G(s+4) = 0, and the
Mellin route gives L_G = t^4 (theta+1)^2 - t^2 theta^2 (4 theta + 1) + 4 theta^2 (theta - 2)^2.
Checks:
 (1) the x-equation on g, symbolically;
 (2) L_G annihilates H_G at four t, derivatives under the integral, 40 digits; 1/200 perturbation;
 (3) singularity type at t = infinity: the top power t^4 carries theta-degree 2 < 4 (irregular by the
     Fuchs criterion); and the monodromy companion J(t) = int I0(tx) K0(x) e^{-x^2} dx, which lies in the
     solution space (K0(e^{i pi} z) = K0(z) - i pi I0(z)), grows like e^{t^2/4}: measured slope of
     log J against t^2;
 (4) control: the same conversion for g = K0(x) cosh(gamma x), gamma = 3/2, gives a Fuchsian operator
     (top theta-degree equals the order at both ends) that the CF record's L_Leg right-divides.
Usage: python debug/ga2_gaussian_avatar.py"""
import json
import sys
from pathlib import Path

import mpmath as mp
import sympy as sp

sys.path.insert(0, str(Path(__file__).parent))
from c2b_feyn_A import leg_op, theta_to_d, right_divide, t, u  # noqa: E402

x = sp.Symbol("x", positive=True)
s = sp.Symbol("s")
DATA = Path(__file__).parent / "data"
OUT = {}


def theta_apply(expr, n):
    for _ in range(n):
        expr = sp.expand(x * sp.diff(expr, x))
    return expr


def check_x_equation():
    g = sp.besselk(0, x) * sp.exp(-x**2)
    e = theta_apply(g, 2) + x**2 * (4 * theta_apply(g, 1) + 3 * g) + 4 * x**4 * g
    e = sp.simplify(e.rewrite(sp.besselk))
    return sp.simplify(e) == 0


def recurrence_to_operator(rec, kmax):
    """rec: dict shift -> polynomial in s with sum_m rec[m](s) G(s + m) = 0 for G = M[g] (shifts 0, 2, ...).
    Returns the t-operator {k: p_k(u)} with sum_k t^(2k) p_k(theta) F = 0 for F = int K0(tx) g(x) dx."""
    # Phi(s) = A(s) G(1 - s), A(s) = 2^(s-2) Gamma(s/2)^2, A(s)/A(s-2) = (s-2)^2.
    # G(1 - s + m) = G(1 - (s - m)) = Phi(s - m)/A(s - m); multiply by A(s): A(s)/A(s-m) = prod_{j=1}^{m/2} (s - 2j)^2.
    terms = {}
    for m, poly in rec.items():
        q = sp.expand(poly.subs(s, 1 - s) * sp.prod([(s - 2 * j) ** 2 for j in range(1, m // 2 + 1)]))
        terms[m] = q                              # q(s) Phi(s - m)
    M = max(terms)
    # q(s) Phi(s - m) = M[t^(-m) q(-theta + m) F]; multiply by t^M on the left: t^(M - m) q(m - theta)
    ops = {}
    for m, q in terms.items():
        k = (M - m) // 2
        ops[k] = sp.expand(q.subs(s, m - u))
    return ops


def Hvals(tv, kind, g=None):
    w = {"gauss": lambda y: mp.exp(-y * y), "cosh": lambda y: mp.cosh(g * y)}[kind]
    f0 = lambda y: mp.besselk(0, tv * y) * mp.besselk(0, y) * w(y)
    f1 = lambda y: -y * mp.besselk(1, tv * y) * mp.besselk(0, y) * w(y)
    f2 = lambda y: y**2 * (mp.besselk(0, tv * y) + mp.besselk(1, tv * y) / (tv * y)) * mp.besselk(0, y) * w(y)
    # K0''' = -(K0'' + ... ) : use d/dz K0 = -K1, K1' = -K0 - K1/z
    def k0d3(z):      # K0''' = -(K1 (1 + 2/z^2) + K0/z)
        return -(mp.besselk(1, z) * (1 + 2 / z**2) + mp.besselk(0, z) / z)
    f3 = lambda y: y**3 * k0d3(tv * y) * mp.besselk(0, y) * w(y)
    def k0d4(z):      # K0'''' = K0 (1 + 3/z^2) + K1 (2/z + 6/z^3)
        return mp.besselk(0, z) * (1 + 3 / z**2) + mp.besselk(1, z) * (2 / z + 6 / z**3)
    f4 = lambda y: y**4 * k0d4(tv * y) * mp.besselk(0, y) * w(y)
    cuts = [0, 0.5, 1, 2, 4, 8, 16, mp.inf]
    return [mp.quad(f, cuts) for f in (f0, f1, f2, f3, f4)]


def residual(cs, vals, tv):
    terms = [sp.lambdify(t, c, "mpmath")(tv) * v for c, v in zip(cs, vals)]
    return abs(sum(terms)) / max(abs(v) for v in terms)


def theta_degrees(ops):
    return {k: sp.Poly(p, u).degree() for k, p in ops.items()}


if __name__ == "__main__":
    mp.mp.dps = 40
    OUT["x_equation"] = bool(check_x_equation())
    print("x-equation on K0 e^{-x^2}:", OUT["x_equation"])
    rec = {0: s**2, 2: -(4 * s + 5), 4: sp.Integer(4)}
    ops = recurrence_to_operator(rec, 2)
    print("L_G theta-form:", {k: sp.factor(p) for k, p in ops.items()})
    OUT["L_G_theta"] = {str(k): str(sp.factor(p)) for k, p in ops.items()}
    OUT["theta_degrees"] = {str(k): v for k, v in theta_degrees(ops).items()}
    L = theta_to_d(ops)
    print("order:", len(L) - 1, " theta-degrees by power of t^2:", theta_degrees(ops))
    # sanity: derivative formulas for K0
    z0 = mp.mpf("1.3")
    OUT["k0_deriv_check"] = [float(abs(mp.diff(lambda z: mp.besselk(0, z), z0, n) - v))
                             for n, v in ((3, -(mp.besselk(1, z0) * (1 + 2 / z0**2) + mp.besselk(0, z0) / z0)),
                                          (4, mp.besselk(0, z0) * (1 + 3 / z0**2) + mp.besselk(1, z0) * (2 / z0 + 6 / z0**3)))]
    print("K0 derivative formulas (abs err):", OUT["k0_deriv_check"])
    pts = [mp.mpf("0.4"), mp.mpf("0.9"), mp.mpf("1.7"), mp.mpf("3.1")]
    res, resp = {}, {}
    Lp = list(L)
    cs = sp.Poly(Lp[0], t).all_coeffs()[::-1]
    kk = max(range(len(cs)), key=lambda i: abs(cs[i]))
    Lp[0] = Lp[0] + sp.Rational(1, 200) * cs[kk] * t**kk
    for tv in pts:
        vals = Hvals(tv, "gauss")
        res[str(tv)] = float(residual(L, vals, tv))
        resp[str(tv)] = float(residual(Lp, vals, tv))
    print("L_G on H_G:", res)
    print("1/200 perturbation:", resp)
    OUT.update({"residual": res, "perturbed": resp})
    # growth of the monodromy companion
    # J lies in the solution space of every annihilator of H_G with rational coefficients: continuing t around 0
    # turns H_G into H_G - i pi J.  Saddle point: log J = (t-1)^2/4 + O(log t).
    J, rem_ = {}, {}
    for tv in (4, 6, 8, 10, 14, 20, 28):
        tv = mp.mpf(tv)
        J[str(tv)] = mp.quad(lambda y: mp.besseli(0, tv * y) * mp.besselk(0, y) * mp.exp(-y * y), [0, 1, 2, 4, 8, 16, 32, mp.inf])
        rem_[str(tv)] = float(mp.log(J[str(tv)]) - (tv - 1) ** 2 / 4)
    ts = [mp.mpf(k) for k in J]
    slopes = [float((mp.log(J[str(ts[i + 1])]) - mp.log(J[str(ts[i])])) / (ts[i + 1] ** 2 - ts[i] ** 2)) for i in range(len(ts) - 1)]
    print("slope of log J against t^2 (expect -> 1/4):", slopes)
    print("log J - (t-1)^2/4 (expect O(log t)):", rem_)
    OUT["logJ_slopes"] = slopes
    OUT["logJ_minus_saddle"] = rem_
    # control: cosh avatar
    g = sp.Rational(3, 2)
    # annihilator of g = K0(x) cosh(gamma x): module basis y w, (theta y) w, y (theta w), (theta y)(theta w)
    # theta^2 y = x^2 y ; theta^2 w = theta w + gamma^2 x^2 w
    def step(vec):
        a, b, c, d = vec           # coefficients of yw, (ty)w, y(tw), (ty)(tw), polynomials in x
        X2 = x**2
        return [sp.expand(theta_apply(a, 1) + X2 * b + g**2 * X2 * c),
                sp.expand(theta_apply(b, 1) + a + g**2 * X2 * d),
                sp.expand(theta_apply(c, 1) + a + c + X2 * d),
                sp.expand(theta_apply(d, 1) + b + c + d)]
    vecs = [[sp.Integer(1), 0, 0, 0]]
    for _ in range(4):
        vecs.append(step(vecs[-1]))
    Mx = sp.Matrix(vecs).T                       # 4 x 5
    ns = Mx.nullspace()
    assert len(ns) == 1
    v = ns[0]
    den = sp.lcm([sp.fraction(sp.cancel(e))[1] for e in v])
    coeffs = [sp.expand(sp.cancel(e * den)) for e in v]
    # sum_k coeffs[k](x) theta^k g = 0 -> Mellin: coeffs[k] = sum_j c_kj x^j ; x^j theta^k g -> (-(s + j))^k G(s + j)
    rec_c = {}
    for k, c in enumerate(coeffs):
        for (j,), cj in sp.Poly(c, x).terms():
            rec_c[j] = sp.expand(rec_c.get(j, 0) + cj * (-(s + j)) ** k)
    assert all(j % 2 == 0 for j in rec_c)
    ops_c = recurrence_to_operator(rec_c, max(rec_c) // 2)
    Lc = theta_to_d(ops_c)
    degs_c = theta_degrees(ops_c)
    print("cosh control: order", len(Lc) - 1, " theta-degrees:", degs_c)
    rem = right_divide(Lc, leg_op(g))
    OUT["cosh_control"] = {"order": len(Lc) - 1, "theta_degrees": {str(k): v for k, v in degs_c.items()},
                           "L_Leg_right_divides": len(rem) == 0}
    print("L_Leg right-divides the cosh operator:", len(rem) == 0)
    json.dump(OUT, open(DATA / "ga2_gaussian_avatar.json", "w"), indent=1)
