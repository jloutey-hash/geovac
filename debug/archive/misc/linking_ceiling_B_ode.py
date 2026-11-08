"""PL2 step B: exact linear ODE in t for F_n(t) = int_0^inf K0(x)^m K0(t x) dx, m = n - 1.
Route: K0^m satisfies the m-th symmetric power of theta^2 - x^2 (theta = x d/dx); its Mellin
transform A(s) obeys a recurrence; Phi(s) = M[F](s) = A(1-s) B(s), B(s) = 2^{s-2} Gamma(s/2)^2,
B(s+2) = s^2 B(s); the recurrence for Phi converts back to an ODE in theta_t.
Then the ODE is checked numerically on F_n by direct quadrature (theta_t^j K0 in closed form)."""
import sympy as sp, mpmath as mp, json, sys
z, s, u = sp.symbols("z s u")

def sym_power_ode(m):
    """Rows: theta^j f in the basis w_i = u^{m-i} (theta u)^i, theta^2 u = z u, theta z = 2z."""
    def theta(vec):  # vec: list of polys in z, coefficients of w_0..w_m
        new = [sp.Integer(0)] * (m + 1)
        for i, cf in enumerate(vec):
            new[i] += z * 2 * sp.diff(cf, z)          # theta(cf(z)) = 2 z cf'
            if i < m: new[i + 1] += cf * (m - i)
            if i > 0: new[i - 1] += cf * i * z
        return [sp.expand(c) for c in new]
    rows = [[sp.Integer(1)] + [sp.Integer(0)] * m]
    for _ in range(m + 1): rows.append(theta(rows[-1]))
    M = sp.Matrix(rows).T                          # (m+1) x (m+2)
    ns = M.nullspace()
    assert len(ns) == 1
    v = ns[0]; den = sp.lcm([sp.fraction(sp.together(e))[1] for e in v])
    v = [sp.factor(sp.cancel(e * den)) for e in v]
    g = sp.gcd_list([e for e in v if e != 0])
    v = [sp.expand(sp.cancel(e / g)) for e in v]
    return v                                       # sum_j v[j](z) theta^j f = 0

def recurrence_to_t_ode(v):
    # sum_j v_j(z) theta^j f = sum_k z^k P_k(theta) f ; P_k(th) = sum_j coeff(v_j, z^k) th^j
    K = max(sp.Poly(e, z).degree() for e in v if e != 0)
    P = [sp.expand(sum(sp.Poly(v[j], z).coeff_monomial(z**k) * u**j for j in range(len(v))))
         for k in range(K + 1)]
    # A-recurrence: sum_k P_k(-(sig+2k)) A(sig+2k) = 0, sig = 1 - s
    # Phi-recurrence: sum_k P_k(s-1-2k) prod_{i=1}^k (s-2i)^2 Phi(s-2k) = 0
    Rt = [sp.expand(P[k].subs(u, s - 1 - 2*k) * sp.prod([(s - 2*i)**2 for i in range(1, k + 1)]))
          for k in range(K + 1)]
    # shift s -> s + 2K, index i = K - k:  sum_i R_i(s) Phi(s + 2i) = 0
    R = [sp.expand(Rt[K - i].subs(s, s + 2*K)) for i in range(K + 1)]
    # ODE: sum_i t^{2i} Q_i(theta_t) F = 0, Q_i(u) = R_i(-u - 2i)
    Q = [sp.factor(R[i].subs(s, -u - 2*i)) for i in range(K + 1)]
    return Q

def theta_K0_coeffs(J):
    """theta_u^j K0(u) = al_j(u) K0(u) + be_j(u) u K1(u); theta K0 = -uK1, theta(uK1) = -u^2 K0."""
    al, be = [sp.Integer(1)], [sp.Integer(0)]
    for _ in range(J):
        a, b = al[-1], be[-1]
        na = sp.expand(u * sp.diff(a, u) - u**2 * b)
        nb = sp.expand(u * sp.diff(b, u) - a)
        al.append(na); be.append(nb)
    return al, be

def check_numeric(m, Q, ts, dps=30):
    mp.mp.dps = dps
    J = max(sp.Poly(q, u).degree() for q in Q)
    al, be = theta_K0_coeffs(J)
    alf = [sp.lambdify(u, a, "mpmath") for a in al]; bef = [sp.lambdify(u, b, "mpmath") for b in be]
    res = []
    for t in ts:
        t = mp.mpf(t)
        th = []
        for j in range(J + 1):
            f = lambda x, j=j: mp.besselk(0, x)**m * (alf[j](t*x) * mp.besselk(0, t*x)
                                                      + bef[j](t*x) * t*x * mp.besselk(1, t*x))
            th.append(mp.quad(f, [0, 1, 4, mp.inf]))
        tot = mp.mpf(0); scale = mp.mpf(0)
        for i, q in enumerate(Q):
            pq = sp.Poly(q, u)
            for (deg,), cf in pq.terms():
                term = t**(2*i) * mp.mpf(sp.Rational(cf).p) / sp.Rational(cf).q * th[deg]
                tot += term; scale += abs(term)
        res.append((float(t), float(abs(tot) / scale)))
    return res

import os
FN = "debug/data/linking_ceiling_B_ode.json"
out = json.load(open(FN)) if os.path.exists(FN) else {}
for n in [int(a) for a in sys.argv[1:]] or [2, 3, 4]:
    m = n - 1
    v = sym_power_ode(m)
    Q = recurrence_to_t_ode(v)
    order = max(sp.Poly(q, u).degree() for q in Q)
    print(f"n={n}: x-ODE order {len(v)-1}; t-ODE in theta: order {order}, t-degree {2*(len(Q)-1)}")
    for i, q in enumerate(Q): print(f"   t^{2*i}:", q)
    chk = check_numeric(m, Q, ["0.3", "0.7", "1.6"])
    print("   numeric relative residual:", chk)
    out[str(n)] = {"Q": [str(q) for q in Q], "order": order, "check": chk}
json.dump(out, open(FN, "w"), indent=1)
