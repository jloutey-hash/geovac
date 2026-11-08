"""Loop test PT2, step B: position-space form of the fingerprint triangle and chain, and their
exact linear ODEs in t by the Mellin route (record docs/predictions/2026-10-17_density_loop.md).

Fingerprint (rho~ = Delta^-3, Delta = sqrt(c k^2 + 1), a = c^-1/2). Fourier back to position:
  triangle  G(t) = (pi/2) t^2 int_0^inf x^2 psi(x)^2 psi(t x) dx,   psi = L_{-1} - I_1
            (3-D inverse FT of rho~/k; psi solves (theta^2 - 1 - x^2) psi = -2/pi, so
             theta (theta^2 - 1 - x^2) psi = 0);
  chain     C(t) = (8 t^3/pi^3) int_0^inf chi(x)^2 K0(t x) dx,  chi(x) = int_0^x s K1(s) ds
            (3-D inverse FT of rho~/k^2 up to 1/(2 pi^2 r); theta chi = x^2 K1, so
             ((theta-2)^2 - 1 - x^2) theta chi = 0);
scales (1, 1, t), t on the density shared by both kernels (chain) / any density (triangle).
Step 1 checks both forms against the 2-D momentum formulas of the record.
Step 2: f^2 obeys the symmetric square of the order-3 operator; its Mellin transform a recurrence;
Phi(s) = M[G](s) = A(alpha - s) E(s + beta) with E(sig + 2) = r(sig) E(sig); back to an ODE in
theta_t. Code self-check: reproduces Paper 59's F_3 operator 4 th^3 - t^2 (th+1)^3.
Step 3: each ODE checked on its function by quadrature of theta_t^j applied under the integral.
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import mpmath as mp
import sympy as sp

OUT = Path(__file__).parent / "data" / "density_loop_B.json"
z, s, u, y = sp.symbols("z s u y")


# ------------------------------------------------------------------ symbolic
def sym_power(rel, N, m):
    """theta^N f = sum_{j<N} rel[j](z) theta^j f  (theta z = 2 z).  Returns v with
    sum_j v[j](z) theta^j (f^m) = 0 (minimal order from this basis)."""
    keys = [k for k in itertools.combinations_with_replacement(range(N), m)]
    idx = {k: i for i, k in enumerate(keys)}

    def theta(vec):
        new = [sp.Integer(0)] * len(keys)
        for k, cf in zip(keys, vec):
            if cf == 0:
                continue
            new[idx[k]] += 2 * z * sp.diff(cf, z)
            for l in range(m):
                kk = list(k)
                kk[l] += 1
                if kk[l] < N:
                    new[idx[tuple(sorted(kk))]] += cf
                else:
                    for j in range(N):
                        k2 = kk.copy()
                        k2[l] = j
                        new[idx[tuple(sorted(k2))]] += cf * rel[j]
        return [sp.expand(c) for c in new]

    rows = [[sp.Integer(1)] + [sp.Integer(0)] * (len(keys) - 1)]
    for order in range(1, len(keys) + 1):
        rows.append(theta(rows[-1]))
        M = sp.Matrix(rows).T
        ns = M.nullspace()
        if ns:
            v = ns[0]
            den = sp.lcm([sp.fraction(sp.together(e))[1] for e in v])
            v = [sp.factor(sp.cancel(e * den)) for e in v]
            g = sp.gcd_list([e for e in v if e != 0])
            return [sp.expand(sp.cancel(e / g)) for e in v]
    raise RuntimeError("no relation")


def mellin_t_ode(v, alpha, beta, r):
    """Phi(s) = A(alpha - s) E(s + beta), A = M[f^m], E(sig+2) = r(sig) E(sig).
    Returns Q with sum_i t^{2i} Q_i(theta_t) F = 0."""
    K = max(sp.Poly(e, z).degree() for e in v if e != 0)
    P = [sp.expand(sum(sp.Poly(v[j], z).coeff_monomial(z**k) * u**j for j in range(len(v))))
         for k in range(K + 1)]
    R = []
    for i in range(K + 1):
        k = K - i
        S = s + 2 * K
        R.append(sp.expand(P[k].subs(u, -(alpha - S + 2 * k)) * sp.prod([r(S + beta - 2 * l) for l in range(1, k + 1)])))
    # a common factor may be removed from the recurrence in s (not from the theta-polynomials:
    # theta does not commute with t^2, g(theta) t^2 = t^2 g(theta + 2))
    g = sp.gcd_list([e for e in R if e != 0])
    R = [sp.expand(sp.cancel(e / g)) for e in R]
    return [sp.factor(sp.expand(R[i].subs(s, -u - 2 * i))) for i in range(K + 1)]


def theta_powers(J, rel2, inhom):
    """theta^j g = al_j g + be_j theta g + ga_j for g with theta^2 g = rel2(y) g + inhom."""
    al, be, ga = [sp.Integer(1)], [sp.Integer(0)], [sp.Integer(0)]
    for _ in range(J):
        a, b, c = al[-1], be[-1], ga[-1]
        al.append(sp.expand(y * sp.diff(a, y) + b * rel2))
        be.append(sp.expand(a + y * sp.diff(b, y)))
        ga.append(sp.expand(y * sp.diff(c, y) + b * inhom))
    return al, be, ga


# ------------------------------------------------------------------ numerics
X0 = 120


def psi(x):
    x = mp.mpf(x)
    if x < X0:
        with mp.workdps(mp.mp.dps + int(x / 2.2) + 15):
            return +(mp.struvel(-1, x) - mp.besseli(1, x))
    tot, m = mp.mpf(0), 0
    while True:
        c = (-1) ** m * mp.factorial(2 * m) * mp.binomial(mp.mpf(-1.5), m) * x ** (-2 * m - 2)
        tot += c
        if abs(c) < mp.eps * abs(tot) or m > x / 2:
            return 2 / mp.pi * tot
        m += 1


def theta_psi(x):
    x = mp.mpf(x)
    if x < X0:
        with mp.workdps(mp.mp.dps + int(x / 2.2) + 15):
            return +(x * mp.struvel(-2, x) + mp.struvel(-1, x) - x * mp.besseli(0, x) + mp.besseli(1, x))
    tot, m = mp.mpf(0), 0
    while True:
        c = (-1) ** m * mp.factorial(2 * m) * mp.binomial(mp.mpf(-1.5), m) * x ** (-2 * m - 2) * (-2 * m - 2)
        tot += c
        if abs(c) < mp.eps * abs(tot) or m > x / 2:
            return 2 / mp.pi * tot
        m += 1


def chi(x):
    x = mp.mpf(x)
    if x == 0:
        return mp.mpf(0)
    with mp.workdps(mp.mp.dps + 10):
        return +(mp.pi * x / 2 * (mp.besselk(1, x) * mp.struvel(0, x) + mp.besselk(0, x) * mp.struvel(1, x)))


def G_pos(t, jt=None):
    """G or theta_t^j G (jt = (al, be, ga) lambdified for (theta+2)^j)."""
    t = mp.mpf(t)

    def h(x):
        yy = t * x
        if jt is None:
            inner = psi(yy)
        else:
            al, be, ga = jt
            inner = al(yy) * psi(yy) + be(yy) * theta_psi(yy) + ga(yy)
        return x * x * psi(x) ** 2 * inner
    pts = [0, mp.mpf(1) / 4, 1, 4, 16, 64, X0, 4 * X0, mp.inf]
    return mp.pi / 2 * t * t * mp.quad(h, pts)


def C_pos(t, jt=None):
    t = mp.mpf(t)

    def h(x):
        yy = t * x
        if jt is None:
            inner = mp.besselk(0, yy)
        else:
            al, be = jt
            inner = al(yy) * mp.besselk(0, yy) - be(yy) * yy * mp.besselk(1, yy)
        return chi(x) ** 2 * inner
    return 8 * t**3 / mp.pi**3 * mp.quad(h, [0, mp.mpf(1) / 4, 1, 4, 16, 64, mp.inf])


def G_mom(t):
    t = mp.mpf(t)
    D = lambda k, c: mp.sqrt(c * k * k + 1)
    ct = 1 / t**2
    hh = lambda w: w / D(w, ct)

    def inner(p):
        return mp.quad(lambda q: D(q, 1) ** -3 * (hh(p + q) - hh(abs(p - q))), [0, p, 2 * p + 1, mp.inf])
    return mp.quad(lambda p: D(p, 1) ** -3 * inner(p), [0, 1, 4, mp.inf]) / mp.pi


def C_mom(t):
    t = mp.mpf(t)
    D = lambda k, c: mp.sqrt(c * k * k + 1)
    ct = 1 / t**2

    def inner(p):
        f = lambda q: D(q, 1) ** -3 / (p * q) * t * t * (1 / D(abs(p - q), ct) - 1 / D(p + q, ct))
        return mp.quad(f, [0, p, 2 * p + 1, mp.inf])
    return 2 / mp.pi**2 * mp.quad(lambda p: D(p, 1) ** -3 * inner(p), [0, 1, 4, mp.inf])


def residual(Q, ts, theta_F):
    """sum_i t^{2i} Q_i(theta) F / max term."""
    out = []
    for t in ts:
        t = mp.mpf(t)
        J = max(sp.Poly(q, u).degree() for q in Q)
        th = [theta_F(t, j) for j in range(J + 1)]
        terms = []
        for i, q in enumerate(Q):
            cf = sp.Poly(q, u).all_coeffs()[::-1]
            terms.append(t ** (2 * i) * sum(mp.mpf(sp.Rational(c).p) / sp.Rational(c).q * th[j] for j, c in enumerate(cf)))
        out.append(abs(sum(terms)) / max(abs(x) for x in terms))
    return out


def main(dps: int) -> None:
    mp.mp.dps = dps
    rec = {"dps": dps}
    # --- step 1: position forms vs the record's momentum formulas
    for t in ["0.7", "1.6"]:
        a, b = G_pos(t), G_mom(t)
        c, d = C_pos(t), C_mom(t)
        print("step1 t", t, "G", mp.nstr(a, 18), mp.nstr(b, 18), "rel", mp.nstr(abs(a / b - 1), 3),
              "| C", mp.nstr(c, 18), mp.nstr(d, 18), "rel", mp.nstr(abs(c / d - 1), 3), flush=True)
        rec[f"step1_t{t}"] = {"G_pos": mp.nstr(a, 25), "G_mom": mp.nstr(b, 25),
                              "C_pos": mp.nstr(c, 25), "C_mom": mp.nstr(d, 25)}
    # --- step 2: operators
    vK0 = sym_power([z, 0], 2, 2)
    QF3 = mellin_t_ode(vK0, 1, 0, lambda sg: sg**2)
    print("self-check F3:", QF3, flush=True)
    vpsi = sym_power([2 * z, 1 + z, 0], 3, 2)     # theta^3 psi = 2z psi + (1+z) theta psi
    vchi = sym_power([0, z - 3, 4], 3, 2)         # theta^3 chi = (z-3) theta chi + 4 theta^2 chi
    QG = mellin_t_ode(vpsi, 1, 2, lambda sg: sg**2 - 1)
    QC = mellin_t_ode(vchi, -2, 3, lambda sg: sg**2)
    for name, Q in [("G", QG), ("C", QC)]:
        print(name, "ODE order", max(sp.Poly(q, u).degree() for q in Q), "t-degree", 2 * (len(Q) - 1), flush=True)
        for i, q in enumerate(Q):
            print(f"   t^{2*i}:", q, flush=True)
    rec["QF3"] = [str(q) for q in QF3]
    rec["QG"] = [str(q) for q in QG]
    rec["QC"] = [str(q) for q in QC]
    # --- step 3: numeric check of each ODE
    JG = max(sp.Poly(q, u).degree() for q in QG)
    al, be, ga = theta_powers(JG + 1, 1 + y**2, -2 / sp.pi)
    lam = {}

    def thG(t, j):
        if j not in lam:
            A = sum(sp.binomial(j, i) * 2 ** (j - i) * al[i] for i in range(j + 1))
            B = sum(sp.binomial(j, i) * 2 ** (j - i) * be[i] for i in range(j + 1))
            Cc = sum(sp.binomial(j, i) * 2 ** (j - i) * ga[i] for i in range(j + 1))
            lam[j] = tuple(sp.lambdify(y, sp.expand(e), "mpmath") for e in (A, B, Cc))
        return G_pos(t, lam[j])

    JC = max(sp.Poly(q, u).degree() for q in QC)
    alK, beK, _ = theta_powers(JC + 1, y**2, 0)
    lamC = {}

    def thC(t, j):
        if j not in lamC:
            A = sum(sp.binomial(j, i) * 3 ** (j - i) * alK[i] for i in range(j + 1))
            B = sum(sp.binomial(j, i) * 3 ** (j - i) * beK[i] for i in range(j + 1))
            lamC[j] = tuple(sp.lambdify(y, sp.expand(e), "mpmath") for e in (A, B))
        return C_pos(t, lamC[j])

    rG = residual(QG, ["0.7", "1.6"], thG)
    rC = residual(QC, ["0.7", "1.6"], thC)
    print("ODE residual G:", [mp.nstr(x, 3) for x in rG], " C:", [mp.nstr(x, 3) for x in rC], flush=True)
    rec["resid_G"] = [mp.nstr(x, 3) for x in rG]
    rec["resid_C"] = [mp.nstr(x, 3) for x in rC]
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 20)
