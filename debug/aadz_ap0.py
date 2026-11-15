"""Record AD (docs/predictions/2026-11-15_aa_bc_d_orbital.md), AD0 self-checks (adapted from aaqx_ap0.py).
(i)  the 'AB' path of aapc_exact still reproduces aabc_exact key for key (Q and mod 2^61 - 1, L <= 30);
(ii) F_L of 'PD' at L = 0, 2, 5, 12 (both cases) against a quadrature: the 3d-sigma density multipole projected
     directly with h2(r - C) = (3 (r t - D1)^2 - s^2)/2 (no T recurrence), against the multipole v_L of the exact
     1s potential V_A; no kernel split.
(iii) F_L for L <= 220 unchanged when the internal L_max goes 240 -> 260 (mod 2^61 - 1).
Log: data/aadz_ap0.log.
"""
import sys
from fractions import Fraction
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import aabc_exact as X  # noqa: E402
import aapc_exact as P  # noqa: E402
import na_exact as E  # noqa: E402


def mf(q):
    q = Fraction(q)
    return mp.mpf(q.numerator) / q.denominator


def legendre_and_deriv(L, t):
    p0, p1 = mp.mpf(1), t
    if L == 0:
        return p0, mp.mpf(0)
    for n in range(1, L):
        p0, p1 = p1, ((2 * n + 1) * t * p1 - n * p0) / (n + 1)
    return p1, L * (t * p1 - p0) / (t * t - 1) if abs(t * t - 1) > 0 else mp.mpf(L * (L + 1)) / 2


def multipole_in_s(func, r, D, nodes):
    """int_{-1}^{1} func(s, t) dt with s = |r - D e_z| as the variable (func s smooth): Gauss-Legendre in s."""
    a, b = abs(r - D), r + D
    h, c = (b - a) / 2, (b + a) / 2
    tot = mp.mpf(0)
    for u, w in nodes:
        s = c + h * u
        t = (r * r + D * D - s * s) / (2 * r * D)
        tot += w * func(s, t) * s
    return tot * h / (r * D)


def V_A(s, za):
    return 1 / s - (za + 1 / s) * mp.exp(-2 * za * s)


def quad_GL(case, L, dps=22):
    """F_L of PD by direct projection: (2L+1)/2 int N_B e^{-zB r} N_d h2 e^{-zC s} P_L(t) dt with
    h2 = (3 (r t - D1)^2 - s^2)/2 taken as is (no T recurrence), against v_L of the exact 1s potential."""
    with mp.workdps(dps):
        D1, D2, zB, zC, zA = (mf(v) for v in P.params(case))
        NB = mp.sqrt(zB ** 3 / mp.pi)
        Nd = mp.sqrt(5 * (2 * zC) ** 7 / (4 * mp.pi * 720))
        nodes = mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(5, mp.mp.prec)

        def rhoL(r):
            def f(s, t):
                return mp.exp(-zC * s) * (3 * (r * t - D1) ** 2 - s * s) / 2 * legendre_and_deriv(L, t)[0]
            return (2 * L + 1) / mp.mpf(2) * NB * Nd * mp.exp(-zB * r) * multipole_in_s(f, r, D1, nodes)

        def vL(r):
            return (2 * L + 1) / mp.mpf(2) * multipole_in_s(
                lambda s, t: V_A(s, zA) * legendre_and_deriv(L, t)[0], r, D2, nodes)
        return mp.quad(lambda r: 4 * mp.pi / (2 * L + 1) * r * r * rhoL(r) * vL(r), [0, D1, D2, D2 + 4, mp.inf])


def main():
    ok = True
    for p in (None, (1 << 61) - 1):
        K = E.Field(p)
        for case in P.CASES:
            a = P.coefficients(case, 30, K, 'AB')
            pc = X.pieces(case, 30, K)
            b = X.combine(pc, P.params(case)[4], K)
            same = a == b
            ok &= same
            print(case, 'Q' if p is None else 'p', '(i) AB reproduced:', same, flush=True)
    K = E.Field(None)
    for case in P.CASES:
        f = P.coefficients(case, 12, K, 'PD')
        mp.mp.dps = 80
        tw = P.twopiN(case, 'PD')
        keys = set(k for q in f for k in q)
        cv = {k: E.value(k, mp) for k in keys}
        for L in (0, 2, 5, 12):
            ex = tw * sum(mp.mpf(int(v.p)) / int(v.q) * cv[k] for k, v in f[L].items())
            q = quad_GL(case, L)
            with mp.workdps(40):
                rel = abs(q / ex - 1)
            good = rel < mp.mpf('1e-18')
            ok &= bool(good)
            print(case, 'L', L, '(ii) exact', mp.nstr(ex, 20), 'quad', mp.nstr(q, 20), 'rel', mp.nstr(rel, 3), good,
                  flush=True)
    # (iii) L_max invariance: F_L for L <= 220 unchanged when the internal L_max goes 240 -> 260 (mod p)
    Kp = E.Field((1 << 61) - 1)
    for case in P.CASES:
        a240 = P.coefficients(case, 240, Kp, 'PD')
        a260 = P.coefficients(case, 260, Kp, 'PD')
        same = a240[:221] == a260[:221]
        ok &= same
        print(case, '(iii) F_L, L <= 220, unchanged 240 -> 260:', same, flush=True)
    print('AD0 self-checks pass:', ok)


if __name__ == '__main__':
    main()
