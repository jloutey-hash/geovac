"""Record AX (docs/predictions/2026-11-15_aa_bc_d_mneq0.md), AX0 self-checks (adapted from aadz_ap0.py, aaqx_ap0.py).
(i)  the 'AB' path of aapc_exact still reproduces aabc_exact key for key (Q and mod 2^61 - 1, L <= 30);
(ii) G_L of 'DX' at L = 1, 3, 7, 12 and K_L of 'DD' at L = 2, 5, 12 (both cases) against a quadrature in which the
     density's m = 1 or m = 2 multipole is projected directly (no T recurrence, no raising identity):
        DX: rho^1_L(r) = (2L+1)/(2L(L+1)) N e^{-zB r} r int e^{-zC s} (r t - D1) (1-t^2) P_L'(t) dt
        DD: rho^2_L(r) = (2L+1)(L-2)!/(2 (L+2)!) N e^{-zB r} r^2 int e^{-zC s} (1-t^2)^2 P_L''(t) dt
     (sqrt(1-t^2)^m P_L^m = (1-t^2)^m d^m P_L, no Condon-Shortley phase), against the multipole v_L of the exact 1s
     potential V_A; F = int 4 pi/(2L+1) r^2 rho_L v_L dr; no kernel split.  Relative 1e-18 at 22 digits.
(iii) every G_L, K_L for L <= 220 unchanged (mod 2^61 - 1) when the internal L_max goes 240 -> 260.
Log: data/aadx_ap0.log.
"""
import sys
import time
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import aabc_exact as X  # noqa: E402
import aapc_exact as P  # noqa: E402
import na_exact as E  # noqa: E402
from aadz_ap0 import mf, legendre_and_deriv, multipole_in_s, V_A  # noqa: E402

LOG = HERE / 'data' / 'aadx_ap0.log'


def say(*a):
    line = ' '.join(str(v) for v in a)
    print(line, flush=True)
    with open(LOG, 'a') as fh:
        fh.write(line + '\n')


def dmP(m, L, t):
    """d^m P_L(t), m = 1, 2, by the fixed-m recurrence (L-m+1) Q_{L+1} = (2L+1) t Q_L - (L+m) Q_{L-1}."""
    if L < m:
        return mp.mpf(0)
    q0, q1 = mp.mpf(0), mp.mpf(1) if m == 1 else mp.mpf(3)
    for n in range(m, L):
        q0, q1 = q1, ((2 * n + 1) * t * q1 - (n + m) * q0) / (n - m + 1)
    return q1


def quad_GL(case, L, integral, dps=22):
    with mp.workdps(dps):
        D1, D2, zB, zC, zA = (mf(v) for v in P.params(case))
        Rn = mp.sqrt((2 * zC) ** 7 / 720)
        ang = mp.sqrt(15 / (4 * mp.pi)) if integral == 'DX' else mp.sqrt(15 / (16 * mp.pi))
        N = mp.sqrt(zB ** 3 / mp.pi) * ang * Rn
        nodes = mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(5, mp.mp.prec)

        if integral == 'DX':
            def rho(r):
                def f(s, t):
                    return mp.exp(-zC * s) * (r * t - D1) * (1 - t * t) * dmP(1, L, t)
                return ((2 * L + 1) / mp.mpf(2 * L * (L + 1)) * N * mp.exp(-zB * r) * r
                        * multipole_in_s(f, r, D1, nodes))
        else:
            def rho(r):
                def f(s, t):
                    return mp.exp(-zC * s) * (1 - t * t) ** 2 * dmP(2, L, t)
                c = (2 * L + 1) * mp.factorial(L - 2) / (2 * mp.factorial(L + 2))
                return c * N * mp.exp(-zB * r) * r * r * multipole_in_s(f, r, D1, nodes)

        def vL(r):
            return (2 * L + 1) / mp.mpf(2) * multipole_in_s(
                lambda s, t: V_A(s, zA) * legendre_and_deriv(L, t)[0], r, D2, nodes)
        return mp.quad(lambda r: 4 * mp.pi / (2 * L + 1) * r * r * rho(r) * vL(r), [0, D1, D2, D2 + 4, mp.inf])


def main(qdps=22, tol='1e-18', extra=(), tol_extra=None):
    """AX: main() (22 digits, 1e-18).  AX2 (amendment): main(32, '1e-25', (20,), '1e-23'), see
    docs/predictions/2026-11-16_aa_bc_d_mneq0_2.md."""
    LOG.write_text('')
    say('quadrature digits', qdps, 'tolerance', tol, 'extra L', extra)
    t0 = time.time()
    ok = True
    for p in (None, (1 << 61) - 1):
        K = E.Field(p)
        for case in P.CASES:
            a = P.coefficients(case, 30, K, 'AB')
            b = X.combine(X.pieces(case, 30, K), P.params(case)[4], K)
            same = a == b
            ok &= same
            say(case, 'Q' if p is None else 'p', '(i) AB reproduced:', same)
    K = E.Field(None)
    for case in P.CASES:
        for integral, Ls in (('DX', (1, 3, 7, 12) + tuple(extra)), ('DD', (2, 5, 12) + tuple(extra))):
            f = P.coefficients(case, max(Ls), K, integral)
            mp.mp.dps = 80
            tw = P.twopiN(case, integral)
            keys = set(k for q in f for k in q)
            cv = {k: E.value(k, mp) for k in keys}
            for L in Ls:
                ex = tw * sum(mp.mpf(int(v.p)) / int(v.q) * cv[k] for k, v in f[L].items())
                q = quad_GL(case, L, integral, dps=qdps)
                with mp.workdps(max(40, qdps + 18)):
                    rel = abs(q / ex - 1)
                good = rel < mp.mpf(tol_extra if (tol_extra and L in extra) else tol)
                ok &= bool(good)
                say(case, integral, 'L', L, '(ii) exact', mp.nstr(ex, 20), 'quad', mp.nstr(q, 20), 'rel',
                    mp.nstr(rel, 3), good)
    say('(i)+(ii) done', round(time.time() - t0), 's')
    Kp = E.Field((1 << 61) - 1)
    for case in P.CASES:
        for integral in ('DX', 'DD'):
            t1 = time.time()
            a240 = P.coefficients(case, 240, Kp, integral)
            a260 = P.coefficients(case, 260, Kp, integral)
            same = a240[:221] == a260[:221]
            ok &= same
            say(case, integral, '(iii) L <= 220 unchanged 240 -> 260:', same, round(time.time() - t1), 's')
    say('AX0 self-checks pass:', ok, ' total', round(time.time() - t0), 's')


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'ax2':
        LOG = HERE / 'data' / 'aadx2_ap0.log'      # keep AX's log
        main(32, '1e-25', (20,), '1e-23')
    else:
        main()
