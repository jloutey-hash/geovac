r"""Fast arbitration: semi-independent reference for X_l^{m,s}(P1,P2).

The inner P-integral is done in CLOSED FORM (monomial partial integrals, exact), and
the OUTER Q-integral is done by DIRECT mpmath quadrature of the reduced d^m Q_l (via
_RQ_mp).  This shares the IBP inner with general_m_engine but replaces the B-moment
recurrence with direct quadrature -- so:
  ref   vs  recurrence  -> isolates the forward-B-recurrence error
  ref   vs  driver-grid -> isolates the driver's q_deriv (differentiation) error
"""
import sys; sys.path.insert(0, 'debug')
import numpy as np
import mpmath as mp
import prolate_ci_general_m as drv
import general_m_engine as eng
from val_stage12 import driver_Xtab

mp.mp.dps = 35


def inner_closed(x2, P1, l, m, s, c, Amono):
    """int_1^{x2} xi^P1 (xi^2-1)^s d^mP_l e^{-c xi} dxi  (closed form)."""
    W = list(eng._W_poly(l, m, s, P1))       # coeffs of xi^P1 (xi^2-1)^s d^mP_l
    tot = mp.mpf(0)
    ex = mp.e ** (-c * x2)
    for j, wj in enumerate(W):
        if wj == 0:
            continue
        # int_1^{x2} xi^j e^{-c xi} = A_j(c) - e^{-c x2} sum_k j!/((j-k)! c^{k+1}) x2^{j-k}
        tail = mp.mpf(0)
        fac = mp.mpf(1)
        for k in range(j + 1):
            if k == 0:
                term = mp.mpf(1) / c
            else:
                fac *= (j - k + 1)
                term = fac / c ** (k + 1)
            tail += term * x2 ** (j - k)
        tot += wj * (Amono[j] - ex * tail)
    return tot


def ref_X(l, m, s, P1, P2, c):
    c = mp.mpf(c)
    n_mono = P1 + P2 + 2 * s + (l - m) + 2
    Amono = eng._mono_moments(c, n_mono)

    def regionI(Pouter, Pinner):
        # outer Q at power Pouter, inner P at power Pinner
        def g(u):
            x = 1 + u
            return (x ** Pouter * (x * x - 1) ** s * mp.e ** (-c * x)
                    * eng._RQ_mp(l, m, x) * inner_closed(x, Pinner, l, m, s, c, Amono))
        return mp.quad(g, [0, mp.mpf('0.05'), mp.mpf('0.3'), 1, 2, 4, 8, mp.inf])

    return regionI(P2, P1) + regionI(P1, P2)


if __name__ == "__main__":
    alpha = 1.0
    c = 2.0 * alpha
    ms_pairs = [(0, 0), (0, 1), (1, 1), (2, 2)]
    m_set = [0, 1, 2]; s_set = [0, 1, 2]
    l_neumann, p_max = 10, 4
    Xd = driver_Xtab(m_set, s_set, l_neumann, p_max, alpha)
    Xr = eng.build_Xtab(ms_pairs, l_neumann, p_max, alpha)

    tests = [(2, 0, 0, 0, 0), (4, 0, 1, 0, 0), (6, 0, 1, 0, 0), (8, 0, 1, 0, 0),
             (10, 0, 1, 0, 0), (6, 0, 1, 0, 2), (2, 2, 2, 0, 0), (6, 2, 2, 0, 0),
             (8, 2, 2, 0, 0), (10, 2, 2, 0, 0), (7, 1, 1, 0, 0), (9, 1, 1, 0, 0)]
    print(" (l,m,s,P1,P2) |      ref(35dps)     | driver-grid rel | recurrence rel")
    for (l, m, s, P1, P2) in tests:
        r = ref_X(l, m, s, P1, P2, c)
        dv = Xd[(l, m, s)][P1, P2]
        rv = Xr[(l, m, s)][P1, P2]
        rd = abs(mp.mpf(dv) - r) / abs(r) if r != 0 else mp.mpf('nan')
        rr = abs(mp.mpf(rv) - r) / abs(r) if r != 0 else mp.mpf('nan')
        print(f" ({l},{m},{s},{P1},{P2}) | {mp.nstr(r,10):>18} | {mp.nstr(rd,4):>13} | {mp.nstr(rr,4)}")
