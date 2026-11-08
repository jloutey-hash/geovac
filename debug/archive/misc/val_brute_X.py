r"""Brute high-precision 2D reference for X_l^{m,s}(P1,P2), to arbitrate the driver
grid method vs the recurrence method at large l.

X = int_1^inf dxi2 xi2^P2 (xi2^2-1)^s e^{-c xi2} RQ(xi2) * int_1^{xi2} dxi1 xi1^P1 (xi1^2-1)^s e^{-c xi1} RP(xi1)
  + (P1<->P2 with the other electron the larger)
where RP = d^m P_l, RQ = d^m Q_l, c = 2*alpha.
"""
import sys; sys.path.insert(0, 'debug')
import numpy as np
import mpmath as mp
import prolate_ci_general_m as drv
import general_m_engine as eng
from numpy.polynomial import polynomial as P

mp.mp.dps = 40


def RP_mp(l, m, x):
    return eng._polyval(list(eng._RP_poly(l, m)), x)


def brute_X(l, m, s, P1, P2, c):
    c = mp.mpf(c)

    def inner(xhi, Ppow):
        # int_1^{xhi} xi^Ppow (xi^2-1)^s e^{-c xi} RP(xi) dxi
        def g(u):
            x = 1 + u
            return x ** Ppow * (x * x - 1) ** s * mp.e ** (-c * x) * RP_mp(l, m, x)
        return mp.quad(g, [0, min(mp.mpf(xhi) - 1, mp.mpf('0.02')),
                           min(mp.mpf(xhi) - 1, mp.mpf('0.3')), mp.mpf(xhi) - 1])

    def outerQ(Pouter, Pinner):
        # int_1^inf xi^Pouter (xi^2-1)^s e^{-c xi} RQ(xi) * inner(xi, Pinner) dxi
        def g(u):
            x = 1 + u
            return (x ** Pouter * (x * x - 1) ** s * mp.e ** (-c * x)
                    * eng._RQ_mp(l, m, x) * inner(x, Pinner))
        return mp.quad(g, [0, mp.mpf('0.05'), mp.mpf('0.3'), 1, 2, 4, 8, mp.inf])

    # region xi1<xi2: outer=xi2(P2,Q), inner=xi1(P1,P)
    A = outerQ(P2, P1)
    # region xi2<xi1: outer=xi1(P1,Q), inner=xi2(P2,P)
    B = outerQ(P1, P2)
    return A + B


if __name__ == "__main__":
    alpha = 1.0
    c = 2.0 * alpha
    # driver grid Xtab and recurrence Xtab, mu<=1
    from val_stage12 import driver_Xtab
    ms_pairs = [(0, 0), (0, 1), (1, 1), (2, 2)]
    m_set = [0, 1, 2]; s_set = [0, 1, 2]
    l_neumann, p_max = 8, 4
    Xd = driver_Xtab(m_set, s_set, l_neumann, p_max, alpha)
    Xr = eng.build_Xtab(ms_pairs, l_neumann, p_max, alpha)

    tests = [(2, 0, 0, 0, 0), (4, 0, 1, 0, 0), (6, 0, 1, 0, 0), (8, 0, 1, 0, 0),
             (6, 0, 1, 0, 2), (2, 2, 2, 0, 0), (6, 2, 2, 0, 0), (8, 2, 2, 0, 0),
             (7, 1, 1, 0, 0)]
    print(" (l,m,s,P1,P2) |    brute(40dps)     |  driver-grid rel  |  recurrence rel")
    for (l, m, s, P1, P2) in tests:
        ref = brute_X(l, m, s, P1, P2, c)
        dv = Xd[(l, m, s)][P1, P2]
        rv = Xr[(l, m, s)][P1, P2]
        rd = abs(mp.mpf(dv) - ref) / abs(ref) if ref != 0 else mp.mpf('nan')
        rr = abs(mp.mpf(rv) - ref) / abs(ref) if ref != 0 else mp.mpf('nan')
        print(f" ({l},{m},{s},{P1},{P2}) | {mp.nstr(ref,10):>18} | {mp.nstr(rd,4):>14} | {mp.nstr(rr,4)}")
