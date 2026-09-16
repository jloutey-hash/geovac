r"""Verify the recurrence X vs the semi-independent reference for the delta-channel
kernel orders m=3,4 (mu=2 regime) at high l -- the regime where the driver blows up.
"""
import sys; sys.path.insert(0, 'debug')
import mpmath as mp
import general_m_engine as eng
from val_arbitrate import ref_X

mp.mp.dps = 35

if __name__ == "__main__":
    alpha = 1.0
    c = 2.0 * alpha
    # (m,s) pairs of the mu<=2 assembly, extended kernel orders
    ms_pairs = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2), (3, 3), (4, 4)]
    l_neumann, p_max = 16, 4
    Xr = eng.build_Xtab(ms_pairs, l_neumann, p_max, alpha)

    tests = [(4, 4, 4, 0, 0), (8, 4, 4, 0, 0), (12, 4, 4, 0, 0), (16, 4, 4, 0, 0),
             (3, 3, 3, 0, 0), (9, 3, 3, 0, 0), (15, 3, 3, 0, 0),
             (12, 2, 2, 0, 0), (16, 2, 2, 0, 0), (14, 0, 2, 0, 0),
             (12, 4, 4, 0, 2), (16, 4, 4, 2, 2)]
    print(" (l,m,s,P1,P2) |      ref(35dps)      |  recurrence rel err")
    worst = mp.mpf(0)
    for (l, m, s, P1, P2) in tests:
        r = ref_X(l, m, s, P1, P2, c)
        rv = Xr[(l, m, s)][P1, P2]
        rr = abs(mp.mpf(rv) - r) / abs(r) if r != 0 else mp.mpf('nan')
        worst = max(worst, rr)
        print(f" ({l},{m},{s},{P1},{P2}) | {mp.nstr(r,10):>19} | {mp.nstr(rr,4)}")
    print(f"\n WORST recurrence rel err (m up to 4, l up to 16) = {mp.nstr(worst,4)}")
