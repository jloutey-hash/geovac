"""WG phase 2: self-test of wg_common on toy values (not scored cases)."""
import mpmath as mp
import wg_common as C
import wg_derive as W
mp.mp.dps = 50
for T0, T1 in (('1.3', '1.1'), ('110', '100')):
    T0, T1 = mp.mpf(T0), mp.mpf(T1)
    Q = dict(pref=mp.mpf(1), T0=T0, T1=T1)
    F = C.F_quartet(Q, 40)
    for l in (0, 7, 40):
        print('toy', mp.nstr(T0, 4), l, mp.nstr(F[l] / W.F_closed(T0, T1, l) - 1, 5))
    M = C.kummer_table(T0, 300)
    for m in (0, 50, 120):
        print('kummer', m, mp.nstr(M[m] / mp.hyp1f1(1, m + 1.5, T0) - 1, 5))
