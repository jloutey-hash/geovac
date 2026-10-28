"""WG phase 2: diagnosis of the WG-0(b) g1 l=30 mismatch (Weber quadrature vs series). Not a scored item."""
import mpmath as mp
import wg_common as C
for dps in (30, 50):
    mp.mp.dps = dps
    Q = C.quartets('g1')[0]
    for l in (20, 30):
        s = C.F_quartet(Q, l)[l]
        w = C.F_weber_quartet(Q, l)
        # third route: direct Legendre projection of pref * F0(T0 - T1 x)
        f = lambda x: Q['pref'] * mp.sqrt(mp.pi / (Q['T0'] - Q['T1'] * x)) / 2 * mp.erf(mp.sqrt(Q['T0'] - Q['T1'] * x)) if x < 1 else Q['pref']
        dd = (2 * l + 1) / mp.mpf(2) * mp.quad(lambda x: f(x) * mp.legendre(l, x), mp.linspace(-1, 1, 41))
        print(dps, l, 'weber/series-1', mp.nstr(w / s - 1, 5), 'direct/series-1', mp.nstr(dd / s - 1, 5))
