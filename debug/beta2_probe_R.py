"""Convergence probe for R(k,w) (graded s-map) and F(k) at high precision."""
import sys; sys.path.insert(0,'debug')
import mpmath as mp, time
import beta2_t2_kw_core as C

mp.mp.dps = 60
for k in (5, 40, 100, 200):
    print(f"k={k}")
    for p in (4, 6, 8):
        vals = []
        for Ns in (80, 140, 220, 320):
            v = C.R_of_w(*C.R_setup(mp.mpf(k), Ns, p), mp.mpf('0.7'))
            vals.append(v)
        ds = [mp.nstr(abs(vals[i+1]-vals[i]), 3) for i in range(len(vals)-1)]
        print(f"   p={p}: R={mp.nstr(vals[-1], 25)}  deltas {ds}")
