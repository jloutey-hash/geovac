"""Validate the Watson series for R against the exact numerical R."""
import sys; sys.path.insert(0,'debug')
import mpmath as mp
import beta2_t2_kw_core as C
import beta2_t2_tail as T
mp.mp.dps = 60
for MX in (6, 12, 20):
    A, B = T.build_AB(MX)
    print(f"MX={MX}: len A={len(A)} len B={len(B)}")
    for k in (60, 100, 200, 300):
        for w in ('0.3','0.85'):
            ww = mp.mpf(w)
            ex = C.R_of_w(*C.R_setup(mp.mpf(k), 320, 6), ww)
            ap = T.R_asym(mp.mpf(k), ww, A, B)
            print(f"   k={k:4d} w={w}: R={mp.nstr(ex,12)} rel err {mp.nstr(abs(ap-ex)/abs(ex),3)}")
