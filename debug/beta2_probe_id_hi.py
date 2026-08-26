import sys, time; sys.path.insert(0,'debug')
import mpmath as mp
import beta2_t2_direct2d as D2
import beta2_t2_kw_core as C
mp.mp.dps = 95
for k in ('2', '50', '150', '300'):
    kk = mp.mpf(k)
    a = D2.F_direct(kk, 260, 6); b = D2.F_direct(kk, 360, 6)
    c = C.F_of_k(kk, 300+int(0.95*float(k)), 60+int(0.7*float(k)), 6)
    print(f"k={k:>4}: F={mp.nstr(b,30)}  2D selfconv {mp.nstr(abs(b-a)/abs(b),3)}"
          f"   |F_2D - F_KW|/F = {mp.nstr(abs(b-c)/abs(b),3)}", flush=True)
