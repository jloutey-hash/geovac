import sys, time; sys.path.insert(0,'debug')
import mpmath as mp
import beta2_t2_kw_core as C
import beta2_t2_tail as T
mp.mp.dps = 60
for MX in (12, 20, 30):
    t0=time.time(); acc = T.F_buckets(MX); tb=time.time()-t0
    print(f"MX={MX}: buckets={len(acc)} build {tb:.1f}s")
    for k in (100, 200, 300):
        ex = C.F_of_k(mp.mpf(k), 360, int(0.8*k)+60, 6)
        ap = T.F_asym(k, acc)
        print(f"   k={k}: F={mp.nstr(ex,14)}  rel err {mp.nstr(abs(ap-ex)/abs(ex),3)}")
