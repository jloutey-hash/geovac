import sys, time; sys.path.insert(0,'debug')
import mpmath as mp
import beta2_t2_kw_core as C
import beta2_t2_tail as T
mp.mp.dps = 90
for MX in (40, 56, 70):
    t0=time.time(); acc = T.F_buckets(MX); tb=time.time()-t0
    print(f"MX={MX}: buckets={len(acc)} build {tb:.1f}s", flush=True)
    for k in (150, 200, 260):
        ex = C.F_of_k(mp.mpf(k), 460, int(1.0*k)+80, 6)
        ap = T.F_asym(k, acc)
        print(f"   k={k}: rel err {mp.nstr(abs(ap-ex)/abs(ex),3)}", flush=True)
