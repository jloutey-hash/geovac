import sys, time; sys.path.insert(0,'debug')
import mpmath as mp
import beta2_t2_kw_core as C
import beta2_t2_tail as T
mp.mp.dps = 100
for MX in (90, 120):
    t0=time.time(); acc = T.F_buckets(MX); tb=time.time()-t0
    print(f"MX={MX}: buckets={len(acc)} build {tb:.0f}s", flush=True)
    for k in (250, 320):
        ex = C.F_of_k(mp.mpf(k), 170+int(0.5*k), 80+int(0.8*k), 6)
        ap = T.F_asym(k, acc)
        print(f"   k={k}: rel err {mp.nstr(abs(ap-ex)/abs(ex),3)}   TAIL({k})={mp.nstr(T.TAIL(k,acc),12)}", flush=True)
