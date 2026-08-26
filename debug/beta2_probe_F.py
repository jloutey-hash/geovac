import sys; sys.path.insert(0,'debug')
import mpmath as mp, time
import beta2_t2_kw_core as C
mp.mp.dps = 60
for k in (5, 40, 100, 200):
    kk = mp.mpf(k); vals=[]; times=[]
    for Nw in (int(0.35*k)+30, int(0.55*k)+40, int(0.8*k)+60, int(1.2*k)+90):
        t0=time.time(); v = C.F_of_k(kk, 240, Nw, 6); times.append(time.time()-t0); vals.append(v)
    ds=[mp.nstr(abs(vals[i+1]-vals[i]),3) for i in range(len(vals)-1)]
    print(f"k={k:4d} F={mp.nstr(vals[-1],25)}  deltas {ds}  t={['%.2f'%t for t in times]}")
