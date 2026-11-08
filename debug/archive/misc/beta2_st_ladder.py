import sys, time; sys.path.insert(0,'debug')
import mpmath as mp
import beta2_t2_st_route as ST
mp.mp.dps = int(sys.argv[1]) if len(sys.argv)>1 else 50
delta = mp.mpf(sys.argv[2]) if len(sys.argv)>2 else mp.mpf('0.08')
M = int(sys.argv[3]) if len(sys.argv)>3 else 10
print(f"dps={mp.mp.dps} delta={delta} M={M}")
prev=None
for N in (20, 28, 36, 44):
    t0=time.time(); v = ST.T2_st(delta, N, N, M); dt=time.time()-t0
    d = '' if prev is None else f"  |dprev| {mp.nstr(abs(v-prev),3)}"
    print(f"  N={N:3d}: {mp.nstr(v, mp.mp.dps-6)}  ({dt:.0f}s){d}", flush=True)
    prev=v
