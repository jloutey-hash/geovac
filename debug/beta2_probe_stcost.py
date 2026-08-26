import sys, time; sys.path.insert(0,'debug')
import mpmath as mp
import beta2_t2_st_route as ST
mp.mp.dps=50
t0=time.time(); v=ST.T2_st(mp.mpf('0.08'), 16, 16, 10); print(f"N=16: {mp.nstr(v,30)} ({time.time()-t0:.0f}s)", flush=True)
t0=time.time(); v=ST.T2_st(mp.mpf('0.08'), 24, 24, 10); print(f"N=24: {mp.nstr(v,30)} ({time.time()-t0:.0f}s)", flush=True)
