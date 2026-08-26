"""Lean high-precision T2 with fast_gl patch. Nc-convergence at delta=0.08 + delta=0.05 cross-check."""
import mpmath as mp, sys, time
sys.path.insert(0,'debug')
import routeC_T2_highprec as H
from _fastgl import fast_gl
# patch the slow gl with fast_gl (order-independent sums => safe; validated)
H.gl = fast_gl
def _fn(Nk):
    if Nk not in H._FIB: H._FIB[Nk]=fast_gl(Nk)
    return H._FIB[Nk]
H._FIB.clear(); H.fiber_nodes=_fn
import _t2_fast as F   # imports H.J etc; uses patched gl via H

mp.mp.dps=50
F.set_sched(50,160,1400)
anchor=mp.mpf('0.3953557659017139641')
def show(tag,delta,Nc):
    t0=time.time(); v=F.T2a(mp.mpf(delta),Nc); dt=time.time()-t0
    print(f'{tag} delta={delta} Nc={Nc}: {mp.nstr(v,44)}  ({dt:.0f}s)  |anch|={mp.nstr(abs(v-anchor),3)}',flush=True)
    return v
a56=show('A','0.08',56)
a72=show('A','0.08',72)
print(f'  outer conv |a72-a56| = {mp.nstr(abs(a72-a56),3)}',flush=True)
b72=show('B','0.05',72)
print(f'  CROSS-DELTA |a72(.08)-b72(.05)| = {mp.nstr(abs(a72-b72),3)}  <-- validated digits',flush=True)
print(f'T2_BEST = {mp.nstr(a72,44)}',flush=True)
print('DONE',flush=True)
