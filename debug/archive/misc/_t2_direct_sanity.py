"""Crude INDEPENDENT assembly check: T2=(8/pi) int int J directly by GL NxN over [0,1]^2
with Jsmart (no domain split; (0,0) rho^{3/2} gives ~4-5 digit error only). Confirms the
scheme assembly + anchor at low precision."""
import sys, time; sys.path.insert(0,'debug')
import mpmath as mp
mp.mp.dps=30
from _fastgl import fast_gl
from routeC_T2_corner_subtraction import Jsmart
def direct(N):
    xs,ws=fast_gl(N); tot=mp.mpf(0)
    for xi,wi in zip(xs,ws):
        s=(xi+1)/2; wsi=wi/2
        for xj,wj in zip(xs,ws):
            t=(xj+1)/2; wtj=wj/2
            tot+=wsi*wtj*Jsmart(s,t)
    return (8/mp.pi)*tot
ref=mp.mpf('0.3953557659017139641')
for N in (24,40):
    t0=time.time(); v=direct(N)
    print(f'N={N}: T2~{mp.nstr(v,14)}  |vs anchor|={mp.nstr(abs(v-ref),3)}  ({time.time()-t0:.0f}s)',flush=True)
