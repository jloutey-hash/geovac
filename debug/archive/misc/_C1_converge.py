import sys,time; sys.path.insert(0,'debug')
import mpmath as mp
import _watson_fibre as W
from _fastgl import fast_gl
mp.mp.dps=34
F0=7/mp.e
def c1_formula(u):
    u=mp.mpf(u); sm=(1-mp.sqrt(1-4*u))/2
    return F0*u*sum(W.m_n(u,b,0) for b in [sm,sm+1,1-sm,2-sm])
def C1(Nu):
    umax=mp.mpf(1)/4; xs,ws=fast_gl(Nu); Hh=mp.pi/2; tot=mp.mpf(0)
    for xg,wg in zip(xs,ws):
        phi=Hh*(xg+1)/2; u=umax*mp.sin(phi)**2; du=umax*mp.sin(2*phi); wj=Hh*wg/2
        if u<mp.mpf('1e-7'): continue
        tot+=wj*du*c1_formula(u)*u/mp.sqrt(1-4*u)
    return tot
prev=None
for Nu in (48,96,160):
    t0=time.time(); v=C1(Nu); d='' if prev is None else f" |dNu|={mp.nstr(abs(v-prev),3)}"
    print(f"Nu={Nu}: C1={mp.nstr(v,22)}{d}  ({time.time()-t0:.0f}s)",flush=True); prev=v
