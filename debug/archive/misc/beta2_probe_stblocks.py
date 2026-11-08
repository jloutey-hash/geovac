import sys, time; sys.path.insert(0,'debug')
import mpmath as mp
import beta2_t2_st_route as ST
from _fastgl import fast_gl
mp.mp.dps=50
H=mp.pi/2; delta=mp.mpf('0.08')
def timed(tag, fn):
    t0=time.time(); v=fn(); print(f"  {tag}: {mp.nstr(v,25)}  ({time.time()-t0:.1f}s)", flush=True)
    return v
N=12
def corner():
    xs,ws=fast_gl(N); xa,wa=fast_gl(N); smax=mp.sqrt(delta); tot=mp.mpf(0)
    for xi,wi in zip(xs,ws):
        sig=smax*(xi+1)/2; wsig=smax*wi/2
        for xj,wj in zip(xa,wa):
            psi=H*(xj+1)/2; wpsi=H*wj/2; al=mp.sin(psi)**2; jal=mp.sin(2*psi)
            tot += wsig*wpsi*jal*2*sig**3*ST.J_acc(sig*sig*al, sig*sig*(1-al))
    return tot
def blk1():
    xs,ws=fast_gl(N); xa,wa=fast_gl(N); tot=mp.mpf(0)
    for xi,wi in zip(xs,ws):
        phis=H*(xi+1)/2; s=delta*mp.sin(phis)**2; js=delta*mp.sin(2*phis); wphis=H*wi/2; lo=delta-s
        for xj,wj in zip(xa,wa):
            phit=H*(xj+1)/2; t=lo+(1-lo)*mp.sin(phit)**2; jt=(1-lo)*mp.sin(2*phit); wphit=H*wj/2
            tot += wphis*js*wphit*jt*ST.J_acc(s,t)
    return tot
def blk2():
    xs,ws=fast_gl(N); xa,wa=fast_gl(N); tot=mp.mpf(0)
    for xi,wi in zip(xs,ws):
        phis=H*(xi+1)/2; s=delta+(1-delta)*mp.sin(phis)**2; js=(1-delta)*mp.sin(2*phis); wphis=H*wi/2
        for xj,wj in zip(xa,wa):
            phit=H*(xj+1)/2; t=mp.sin(phit)**2; jt=mp.sin(2*phit); wphit=H*wj/2
            tot += wphis*js*wphit*jt*ST.J_acc(s,t)
    return tot
timed('corner', corner); timed('blk1', blk1); timed('blk2', blk2)
