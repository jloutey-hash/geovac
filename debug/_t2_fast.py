"""Fast T2 via ADAPTIVE fibre Nk: interior spectral at Nk~200, corner needs ~1000.
Nk_eff(s,t) = clamp(round(A/sqrt(c_min)), Nk_lo, Nk_hi). Validated vs uniform-Nk anchor."""
import mpmath as mp, sys, time
sys.path.insert(0,'debug'); import routeC_T2_highprec as H

_A=[35]; _LO=[160]; _HI=[1400]
def set_sched(A,lo,hi): _A[0]=A; _LO[0]=lo; _HI[0]=hi

_GRID=[160,220,300,400,520,660,820,1000,1200,1400]  # small cached set of Nk values
def nk_of(cs,ct):
    cmin=min(cs,ct)
    if cmin<=0: return _GRID[-1]
    want=int(_A[0]/mp.sqrt(cmin))
    for g in _GRID:
        if g>=want: return g
    return _GRID[-1]

def Ja(s,t):
    cs=s*(1-s); ct=t*(1-t)
    if cs==0 or ct==0: return mp.mpf(0)
    Nk=nk_of(cs,ct)
    return H.J(s,t,Nk)

def cornera(delta,Ns,Na):
    Hh=mp.pi/2; xs,ws=H.gl(Ns); xa,wa=H.gl(Na); smax=mp.sqrt(delta); tot=mp.mpf(0)
    for xi,wi in zip(xs,ws):
        sig=smax*(xi+1)/2; wsig=smax*wi/2
        for xj,wj in zip(xa,wa):
            psi=Hh*(xj+1)/2; wpsi=Hh*wj/2; al=mp.sin(psi)**2; jal=mp.sin(2*psi)
            s=sig*sig*al; t=sig*sig*(1-al)
            tot+=wsig*wpsi*jal*2*sig**3*Ja(s,t)
    return tot

def trapa(delta,No,Ni):
    Hh=mp.pi/2; xs,ws=H.gl(No); xa,wa=H.gl(Ni); tot=mp.mpf(0)
    for xi,wi in zip(xs,ws):
        phis=Hh*(xi+1)/2; s=delta*mp.sin(phis)**2; js=delta*mp.sin(2*phis); wphis=Hh*wi/2; lo=delta-s
        for xj,wj in zip(xa,wa):
            phit=Hh*(xj+1)/2; t=lo+(1-lo)*mp.sin(phit)**2; jt=(1-lo)*mp.sin(2*phit); wphit=Hh*wj/2
            tot+=wphis*js*wphit*jt*Ja(s,t)
    for xi,wi in zip(xs,ws):
        phis=Hh*(xi+1)/2; s=delta+(1-delta)*mp.sin(phis)**2; js=(1-delta)*mp.sin(2*phis); wphis=Hh*wi/2
        for xj,wj in zip(xa,wa):
            phit=Hh*(xj+1)/2; t=mp.sin(phit)**2; jt=mp.sin(2*phit); wphit=Hh*wj/2
            tot+=wphis*js*wphit*jt*Ja(s,t)
    return tot

def T2a(delta,Nc):
    return (8/mp.pi)*(cornera(delta,Nc,Nc)+trapa(delta,Nc,Nc))

if __name__=='__main__':
    mp.mp.dps=45; ref=mp.mpf('0.3953557659017139641')
    # validate against 19-digit anchor at growing Nc, adaptive Nk
    set_sched(40,160,1400)
    prev=None
    for Nc in [40,52,64,76]:
        t0=time.time(); v=T2a(mp.mpf('0.08'),Nc); dt=time.time()-t0
        d='' if prev is None else f'  |dprev|={mp.nstr(abs(v-prev),3)}'
        print(f'Nc={Nc} adaptiveNk: {mp.nstr(v,30)}  ({dt:.0f}s){d}  |anch|={mp.nstr(abs(v-ref),3)}',flush=True)
        prev=v
    print('DONE',flush=True)
