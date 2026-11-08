"""Twist-free sibling of T2: T2_0 = (8/pi) int int [int_0^inf P(s,k)P(t,k) dk] ds dt (j0->1).
The pure length-2 X(2) modular integral WITHOUT the b=s+t character twist. Faster (no
oscillation); if it closes in {pi, varpi, G} it validates the ring and isolates the twist."""
import mpmath as mp, sys, time
sys.path.insert(0,'debug'); import routeC_T2_highprec as H

def Jnt(s, t, Nk):
    cs=s*(1-s); ct=t*(1-t)
    if cs==0 or ct==0: return mp.mpf(0)
    L=1/(mp.sqrt(cs)+mp.sqrt(ct))
    xs,ws=H.fiber_nodes(Nk); tot=mp.mpf(0)
    for x,w in zip(xs,ws):
        u=(x+1)/2; k=L*u/(1-u); dk=L/(1-u)**2
        tot+=(w/2)*H.P(s,k)*H.P(t,k)*dk
    return tot

def cornernt(delta,Ns,Na,Nk):
    Hh=mp.pi/2; xs,ws=H.gl(Ns); xa,wa=H.gl(Na); smax=mp.sqrt(delta); tot=mp.mpf(0)
    for xi,wi in zip(xs,ws):
        sig=smax*(xi+1)/2; wsig=smax*wi/2
        for xj,wj in zip(xa,wa):
            psi=Hh*(xj+1)/2; wpsi=Hh*wj/2; al=mp.sin(psi)**2; jal=mp.sin(2*psi)
            s=sig*sig*al; t=sig*sig*(1-al)
            tot+=wsig*wpsi*jal*2*sig**3*Jnt(s,t,Nk)
    return tot

def trapnt(delta,No,Ni,Nk):
    Hh=mp.pi/2; xs,ws=H.gl(No); xa,wa=H.gl(Ni); tot=mp.mpf(0)
    for xi,wi in zip(xs,ws):
        phis=Hh*(xi+1)/2; s=delta*mp.sin(phis)**2; js=delta*mp.sin(2*phis); wphis=Hh*wi/2; lo=delta-s
        for xj,wj in zip(xa,wa):
            phit=Hh*(xj+1)/2; t=lo+(1-lo)*mp.sin(phit)**2; jt=(1-lo)*mp.sin(2*phit); wphit=Hh*wj/2
            tot+=wphis*js*wphit*jt*Jnt(s,t,Nk)
    for xi,wi in zip(xs,ws):
        phis=Hh*(xi+1)/2; s=delta+(1-delta)*mp.sin(phis)**2; js=(1-delta)*mp.sin(2*phis); wphis=Hh*wi/2
        for xj,wj in zip(xa,wa):
            phit=Hh*(xj+1)/2; t=mp.sin(phit)**2; jt=mp.sin(2*phit); wphit=Hh*wj/2
            tot+=wphis*js*wphit*jt*Jnt(s,t,Nk)
    return tot

def T2nt(delta,Nc,Nk):
    return (8/mp.pi)*(cornernt(delta,Nc,Nc,Nk)+trapnt(delta,Nc,Nc,Nk))

if __name__=='__main__':
    mp.mp.dps=40; d=mp.mpf('0.08'); prev=None
    for Nc,Nk in [(40,300),(56,400),(72,500)]:
        t0=time.time(); v=T2nt(d,Nc,Nk); dt=time.time()-t0
        dd='' if prev is None else f'  |dprev|={mp.nstr(abs(v-prev),3)}'
        print(f'Nc={Nc} Nk={Nk}: {mp.nstr(v,32)}  ({dt:.0f}s){dd}',flush=True); prev=v
    print('DONE',flush=True)
