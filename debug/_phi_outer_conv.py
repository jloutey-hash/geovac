"""Decisive feasibility test: does the co-area u-integration for Phi(rho) at MODERATE rho converge
past ~14 digits (with exact analytic-tail fibre + a u->0 guard), or does it cap like the 2D (s,t)?"""
import sys,time; sys.path.insert(0,'debug')
import mpmath as mp
import routeC_T2_coarea_precision as CP
from _fastgl import fast_gl
mp.mp.dps=44
def Phi_guarded(rho,Nu,guard):
    rho=mp.mpf(rho); umax=min(mp.mpf(1)/4,1/(4*rho))
    xs,ws=fast_gl(Nu); Hh=mp.pi/2; tot=mp.mpf(0)
    for xg,wg in zip(xs,ws):
        phi=Hh*(xg+1)/2; u=umax*mp.sin(phi)**2; du=umax*mp.sin(2*phi); wj=Hh*wg/2
        bf=CP.branch_fibre(u,rho,M=8,guard=mp.mpf(guard))
        tot+=wj*du*bf*u/(mp.sqrt(1-4*u)*mp.sqrt(1-4*rho*u))
    return tot
for rho in ['0.3']:
    print(f"rho={rho}, guard=1e-4:",flush=True); prev=None
    for Nu in (40,80,160,320):
        t0=time.time(); v=Phi_guarded(rho,Nu,'1e-4'); dt=time.time()-t0
        d='' if prev is None else f"  |dNu|={mp.nstr(abs(v-prev),3)}"
        print(f"  Nu={Nu:4d}: {mp.nstr(v,36)}{d}  ({dt:.0f}s)",flush=True); prev=v
