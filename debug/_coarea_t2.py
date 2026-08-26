"""End-to-end co-area 1D T2 at moderate precision: validate the pipeline vs the 19-digit anchor.
Guard deep-u (negligible, contributes ~ guard^1) and cap Nq for speed; log-subtract the rho=1 cusp."""
import sys,time; sys.path.insert(0,'debug')
import mpmath as mp
import routeC_T2_coarea_precision as CP
from _fastgl import fast_gl
mp.mp.dps=30
ANCHOR='0.3953557659017139641'

# moderate-precision branch fibre: guard deep-u, cap Nq
def bf(u,rho): return CP.branch_fibre(u,rho,M=6,guard=mp.mpf('1e-4'))

def Phi(rho,Nu):
    rho=mp.mpf(rho); umax=min(mp.mpf(1)/4,1/(4*rho))
    xs,ws=fast_gl(Nu); Hh=mp.pi/2; tot=mp.mpf(0)
    for xg,wg in zip(xs,ws):
        phi=Hh*(xg+1)/2; u=umax*mp.sin(phi)**2; du=umax*mp.sin(2*phi); wj=Hh*wg/2
        b=bf(u,rho)
        tot+=wj*du*b*u/(mp.sqrt(1-4*u)*mp.sqrt(1-4*rho*u))
    return tot

# A for log-subtraction at rho=1 (Phi ~ -A ln(1-rho)+C)
A=mp.mpf('0.07905403211681687671461246470412103')
def T2_coarea(Nu,Nrho):
    # int_0^1 Phi drho = A + int_0^1 [Phi + A ln(1-rho)] drho ; sub rho=1-v^2 to smooth the log
    xs,ws=fast_gl(Nrho); tot=mp.mpf(0)
    for xg,wg in zip(xs,ws):
        v=(xg+1)/2                     # v in [0,1]; rho=1-v^2, drho=2v dv, ln(1-rho)=2 ln v
        rho=1-v*v
        if rho<=0 or rho>=1: continue
        integ=(Phi(rho,Nu)+A*2*mp.log(v))*2*v
        tot+=(wg/2)*integ
    intPhi=A+tot
    return (16/mp.pi)*intPhi

for Nu,Nrho in [(40,32),(60,48)]:
    t0=time.time(); v=T2_coarea(Nu,Nrho); dt=time.time()-t0
    print(f"Nu={Nu} Nrho={Nrho}: T2={mp.nstr(v,22)}  |v-anchor|={mp.nstr(abs(v-mp.mpf(ANCHOR)),3)}  ({dt:.0f}s)",flush=True)
