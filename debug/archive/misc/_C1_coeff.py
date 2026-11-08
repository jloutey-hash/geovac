"""Leading rho->0 coefficient of the twist: Phi(rho) = C1 rho + C2 rho^2 + ... ; compute C1.
c1(u) = lim_{rho->0} S(u,rho)/rho = F0 * u * sum_i m0(u, b_i^(0)),  F0=7/e,
  b_i^(0) = {sm, sm+1, 1-sm, 2-sm}, sm=(1-sqrt(1-4u))/2  (tm->0 limit of the 4 branches),
  m0(u,b)=int_0^inf j0(kb) P(u,k) dk.   C1 = int_0^{1/4} c1(u) u/sqrt(1-4u) du."""
import sys,time; sys.path.insert(0,'debug')
import mpmath as mp
import routeC_T2_coarea_precision as CP
import _watson_fibre as W
mp.mp.dps=34
F0=7/mp.e

def c1_formula(u):
    u=mp.mpf(u); sm=(1-mp.sqrt(1-4*u))/2
    bs=[sm, sm+1, 1-sm, 2-sm]
    return F0*u*sum(W.m_n(u,b,0) for b in bs)

# validate c1(u) formula vs Richardson limit of S(u,rho)/rho (from the data trend)
print("validate c1(u) = F0 u sum m0  vs  lim S/rho:",flush=True)
for u in ['0.15','0.05']:
    cf=c1_formula(u)
    # quick Richardson: S/rho at rho=0.0125, 0.00625 -> extrapolate (linear in rho)
    s1=CP.branch_fibre(mp.mpf(u),mp.mpf('0.0125'),M=8)/mp.mpf('0.0125')
    s2=CP.branch_fibre(mp.mpf(u),mp.mpf('0.00625'),M=8)/mp.mpf('0.00625')
    extrap=2*s2-s1   # Richardson for ~linear-in-rho correction
    print(f"  u={u}: c1_formula={mp.nstr(cf,14)}  S/rho-extrap={mp.nstr(extrap,14)}  |d|={mp.nstr(abs(cf-extrap),3)}",flush=True)

# C1 = int_0^{1/4} c1(u) u/sqrt(1-4u) du  (sin^2 map handles the u=1/4 endpoint)
def C1(Nu):
    umax=mp.mpf(1)/4; xs,ws=W.__dict__.get('_gl',None), None
    from _fastgl import fast_gl
    xs,ws=fast_gl(Nu); Hh=mp.pi/2; tot=mp.mpf(0)
    for xg,wg in zip(xs,ws):
        phi=Hh*(xg+1)/2; u=umax*mp.sin(phi)**2; du=umax*mp.sin(2*phi); wj=Hh*wg/2
        if u<mp.mpf('1e-7'): continue
        tot+=wj*du*c1_formula(u)*u/mp.sqrt(1-4*u)
    return tot
print("C1 = lim Phi(rho)/rho :",flush=True)
prev=None
for Nu in (24,48,96):
    t0=time.time(); v=C1(Nu); d='' if prev is None else f" |dNu|={mp.nstr(abs(v-prev),3)}"
    print(f"  Nu={Nu}: C1={mp.nstr(v,20)}{d}  ({time.time()-t0:.0f}s)",flush=True); prev=v
