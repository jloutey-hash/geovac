import numpy as np, tc_flex_variance_he as F
from scipy.optimize import minimize
P=F.setup(); M=len(F.GAMMAS); rg=P["r"]
amp=P["amp"]
def f2c(free): return np.append(free,1-np.sum(free))
def u_of(cs):
    u=np.zeros_like(rg)
    for c,a,g in zip(cs,amp,F.GAMMAS): u+= c*a*np.exp(-g*rg)
    return u
# physical penalty: u(r) <= 0 everywhere (correlation hole); weight by r^2 dr
w=rg*rg*P["wr"]
def phys_pen(cs):
    u=u_of(cs); return np.sum(w*np.maximum(0.0,u)**2)
print("beta   max|c|   Var        dEw_mHa   dTC_mHa   kV     u(0)  u(1)  u(2)  physpen")
for beta in [0.0, 1.0, 10.0, 100.0]:
    def obj(free):
        cs=f2c(free); v,_=F.variance(P,cs)
        return v + beta*phys_pen(cs) + 1e3*np.sum(np.maximum(0,np.abs(cs)-6)**2)
    best=None
    for s0 in [np.array([0.,0.,1.,0.]), np.array([0.,1.,0.,0.]), np.zeros(M-1)]:
        r=minimize(obj,s0,method='Nelder-Mead',options=dict(xatol=1e-6,fatol=1e-13,maxiter=8000,maxfev=12000))
        if best is None or r.fun<best.fun: best=r
    cs=f2c(best.x); d=F.evaluate(P,cs); u=u_of(cs)
    ur=lambda x: float(np.interp(x,rg,u))
    print("%5.0f  %6.2f  %.4e  %+7.2f  %+8.3f  %5.2f  %+5.2f %+5.2f %+5.2f  %.2e"%(
        beta,np.max(np.abs(cs)),d['variance'],d['dEw_mHa'],d['dTC_mHa'],d['kV'],ur(0),ur(1),ur(2),phys_pen(cs)))
