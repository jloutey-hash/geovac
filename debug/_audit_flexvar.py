import numpy as np, tc_flex_variance_he as F
from scipy.optimize import minimize
P=F.setup(); M=len(F.GAMMAS)
def f2c(free): return np.append(free,1-np.sum(free))
warm=np.array([-3.577,4.0,2.248,-0.463])   # the CMAX=4 optimum (free part)
print("\nCMAX sweep (does it run away toward spurious large-J, or stabilize?):")
print("%5s %8s %11s %9s %10s %6s"%("CMAX","max|c|","Var","dEw_mHa","dTC_mHa","kV"))
prev=warm
for CMAX in [4,6,10,20]:
    def obj(free):
        cs=f2c(free); v,_=F.variance(P,cs)
        return v+1e3*np.sum(np.maximum(0,np.abs(cs)-CMAX)**2)
    best=None
    for s0 in [prev, np.zeros(M-1)+0.0]:
        r=minimize(obj,s0,method='Nelder-Mead',options=dict(xatol=1e-6,fatol=1e-13,maxiter=8000,maxfev=12000))
        if best is None or r.fun<best.fun: best=r
    cs=f2c(best.x); d=F.evaluate(P,cs); prev=best.x
    print("%5d %8.2f %11.4e %+9.2f %+10.3f %6.2f"%(CMAX,np.max(np.abs(cs)),d['variance'],d['dEw_mHa'],d['dTC_mHa'],d['kV']))
    if CMAX==6:
        # u(r) shape for this optimum: u(r)=sum c_mu (-1/2g)e^{-g r}
        rr=np.array([0.0,0.5,1.0,2.0,4.0])
        u=np.zeros_like(rr)
        for c,g in zip(cs,F.GAMMAS): u+= c*(-1/(2*g))*np.exp(-g*rr)
        up0=0.5*np.sum(cs)  # u'(0)=1/2*sum c = cusp
        print("   u(r) at r=[0,.5,1,2,4] =",np.round(u,3),"  u'(0)=%.3f (cusp=0.5)"%up0)
