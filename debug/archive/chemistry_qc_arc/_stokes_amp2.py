"""Confirm the closed-form amplitude 2|A|=(c1+b^2)^{3/2}/(4 sqrt(pi) sqrt(c1) b) at a 2nd point."""
import sys; sys.path.insert(0,'debug')
import mpmath as mp
from _fastgl import fast_gl
import _watson_fibre as W
mp.mp.dps=45
def Pcs(cs,k):
    d=mp.sqrt(cs*k*k+1); return cs*mp.e**(-d)*(d**-3+3*d**-4+3*d**-5)
def dcoeffs(cs,A,M):
    k0=max(mp.mpf(60),35/mp.sqrt(cs)); ks=[mp.mpf(k0)*(i+1) for i in range(M+1)]
    gs=[Pcs(cs,k)*mp.e**(A*k)*k**3 for k in ks]; V=mp.matrix(M+1,M+1)
    for i,k in enumerate(ks):
        for m in range(M+1): V[i,m]=1/k**m
    return list(mp.lu_solve(V,mp.matrix(gs)))
def m_hi(cs,b,n,d,A,K=mp.mpf(120),Nq=2600,M=12):
    xs,ws=fast_gl(Nq); bnd=mp.mpf(0)
    for x,w in zip(xs,ws):
        k=K*(x+1)/2; wk=K*w/2; j0=mp.sin(k*b)/(k*b) if k*b>1e-40 else mp.mpf(1)
        bnd+=wk*(k**(2*n))*j0*Pcs(cs,k)
    z=A-1j*b; t=mp.mpc(0)
    for j in range(M+1):
        o=2*n-3-j; t+=d[j]*z**(-o)*mp.gammainc(o,z*K)
    return bnd+(t/b).imag
def extract(cs,b,N=52):
    cs=mp.mpf(cs); b=mp.mpf(b); A=mp.sqrt(cs); d=dcoeffs(cs,A,12); Fn=W.F_taylor(N)
    a=[Fn[n]*m_hi(cs,b,n,d,A) for n in range(N+1)]
    zst=-(A-1j*b)**2; absz=abs(zst); th=mp.arg(zst)
    rows=[]; rhs=[]
    for n in range(38,N+1):
        y=a[n]*absz**n/mp.factorial(2*n)*mp.mpf(n)**(mp.mpf(5)/2)
        c=mp.cos(n*th); s=mp.sin(n*th); rows.append([c,s,c/n,s/n]); rhs.append(y)
    Mx=mp.matrix(rows); sol=mp.lu_solve(Mx.T*Mx, Mx.T*mp.matrix(rhs))
    amp=mp.sqrt(sol[0]**2+sol[1]**2)
    pred=(cs+b*b)**(mp.mpf(3)/2)/(4*mp.sqrt(mp.pi)*mp.sqrt(cs)*b)
    return amp,pred
for cs,b in [('0.1','0.5'),('0.3','0.5'),('0.2','1.0')]:
    amp,pred=extract(cs,b)
    print(f"(c1,b)=({cs},{b}): 2|A|_fit={mp.nstr(amp,10)}  closed={mp.nstr(pred,10)}  ratio={mp.nstr(amp/pred,8)}",flush=True)
