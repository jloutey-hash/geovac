"""Exact algebraic Stokes amplitude: with z* and p=-5/2 known, s_n*n^{5/2}=P cos(n th)+Q sin(n th)+
(P1 cos+Q1 sin)/n (subleading), th=arg(z*). Linear LSQ over a high-n window -> 2|A|=sqrt(P^2+Q^2).
Then the correction factor 2|A|/leading (leading=[3/Gamma(5/2)][1/(sqrt(c1) b)](c1+b^2)^{3/2}/16)."""
import sys,time; sys.path.insert(0,'debug')
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
    dd=mp.lu_solve(V,mp.matrix(gs)); return [dd[m] for m in range(M+1)]
def m_hi(cs,b,n,d,A,K=mp.mpf(120),Nq=2600,M=12):
    xs,ws=fast_gl(Nq); bnd=mp.mpf(0)
    for x,w in zip(xs,ws):
        k=K*(x+1)/2; wk=K*w/2; j0=mp.sin(k*b)/(k*b) if k*b>1e-40 else mp.mpf(1)
        bnd+=wk*(k**(2*n))*j0*Pcs(cs,k)
    z=A-1j*b; t=mp.mpc(0)
    for j in range(M+1):
        o=2*n-3-j; t+=d[j]*z**(-o)*mp.gammainc(o,z*K)
    return bnd+(t/b).imag
cs=mp.mpf('0.2'); b=mp.mpf('0.5'); N=54
A=mp.sqrt(cs); d=dcoeffs(cs,A,12); Fn=W.F_taylor(N)
a=[Fn[n]*m_hi(cs,b,n,d,A) for n in range(N+1)]
zst=-(A-1j*b)**2; absz=abs(zst); th=mp.arg(zst)
# s_n n^{5/2}, fit P cos + Q sin + (P1 cos+Q1 sin)/n over window
def fit(nlo,nhi):
    rows=[]; rhs=[]
    for n in range(nlo,nhi+1):
        sn=a[n]*absz**n/mp.factorial(2*n); y=sn*mp.mpf(n)**(mp.mpf(5)/2)
        c=mp.cos(n*th); s=mp.sin(n*th)
        rows.append([c,s,c/n,s/n]); rhs.append(y)
    M_=mp.matrix(rows); r=mp.matrix(rhs)
    # normal equations
    sol=mp.lu_solve(M_.T*M_, M_.T*r)
    P,Q=sol[0],sol[1]; return mp.sqrt(P*P+Q*Q)
for w in [(30,42),(36,48),(40,54)]:
    amp=fit(*w); print(f"window {w}: 2|A|=sqrt(P^2+Q^2)={mp.nstr(amp,14)}",flush=True)
lead=(3/mp.gamma(mp.mpf(5)/2))*(1/(mp.sqrt(cs)*b))*(cs+b*b)**(mp.mpf(3)/2)/16
amp=fit(40,54)
print(f"\nleading 2|A|_lead = {mp.nstr(lead,12)}",flush=True)
print(f"exact/leading factor = {mp.nstr(amp/lead,12)}",flush=True)
# PSLQ the factor against small algebraic candidates
fac=amp/lead
for name,val in [('1',mp.mpf(1)),('c1+b2',cs+b*b),('b2/(c1+b2)',b*b/(cs+b*b)),('c1/(c1+b2)',cs/(cs+b*b))]:
    print(f"  factor vs {name}={mp.nstr(val,8)}: ratio {mp.nstr(fac/val,10)}",flush=True)
