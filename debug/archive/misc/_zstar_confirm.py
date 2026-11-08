"""Confirm z* = -(sqrt(c_s)-i b)^2 (|z*|=c_s+b^2): Pade trend in N + second (c_s,b) point.
Re(z*) is the clean coordinate (branch-cut slows |z*|). Show Re(z*)-> b^2-c_s and |z*|-> c_s+b^2."""
import sys,time; sys.path.insert(0,'debug')
import mpmath as mp
from _fastgl import fast_gl
import _watson_fibre as W
import importlib.util
spec=importlib.util.spec_from_file_location("za","debug/_zstar_analytic.py")
# reuse the moment engine by re-defining inline (avoid running its main)
mp.mp.dps=50
def Pcs(cs,k):
    d=mp.sqrt(cs*k*k+1); return cs*mp.e**(-d)*(d**-3+3*d**-4+3*d**-5)
def dcoeffs(cs,A,M):
    k0=max(mp.mpf(60),35/mp.sqrt(cs)); ks=[mp.mpf(k0)*(i+1) for i in range(M+1)]
    gs=[Pcs(cs,k)*mp.e**(A*k)*k**3 for k in ks]; V=mp.matrix(M+1,M+1)
    for i,k in enumerate(ks):
        for m in range(M+1): V[i,m]=1/k**m
    dd=mp.lu_solve(V,mp.matrix(gs)); return [dd[m] for m in range(M+1)]
def m_hi(cs,b,n,K,Nq,d,A,M):
    xs,ws=fast_gl(Nq); bnd=mp.mpf(0)
    for x,w in zip(xs,ws):
        k=K*(x+1)/2; wk=K*w/2; j0=mp.sin(k*b)/(k*b) if k*b>1e-40 else mp.mpf(1)
        bnd+=wk*(k**(2*n))*j0*Pcs(cs,k)
    z=A-1j*b; t=mp.mpc(0)
    for j in range(M+1):
        o=2*n-3-j; t+=d[j]*z**(-o)*mp.gammainc(o,z*K)
    return bnd+(t/b).imag
def zstar_pade(cs,b,N,L,M):
    cs=mp.mpf(cs); b=mp.mpf(b); A=mp.sqrt(cs); d=dcoeffs(cs,A,12)
    Fn=W.F_taylor(N); ms=[m_hi(cs,b,n,mp.mpf(120),2600,d,A,12) for n in range(N+1)]
    a=[Fn[n]*ms[n] for n in range(N+1)]; bc=[a[n]/mp.factorial(2*n) for n in range(N+1)]
    p,q=mp.pade(bc,L,M); r=sorted(mp.polyroots(q[::-1],maxsteps=400,extraprec=300),key=lambda z:abs(z.imag) if abs(z.imag)>1e-6 else 1e9)
    # pick nearest genuinely-complex root
    cr=[z for z in r if abs(z.imag)>1e-4]
    return cr[0] if cr else r[0]
for cs,b in [('0.2','0.5'),('0.1','0.5')]:
    cf=-(mp.sqrt(mp.mpf(cs))-1j*mp.mpf(b))**2
    print(f"\n(c_s,b)=({cs},{b}): closed form z*={mp.nstr(cf,10)}  Re={mp.nstr(cf.real,8)} |z*|=c_s+b^2={mp.nstr(mp.mpf(cs)+mp.mpf(b)**2,8)}",flush=True)
    for N,L,M in [(30,15,15),(48,24,24),(60,28,32)]:
        t0=time.time(); z=zstar_pade(cs,b,N,L,M)
        print(f"  N={N} Pade[{L}/{M}]: z*={mp.nstr(z,10)}  Re={mp.nstr(z.real,7)} |z*|={mp.nstr(abs(z),8)}  ({time.time()-t0:.0f}s)",flush=True)
