"""Validate the CLOSED-FORM Borel singularity z* = -(sqrt(c_s) -+ i b)^2 (|z*|=c_s+b^2).
High-n moments via analytic tail: m_n=int_0^K + int_K^inf, tail =
  (1/b) sum_j d_j Im[ z^{-(2n-3-j)} Gamma(2n-3-j, zK) ],  z=A-ib, A=sqrt(c_s),
  d_j = 1/k-series coeffs of G(k)=P(c_s,k) e^{Ak} k^3 (-> c_s^{-1/2} at k->inf).
Compare converged Borel-Pade z* to the closed form."""
import sys,time; sys.path.insert(0,'debug')
import mpmath as mp
from _fastgl import fast_gl
import _watson_fibre as W
mp.mp.dps=50

def Pcs(cs,k):
    d=mp.sqrt(cs*k*k+1); return cs*mp.e**(-d)*(d**-3+3*d**-4+3*d**-5)

def dcoeffs(cs,A,M):
    k0=max(mp.mpf(60),35/mp.sqrt(cs))
    ks=[mp.mpf(k0)*(i+1) for i in range(M+1)]
    gs=[Pcs(cs,k)*mp.e**(A*k)*k**3 for k in ks]
    V=mp.matrix(M+1,M+1)
    for i,k in enumerate(ks):
        for m in range(M+1): V[i,m]=1/k**m
    dd=mp.lu_solve(V,mp.matrix(gs)); return [dd[m] for m in range(M+1)]

def m_tail(cs,b,n,K,d,A,M):
    z=A-1j*b; tot=mp.mpc(0)
    for j in range(M+1):
        order=2*n-3-j
        tot+=d[j]*z**(-order)*mp.gammainc(order, z*K)
    return (tot/b).imag

def m_bounded(cs,b,n,K,Nq):
    xs,ws=fast_gl(Nq); tot=mp.mpf(0)
    for x,w in zip(xs,ws):
        k=K*(x+1)/2; wk=K*w/2
        j0=mp.sin(k*b)/(k*b) if k*b>1e-40 else mp.mpf(1)
        tot+=wk*(k**(2*n))*j0*Pcs(cs,k)
    return tot

def m_hi(cs,b,n,K=mp.mpf(120),Nq=2600,M=10,d=None,A=None):
    if d is None: A=mp.sqrt(cs); d=dcoeffs(cs,A,M)
    return m_bounded(cs,b,n,K,Nq)+m_tail(cs,b,n,K,d,A,M)

cs=mp.mpf('0.2'); b=mp.mpf('0.5'); N=36
A=mp.sqrt(cs); d=dcoeffs(cs,A,12)
# validate m_hi vs quadosc at low n
print("validate m_hi vs quadosc (low n):",flush=True)
for n in (0,2,5):
    mh=m_hi(cs,b,n,d=d,A=A,M=12); mq=W.m_n(cs,b,n)
    print(f"  n={n}: m_hi={mp.nstr(mh,18)} |m_hi-quadosc|={mp.nstr(abs(mh-mq),3)}",flush=True)
t0=time.time()
Fn=W.F_taylor(N); ms=[m_hi(cs,b,n,d=d,A=A,M=12) for n in range(N+1)]
a=[Fn[n]*ms[n] for n in range(N+1)]
bc=[a[n]/mp.factorial(2*n) for n in range(N+1)]
print(f"\nbuilt a_0..a_{N} in {time.time()-t0:.0f}s; Borel-Pade z*:",flush=True)
for L,M in [(16,16),(18,18),(14,18)]:
    try:
        p,q=mp.pade(bc,L,M)
        r=sorted(mp.polyroots(q[::-1],maxsteps=400,extraprec=300),key=lambda z:abs(z))
        print(f"  Pade[{L}/{M}]: z*={mp.nstr(r[0],12)}  |z*|={mp.nstr(abs(r[0]),10)}",flush=True)
    except Exception as e: print(f"  Pade[{L}/{M}] fail {e}",flush=True)
zcf=-(mp.sqrt(cs)-1j*b)**2
print(f"\nCLOSED FORM  -(sqrt(c_s)-i b)^2 = {mp.nstr(zcf,12)}  |.|=c_s+b^2={mp.nstr(cs+b*b,10)}",flush=True)
