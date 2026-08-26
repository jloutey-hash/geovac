"""Next ingredient: the Stokes DATA at z*=-(sqrt(c_s)-ib)^2. With z* known, remove the exponential:
s_n = a_n * |z*|^n / (2n)! = 2|A|n^p cos(n*theta+phi), theta=arg(z*). Extract exponent p (envelope)
and amplitude A. p=singularity type (N(D) was sqrt => p=-1/2). Tests algebraic vs nested."""
import sys,time; sys.path.insert(0,'debug')
import mpmath as mp
from _fastgl import fast_gl
import _watson_fibre as W
mp.mp.dps=50
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
cs=mp.mpf('0.2'); b=mp.mpf('0.5'); N=40
A=mp.sqrt(cs); d=dcoeffs(cs,A,12); Fn=W.F_taylor(N)
a=[Fn[n]*m_hi(cs,b,n,d,A) for n in range(N+1)]
zst=-(A-1j*b)**2; absz=abs(zst); th=mp.arg(zst)
print(f"z*={mp.nstr(zst,10)} |z*|={mp.nstr(absz,10)} theta=arg={mp.nstr(th,8)} rad",flush=True)
# s_n = a_n |z*|^n/(2n)!  ; find p from |s_n| envelope: |s_n| ~ 2|A| n^p |cos(...)|
print("\nn   s_n=a_n|z*|^n/(2n)!     s_n/n^p for p=-1/2, 1/2, 3/2:",flush=True)
for n in range(12,N+1):
    sn=a[n]*absz**n/mp.factorial(2*n)
    row=" ".join(mp.nstr(sn/ (mp.mpf(n)**p),9) for p in (mp.mpf(-1)/2,mp.mpf(1)/2,mp.mpf(3)/2))
    print(f"{n:2d}  {mp.nstr(sn,12):>16}   {row}",flush=True)
# oscillation check: consecutive n differ by phase theta => s_n ~ cos(n th+phi). Print n*theta mod 2pi
print("\nphase n*theta mod 2pi (compare to sign flips of s_n):",flush=True)
for n in range(30,N+1):
    print(f"  n={n}: n*theta mod 2pi = {mp.nstr(mp.fmod(n*th,2*mp.pi),6)}  sign(s_n)={mp.sign(a[n])}",flush=True)
