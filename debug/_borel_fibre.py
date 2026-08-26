"""Target B, step 1: characterize the fibre's c_t-series resurgence (two-scale analog of N(D)).
J(c_s,c_t,b)/c_t = sum_n a_n c_t^n, a_n = F_n * m_n(c_s,b), m_n=int k^{2n} j0(kb) P(c_s,k) dk.
Study large-n growth (|a_{n+1}/a_n| ~ linear => Gevrey-1; slope=1/Borel-radius) and Borel-Pade the
singularity locations.  Compare to N(D): Gevrey-1, Borel radius 2, singularities at branch points."""
import sys,time; sys.path.insert(0,'debug')
import mpmath as mp
import _watson_fibre as W
mp.mp.dps=50
cs=mp.mpf('0.2'); b=mp.mpf('0.5'); N=24
t0=time.time()
Fn=W.F_taylor(N)
ms=[W.m_n(cs,b,n) for n in range(N+1)]
a=[Fn[n]*ms[n] for n in range(N+1)]
print(f"built a_0..a_{N} in {time.time()-t0:.0f}s",flush=True)
print("n   a_n                      |a_n|^(1/n)      |a_{n+1}/a_n|",flush=True)
for n in range(N+1):
    root = mp.nstr(abs(a[n])**(mp.mpf(1)/n),8) if n>0 else '-'
    rat = mp.nstr(abs(a[n]/a[n-1]),8) if n>0 else '-'
    print(f"{n:2d}  {mp.nstr(a[n],14):>22}  {root:>14}  {rat:>12}",flush=True)
# ratio/n -> 1/Borel-radius if Gevrey-1 (a_n ~ n! / R^n)
print("\n|a_{n+1}/a_n| / n  (-> 1/R if a_n ~ n!/R^n Gevrey-1):",flush=True)
for n in range(6,N):
    print(f"  n={n}: {mp.nstr(abs(a[n+1]/a[n])/n,8)}",flush=True)
# also test double-factorial (2n)! growth like N(D): |a_{n+1}/a_n|/(n^2)
print("\n|a_{n+1}/a_n| / n^2  (-> const if a_n ~ (2n)!-type):",flush=True)
for n in range(6,N):
    print(f"  n={n}: {mp.nstr(abs(a[n+1]/a[n])/(n*n),8)}",flush=True)
