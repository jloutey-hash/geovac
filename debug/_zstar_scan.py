"""Closed-form route: locate the dominant Borel singularity z*(c_s,b) of the fibre c_t-series over a
(c_s,b) grid, and hunt for its analytic form (relate to curve branch points: candidates c_s, b^2,
c_s+b^2, (sqrt(c_s)+-i b)^2, etc.)."""
import sys,time; sys.path.insert(0,'debug')
import mpmath as mp
import _watson_fibre as W
mp.mp.dps=50
def zstar(cs,b,N=18):
    cs=mp.mpf(cs); b=mp.mpf(b)
    Fn=W.F_taylor(N); ms=[W.m_n(cs,b,n) for n in range(N+1)]
    a=[Fn[n]*ms[n] for n in range(N+1)]
    bc=[a[n]/mp.factorial(2*n) for n in range(N+1)]
    p,q=mp.pade(bc,9,9)
    roots=sorted(mp.polyroots(q[::-1],maxsteps=300,extraprec=200),key=lambda z:abs(z))
    return roots[0]
print("c_s    b     z*=Re+Imj                |z*|        candidates",flush=True)
for cs,b in [('0.2','0.5'),('0.2','1.0'),('0.2','0.3'),('0.1','0.5'),('0.15','0.5'),('0.3','0.5')]:
    t0=time.time(); z=zstar(cs,b)
    cf=mp.mpf(cs); bf=mp.mpf(b)
    cands=f"c_s={mp.nstr(cf,4)} b^2={mp.nstr(bf*bf,4)} c_s+b^2={mp.nstr(cf+bf*bf,4)} (sqrtcs+ib)^2={mp.nstr((mp.sqrt(cf)+1j*bf)**2,6)}"
    print(f"  {cs:>4} {b:>4}  {mp.nstr(z,10):>26}  {mp.nstr(abs(z),8):>10}  {cands}",flush=True)
