"""Stress-test fibre evaluators at oscillatory (b->2) corner points, and develop an
oscillation-free fibre by CONTOUR ROTATION.

J = int_0^inf j0(bk) E(k) dk,  E(k)=P(s,k)P(t,k),  b=s+t.
Contour idea: j0(bk)=sin(bk)/(bk); with E even & analytic in |Im k|<1/sqrt(cmin),
   int_0^inf sin(bk)/(bk) E(k) dk .
Rotate the e^{ibk} half to k=iy (decays e^{-by}); the branch point of sqrt(c k^2+1)
at k=i/sqrt(c) forces a cut contribution for y>1/sqrt(cmin).
Here we just BENCHMARK existing evaluators + mpmath.quadosc to find a 32-digit reference.
"""
import sys, time
import mpmath as mp
sys.path.insert(0, 'debug')
import routeC_T2_highprec as H
from _t2_decayfibre import Jdecay
from _t2_tailfibre import Jtail

def J_quadosc(s, t):
    """mpmath oscillatory quadrature: int_0^inf sin(bk)/(bk) E(k) dk. Period 2pi/b."""
    cs = s*(1-s); ct = t*(1-t); b = s+t
    def integ(k):
        kb = k*b
        j0 = mp.sin(kb)/kb if kb > mp.mpf('1e-30') else mp.mpf(1)
        return j0*H.P(s, k)*H.P(t, k)
    # zeros of sin(bk) at k=n*pi/b -> quadosc period
    return mp.quadosc(integ, [0, mp.inf], period=2*mp.pi/b)

if __name__ == '__main__':
    mp.mp.dps = 55
    pts = [(mp.mpf('0.99'), mp.mpf('0.99')),      # (1,1) rho=0.02
           (mp.mpf('0.995'), mp.mpf('0.995')),    # rho=0.01
           (mp.mpf('0.9975'), mp.mpf('0.9975')),  # rho=0.005
           (mp.mpf('0.999'), mp.mpf('0.999'))]    # rho=0.002
    for (s, t) in pts:
        cmin = float(min(s*(1-s), t*(1-t)))
        print(f"\n-- s=t={mp.nstr(s,6)}  cmin={cmin:.2e}  b={float(s+t):.4f} --", flush=True)
        t0 = time.time(); jq = J_quadosc(s, t); print(f"  quadosc      = {mp.nstr(jq,40)}  ({time.time()-t0:.1f}s)", flush=True)
        for Nk in (800, 1400, 2200, 3200):
            t0 = time.time(); jd = Jdecay(s, t, Nk)
            print(f"  Jdecay Nk={Nk:<4} = {mp.nstr(jd,40)}  d(quadosc)={mp.nstr(abs(jd-jq),3)}  ({time.time()-t0:.1f}s)", flush=True)
        t0 = time.time(); jt = Jtail(s, t); print(f"  Jtail        = {mp.nstr(jt,40)}  d(quadosc)={mp.nstr(abs(jt-jq),3)}  ({time.time()-t0:.1f}s)", flush=True)
