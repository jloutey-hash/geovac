"""INDEPENDENT cross-check evaluator for T2. Different OUTER partition from
routeC_T2_corner_subtraction (rectangular split): here the memo's triangular scheme
  - (0,0) corner region {s+t<delta}: Duffy rho=sigma^2 + angular sin^2
  - trap {s+t>delta}: split at s=delta, sin^2 edge maps
but with the VALIDATED Jsmart fibre (auto-routes quadosc at oscillatory corners).
Shares only the fibre (itself cross-validated to 45 dig by 3 independent methods);
the outer quadrature is structurally different -> a real cross-check of the outer.
"""
import sys, time
import mpmath as mp
sys.path.insert(0, 'debug')
from _fastgl import fast_gl
from routeC_T2_corner_subtraction import Jsmart

def corner(delta, Nsig, Nal):
    H = mp.pi/2
    xs, ws = fast_gl(Nsig); xa, wa = fast_gl(Nal)
    smax = mp.sqrt(delta); tot = mp.mpf(0)
    for xi, wi in zip(xs, ws):
        sig = smax*(xi+1)/2; wsig = smax*wi/2
        for xj, wj in zip(xa, wa):
            psi = H*(xj+1)/2; wpsi = H*wj/2
            al = mp.sin(psi)**2; jal = mp.sin(2*psi)
            s = sig*sig*al; t = sig*sig*(1-al)
            tot += wsig*wpsi*jal*2*sig**3*Jsmart(s, t)
    return tot

def trap(delta, Nout, Nin):
    H = mp.pi/2
    xs, ws = fast_gl(Nout); xa, wa = fast_gl(Nin)
    tot = mp.mpf(0)
    # piece 1: s in [0,delta], t in [delta-s,1]
    for xi, wi in zip(xs, ws):
        phis = H*(xi+1)/2; s = delta*mp.sin(phis)**2; js = delta*mp.sin(2*phis); wphis = H*wi/2
        lo = delta - s
        for xj, wj in zip(xa, wa):
            phit = H*(xj+1)/2; t = lo+(1-lo)*mp.sin(phit)**2; jt = (1-lo)*mp.sin(2*phit); wphit = H*wj/2
            tot += wphis*js*wphit*jt*Jsmart(s, t)
    # piece 2: s in [delta,1], t in [0,1]
    for xi, wi in zip(xs, ws):
        phis = H*(xi+1)/2; s = delta+(1-delta)*mp.sin(phis)**2; js = (1-delta)*mp.sin(2*phis); wphis = H*wi/2
        for xj, wj in zip(xa, wa):
            phit = H*(xj+1)/2; t = mp.sin(phit)**2; jt = mp.sin(2*phit); wphit = H*wj/2
            tot += wphis*js*wphit*jt*Jsmart(s, t)
    return tot

def T2_cross(delta, Nc, Nt):
    ci = corner(delta, Nc, Nc); tr = trap(delta, Nt, Nt)
    return (8/mp.pi)*(ci + tr), ci, tr

if __name__ == '__main__':
    mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    delta = mp.mpf(sys.argv[2]) if len(sys.argv) > 2 else mp.mpf('0.08')
    Nc = int(sys.argv[3]) if len(sys.argv) > 3 else 36
    Nt = int(sys.argv[4]) if len(sys.argv) > 4 else 40
    ref = mp.mpf('0.3953557659017139641')
    print(f"CROSS dps={mp.mp.dps} delta={delta} Nc={Nc} Nt={Nt}", flush=True)
    t0 = time.time(); v, ci, tr = T2_cross(delta, Nc, Nt)
    print(f"  corner={mp.nstr(ci,26)} trap={mp.nstr(tr,26)}", flush=True)
    print(f"  T2 = {mp.nstr(v, mp.mp.dps-4)}  ({time.time()-t0:.0f}s)", flush=True)
    print(f"  |T2 - anchor19| = {mp.nstr(abs(v-ref),4)}", flush=True)
