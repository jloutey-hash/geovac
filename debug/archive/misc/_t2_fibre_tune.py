"""Find an AFFORDABLE accurate fibre at oscillatory-corner depths the outer grid samples.
Reference = quadosc (proven 40 dig). Compare Jdecay(Nk) and Jtail(K,Nq,M) at rho=0.005,0.002,0.001
near (1,1) (b~2, worst) and (0,1) (b~1)."""
import sys, time
import mpmath as mp
sys.path.insert(0, 'debug')
import routeC_T2_highprec as H
from _t2_decayfibre import Jdecay
from _t2_tailfibre import Jtail, bounded, tail, A_of
from _t2_fibre_stress import J_quadosc

def Jtail_tuned(s, t, K, Nq, M):
    b = s+t
    return bounded(s, t, b, K, Nq) + tail(s, t, b, K, M)

def probe(s, t, label):
    cmin = float(min(s*(1-s), t*(1-t))); b = float(s+t)
    print(f"\n-- {label}: s={mp.nstr(s,6)} t={mp.nstr(t,6)} cmin={cmin:.2e} b={b:.4f} --", flush=True)
    t0 = time.time(); ref = J_quadosc(s, t); tq = time.time()-t0
    print(f"  quadosc REF = {mp.nstr(ref,38)}  ({tq:.1f}s)", flush=True)
    for Nk in (900, 1500, 2500):
        t0 = time.time(); v = Jdecay(s, t, Nk); dt = time.time()-t0
        ag = float(-mp.log10(abs(v-ref)/abs(ref))) if v != ref else 99
        print(f"  Jdecay Nk={Nk:<5} agree={ag:5.1f}dig ({dt:.1f}s)", flush=True)
    Ksc = max(mp.mpf(95), mp.mpf('3.2')/mp.sqrt(cmin))
    for (K, Nqf, M) in [(Ksc, None, 8), (Ksc, None, 12), (1.6*Ksc, None, 12)]:
        Kb = float(K*b); Nq = Nqf or max(400, int(2.4*Kb)+300)
        t0 = time.time(); v = Jtail_tuned(s, t, K, Nq, M); dt = time.time()-t0
        ag = float(-mp.log10(abs(v-ref)/abs(ref))) if v != ref else 99
        print(f"  Jtail K={float(K):.0f} Nq={Nq} M={M}: agree={ag:5.1f}dig ({dt:.1f}s)", flush=True)

if __name__ == '__main__':
    mp.mp.dps = 50
    probe(mp.mpf('0.9975'), mp.mpf('0.9975'), '(1,1) rho=0.005')
    probe(mp.mpf('0.999'), mp.mpf('0.999'), '(1,1) rho=0.002')
    probe(mp.mpf('0.0025'), mp.mpf('0.9975'), '(0,1) rho=0.005')
