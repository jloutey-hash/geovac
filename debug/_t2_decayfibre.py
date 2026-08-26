"""Fast decay-map fibre (same math as routeC_T2_highprec.J, but fast_gl nodes).
J(s,t)=int_0^inf j0(kb)P(s,k)P(t,k)dk via k=L u/(1-u), L=1/(sqrt cs+sqrt ct).
Independent family from the tail-analytic fibre (Jtail): shares NO tail machinery."""
import mpmath as mp
import sys
sys.path.insert(0, 'debug')
from _fastgl import fast_gl
import routeC_T2_highprec as H   # for P()

def Jdecay(s, t, Nk=600):
    cs = s*(1-s); ct = t*(1-t)
    if cs == 0 or ct == 0:
        return mp.mpf(0)
    L = 1/(mp.sqrt(cs)+mp.sqrt(ct)); b = s+t
    xs, ws = fast_gl(Nk)
    tot = mp.mpf(0)
    for x, w in zip(xs, ws):
        u = (x+1)/2
        k = L*u/(1-u); dk = L/(1-u)**2
        kb = k*b
        j0 = mp.sin(kb)/kb if kb > mp.mpf('1e-40') else mp.mpf(1)
        tot += (w/2)*j0*H.P(s, k)*H.P(t, k)*dk
    return tot

if __name__ == '__main__':
    mp.mp.dps = 50
    import time
    from _t2_tailfibre import Jtail
    for (s, t) in [(mp.mpf('0.3'), mp.mpf('0.17')), (mp.mpf('0.01'), mp.mpf('0.01')),
                   (mp.mpf('0.005'), mp.mpf('0.02'))]:
        t0 = time.time(); jd = Jdecay(s, t, 600); td = time.time()-t0
        t0 = time.time(); jt = Jtail(s, t); tt = time.time()-t0
        jd2 = Jdecay(s, t, 1000)
        print(f"s={mp.nstr(s,3)} t={mp.nstr(t,3)}: Jdecay600={mp.nstr(jd,20)} ({td:.2f}s) "
              f"Jtail={mp.nstr(jt,20)} ({tt:.2f}s)  |tail-decay|={mp.nstr(abs(jt-jd2),3)} "
              f"|Nk600-1000|={mp.nstr(abs(jd-jd2),3)}", flush=True)
