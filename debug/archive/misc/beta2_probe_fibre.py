import sys, time; sys.path.insert(0,'debug')
import mpmath as mp
import beta2_t2_st_route as ST
mp.mp.dps = 50
print("fibre J(s,t): analytic-tail vs brute quadosc, and self-convergence in (M,Kf,Nq)")
for (s,t) in (('0.3','0.4'), ('0.02','0.5'), ('0.001','0.9'), ('0.45','0.48')):
    ss, tt = mp.mpf(s), mp.mpf(t)
    a = ST.J_acc(ss, tt, 8)
    b = ST.J_acc(ss, tt, 14, mp.mpf('5.0'), mp.mpf('3.2'))
    bq = ss+tt
    f = lambda k: (mp.sin(k*bq)/(k*bq))*ST.P(ss,k)*ST.P(tt,k)
    t0=time.time(); q = mp.quadosc(f,[0,mp.inf],period=2*mp.pi/bq); dt=time.time()-t0
    print(f" s={s} t={t}: J={mp.nstr(b,30)}  selfconv {mp.nstr(abs(b-a)/abs(b),3)}"
          f"  vs quadosc {mp.nstr(abs(b-q)/abs(b),3)} ({dt:.0f}s)", flush=True)
