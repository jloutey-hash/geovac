import sys, time; sys.path.insert(0,'debug')
import mpmath as mp
import beta2_t2_st_route as ST
mp.mp.dps = 50
print("J(s,t): analytic-tail/cut fibre vs brute quadosc")
for (s,t) in (('0.3','0.4'), ('0.02','0.5'), ('0.001','0.9'), ('0.45','0.48'),
              ('0.00001','0.5'), ('0.0003','0.0004')):
    ss, tt = mp.mpf(s), mp.mpf(t)
    t0=time.time(); a = ST.J_acc(ss, tt); ta=time.time()-t0
    b2 = ST.J_acc(ss, tt, 30, mp.mpf(35), mp.mpf(160))
    bq = ss+tt
    f = lambda k: (mp.sin(k*bq)/(k*bq))*ST.P(ss,k)*ST.P(tt,k)
    q = mp.quadosc(f,[0,mp.inf],period=2*mp.pi/bq)
    print(f" s={s:>9} t={t:>7}: J={mp.nstr(a,26)}  selfconv {mp.nstr(abs(b2-a)/abs(a),3)}"
          f"  vs quadosc {mp.nstr(abs(a-q)/abs(a),3)}  ({ta*1000:.0f}ms)", flush=True)
