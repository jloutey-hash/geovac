import sys, time; sys.path.insert(0,'debug')
import mpmath as mp
import beta2_t2_kw_core as C
import beta2_t2_tail as T
mp.mp.dps = 70
acc = T.F_buckets(56)
# sanity on gammainc: compare _E(P,a,K) against direct quadosc for a few P
for P in (10, 40, 90):
    for a in (1, 2):
        E = T._E(P, a, mp.mpf(200))
        f = lambda k: k**(-P)*mp.cos(a*k)
        q = mp.quadosc(f, [200, mp.inf], period=2*mp.pi/a)
        print(f"  P={P} a={a}: Re E={mp.nstr(E.real,20)} quadosc={mp.nstr(q,20)} rel {mp.nstr(abs(E.real-q)/abs(q),3)}")
print(flush=True)
t0=time.time()
num = C.int_F(200, 260, 20, 24, lambda k: 340, 6, nw=lambda k: int(0.75*float(k))+70)
print(f"numeric int_200^260 F dk = {mp.nstr(num, 40)}   ({time.time()-t0:.0f}s)")
ana = T.TAIL(200, acc) - T.TAIL(260, acc)
print(f"analytic TAIL(200)-TAIL(260) = {mp.nstr(ana, 40)}")
print(f"rel diff = {mp.nstr(abs(num-ana)/abs(num), 4)}")
print(f"TAIL(200) = {mp.nstr(T.TAIL(200,acc),30)}   TAIL(150) = {mp.nstr(T.TAIL(150,acc),30)}")
