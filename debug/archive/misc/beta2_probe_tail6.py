import sys, time; sys.path.insert(0,'debug')
import mpmath as mp
import beta2_t2_kw_core as C
import beta2_t2_tail as T
mp.mp.dps = 70
K = mp.mpf(200)
print("--- (a) recurrence check on _E(P,a,K):  E_{P+1} = (i a E_P + K^-P e^{iaK})/P")
for a in (1,2):
    for P in (9, 30, 60, 100):
        E1 = T._E(P, a, K); E2 = T._E(P+1, a, K)
        rhs = (mp.mpc(0,a)*E1 + K**(-P)*mp.e**(mp.mpc(0,a*K)))/P
        print(f"   a={a} P={P}: rel {mp.nstr(abs(E2-rhs)/abs(E2),3)}")
print("--- (b) TAIL analytic vs numeric integral of the SAME asymptotic F")
acc = T.F_buckets(56)
def num_asym(a,b,npan,nn):
    a=mp.mpf(a); b=mp.mpf(b); xs,ws=C.fast_gl(nn); h=(b-a)/npan; tot=mp.mpf(0)
    for ip in range(npan):
        lo=a+ip*h
        for xg,wg in zip(xs,ws):
            tot += (h*wg/2)*T.F_asym(lo+h*(xg+1)/2, acc)
    return tot
for (a,b) in ((200,260),(200,400),(150,200)):
    q = num_asym(a,b,int((b-a)/1.5)+1,22)
    an = T.TAIL(a,acc)-T.TAIL(b,acc)
    print(f"   [{a},{b}]: num {mp.nstr(q,30)}  ana {mp.nstr(an,30)}  rel {mp.nstr(abs(q-an)/abs(q),3)}")
print("--- (c) coarse direct check: exact F integrated vs TAIL difference")
t0=time.time()
num = C.int_F(200, 215, 10, 16, 200, 6, nw=lambda k: int(0.55*float(k))+45)
ana = T.TAIL(200,acc)-T.TAIL(215,acc)
print(f"   num {mp.nstr(num,30)}  ana {mp.nstr(ana,30)}  rel {mp.nstr(abs(num-ana)/abs(num),3)}  ({time.time()-t0:.0f}s)")
print(f"--- TAIL(150)={mp.nstr(T.TAIL(150,acc),30)}  TAIL(200)={mp.nstr(T.TAIL(200,acc),30)}")
