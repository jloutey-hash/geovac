import sys; sys.path.insert(0,'debug')
import mpmath as mp
import routeC_T2_eichler_lambert as EL
import _t2_tailfibre as TF
mp.mp.dps = 40
for s,t in [(mp.mpf('0.3'),mp.mpf('0.2')),(mp.mpf('0.03'),mp.mpf('0.45'))]:
    cs=s*(1-s); ct=t*(1-t); b=s+t
    # Jtail self-consistency: vary internal (K, Nq, M)
    j_a = TF.Jtail(s,t)                    # defaults
    j_b = TF.Jtail(s,t, K=mp.mpf(60), Nq=700, M=8)
    j_c = TF.Jtail(s,t, K=mp.mpf(90), Nq=900, M=10)
    j200 = EL.Jsum(cs,ct,[b],200)
    print(f"s={mp.nstr(s,3)} t={mp.nstr(t,3)} b={mp.nstr(b,3)}:",flush=True)
    print(f"   Jtail_default = {mp.nstr(j_a,34)}",flush=True)
    print(f"   selfconv |a-b|={mp.nstr(abs(j_a-j_b),3)}  |b-c|={mp.nstr(abs(j_b-j_c),3)}  vs Jsum200 |a-200|={mp.nstr(abs(j_a-j200),3)}",flush=True)
