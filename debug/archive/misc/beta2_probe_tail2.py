import sys, time; sys.path.insert(0,'debug')
import mpmath as mp
import beta2_t2_kw_core as C
import beta2_t2_tail as T
mp.mp.dps = 80
env = {}
for MX in (20, 30, 40, 56):
    t0=time.time(); A, B = T.build_AB(MX); tb=time.time()-t0
    print(f"MX={MX}: |A|={len(A)} |B|={len(B)} build {tb:.1f}s")
    for k in (100, 150, 200, 300):
        row=[]
        for w in ('0.3','0.85'):
            ww = mp.mpf(w); kk=mp.mpf(k)
            if (k,w) not in env:
                env[(k,w)] = C.R_of_w(*C.R_setup(kk, 420, 6), ww)
            ex = env[(k,w)]
            ap = T.R_asym(kk, ww, A, B)
            scale = 4*kk**-4*mp.mpf('0.9')   # rough envelope
            row.append(mp.nstr(abs(ap-ex)/scale, 3))
        print(f"   k={k:4d}: abs err / envelope = {row}")
