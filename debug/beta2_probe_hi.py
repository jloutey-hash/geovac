import sys, time; sys.path.insert(0,'debug')
import mpmath as mp
import beta2_t2_kw_core as C
mp.mp.dps = 90
for k in (5, 80, 200, 320):
    kk=mp.mpf(k); vals=[]
    for Ns in (200+int(0.5*k), 260+int(0.8*k), 340+int(1.1*k)):
        vals.append(C.R_of_w(*C.R_setup(kk, Ns, 6), mp.mpf('0.7')))
    print(f"k={k:4d} R deltas: {[mp.nstr(abs(vals[i+1]-vals[i]),3) for i in range(2)]}  R={mp.nstr(vals[-1],8)}", flush=True)
for k in (80, 200, 320):
    kk=mp.mpf(k); vals=[]; 
    for f in (0.75, 1.0, 1.3):
        vals.append(C.F_of_k(kk, 260+int(0.8*k), int(f*k)+70, 6))
    print(f"k={k:4d} F deltas: {[mp.nstr(abs(vals[i+1]-vals[i]),3) for i in range(2)]}  F={mp.nstr(vals[-1],8)}", flush=True)
