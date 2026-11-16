import sys, time; sys.path.insert(0,"debug"); sys.path.insert(0,"."); sys.path.insert(0,"tests")
from exch_scope_series import *
p = 1000003
f1 = Q.exchange_closed_form(Fraction(1), S1, S1, Fraction(1), S1, S1, R, tau_max=6)
f2 = Q.exchange_closed_form(Fraction(2), S1, S1, Fraction(2), S1, S1, R, tau_max=6)
t=time.time(); ser = series_of(f1 + f2, 120); print("series", time.time()-t)
for r in range(4, 10):
    for d in range(4, 40, 2):
        ns = guess_mod_p(ser, r, d, p, mmax=110)
        if ns:
            print(f"LCLM control: first hit r={r} d={d} nullity {len(ns)}"); break
    else:
        continue
    break
# planted non-holonomic: tan(R) (+ R^9 lnR tan R), exact Taylor
tser = sp.series(sp.tan(R), R, 0, 122).removeO()
T = {k: Fraction(int(c.p), int(c.q)) for (k,), c in sp.Poly(tser, R).terms()}
ser_t = {(sp.Integer(1), 0): T, (sp.Integer(1), 1): {k+9: v for k, v in T.items()}}
hits = [(r, d) for r in range(1, 8) for d in range(0, 14) if guess_mod_p(ser_t, r, d, p, mmax=110)]
print("tan control hits within (r<=7, d<=13), mmax 110:", hits[:5])
