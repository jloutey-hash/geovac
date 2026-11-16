import sys; sys.path.insert(0,"debug"); sys.path.insert(0,"."); sys.path.insert(0,"tests")
from exch_scope_series import *
from exch_scope_series import _poly_coeffs
import mpmath as mp
mp.mp.dps=40
f = Q.exchange_closed_form(Fraction(1), S1, S1, Fraction(1), S1, S1, R, tau_max=6)
atoms, reps = decompose(f)
for a,v in atoms.items():
    num,den=sp.fraction(sp.together(v)); print(a, "den val", min(_poly_coeffs(den)), "deg num", sp.Poly(num,R).degree() if not num.free_symbols-{R} else "c")
ser = series_of(f, 70)
for x in ["0.3","1.0"]:
    xv=mp.mpf(x); tot=0
    for (ce,lp),d in ser.items():
        cv=mp.mpf(sp.N(ce.subs(G,sp.EulerGamma),50))
        for k,v in d.items(): tot+=cv*mp.mpf(v.numerator)/v.denominator*xv**k*mp.log(xv)**lp
    print(x, tot, sp.N(f.subs(R,sp.Rational(x)),30))
for key,d in ser.items(): print(key, min(d), max(d))
mp.mp.dps=80
fl = sp.lambdify(R, f, modules="mpmath")
for x in ["0.05","0.3","0.5","1.0","2.0"]:
    xv=mp.mpf(x); tot=0
    for (ce,lp),d in ser.items():
        cv=mp.mpf(sp.N(ce.subs(G,sp.EulerGamma),90))
        for k,v in d.items(): tot+=cv*mp.mpf(v.numerator)/v.denominator*xv**k*mp.log(xv)**lp
    print(x, mp.nstr(tot,25), mp.nstr(fl(xv),25))
