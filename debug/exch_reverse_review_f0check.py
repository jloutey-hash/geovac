"""exch_reverse_review_f0check.py -- REVIEW: is the t = 0.99 closed-form/hp gap an f0-quadrature artifact?
f0(R) at t = 99/100 by mp.quad with error estimate, at increasing maxdegree, and with a finer subdivision."""
import sys
import mpmath as mp
sys.path.insert(0, "debug")
import exch_reverse_review_closed as C
mp.mp.dps = 60
t, R = mp.mpf(99) / 100, mp.mpf(3) / 2
f = lambda s: C.funcs(t, s, "real")[1] / s
for md in (6, 8, 10):
    v, e = mp.quad(f, [R, R + 1, R + 2, R + 4, R + 10, R + 30, mp.inf], error=True, maxdegree=md)
    print("maxdegree", md, mp.nstr(v, 45), "err", mp.nstr(e, 3), flush=True)
pts = [R + mp.mpf(k) / 4 for k in range(0, 41)] + [R + 20, R + 40, R + 80, mp.inf]
v2, e2 = mp.quad(f, pts, error=True, maxdegree=10)
print("fine", mp.nstr(v2, 45), "err", mp.nstr(e2, 3))
# integrand values near s = R: is there cancellation?
g, fm1 = C.funcs(t, R, "real")
print("FM1(R)", mp.nstr(fm1, 20), "psi_a", mp.nstr(g[7], 20), "psi_b", mp.nstr(g[8], 20))
