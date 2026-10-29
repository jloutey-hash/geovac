"""Record B2X: large-L behaviour of the true Legendre coefficients F_L (exact rationals times the constants at high
precision), for the two-electron case W1 or the one-electron control N1.  Reports r_L = F_{L+1}/F_L and the local
exponent beta_eff(L) = ln(G_{L+1}/G_L)/ln(1 + 1/L), G_L = F_L/lam^L, with lam = rho; and the digits that agree between
two working precisions.        python debug/b2_large_l.py <W1|N1> <LMAX>  ->  data/b2x_<case>.json
"""
import json
import sys
import time

import mpmath as mp

import ao_exact
import na_exact

CASES = {'W1': (ao_exact, ('1.7', '2.5', '1.4', '1.1'), (17, 25)),
         'N1': (na_exact, ('1.8', '2.6', '1.3', '1.0'), (9, 13))}
case, LMAX = sys.argv[1], int(sys.argv[2])
X, args, (a, b) = CASES[case]
t0 = time.time()
Q = X.coefficients(*args, LMAX, na_exact.Field(None))
print(case, 'exact coefficients', round(time.time() - t0), 's', flush=True)


def values(dps):
    mp.mp.dps = dps
    cv = {}
    out = []
    for q in Q:
        s = mp.mpf(0)
        for k, v in q.items():
            if k not in cv:
                cv[k] = X.value(k, mp)
            s += mp.mpf(int(v.p)) / int(v.q) * cv[k]
        out.append(s)
    return out


D = 6 * LMAX + 200
lo, hi = values(D), values(D + D // 2)
agree = [float(-mp.log10(abs(x / y - 1))) if x != y else float(D) for x, y in zip(lo, hi)]
mp.mp.dps = D + D // 2
lam = mp.mpf(a) / b
rows = []
for L in range(1, len(hi) - 1):
    r = hi[L + 1] / hi[L]
    be = mp.log(r / lam) / mp.log(1 + mp.mpf(1) / L)
    rows.append(dict(L=L, F=mp.nstr(hi[L], 12), r=mp.nstr(r, 12), beta_eff=mp.nstr(be, 8), digits=round(agree[L], 1)))
for row in rows:
    if row['L'] % 25 == 0 or row['L'] < 5:
        print(row, flush=True)
json.dump(dict(case=case, lam=f'{a}/{b}', rows=rows), open(f'data/b2x_{case}.json', 'w'), indent=0)
print('done', round(time.time() - t0), 's')
