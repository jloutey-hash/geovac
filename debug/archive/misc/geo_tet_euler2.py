"""Euler characteristics of two sector cuts at a split value of x (geometry only).

chi is independent of x for generic x, so take x = -16/9 (sqrt(-x) rational) and primes
p = 1 mod 4 (i in F_p): every point stratum below is then defined over F_p and the
point count is a polynomial in p; chi = its value at p = 1.

  sunrise (no integer line): U = (C*)^3 in (A,B,T) minus {lambda(d_a,d_b,d_t) = 0}.
      Predicted chi = -4 (the cone through (1,4,x) minus three curves meeting in 6 points),
      i.e. 4 masters = the rank of the K3 piece.
  sector {g1,C} (Coulomb + geminal 1, collapsed to a Yukawa in D = 3): U = (C*)^2 in (B,T)
      minus the conic lambda(-g1^2, d_b, d_t) = 0.  Predicted chi = -3: E_B (rank 2) + 1.
"""
import json
from fractions import Fraction as Fr
import numpy as np

X = Fr(-16, 9)


def primes(lo, hi):
    return [n for n in range(lo, hi) if n > 1 and all(n % d for d in range(2, int(n**0.5) + 1))]


def red(q, p):
    q = Fr(q)
    return (q.numerator * pow(q.denominator, -1, p)) % p


def lam(u, v, w, p):
    return (u * u + v * v + w * w - 2 * u * v - 2 * v * w - 2 * w * u) % p


out = {}
pts_s, pts_e = [], []
for p in primes(13, 400):
    if p % 4 != 1 or p in (2, 3, 5, 7):
        continue
    r = np.arange(p, dtype=np.int64)
    xa = red(X, p)
    # sunrise
    tot = 0
    for Av in range(1, p):
        da = (Av - 1) % p
        B = r[:, None]
        T = r[None, :]
        L = lam(da, (B - 4) % p, (T - xa) % p, p)
        tot += int(((B != 0) & (T != 0) & (L != 0)).sum())
    pts_s.append((p, tot))
    # sector {g1, C}
    g1 = red(Fr(9, 4), p)
    B = r[:, None]
    T = r[None, :]
    L = lam((-g1) % p, (B - 4) % p, (T - xa) % p, p)
    pts_e.append((p, int(((B != 0) & (T != 0) & (L != 0)).sum())))


def fit(pts, deg):
    ps = np.array([q for q, _ in pts], dtype=float)
    Ns = np.array([n for _, n in pts], dtype=float)
    ic = [int(round(c)) for c in np.polyfit(ps, Ns, deg)]
    ok = sum(int(np.polyval(ic, q)) == n for q, n in pts)
    return {'coeffs': ic, 'exact_at': f'{ok}/{len(pts)}', 'chi': int(np.polyval(ic, 1))}


out['sunrise'] = fit(pts_s, 3)
out['sector_g1C'] = fit(pts_e, 2)
out['counts_sunrise'] = pts_s
out['counts_g1C'] = pts_e
with open('debug/data/geo_tet_euler2.json', 'w') as fh:
    json.dump(out, fh, indent=1)
print(out['sunrise'], out['sector_g1C'])
