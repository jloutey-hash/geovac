"""HG3 addendum: the k integral of debug/hg3_master.py with an analytic tail (the T2 tail-analytic idea).
F(k) = k^2/(k^2 + w1^2) Int Int dt ds j0(k R (t - s)) E1(t, k) E2(s, k).  For large k only the bond ends contribute,
so F(k) ~ sum_n a_n k^-n (even n >= 4) + sin(R k) sum_n b_n k^-n (odd n >= 5) + cos(R k) sum_n c_n k^-n (even n >= 6).
Int_0^K F dk by Gauss-Legendre panels; the a, b, c are solved for (linear least squares) from F on [K, 4K], and the
tail Int_K^inf is done exactly: K^(1-n)/(n-1) and K^(1-n) E_n(-i R K).
    python debug/hg3_tail.py
"""
import math

import mpmath as mp
import numpy as np

import hg3_master as H

U2, U3, W2, W3, R = 1.0, 2.0, 1.0, 2.0, 2.0


def F(kk, w1, t, tc, wt):
    T, S = t[:, None], t[None, :]
    D1 = t * tc * kk ** 2 + t * U3 ** 2 + tc * U2 ** 2
    D2 = t * tc * kk ** 2 + t * W2 ** 2 + tc * W3 ** 2
    E1 = wt * np.exp(-R * np.sqrt(D1)) / np.sqrt(D1)
    E2 = wt * np.exp(-R * np.sqrt(D2)) / np.sqrt(D2)
    return kk ** 2 / (kk ** 2 + w1 ** 2) * (E1 @ H.j0(kk * R * (T - S)) @ E2)


def basis(k, nord):
    cols, kinds = [], []
    for j in range(nord):
        cols.append(k ** -(4 + 2 * j)); kinds.append(('a', 4 + 2 * j))
        cols.append(np.sin(R * k) * k ** -(5 + 2 * j)); kinds.append(('s', 5 + 2 * j))
        cols.append(np.cos(R * k) * k ** -(6 + 2 * j)); kinds.append(('c', 6 + 2 * j))
    return np.array(cols).T, kinds


def tail(K, coef, kinds):
    out = mp.mpf(0)
    for c, (kind, n) in zip(coef, kinds):
        if kind == 'a':
            out += c * mp.mpf(K) ** (1 - n) / (n - 1)
        else:
            z = mp.mpf(K) ** (1 - n) * mp.expint(n, -1j * R * K)   # Int_K^inf e^{i R k} k^-n dk
            out += c * (z.imag if kind == 's' else z.real)
    return float(out)


def master_tail(w1, lt, K, nord, npan=None, ngl=24):
    t, tc, wt = H.ts(int(4.5 * 2 ** lt), 2.0 ** -lt)
    npan = npan or int(2 * K)
    x, w = np.polynomial.legendre.leggauss(ngl)
    parts = []
    for p in range(npan):
        a, b = K * p / npan, K * (p + 1) / npan
        for xi, wi in zip(x, w):
            parts.append(0.5 * (b - a) * wi * F(0.5 * (b - a) * xi + 0.5 * (a + b), w1, t, tc, wt))
    ks = np.linspace(K, 4 * K, 240)
    Fs = np.array([F(kk, w1, t, tc, wt) for kk in ks])
    A, kinds = basis(ks, nord)
    sc = np.abs(A).max(axis=0)
    coef, *_ = np.linalg.lstsq(A / sc, Fs, rcond=None)
    coef = coef / sc
    resid = np.abs(A @ coef - Fs).max() / np.abs(Fs).max()
    tl = tail(K, coef, kinds)
    return R / (2 * math.pi) * (math.fsum(parts) + tl), resid, R / (2 * math.pi) * tl


if __name__ == '__main__':
    mp.mp.dps = 30
    LM1 = 2.39121170158262983284e-4
    for w1 in (1.0, 2.0):
        for lt, K, nord in ((5, 20, 4), (5, 30, 4), (5, 30, 5), (5, 40, 5), (6, 30, 5), (6, 40, 6)):
            v, res, tl = master_tail(w1, lt, K, nord)
            print('w1', w1, 'tlev', lt, 'K', K, 'orders', nord, repr(v),
                  'rel LM I' if w1 == 1 else '', (v - LM1) / LM1 if w1 == 1 else '',
                  'fit resid', '%.1e' % res, 'tail', '%.2e' % tl, flush=True)
