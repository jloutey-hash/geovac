"""HG3 addendum: refine the k rule of debug/hg3_master.py at fixed t level 5, at LM's two test points
(u2 = w2 = 1, u3 = w3 = 2, r = 2; w1 = 1 is LM Table I, w1 = 2 is LM Table III).
    python debug/hg3_krefine.py
"""
import hg3_master as H

LM1 = 2.39121170158262983284e-4
for w1 in (1.0, 2.0):
    for lk in (6, 7, 8, 9):
        import math
        import numpy as np
        ht, hk = 2.0 ** -5, 2.0 ** -lk
        t, tc, wt = H.ts(int(4.5 / ht), ht)
        u, uc, wu = H.ts(int(4.5 / hk), hk)
        k = 4.0 * u / uc
        wk = wu * 4.0 / uc ** 2
        T, S = t[:, None], t[None, :]
        tot = []
        for kk, ww in zip(k, wk):
            D1 = t * tc * kk ** 2 + t * 4 + tc * 1
            D2 = t * tc * kk ** 2 + t * 1 + tc * 4
            E1 = wt * np.exp(-2 * np.sqrt(D1)) / np.sqrt(D1)
            E2 = wt * np.exp(-2 * np.sqrt(D2)) / np.sqrt(D2)
            tot.append(ww * kk ** 2 / (kk ** 2 + w1 ** 2) * (E1 @ H.j0(kk * 2 * (T - S)) @ E2))
        v = 2 / (2 * math.pi) * math.fsum(tot)
        rel = (v - LM1) / LM1 if w1 == 1.0 else None
        print('w1', w1, 'k level', lk, repr(v), 'rel to LM Table I' if rel is not None else '', rel, flush=True)
