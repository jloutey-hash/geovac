"""Record FC control K4ext (added after the cold review): the engine's large-L tail against record AB's exact route.

K4A1 (b = a, Z = W) is record AB's case A1. Record AB stored beta_eff of its exact coefficients on F / rho^L,
rho = 9/13, at L = 20, 40, ..., 220 (debug/data/aabc_p1p3b_2600.json, key A1.beta_f; AB is scored and public). This
script computes the same beta_eff from fc_engine at levels A and B to L = 222 and writes the differences to
debug/data/fc_k4ext.json. An L-dependent normalisation or recursion error in the engine shifts beta_eff by O(1).

    python debug/fc_k4ext.py
"""
import json
import sys
import time
from pathlib import Path

import mpmath as mp

sys.path.insert(0, str(Path(__file__).parent))
import fc_engine as FE  # noqa: E402

mp.mp.dps = 120
AB = json.load(open(Path(__file__).parent / 'data' / 'aabc_p1p3b_2600.json'))['A1']['beta_f']
lan = FE.landau('K4A1')
rho = lan['rhoZ']
assert rho == lan['rhoW']
rm = mp.mpf(rho.numerator) / rho.denominator
out = {'rho': str(rho), 'AB_rho': '9/13'}
for lev in ('A', 'B'):
    t0 = time.time()
    F, _ = FE.coefficients('K4A1', 222, lev)
    G = [mp.mpf(f.mid().str(110, radius=False)) / rm ** L for L, f in enumerate(F)]
    rows = {}
    for Ls, b_ab in AB.items():
        L = int(Ls)
        b = mp.log(abs(G[L + 1] / G[L - 1])) / mp.log(mp.mpf(L + 1) / (L - 1))
        rows[Ls] = {'fc': float(b), 'AB': b_ab, 'diff': float(b - b_ab)}
    out[lev] = {'rows': rows, 'max_abs_diff': max(abs(r['diff']) for r in rows.values()),
                'all_positive': all(f > 0 for f in F), 'time_s': round(time.time() - t0, 1)}
    print(lev, out[lev]['max_abs_diff'], out[lev]['all_positive'], out[lev]['time_s'], flush=True)
(Path(__file__).parent / 'data' / 'fc_k4ext.json').write_text(json.dumps(out, indent=1))
