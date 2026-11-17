"""WX phase 0 (before registering; no scored case): measure the grid-level noise of the channel sum F_l on the OLD,
already-voided WL cases A, C, D, so the fresh record's control bound is measured, not guessed.
Level 3 vs level 4 (the level the fresh record will score on, against the next finer one), l = 16..40.
Also A from level 3 and level 4, to see what the noise does to the extracted size.
Usage: python debug/wx_measure_L.py <case>   -> debug/data/wx_measure_L_<case>.json"""
import json
import sys
import time
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent / "tests"))
sys.path.insert(0, str(HERE))
from test_paper59_water_bond_angle import F_coeffs   # noqa: E402
from wl_local import phi                             # noqa: E402
from wl_extrap import leg                            # noqa: E402

mp.mp.dps = 40
LS = list(range(16, 41, 2))
OLD = {'A': (4, 1.24, 1.809, 1.4), 'B': (4, 1.24, 1.809, 2.2), 'C': (4, 1.24, 1.809, 2.6),
       'D': (2.25, 1.24, 1.809, 2.6)}


def basis(xs, nb):
    xs = mp.mpf(xs)
    return {l: [leg(lambda x, a=mp.mpf(7) / 2 + j: (xs - x) ** a, l) for j in range(nb)] for l in LS}


def extract(F, B, nb):
    A = mp.matrix([[B[l][j] / B[l][0] for j in range(nb)] for l in LS])
    y = mp.matrix([mp.mpf(F[l]) / B[l][0] for l in LS])
    return mp.lu_solve(A.T * A, A.T * y)[0]


name = sys.argv[1]
zX, zH, D1, D2 = OLD[name]
xs = (D1 ** 2 + D2 ** 2) / (2 * D1 * D2)
out = {}
F = {}
for lev in (2, 3, 4):
    t = time.time()
    F[lev] = F_coeffs(D1, D2, zX, zH, 40, level=lev)
    out[f"time_l{lev}"] = time.time() - t
B5 = basis(xs, 5)
sc = phi(zX, zH, D1, D2) * mp.mpf(2 * D1 * D2) ** mp.mpf(3.5)
for lev in (2, 3, 4):
    out[f"A4_l{lev}"] = float(extract(F[lev], {l: v[:4] for l, v in B5.items()}, 4) / sc)
    out[f"A5_l{lev}"] = float(extract(F[lev], B5, 5) / sc)
out["noise23"] = max(abs(F[3][l] / F[2][l] - 1) for l in LS)
out["noise34"] = max(abs(F[4][l] / F[3][l] - 1) for l in LS)
out["noise34_by_l"] = {l: abs(F[4][l] / F[3][l] - 1) for l in LS}
print(name, json.dumps(out, indent=1), flush=True)
(HERE / "data" / f"wx_measure_L_{name}.json").write_text(json.dumps(out, indent=1))
