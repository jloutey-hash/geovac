"""WL diagnostic after scoring (not registered): level 2 vs level 3 of the channel sum, and A from level 3."""
import json, sys
from pathlib import Path
import mpmath as mp
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent / "tests")); sys.path.insert(0, str(HERE))
from test_paper59_water_bond_angle import F_coeffs
from wl_local import phi
from wl_extrap import leg
mp.mp.dps = 40
LS = list(range(16, 41, 2))
CASES = {'A': (4, 1.24, 1.809, 1.4), 'C': (4, 1.24, 1.809, 2.6), 'D': (2.25, 1.24, 1.809, 2.6)}


def basis(xs, nb):
    xs = mp.mpf(xs)
    return {l: [leg(lambda x, a=mp.mpf(7) / 2 + j: (xs - x) ** a, l) for j in range(nb)] for l in LS}


def extract(F, B, nb):
    A = mp.matrix([[B[l][j] / B[l][0] for j in range(nb)] for l in LS])
    y = mp.matrix([mp.mpf(F[l]) / B[l][0] for l in LS])
    return mp.lu_solve(A.T * A, A.T * y)[0]
out = {}
for name in ('A', 'C', 'D'):
    zX, zH, D1, D2 = CASES[name]
    xs = (D1 ** 2 + D2 ** 2) / (2 * D1 * D2)
    F2 = F_coeffs(D1, D2, zX, zH, 40, level=2)
    F3 = F_coeffs(D1, D2, zX, zH, 40, level=3)
    B4 = {l: v[:4] for l, v in basis(xs, 4).items()}
    sc = phi(zX, zH, D1, D2) * mp.mpf(2 * D1 * D2) ** mp.mpf(3.5)
    out[name] = dict(noise23=max(abs(F3[l] / F2[l] - 1) for l in LS), A4_l3=float(extract(F3, B4, 4) / sc))
    print(name, out[name], flush=True)
(HERE / "data" / "wl_diag_level3.json").write_text(json.dumps(out, indent=1))
