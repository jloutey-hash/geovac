import json
from pathlib import Path
from flint import nmod_mat, nmod_poly
import ib_kite_de as K
from ib_reduce import P
res = json.loads((Path(__file__).parent / "data/ib_tetra_pf_samples.json").read_text())
ms = res[0][1]
fitp = res[:-10]
xs = [r[0] for r in fitp]
print("len _NORM", len(K._NORM))
for m in ms[:3]:
    ys = [r[2][m] for r in fitp]
    f = K.fit_rational(xs, ys, len(xs) - 3)
    print(m, None if f is None else (f[0].degree(), f[1].degree()))
