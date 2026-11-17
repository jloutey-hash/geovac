"""Record P63L2 controls D1 and D6 (no branch at eps < 0.1 of the scored integrands is evaluated).
D1: the referee's 30-term tail series vs the float Laguerre route at x = 32, 40, 60 (L = 0, 2).
D6: the referee's tail quadrature (breakpoints of p63l2_run.referee) on the toy x e^{-2x} weight with the B1 kernel,
    eps = 1e-2, 1e-3, 1e-4, against a dps-45 quadrature with fine breakpoints (toy, not the scored integrand).
Writes debug/data/p63l2_d1.json."""
import json
import sys
from pathlib import Path

import mpmath as mp
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import p63l2_run as R  # noqa: E402
import p63odd_engine as EN  # noqa: E402

mp.mp.dps = 30
ser = R.series_coeffs()
out = {"D1": {}, "D6": {}}
for L in (0, 2):
    xs = np.array([32.0, 40.0, 60.0])
    f = EN.phix_f(xs, L)
    out["D1"][L] = [float(abs(fv / R.phi_series(x, ser[L]) - 1)) for x, fv in zip(xs, f)]
pts = [R.XCUT, 64, 128, 256, 512, 1024, 4096, 16384, 65536, mp.inf]
for eps in ("1e-2", "1e-3", "1e-4"):
    e = mp.mpf(eps)
    toy = lambda x: R.kern((2 - e) * x, 2) * mp.exp(-e * x) * x  # noqa: E731
    a = mp.quad(toy, pts)
    with mp.workdps(45):
        fine = [R.XCUT] + [R.XCUT * mp.mpf(2) ** k for k in range(1, 40)] + [mp.inf]
        b = mp.quad(toy, fine)
    out["D6"][eps] = float(abs(a / b - 1))
print(json.dumps(out, indent=1))
(HERE / "data" / "p63l2_d1.json").write_text(json.dumps(out, indent=1))
