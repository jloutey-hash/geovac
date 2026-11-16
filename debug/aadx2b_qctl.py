"""Record AX2b bar control (measured before registering): self-differences of BOTH sides of the AX0 (ii)
comparison at the fresh points L = 16, 25 -- the direct-projection quadrature (44 vs 56 digits) and the exact
split's evaluation (300 vs 400 digits).  No engine-vs-quadrature comparison is printed."""
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import mpmath as mp  # noqa: E402
import aadx_ap0 as A  # noqa: E402
import aapc_exact as P  # noqa: E402
import na_exact as E  # noqa: E402

K = E.Field(None)
for case in P.CASES:
    for integral in ('DX', 'DD'):
        f = P.coefficients(case, 25, K, integral)
        keys = set(k for q in f for k in q)
        for L in (16, 25):
            t = time.time()
            q44 = A.quad_GL(case, L, integral, dps=44)
            q56 = A.quad_GL(case, L, integral, dps=56)
            ev = {}
            for dps in (300, 400):
                mp.mp.dps = dps
                tw = P.twopiN(case, integral)
                cv = {k: E.value(k, mp) for k in keys}
                ev[dps] = tw * sum(mp.mpf(int(v.p)) / int(v.q) * cv[k] for k, v in f[L].items())
            mp.mp.dps = 80
            print(case, integral, 'L', L, 'quad |q44/q56 - 1|', mp.nstr(abs(q44 / q56 - 1), 3),
                  'exact |e300/e400 - 1|', mp.nstr(abs(ev[300] / ev[400] - 1), 3), round(time.time() - t), 's',
                  flush=True)
