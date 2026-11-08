"""FD1-FD3 step D: theta_t^j S(t), j = 0..19, by quadrature with theta applied under the integral,
cached to debug/data/four_density_star_D.json.  S = (16/pi^4) t^3 int chi^3 K0(t x) x^-1 dx;
theta_t^j [t^3 K0(t x)] = t^3 [(theta_y + 3)^j K0](t x), theta_y^2 K0 = y^2 K0, theta_y K0 = -y K1.
Usage: python four_density_star_D_thetas.py dps t1 t2 ...
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import mpmath as mp
import sympy as sp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import density_loop_B_ode as B  # noqa: E402

y = B.y
OUT = HERE / "data" / "four_density_star_D.json"
JMAX = 19
_CHI = {}


def chi3(xx):
    key = (mp.mp.dps, xx)
    if key not in _CHI:
        _CHI[key] = B.chi(xx) ** 3
    return _CHI[key]


def theta_funcs(J):
    al, be, _ = B.theta_powers(J, y**2, 0)
    out = []
    for j in range(J + 1):
        A = sp.expand(sum(sp.binomial(j, i) * 3 ** (j - i) * al[i] for i in range(j + 1)))
        Bb = sp.expand(sum(sp.binomial(j, i) * 3 ** (j - i) * be[i] for i in range(j + 1)))
        out.append((sp.lambdify(y, A, "mpmath"), sp.lambdify(y, Bb, "mpmath")))
    return out


def theta_S(tv, fns):
    tv = mp.mpf(tv)
    vals = []
    for A, Bb in fns:
        def h(xx):
            yy = tv * xx
            return chi3(xx) / xx * (A(yy) * mp.besselk(0, yy) - Bb(yy) * yy * mp.besselk(1, yy))
        vals.append(16 / mp.pi**4 * tv**3 * mp.quad(h, [0, mp.mpf(1) / 8, 1, 4, 16, 64, mp.inf]))
    return vals


def main() -> None:
    dps = int(sys.argv[1])
    mp.mp.dps = dps
    fns = theta_funcs(JMAX)
    rec = json.loads(OUT.read_text()) if OUT.exists() else {}
    for tv in sys.argv[2:]:
        key = f"{tv}@{dps}"
        if key in rec:
            continue
        t0 = time.time()
        vals = theta_S(tv, fns)
        rec[key] = [mp.nstr(v, dps) for v in vals]
        OUT.write_text(json.dumps(rec, indent=1))
        print(f"t={tv} dps={dps}: S = {mp.nstr(vals[0], 20)}  ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
