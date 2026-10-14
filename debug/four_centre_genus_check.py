"""Four-centre genus check (pre-registered: docs/predictions/2026-10-13_four_centre_genus.md).

Claim under test (P1): the general four-centre 1s Slater ERI (XY|ZW) reduces to
    (8/pi) int_0^1 ds int_0^1 dt int_0^inf dk  j0(k |W(s,t)|) G(k,s;D1) G(k,t;D2)
with  G(k,t;D) = c e^{-D Delta} (D^2/Delta^3 + 3D/Delta^4 + 3/Delta^5),
      c = t(1-t),  Delta = sqrt(c k^2 + 1)          (all exponents zeta = 1)
      W(s,t) = P1(s) - P2(t),  P(t) = t B + (1-t) A for the pair (A, B).
The ONLY square roots are Delta1(s), Delta2(t): the Paper 59 fibre curve.
G is the analytic d^2/dzA dzB of the Yukawa-product kernel e^{-D Delta}/Delta.

Routes compared:
  (a) McMurchie-Davidson on n-Gaussian STO fits, n = 6..16 (converges to Slater)
  (b) the reduced (s,t,k) formula above (true Slater, no Gaussians)
Sanity anchor: Paper 59's three-centre (XY|XZ) Slater value 0.2049416525.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from geovac import noci_engine as GE  # noqa: E402


def G(k: float, t: float, D: float) -> float:
    c = t * (1.0 - t)
    De = np.sqrt(c * k * k + 1.0)
    return c * np.exp(-D * De) * (D * D / De**3 + 3.0 * D / De**4 + 3.0 / De**5)


def j0(x: float) -> float:
    return 1.0 if x < 1e-12 else np.sin(x) / x


def reduced_eri(A, B, C, Dd, n_st: int = 40) -> float:
    """(AB|CD) via the reduced two-square-root formula. sin^2 maps on s, t."""
    A, B, C, Dd = map(np.asarray, (A, B, C, Dd))
    D1 = float(np.linalg.norm(A - B)); D2 = float(np.linalg.norm(C - Dd))
    x, wx = leggauss(n_st)
    th = 0.25 * np.pi * (x + 1.0)                 # theta in (0, pi/2)
    u = np.sin(th) ** 2                           # u in (0, 1)
    wu = wx * 0.25 * np.pi * 2.0 * np.sin(th) * np.cos(th)
    total = 0.0
    for si, ws in zip(u, wu):
        P1 = si * B + (1.0 - si) * A
        for ti, wt in zip(u, wu):
            P2 = ti * Dd + (1.0 - ti) * C
            Wn = float(np.linalg.norm(P1 - P2))
            val, _ = quad(lambda k: j0(k * Wn) * G(k, si, D1) * G(k, ti, D2),
                          0.0, np.inf, limit=400, epsabs=1e-13, epsrel=1e-11)
            total += ws * wt * val
    return 8.0 / np.pi * total


def md_eri(A, B, C, Dd, n_gauss: int) -> float:
    sh_a, sh_d, _ = GE.fit_sto_shape(0, 1, n_gauss=n_gauss)
    shapes = {"1s": (sh_a, sh_d)}
    f = [GE.sto_shape_basis(np.asarray(P, float), "1s", 1.0, shapes, (0, 0, 0))
         for P in (A, B, C, Dd)]
    return float(GE.eri_md(*f))


def min_W(A, B, C, Dd, n: int = 401) -> float:
    s = np.linspace(0, 1, n)
    P1 = s[:, None] * np.asarray(B) + (1 - s[:, None]) * np.asarray(A)
    P2 = s[:, None] * np.asarray(Dd) + (1 - s[:, None]) * np.asarray(C)
    return float(np.min(np.linalg.norm(P1[:, None, :] - P2[None, :, :], axis=-1)))


def main() -> None:
    X = [0.0, 0.0, 0.0]; Y = [0.0, 0.0, 2.0]; Z = [1.5, 0.0, 0.5]
    four = ([0.0, 0.0, 0.0], [0.0, 0.0, 1.4], [1.9, 0.4, -0.2], [2.2, 1.1, 1.0])
    out = {}
    cases = {"three_centre_XYXZ": (X, Y, X, Z), "four_centre_XYZW": four}
    for name, geo in cases.items():
        md = {n: md_eri(*geo, n) for n in (6, 8, 10, 12, 14, 16)}
        red = {n: reduced_eri(*geo, n_st=n) for n in (24, 40, 56)}
        out[name] = {"md_by_ngauss": md, "reduced_by_nst": red,
                     "min_W": min_W(*geo)}
        print(name)
        for n, v in md.items():
            print(f"  MD  n_gauss={n:2d}  {v:.12f}")
        for n, v in red.items():
            print(f"  RED n_st={n:2d}     {v:.12f}")
        print(f"  min |W(s,t)| = {out[name]['min_W']:.6f}")
    Path("debug/data").mkdir(parents=True, exist_ok=True)
    Path("debug/data/four_centre_genus_check.json").write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
