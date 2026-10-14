"""P3d test (pre-registered: docs/predictions/2026-10-13_three_electron_k3.md, commit 5542244e).

Holomorphic K3 period of the three-electron chain fingerprint (the analogue of
Paper 59 eq:period):
    P(ca,cb,cc) = int_{R^2} dp dq / sqrt((cb p^2+1)(cc q^2+1)(ca (p+q)^2+1)).
1-D Fourier (Paper 59 eq:K0 at D=0: int cos(kb)/sqrt(c k^2+1) dk = K0(|b|/sqrt c)/sqrt c)
turns it into a three-scale Bessel moment
    P = 8/(pi sqrt(ca cb cc)) * V3(1/sqrt ca, 1/sqrt cb, 1/sqrt cc),
    V3(a,b,c) = int_0^inf K0(at) K0(bt) K0(ct) dt.
Bailey (1936) / Bailey-Borwein-Broadhurst-Glasser (2008) eq. (33), c >= b >= a:
    (2c/pi) V3 = K(k-) K(k+') + K(k+) K(k-'),
    k+- = ( sqrt((c+a)^2-b^2) +- sqrt((c-a)^2-b^2) ) / (2c),   k' = sqrt(1-k^2).
Checks: route (i) direct 2-D quadrature of P; route (ii) Bessel moment (mpmath,
30 digits); route (iii) the product-of-elliptic closed form. Equal-mass anchor:
V3(1,1,1) = (pi/2) K3 K3' = 3 Gamma(1/3)^6 / (32 pi 2^(2/3)), BBBG eq. (47).
"""
from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.integrate import dblquad

mp.mp.dps = 34


def P_direct(ca: float, cb: float, cc: float) -> float:
    def f(tq: float, tp: float) -> float:
        p = np.tan(tp); q = np.tan(tq)
        jac = (1 + p * p) * (1 + q * q)
        return jac / np.sqrt((cb * p * p + 1) * (cc * q * q + 1) * (ca * (p + q) ** 2 + 1))
    h = np.pi / 2
    val, _ = dblquad(f, -h, h, -h, h, epsabs=1e-11, epsrel=1e-10)
    return val


def V3_quad(a, b, c):
    g = lambda t: mp.besselk(0, a * t) * mp.besselk(0, b * t) * mp.besselk(0, c * t)
    return mp.quad(g, [0, mp.mpf(1) / 8, 1, 4, mp.inf])


def V3_closed(a, b, c):
    a, b, c = sorted([mp.mpf(a), mp.mpf(b), mp.mpf(c)])        # c largest
    s1 = mp.sqrt((c + a) ** 2 - b ** 2); s2 = mp.sqrt((c - a) ** 2 - b ** 2 + 0j)
    kp = (s1 + s2) / (2 * c); km = (s1 - s2) / (2 * c)
    K = lambda k: mp.ellipk(k * k)                             # mpmath takes m = k^2
    Kc = lambda k: mp.ellipk(1 - k * k)                        # K(k')
    return mp.pi / (2 * c) * (K(km) * Kc(kp) + K(kp) * Kc(km))


def main() -> None:
    out = {}
    # equal-mass anchor
    v111 = V3_quad(1, 1, 1)
    k3 = mp.sin(mp.pi / 12)
    anchor = {"quad": mp.nstr(v111, 30),
              "gamma_form": mp.nstr(3 * mp.gamma(mp.mpf(1) / 3) ** 6 / (32 * mp.pi * mp.cbrt(4)), 30),
              "K3K3p": mp.nstr(mp.pi / 2 * mp.ellipk(k3**2) * mp.ellipk(1 - k3**2), 30),
              "closed33": mp.nstr(mp.re(V3_closed(1, 1, 1)), 30)}
    print("equal mass:", anchor)
    out["equal_mass"] = anchor
    rows = []
    for cs in [(0.21, 0.09, 0.16), (0.05, 0.24, 0.13), (0.25, 0.25, 0.01), (0.2, 0.2, 0.2)]:
        ca, cb, cc = (mp.mpf(x) for x in cs)
        m = [1 / mp.sqrt(x) for x in (ca, cb, cc)]
        vq = V3_quad(*m); vc = V3_closed(*m)
        pb = 8 / (mp.pi * mp.sqrt(ca * cb * cc)) * vq
        pd = P_direct(*cs)
        srt = sorted(m)
        regime = "real k" if (srt[2] - srt[0]) ** 2 >= srt[1] ** 2 else "complex-conjugate k"
        row = {"c": cs, "masses": [float(x) for x in m], "regime": regime,
               "V3_quad": mp.nstr(vq, 28), "V3_closed": mp.nstr(vc, 28),
               "quad_vs_closed": float(abs(vq - vc)),
               "P_direct": pd, "P_bessel": float(pb), "direct_vs_bessel": abs(pd - float(pb))}
        print(f"c={cs} [{regime}]  V3 quad={mp.nstr(vq,20)}  closed={mp.nstr(vc,20)}"
              f"  |diff|={row['quad_vs_closed']:.1e}   P direct={pd:.12f}"
              f"  bessel={float(pb):.12f}  |diff|={row['direct_vs_bessel']:.1e}")
        rows.append(row)
    out["rows"] = rows
    Path("debug/data/three_electron_p3d_check.json").write_text(json.dumps(out, indent=2, default=str))


if __name__ == "__main__":
    main()
