"""Paper 63, l > 0 rung bound: end-to-end check of the multipole/Gaunt reduction on a shape with an internal l > 0
density and an L = 2 channel (the cases the proof of the l > 0 extension adds).

Chain v1 - v2 - v3 at the fibre (all centres at the origin), Coulomb kernels, fingerprint weights:
    rho1 = d/dz [K0(r)/(2 pi^2)]          = f(r) cos th,  f = -K1/(2 pi^2)      (p_z leaf)
    rho2 = d/dz [K0(r)/(2 pi^2)]          = f(r) cos th                           (p_z, internal)
    rho3 = h(d) [t^3 K0(t r)/(2 pi^2)],  h = z^2 - (x^2 + y^2)/2 (harmonic, r^2 P2)
         = t^5/(2 pi^2) b(t r) P2(cos th),  b(rho) = K0 + 2 K1/rho                (d root, scale t)
Route R3 (rooted at v3): leaf potential g1 = Phi_1 f; W2 = f g1 cos^2 = f g1 (1/3 + (2/3) P2);
    I = (4 pi/5) int r^2 R3(r) Phi_2[(2/3) f g1] dr.
Route R1 (rooted at v1): Phi_2 R3 at v2; W2 = f cos * Phi_2 R3 P2 = f Phi_2R3 ((2/5) P1 + (3/5) P3);
    I = (4 pi/3) int r^2 f Phi_1[(2/5) f Phi_2 R3] dr.
Phi_L u(r) = 4 pi/(2L+1) [ r^-(L+1) int_0^r s^(L+2) u + r^L int_r^inf s^(1-L) u ].
Route MC: the full 9-dimensional integral, importance-sampled, no multipoles.

Output debug/data/p63rb_chain_check.json.  Run: python debug/p63rb_chain_check.py [n_mc_millions]
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_simpson, simpson
from scipy.special import k0, k1, k0e, k1e

OUT = Path(__file__).parent / 'data' / 'p63rb_chain_check.json'
TS = (0.5, 1.3)

u = np.linspace(np.log(1e-9), np.log(80.0), 400001)      # r = e^u
r = np.exp(u)


def cum0(y):
    """int_0^r y(s) ds on the grid, y given in s (ds = s du)."""
    return cumulative_simpson(y * r, x=u, initial=0.0)


def tail(y):
    """int_r^inf y(s) ds, accumulated from the far end (no cancellation against a large total)."""
    z = (y * r)[::-1]
    return cumulative_simpson(z, x=-u[::-1], initial=0.0)[::-1]


def phi(Lq, w):
    return 4 * np.pi / (2 * Lq + 1) * (r ** (-(Lq + 1)) * cum0(r ** (Lq + 2) * w) + r ** Lq * tail(r ** (1 - Lq) * w))


f = -k1(r) / (2 * np.pi ** 2)


def radial(t):
    R3 = t ** 5 / (2 * np.pi ** 2) * (k0(t * r) + 2 * k1(t * r) / (t * r))
    g1 = phi(1, f)
    I3 = 4 * np.pi / 5 * simpson(r ** 2 * R3 * phi(2, (2 / 3) * f * g1) * r, x=u)
    I1 = 4 * np.pi / 3 * simpson(r ** 2 * f * phi(1, (2 / 5) * f * phi(2, R3)) * r, x=u)
    return I3, I1


def mc(t, n_total, chunk=2_000_000, seed=20261117):
    """Importance sampling: v1, v2 radius ~ Gamma(2) (3-D pdf e^-r/(4 pi r)), v3 radius ~ Exp(t) (3-D pdf
    t e^{-t r}/(4 pi r^2)); directions uniform.  Antithetic in the v3 direction (P2 even, the rest odd x odd)."""
    rng = np.random.default_rng(seed)
    acc = []
    done = 0
    while done < n_total:
        m = min(chunk, n_total - done)

        def dirs(k):
            v = rng.normal(size=(k, 3))
            return v / np.linalg.norm(v, axis=1)[:, None]
        a1, a2, a3 = rng.gamma(2.0, size=m), rng.gamma(2.0, size=m), rng.exponential(1.0 / t, size=m)
        d1, d2, d3 = dirs(m), dirs(m), dirs(m)
        x1, x2, x3 = a1[:, None] * d1, a2[:, None] * d2, a3[:, None] * d3
        p1 = np.exp(-a1) / (4 * np.pi * a1)
        p2 = np.exp(-a2) / (4 * np.pi * a2)
        p3 = t * np.exp(-t * a3) / (4 * np.pi * a3 ** 2)
        # densities with exponentially scaled Bessels to avoid underflow: K(x) = Ke(x) e^{-x}
        rho1 = -k1e(a1) * np.exp(-a1) * d1[:, 2] / (2 * np.pi ** 2)
        rho2 = -k1e(a2) * np.exp(-a2) * d2[:, 2] / (2 * np.pi ** 2)
        tr = t * a3
        P2 = 1.5 * d3[:, 2] ** 2 - 0.5
        rho3 = t ** 5 / (2 * np.pi ** 2) * (k0e(tr) + 2 * k1e(tr) / tr) * np.exp(-tr) * P2
        w = rho1 * rho2 * rho3 / (p1 * p2 * p3)
        k12 = 1 / np.linalg.norm(x1 - x2, axis=1)
        val = w * k12 * (1 / np.linalg.norm(x2 - x3, axis=1) + 1 / np.linalg.norm(x2 + x3, axis=1)) / 2
        acc.append(val)
        done += m
    v = np.concatenate(acc)
    return float(v.mean()), float(v.std(ddof=1) / np.sqrt(v.size))


def main():
    n_mc = int(float(sys.argv[1]) * 1e6) if len(sys.argv) > 1 else 20_000_000
    res = {}
    for t in TS:
        t0 = time.time()
        I3, I1 = radial(t)
        m, se = mc(t, n_mc)
        res[str(t)] = dict(root_v3=I3, root_v1=I1, root_rel=abs(I3 / I1 - 1), mc=m, mc_se=se,
                           z=(m - I3) / se, n_mc=n_mc, seconds=time.time() - t0)
        print(t, res[str(t)], flush=True)
    OUT.write_text(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
