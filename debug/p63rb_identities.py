"""Paper 63, the l > 0 extension of prop:tree_rung: measured checks of the three identities the proof adds.

(a) Multipole Poisson: Phi_L f(r) = 4 pi/(2L+1) int_0^inf r<^L / r>^(L+1) f(s) s^2 ds satisfies
    (theta - L)(theta + L + 1) Phi_L f = -4 pi r^2 f,   theta = r d/dr,   for L = 1, 2 and f = K0 K1, r^-1 K1^2.
(b) Hobson: for the harmonic polynomial z^2 - (x^2 + y^2)/2 = r^2 P_2(cos th), applying it as a differential
    operator to K0(r) gives r^2 P_2(cos th) (r^-1 d/dr)^2 K0, and (r^-1 d/dr)^2 K0 = r^-4 (theta - 2) theta K0 lies in
    r^-4 B (B = the Bessel module of K0).  Checked at an off-axis point.
(c) Transfer for a d-type root: with b(rho) = rho^2 (rho^-1 d/drho)^2 K0(rho) = rho^-2 (theta - 2) theta K0,
    int_0^inf b(t r) g(r) dr = t^-2 (theta_t - 2) theta_t T(r^-2 g),   T(h) = int_0^inf K0(t r) h(r) dr.

Output: debug/data/p63rb_identities.json.  Run: python debug/p63rb_identities.py
"""
from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp

mp.mp.dps = 30
H = mp.mpf('1e-6')


def d12(F, x):
    """First and second derivatives by the 5-point stencil (error ~h^4 ~ 1e-24)."""
    v = {k: F(x + k * H) for k in (-2, -1, 0, 1, 2)}
    d1 = (v[-2] - 8 * v[-1] + 8 * v[1] - v[2]) / (12 * H)
    d2 = (-v[-2] + 16 * v[-1] - 30 * v[0] + 16 * v[1] - v[2]) / (12 * H * H)
    return v[0], d1, d2
OUT = Path(__file__).parent / 'data' / 'p63rb_identities.json'


def K0(x):
    return mp.besselk(0, x)


def K1(x):
    return mp.besselk(1, x)


def phi_L(f, L, r):
    inner = mp.quad(lambda s: s ** L / r ** (L + 1) * f(s) * s ** 2, [0, r])
    outer = mp.quad(lambda s: r ** L / s ** (L + 1) * f(s) * s ** 2, [r, 2 * r, mp.inf])
    return 4 * mp.pi / (2 * L + 1) * (inner + outer)


def check_a():
    fs = {'K0K1': lambda s: K0(s) * K1(s), 'K1^2/r': lambda s: K1(s) ** 2 / s}
    rows = []
    for name, f in fs.items():
        for L in (1, 2):
            for r in (mp.mpf('0.7'), mp.mpf(2), mp.mpf('4.5')):
                P0, d1, d2 = d12(lambda u: phi_L(f, L, u), r)
                th = r * d1
                th2 = r * d1 + r ** 2 * d2
                lhs = th2 + th - L * (L + 1) * P0        # (theta - L)(theta + L + 1) = theta^2 + theta - L(L+1)
                rhs = -4 * mp.pi * r ** 2 * f(r)
                rows.append(dict(f=name, L=L, r=float(r), rel=float(abs(lhs / rhs - 1))))
    return rows


def check_b():
    x, y, z = mp.mpf('0.4'), mp.mpf('-0.9'), mp.mpf('1.1')
    F = lambda X, Y, Z: K0(mp.sqrt(X * X + Y * Y + Z * Z))
    lhs = (mp.diff(F, (x, y, z), (0, 0, 2))
           - (mp.diff(F, (x, y, z), (2, 0, 0)) + mp.diff(F, (x, y, z), (0, 2, 0))) / 2)
    r = mp.sqrt(x * x + y * y + z * z)
    c = z / r
    P2 = (3 * c * c - 1) / 2
    G = lambda u: K0(u)
    th = r * mp.diff(G, r, 1)
    th2 = th + r ** 2 * mp.diff(G, r, 2)
    rad = r ** -4 * (th2 - 2 * th)                     # r^-4 (theta - 2) theta K0
    rhs = r ** 2 * P2 * rad
    return dict(point=[float(x), float(y), float(z)], rel=float(abs(lhs / rhs - 1)))


def check_c():
    g = lambda r: r ** 3 * K1(r) * K0(r)
    rows = []
    for t in (mp.mpf('0.6'), mp.mpf('1.3'), mp.mpf('2.2')):
        def b(rho):
            th = rho * mp.diff(K0, rho, 1)
            th2 = th + rho ** 2 * mp.diff(K0, rho, 2)
            return rho ** -2 * (th2 - 2 * th)
        # b in closed form: theta K0 = -rho K1, theta^2 K0 = rho^2 K0, so b = rho^-2 (rho^2 K0 + 2 rho K1) = K0 + 2 K1/rho
        bc = lambda rho: K0(rho) + 2 * K1(rho) / rho
        lhs = mp.quad(lambda r: bc(t * r) * g(r), [0, 1, mp.inf])
        T = lambda u: mp.quad(lambda r: K0(u * r) * g(r) / r ** 2, [0, 1, mp.inf])
        _, d1, d2 = d12(T, t)
        th = t * d1
        th2 = th + t ** 2 * d2
        rhs = t ** -2 * (th2 - 2 * th)
        rows.append(dict(t=float(t), rel=float(abs(lhs / rhs - 1)),
                         b_closed_vs_operator=float(abs(bc(mp.mpf('1.7')) / b(mp.mpf('1.7')) - 1))))
    return rows


def main():
    res = dict(a=check_a(), b=check_b(), c=check_c(), dps=mp.mp.dps)
    res['worst'] = dict(a=max(r['rel'] for r in res['a']), b=res['b']['rel'], c=max(r['rel'] for r in res['c']))
    OUT.write_text(json.dumps(res, indent=1))
    print(json.dumps(res['worst'], indent=1))


if __name__ == '__main__':
    main()
