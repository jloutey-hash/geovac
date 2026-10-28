"""WA7 derivation: the bond-angle singular part of (XY|XZ) when an H function has zero contact.

The H function is a cusp-free pair of 1s Slater functions, chi = phi_a - (a/b)^(5/2) phi_b (normalized phi_z), whose
first-order cusp cancels (Kato coefficient zero), so the density's non-smooth part at that H starts at |a|^3.
Uses the exact local machinery of debug/wa6_local.py.  Cases at D1 = D2 = D:
  'both': cusp-free on Y and Z;  'one': 1s (zeta_H) on Y, cusp-free on Z.
Usage: python debug/wa7_local.py
"""
import sympy as sp

from wa6_local import singular_part, taylor_poly, ax, ay, az, wx, wy, wz, Rw


def G_poly(kind, zX, zH, pair, D, deg):
    """Taylor polynomial of G(a), with the density's non-smooth part at Y = (0,0,D) equal to |a| G(a)."""
    NX = sp.sqrt(zX ** 3 / sp.pi)
    u = sp.sqrt(ax ** 2 + ay ** 2 + az ** 2)
    rr = sp.sqrt(ax ** 2 + ay ** 2 + (D + az) ** 2)
    if kind == '1s':
        odd = sp.sqrt(zH ** 3 / sp.pi) * sp.sinh(zH * u)
    else:
        a, b = pair
        c2 = -(a / b) ** sp.Rational(5, 2)
        odd = sp.sqrt(a ** 3 / sp.pi) * sp.sinh(a * u) + c2 * sp.sqrt(b ** 3 / sp.pi) * sp.sinh(b * u)
    return taylor_poly(-NX * sp.exp(-zX * rr) * odd / u, deg)


def water(case, zX, zH, pair, D):
    zX, zH, D = sp.nsimplify(zX), sp.nsimplify(zH), sp.nsimplify(D)
    pair = tuple(sp.nsimplify(p) for p in pair)
    c, s = sp.symbols('c s')
    kinds = ('free', 'free') if case == 'both' else ('1s', 'free')
    lead = 4 if case == 'both' else 2
    degs = {'1s': 2, 'free': 4}
    PY = G_poly(kinds[0], zX, zH, pair, D, degs[kinds[0]])
    PZ0 = G_poly(kinds[1], zX, zH, pair, D, degs[kinds[1]])
    bx, bz = c * ax - s * az, s * ax + c * az
    PZ = sp.expand(PZ0.subs({ax: bx, az: bz}, simultaneous=True))
    S = singular_part(PY, PZ, maxdeg=lead + 2)
    Y = sp.Matrix([0, 0, D])
    Z = sp.Matrix([D * s, 0, D * c])
    w = Y - Z
    S = sp.expand(S.subs({wx: w[0], wy: w[1], wz: w[2]}))
    S = sp.expand(S.subs(s ** 2, 1 - c ** 2))
    assert not S.has(s)
    v = sp.symbols('v', positive=True)
    S = sp.expand(S.subs(c, 1 - v ** 2).subs(Rw, D * sp.sqrt(2) * v))
    p0 = 7 + lead
    return {k: sp.N(S.coeff(v, k), 16) for k in range(1, p0 + 3, 2)}


if __name__ == "__main__":
    for case in ('both', 'one'):
        for zX in (1, '2.25'):
            co = water(case, zX, 1, (1, 2), '1.809')
            print(case, 'zX =', zX, {k: float(x) for k, x in co.items()})
