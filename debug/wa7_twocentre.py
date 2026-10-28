"""Degree-4 two-centre check for WA7: odd terms of J(R) for two e^-2r/pi densities through R^11."""
import sympy as sp

from wa6_local import singular_part, taylor_poly, ax, ay, az, wx, wy, wz, Rw

u = sp.sqrt(ax ** 2 + ay ** 2 + az ** 2)
G = taylor_poly(-(1 / sp.pi) * sp.sinh(2 * u) / u, 4)
R = sp.symbols('R', positive=True)
print(sp.expand(singular_part(G, G, maxdeg=4).subs({wx: 0, wy: 0, wz: R, Rw: R})))
