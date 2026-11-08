"""self-test of four_density_path_gauge: the known map t^3 (1 + theta) from Legendre into E2 (star record)."""
import sys
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import four_density_path_gauge as G  # noqa: E402
u, t = G.u, G.t
E2 = [(u - 3) ** 2, -u * (u - 2)]
Leg = [u**2, -(u + 1) ** 2]
print("nullity box den=1, t^1..t^5 odd:", G.nullity(E2, Leg, 0, 2, {}, 1, 5, NT=200))
a = G.reconstruct(E2, Leg, 0, 2, {}, 1, 5, NT=200)
print("A =", a)
ok, R, odd = G.exact_checks(E2, Leg, a)
print("E2 o A = 0 mod Leg:", ok, " odd:", odd, " image annihilator (delta form):", [sp.factor(c.as_expr()) for c in R] if R else None)
