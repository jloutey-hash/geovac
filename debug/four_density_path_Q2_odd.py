"""FP step Q2: algebraic solutions of the odd-weight piece Qodd = 4 th (th-3)^2 - t^2 (th+1)^3 and of its adjoint
(reducibility), as a reference for R3."""
from __future__ import annotations
import sys
from fractions import Fraction
from pathlib import Path
import sympy as sp
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_star_F_hyperexp as F  # noqa: E402
u = F.u
Q = [4 * u * (u - 3) ** 2, -(u + 1) ** 3]
Qs = [sp.expand(q.subs(u, -u - 2 * i)) for i, q in enumerate(Q)]
for name, QQ in [("Qodd", Q), ("Qodd*", Qs)]:
    mu = max(r for r in sp.roots(sp.Poly(QQ[-1], u).as_expr(), u) if r.is_rational)
    print(name, F.algebraic_solutions(QQ, [0], [0, Fraction(-1, 2), Fraction(-3, 2), -1, Fraction(1, 2), Fraction(3, 2), 1], [0], mu + 8, N=200, verbose=True))
