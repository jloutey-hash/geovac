r"""Independent verification of Paper 59's irregularity + irreducibility wall.

Checks, from scratch (sympy), the four CONCRETE pillars of Sec.~obstruction:
  (1) eq:pf is IRREGULAR at D=oo with four exponential rates
      lambda in {+-1, +-i sqrt((1-rho)/rho)} and a slope-one Newton edge
      (Poincare rank one) -- the reason the Fuchsian sunrise machinery cannot reach it.
  (2) eq:pf is FORMALLY SELF-ADJOINT (L* = L) -> differential Galois group in Sp_4.
  (3) indicial exponents at the regular singular point D=0 are {0,1,1,2}
      (the double 1 is the D ln D log source).
  (4) the explicit monodromy M0 is unipotent with a SINGLE 2x2 Jordan block and
      NO invariant coordinate subspace (the explicit half of "irreducible").

The D-module "in principle" argument (FL(simple)=simple, Katz 1990) is standard
and is audited in prose, not recomputed here.

eq:pf:  D rho y'''' + 2 rho y''' + D(1-2rho) y'' + (1-2rho) y' - D(1-rho) y = 0
"""
from __future__ import annotations

import sympy as sp

P = lambda *a: print(*a, flush=True)
D, rho, s, lam = sp.symbols('D rho s lambda')
y = sp.Function('y')

# coefficients a4..a0 (of y'''', y''', y'', y', y)
a = [D * rho, 2 * rho, D * (1 - 2 * rho), (1 - 2 * rho), -D * (1 - rho)]  # a4,a3,a2,a1,a0
ok = True


def check(name, cond):
    global ok
    ok = ok and bool(cond)
    P(f"  [{'PASS' if cond else 'FAIL'}] {name}")


P("=== (1) irregularity at D=oo: exponential rates + Newton slope ===")
# Seek y ~ exp(lambda D): the coefficient of the highest power of D in the balance.
# Terms carrying an explicit factor D are a4,a2,a0; their lambda-polynomial is the symbol.
symbol = sp.expand(a[0] / D * lam**4 + a[2] / D * lam**2 + a[4] / D)  # divide the common D
symbol = sp.simplify(symbol)
claimed = sp.expand(rho * (lam**2 - 1) * (lam**2 + (1 - rho) / rho))
P(f"    symbol(lambda) = {symbol}")
check("symbol == rho(l^2-1)(l^2+(1-rho)/rho)", sp.simplify(symbol - claimed) == 0)
roots = sp.solve(sp.Eq(symbol, 0), lam)
P(f"    lambda roots (symbolic) = {roots}")
# compare NUMERICALLY at rho=1/5 (0<rho<1): sqrt(1-1/rho)=i*sqrt((1-rho)/rho).
r0 = sp.Rational(1, 5)
got = sorted((complex(sp.N(rt.subs(rho, r0))) for rt in roots), key=lambda z: (z.real, z.imag))
want = sorted((complex(sp.N(w.subs(rho, r0))) for w in
               (sp.Integer(1), sp.Integer(-1),
                sp.I * sp.sqrt((1 - rho) / rho), -sp.I * sp.sqrt((1 - rho) / rho))),
              key=lambda z: (z.real, z.imag))
P(f"    at rho=1/5: roots={got}  vs paper's {want}")
check("four rates {+-1, +-i sqrt((1-rho)/rho)} (numeric, rho=1/5)",
      all(abs(g - w) < 1e-12 for g, w in zip(got, want)))
# slope-one Newton edge / Poincare rank one: rates are finite & nonzero (not lambda=0),
# i.e. the operator is irregular (an exp(lambda D), not a power, dominates).
check("all four rates nonzero (irregular, not regular-singular at oo)",
      all(sp.simplify(r) != 0 for r in roots))

P("\n=== (2) formal self-adjointness L* == L ===")
f = sp.Function('f')(D)
L = sum(a[i] * sp.diff(f, D, 4 - i) for i in range(5))
# formal adjoint: L*[f] = sum_k (-1)^k d^k/dD^k (a_k f), k = order = 4-i
Lstar = sum((-1)**(4 - i) * sp.diff(a[i] * f, D, 4 - i) for i in range(5))
check("L*[f] - L[f] == 0 identically", sp.simplify(sp.expand(Lstar - L)) == 0)

P("\n=== (3) indicial exponents at D=0 ===")
# substitute y = D**s, take the lowest power of D; its coefficient is the indicial poly.
expr = sum(a[i] * sp.diff(D**s, D, 4 - i) for i in range(5))
expr = sp.expand(expr)
# lowest power of D present:
powers = {sp.simplify(t.as_independent(D)[1]).as_base_exp()[1]
          if t.has(D) else sp.oo for t in expr.as_ordered_terms()}
# robustly: factor D**(s-3) and evaluate the bracket at leading order
lead = sp.simplify(expr / D**(s - 3))
lead0 = sp.simplify(sp.limit(lead, D, 0))
indicial = sp.factor(lead0)
P(f"    indicial polynomial (coeff of D^(s-3)) = {indicial}")
sol = sp.roots(sp.Poly(sp.expand(lead0), s))  # {root: multiplicity}
P(f"    exponents (with multiplicity) = {sol}")
check("indicial roots {0, 1(double), 2}",
      sol == {sp.Integer(0): 1, sp.Integer(1): 2, sp.Integer(2): 1})

P("\n=== (4) explicit monodromy M0: unipotent, single 2x2 Jordan block, no invariant coord subspace ===")
M0 = sp.Matrix([[-1, 2, 2, 2], [-2, 3, 2, 2], [-2, 2, 3, 2], [2, -2, -2, -1]])
I4 = sp.eye(4)
check("all eigenvalues == 1 (unipotent)", M0.eigenvals() == {sp.Integer(1): 4})
N = M0 - I4
check("(M0-I)^2 == 0 (Jordan blocks size <= 2)", (N * N) == sp.zeros(4, 4))
check("rank(M0-I) == 1 (exactly one 2x2 Jordan block)", N.rank() == 1)
# no invariant COORDINATE subspace: no nonempty proper subset of basis vectors
# spans an M0-invariant subspace.
import itertools
invariant_coord = False
for k in (1, 2, 3):
    for combo in itertools.combinations(range(4), k):
        cols = M0[:, list(combo)]
        # invariant iff every image column lies in span of the chosen coords,
        # i.e. rows outside `combo` are all zero in these columns
        outside = [r for r in range(4) if r not in combo]
        if all(cols[r, c] == 0 for r in outside for c in range(k)):
            invariant_coord = True
check("no invariant coordinate subspace", not invariant_coord)

P("\n=== VERDICT ===")
P("  ALL PILLARS VERIFIED" if ok else "  *** A PILLAR FAILED ***")
P("  (remaining leg = the standard Fourier-Laplace-preserves-simplicity theorem,")
P("   Katz 1990 -> FL(j!* L) simple of generic rank 4 -> irreducible; audited in prose.)")
raise SystemExit(0 if ok else 1)
