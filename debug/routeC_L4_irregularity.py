"""Operator-level confirmation that Paper 59's L4 (eq:pf) is IRREGULAR at D=infinity,
head-to-head against the (Fuchsian) sunrise/Legendre elliptic-period operator.

This is the operator-level backing for the paper's structural claim that N(D) is the
Fourier-Laplace / Laplace-dual of the sunrise's elliptic differential -- a rank-4
IRREGULAR connection "one transcendence level above the sunrise" -- and therefore that
the 2019 all-orders elliptic-polylog closed form of the unequal-mass sunrise
(Bogner-Mueller-Stach-Weinzierl, arXiv:1907.01251, iterated integrals on Mbar_{1,3},
a Fuchsian object) does NOT transfer to L4.

Method: a linear ODE has an irregular singularity of Poincare rank 1 at infinity iff it
admits formal solutions e^{lambda x} with constant lambda != 0 (the Bessel equation is
the archetype). We extract the leading exponential balance (the top Newton-polygon edge)
of each operator and read off the lambda-spectrum.

Memo: this file. Backing test: tests/test_routeC_momentum.py::test_L4_irregular_at_infinity.
"""
import sympy as sp

rho = sp.symbols('rho', positive=True)
lam, X = sp.symbols('lambda X')   # X = the working variable (D for L4, rho-hat for Legendre)


def exponential_balance(coeffs, var):
    """coeffs[k] multiplies d^k/dvar^k. Substitute y=e^{lam*var}, collect in `var`,
    return (top_degree, leading_polynomial_in_lambda). The leading poly's nonzero roots
    are the rank-1 irregular exponents; if the only root is lam=0 the point is regular."""
    expr = sum(c * lam**k for k, c in coeffs.items())     # e^{lam var} factored out
    poly_in_var = sp.Poly(sp.expand(expr), var)
    top = poly_in_var.degree()
    lead = sp.factor(poly_in_var.nth(top))
    return top, sp.expand(lead)


# ---------------------------------------------------------------- L4  (Paper 59 eq:pf)
# D rho L'''' + 2 rho L''' + D(1-2rho) L'' + (1-2rho) L' - D(1-rho) L = 0
D = X
L4 = {4: D*rho, 3: 2*rho, 2: D*(1-2*rho), 1: (1-2*rho), 0: -D*(1-rho)}
top4, P_lead = exponential_balance(L4, D)
lam_roots = sp.roots(sp.Poly(P_lead, lam))

print("=" * 72)
print("L4  (rank-4, Laplace variable D) -- Paper 59 eq:pf")
print(f"  top Newton edge sits at D^{top4}  (operator order 4)  => slope-1 edge of full width")
print(f"  leading exponential-balance polynomial P(lambda) = {P_lead}")
print(f"  formal exponential exponents lambda (roots of P): {dict(lam_roots)}")
allnonzero = all(r != 0 for r in lam_roots)
print(f"  ALL lambda nonzero? {allnonzero}   (=> every solution ~ e^{{lambda D}}: fully irregular, no regular soln at inf)")

# cross-check: the four lambda are exactly the four branch points s = +-1, +-i sqrt((1-rho)/rho)
expected = {sp.Integer(1): 1, sp.Integer(-1): 1,
            sp.I*sp.sqrt((1-rho)/rho): 1, -sp.I*sp.sqrt((1-rho)/rho): 1}
match = sp.simplify(P_lead - rho*(lam**2-1)*(lam**2+(1-rho)/rho)) == 0
print(f"  P(lambda) == rho*(lam^2-1)*(lam^2+(1-rho)/rho) ? {match}   (lambda-spectrum = the 4 branch pts: K0/I0 & J0/Y0 Bessel sectors)")

# ------------------------------------------------ Legendre / elliptic-period operator
# M_rho = rho(1-rho) d^2 + (1-2rho) d - 1/4   (periods K(rho), K(1-rho); the sunrise CURVE operator)
r = X
Leg = {2: r*(1-r), 1: (1-2*r), 0: sp.Rational(-1, 4)}
top2, Q_lead = exponential_balance(Leg, r)
lam_roots_leg = sp.roots(sp.Poly(Q_lead, lam))
print("=" * 72)
print("Legendre elliptic-period operator M_rho -- the sunrise-class (Fuchsian) object")
print(f"  top edge at rho^{top2}; leading balance Q(lambda) = {Q_lead}")
print(f"  formal exponents lambda: {dict(lam_roots_leg)}   (only lambda=0 => REGULAR singular at inf: Fuchsian)")
print("=" * 72)
print("VERDICT:")
print("  L4:  irregular singularity at D=inf (Poincare rank 1, 4 nonzero Bessel exponents).")
print("  Sunrise/Legendre operator: Fuchsian (regular singular at inf, only lambda=0).")
print("  => distinct D-module classes across the Laplace transform. The Fuchsian")
print("     elliptic-polylog technology that closes the unequal-mass sunrise (1907.01251)")
print("     does not apply to the irregular L4. Paper 59's 'one level above the sunrise' holds.")
