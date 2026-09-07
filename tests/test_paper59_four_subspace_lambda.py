"""Backing test for Paper 61 (companion of Paper 59) sec:modular remark "The modulus is a four-point cross-ratio":
the Legendre modulus lambda = 1-rho of the three-center elliptic Bessel moment is a member of
the cross-ratio orbit of the four branch points {0, 1, oo, (rho-1)/rho} of the reduced elliptic
curve v^2 = u(u-1)(rho u + 1-rho) (u = x^2), and that orbit contains BOTH period arguments
rho and 1-rho.  Exact / symbolic (sympy), so it verifies the claim, not a numerical coincidence.

This places Paper 61's period in the classical four-singular-point / tame four-subspace (D~4)
lineage; the four points are the two spectral masses' branch data (rho = c2/c1), NOT four
scattering centers -- the modulus counts densities, not nuclei (the honest disambiguation of
the four-subspace <-> elliptic-lambda cross-corpus lead, 2026-08-25).
"""
import sympy as sp


def _cross_ratio_orbit(mu):
    """S3 orbit of the modulus for four branch points {0, 1, oo, mu}: {mu, 1-mu, 1/mu, ...}."""
    return [sp.simplify(v) for v in (mu, 1 - mu, 1 / mu, 1 / (1 - mu), mu / (mu - 1), (mu - 1) / mu)]


def test_lambda_is_cross_ratio_of_branch_points():
    rho = sp.symbols('rho', positive=True)
    mu = (rho - 1) / rho                                   # 4th branch point of the reduced curve
    orbit = _cross_ratio_orbit(mu)

    # Paper 61's Legendre lambda = 1 - rho is in the orbit ...
    assert any(sp.simplify(o - (1 - rho)) == 0 for o in orbit)
    # ... and so is rho itself (the other Legendre period's argument K(rho))
    assert any(sp.simplify(o - rho) == 0 for o in orbit)
    # the orbit is exactly the six expected rational functions (no extra/missing members)
    expected = {sp.simplify(e) for e in
                ((rho - 1) / rho, 1 / rho, rho / (rho - 1), rho, 1 - rho, -1 / (rho - 1))}
    assert {sp.simplify(o) for o in orbit} == expected


def test_reduced_curve_branch_points_and_period_reduction():
    """The reduced elliptic curve v^2 = u(u-1)(rho u + 1-rho) has finite branch points
    {0, 1, (rho-1)/rho} (fourth at infinity); and u = x^2 reduces the quartic period 2:1 --
    dx/sqrt(Q) = (1/2) du/sqrt(cubic) -- verified as the differential identity
    x/sqrt(cubic(x^2)) = 1/sqrt(Q(x)) for x>0 (the du/(2 sqrt u) Jacobian, not quartic=cubic(x^2))."""
    x, rho = sp.symbols('x rho', positive=True)
    u = sp.symbols('u')                                    # generic (roots include 0 and a negative)
    cubic = u * (u - 1) * (rho * u + (1 - rho))
    # finite branch points = roots of the cubic
    roots = set(sp.solve(cubic, u))
    assert roots == {sp.Integer(0), sp.Integer(1), sp.simplify((rho - 1) / rho)}
    # period reduction (differential level): x / sqrt(cubic(x^2)) == 1 / sqrt(Q(x)) for x>0
    Q = (x**2 - 1) * (rho * x**2 + (1 - rho))
    assert sp.simplify(x / sp.sqrt(cubic.subs(u, x**2)) - 1 / sp.sqrt(Q)) == 0
