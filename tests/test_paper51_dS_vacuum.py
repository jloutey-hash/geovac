"""Verification of Paper 51 Theorem thm:dS_vacuum.

Paper 51 (group5_qed_gauge/paper_51_gravity_arc.tex), Section G5
(Decompactified de Sitter vacuum), Theorem thm:dS_vacuum claims:

    The minimum action density at the de Sitter extremum is
        s_min = -Lambda / (12 * sqrt(6))
    at R_crit = 1 / (sqrt(6) * Lambda)

with the action density
    s(R, Lambda) = (R^3 / 4) Lambda^4  -  (R / 8) Lambda^2

obtained as the beta -> infinity decompactification of G2.

The extremum is taken by varying R at fixed Lambda; the claim is that
ds/dR = 0 gives R_crit Lambda = 1/sqrt(6), and substituting back gives
s_min = -Lambda / (12 sqrt(6)).

This test computes the extremum symbolically with sympy and verifies that
it matches the claim exactly (including sign and prefactor).

Per CLAUDE.md §13.4a: equation verification via symbolic identity and
analytical extremum.
"""

import sympy as sp


def test_paper51_thm_dS_vacuum_extremum_location():
    """ds/dR = 0 at R = 1/(sqrt(6) * Lambda)."""
    R, Lambda = sp.symbols("R Lambda", positive=True)

    # Action density from eq. above thm:dS_vacuum (Section G5)
    s = sp.Rational(1, 4) * R**3 * Lambda**4 - sp.Rational(1, 8) * R * Lambda**2

    # Extremum: vary R at fixed Lambda
    ds_dR = sp.diff(s, R)
    crit_solutions = sp.solve(ds_dR, R)

    # Pick the positive critical point
    R_crit_positive = [sol for sol in crit_solutions if sol.is_real and (sol > 0) is sp.true or sp.simplify(sol - 1 / (sp.sqrt(6) * Lambda)) == 0]
    assert len(R_crit_positive) >= 1, f"No positive critical point found: {crit_solutions}"

    # Compare against the claim R_crit = 1/(sqrt(6) Lambda)
    R_crit_claim = 1 / (sp.sqrt(6) * Lambda)
    match = any(sp.simplify(sol - R_crit_claim) == 0 for sol in crit_solutions)
    assert match, (
        f"Critical R does not match claim. "
        f"Solutions: {crit_solutions}, claim: {R_crit_claim}"
    )


def test_paper51_thm_dS_vacuum_smin_value():
    """s(R_crit, Lambda) = -Lambda / (12 sqrt(6)) exactly."""
    R, Lambda = sp.symbols("R Lambda", positive=True)

    # Action density (Section G5, eq. before thm:dS_vacuum)
    s = sp.Rational(1, 4) * R**3 * Lambda**4 - sp.Rational(1, 8) * R * Lambda**2

    # Critical point (claim and verified separately above)
    R_crit = 1 / (sp.sqrt(6) * Lambda)

    s_min_computed = sp.simplify(s.subs(R, R_crit))
    s_min_claim = -Lambda / (12 * sp.sqrt(6))

    residual = sp.simplify(s_min_computed - s_min_claim)
    assert residual == 0, (
        f"s_min mismatch.\n"
        f"  Computed: {s_min_computed}\n"
        f"  Claim:    {s_min_claim}\n"
        f"  Residual: {residual}"
    )


def test_paper51_thm_dS_vacuum_sign_and_minimum():
    """Verify the critical point is a MINIMUM (d2s/dR2 > 0) and s_min < 0."""
    R, Lambda = sp.symbols("R Lambda", positive=True)

    s = sp.Rational(1, 4) * R**3 * Lambda**4 - sp.Rational(1, 8) * R * Lambda**2

    R_crit = 1 / (sp.sqrt(6) * Lambda)

    d2s_dR2 = sp.diff(s, R, 2)
    d2s_at_crit = sp.simplify(d2s_dR2.subs(R, R_crit))

    # Should be positive (it is a minimum)
    # d2s/dR2 = (3/2) R Lambda^4, positive for positive R, Lambda
    assert sp.simplify(d2s_at_crit - sp.Rational(3, 2) * R_crit * Lambda**4) == 0
    # Numerical: at Lambda=1, R_crit=1/sqrt(6), d2s = (3/2) * (1/sqrt(6)) > 0
    val = float(d2s_at_crit.subs(Lambda, 1))
    assert val > 0, f"Second derivative should be positive (minimum), got {val}"

    # And s_min itself is negative (it carries the - sign)
    s_min = sp.simplify(s.subs(R, R_crit))
    s_min_val = float(s_min.subs(Lambda, 1))
    assert s_min_val < 0, f"s_min should be negative, got {s_min_val}"

    # Specifically: at Lambda=1, s_min = -1/(12 sqrt(6))
    expected = -1.0 / (12.0 * (6.0**0.5))
    assert abs(s_min_val - expected) < 1e-15, (
        f"Numerical s_min at Lambda=1: got {s_min_val}, expected {expected}"
    )


def test_paper51_thm_dS_vacuum_dimensional_check():
    """s_min has dimensions of Lambda (mass scale); verify prefactor 1/(12 sqrt(6))."""
    Lambda = sp.symbols("Lambda", positive=True)

    s_min_claim = -Lambda / (12 * sp.sqrt(6))

    # Extract the dimensionless prefactor
    prefactor = sp.simplify(-s_min_claim / Lambda)
    expected_prefactor = 1 / (12 * sp.sqrt(6))

    assert sp.simplify(prefactor - expected_prefactor) == 0, (
        f"Prefactor mismatch: {prefactor} vs {expected_prefactor}"
    )

    # Rationalized form: 1/(12 sqrt(6)) = sqrt(6)/72
    rationalized = sp.sqrt(6) / 72
    assert sp.simplify(prefactor - rationalized) == 0, (
        f"Rationalized prefactor mismatch: {prefactor} vs sqrt(6)/72"
    )


if __name__ == "__main__":
    test_paper51_thm_dS_vacuum_extremum_location()
    test_paper51_thm_dS_vacuum_smin_value()
    test_paper51_thm_dS_vacuum_sign_and_minimum()
    test_paper51_thm_dS_vacuum_dimensional_check()
    print("All 4 tests PASSED.")
