"""Guard: Paper 13's He partial-wave (Schwartz) coefficient, derived independently.

Claim (Paper 13; geovac/cusp_correction.py): the leading coefficient of the He
partial-wave increments Delta E_l ~ -A/(l+1/2)^4 is A = 0.074226 Ha (Bromley &
Mitroy, arXiv:physics/0601049), NOT (10/pi)<delta^3(r12)> = 0.339 Ha, the value
the module used before the group2 /qa remediation (v6.2.0).

This guard does not read the expected value from the module.  It derives the
coefficient from the repo's own 2D variational He solver
(level3_variational.solve_he_variational_2d), which never touches
cusp_correction.py:

  A_eff(l) = [E(l_max = l-1) - E(l_max = l)] * (l + 1/2)^4,   l = 3..7

falls 0.083 -> 0.078 as the (l+1/2)^-5 term dies away (the same approach from
above seen in Bromley-Mitroy Table I).  Extrapolating A_eff(l) linearly in
1/(l+1/2) from the last three points gives the leading coefficient; the module
value must agree with that extrapolation and lie below every A_eff(l).

Wrong answers this rejects: A = 0.339 (the retired (10/pi)<delta> relation);
A = A_eff(7) or A_eff(3) (an un-extrapolated increment).
Fire-tested: restoring `_SCHWARTZ_COEFF_HE = (10/pi)*_COALESCENCE_DENSITY_HE`
in geovac/cusp_correction.py makes it fail.

Basis: n_basis_R = 15, n_basis_alpha = 20, n_quad_alpha = 60 (a few seconds in
total).  The increments are basis-converged to 1e-5 in A_eff: n_basis_R = 20 /
n_quad = 80 and n_basis_alpha = 40 or 60 give the same A_eff to the printed
digits (measured 2026-11-04).
"""
from __future__ import annotations

import numpy as np
import pytest

from geovac.level3_variational import solve_he_variational_2d
from geovac.cusp_correction import schwartz_coefficient, cusp_correction_he


_L_VALUES = (3, 4, 5, 6, 7)


@pytest.fixture(scope="module")
def a_eff() -> dict:
    """A_eff(l) from the 2D variational partial-wave increments."""
    energies = {}
    for lm in range(_L_VALUES[0] - 1, _L_VALUES[-1] + 1):
        r = solve_he_variational_2d(
            Z=2.0, n_basis_R=15, n_basis_alpha=20, l_max=lm,
            alpha_R=2.0, n_quad_alpha=60,
        )
        energies[lm] = float(r["energies"][0])
    out = {}
    for l in _L_VALUES:
        dE = energies[l - 1] - energies[l]
        assert dE > 0.0, f"increment at l={l} must lower the energy, got {dE}"
        out[l] = dE * (l + 0.5) ** 4
    return out


def _extrapolated_A(a: dict) -> float:
    """Leading coefficient: fit A_eff(l) = A + B/(l+1/2) on l = 5, 6, 7."""
    ls = np.array([5, 6, 7], dtype=float)
    x = 1.0 / (ls + 0.5)
    y = np.array([a[5], a[6], a[7]])
    slope, intercept = np.polyfit(x, y, 1)
    return float(intercept)


def test_effective_coefficients_fall_toward_leading_value(a_eff):
    """A_eff(l) decreases monotonically, from ~0.083 at l=3 to ~0.078 at l=7."""
    vals = [a_eff[l] for l in _L_VALUES]
    assert all(vals[i] > vals[i + 1] for i in range(len(vals) - 1)), vals
    assert 0.080 < a_eff[3] < 0.086, a_eff[3]
    assert 0.076 < a_eff[7] < 0.080, a_eff[7]


def test_module_coefficient_matches_independent_extrapolation(a_eff):
    """The module's A agrees with the 2D-solver extrapolation to 2%."""
    A_mod = schwartz_coefficient(2.0)
    A_ext = _extrapolated_A(a_eff)
    assert 0.070 < A_ext < 0.078, f"extrapolated A = {A_ext:.5f}"
    assert abs(A_mod - A_ext) / A_ext < 0.02, (
        f"module A = {A_mod:.6f}, 2D-solver extrapolation = {A_ext:.6f}"
    )
    # Approach from above: the leading coefficient lies below every A_eff(l).
    assert A_mod < min(a_eff.values()), (A_mod, a_eff)


def test_tail_correction_matches_measured_increment(a_eff):
    """The l_max = 6 -> 7 piece of the module's tail matches the measured step.

    cusp_correction_he(l_max=6) - cusp_correction_he(l_max=7) is the module's
    predicted l = 7 increment (with sign); it must match the 2D solver's
    measured l = 7 increment to within the higher-order (l+1/2)^-5 remainder,
    which is ~5% at l = 7.
    """
    pred = cusp_correction_he(Z=2.0, l_max=7) - cusp_correction_he(Z=2.0, l_max=6)
    meas = a_eff[7] / 7.5 ** 4
    assert pred > 0.0
    assert abs(pred - meas) / meas < 0.08, (
        f"module l=7 increment {pred:.3e} vs 2D solver {meas:.3e}"
    )
