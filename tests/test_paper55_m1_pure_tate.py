"""Tests for Paper 55 Theorem thm:m1_pure_tate (M1 period-ring).

Paper 55 (papers/group3_foundations/paper_55_periods_of_geovac.tex),
Theorem at line 344-356, states:

    The M1 sub-mechanism on the discrete Camporesi-Higuchi S^3 spectral
    triple produces period values in the localised pure-Tate sub-ring

        M_1 \\subset Q[pi, pi^{-1}]

    of mixed-Tate periods over Q (Kontsevich-Zagier 2001), at depth 0
    and Tate weight in Z.

The witness table at Paper 55 lines 378-389 lists five M1 observables
across the corpus:

    W1: Vol(S^2)/4 = pi                  (Paper 25, Hopf-base)        pi^{+1}
    W2: 4/pi (asymptotic propinquity)    (Paper 38 L2 rate)           pi^{-1}
    W3: 4/pi (universal compact Lie)     (Paper 40 Thm 1)             pi^{-1}
    W4: 1/pi^2 (Pythagorean orthog.)     (Paper 43 sec 10.2)          pi^{-2}
    W5: 1/pi^2 (F-theorem cross-prod.)   (Paper 50 Thm 3.4)           pi^{-2}

W3 numerically equals W2 (4/pi) — same value, different source. The
sub-task requests verifying 3-5 distinct numerical values; the three
distinct values are W1=pi, W2=4/pi, W4=1/pi^2. This test verifies all
three at 80 dps via PSLQ against the basis {pi^n}, n = -4..4, and
confirms each is a rational linear combination of pi-powers.

In addition, this test cross-verifies the W5 source (Paper 50 Thm 3.4
F-theorem coefficients F_s = log(2)/8 - 3 zeta(3)/(16 pi^2) and
F_D = log(2)/4 + 3 zeta(3)/(8 pi^2); Paper 55 line 1175 claims the
F-theorem cross-product F_D - 2 F_s = 3 zeta(3)/(4 pi^2), and the
1/pi^2 prefactor is the M1 component of this M1 x M3 cross-product per
Paper 55 lines 1140-1180). The 1/pi^2 factor is the M1 witness.

Per CLAUDE.md Section 13.4a equation verification protocol, the test
uses both "Numerical cross-check" (PSLQ at 80 dps) and "Symbolic
identity" (sympy rational extraction) verification types.
"""
from __future__ import annotations

import sympy as sp
import mpmath as mp
import pytest


# ---------------------------------------------------------------------------
# Test parameters
# ---------------------------------------------------------------------------
WORKING_PRECISION_DPS = 80
PSLQ_MAX_COEFF = 10**10
TOLERANCE_DPS = 70  # require agreement to at least 70 dps after PSLQ


# ---------------------------------------------------------------------------
# Witness computation (high-precision values for the M1 observables)
# ---------------------------------------------------------------------------
def _compute_witnesses(dps: int):
    """Compute the three distinct M1 numerical witnesses at given precision.

    Returns a dict mapping witness label -> (value, expected_tate_weight,
    source_description).
    """
    mp.mp.dps = dps
    pi = mp.pi
    return {
        "W1_Vol_S2_over_4": (
            mp.mpf(4) * pi / mp.mpf(4),   # Vol(S^2) = 4 pi, divided by 4 = pi
            +1,
            "Paper 25 sec II.1 Hopf base measure (= pi)",
        ),
        "W2_propinquity_4_over_pi": (
            mp.mpf(4) / pi,                # Paper 38 L2 asymptotic rate
            -1,
            "Paper 38 L2 rate; Paper 40 universal Lie group constant",
        ),
        "W4_pythagorean_1_over_pi2": (
            mp.mpf(1) / pi**2,             # Paper 43 sec 10.2 prefactor
            -2,
            "Paper 43 sec 10.2 Pythagorean orthogonality prefactor; "
            "Paper 50 F-theorem cross-product factor",
        ),
    }


def _pslq_against_pi_powers(value, dps: int):
    """Run PSLQ on [value, pi^-4, pi^-3, ..., pi^4] at given precision.

    Returns the integer relation found, or None if PSLQ fails.
    """
    mp.mp.dps = dps
    pi = mp.pi
    basis = [value] + [pi**k for k in range(-4, 5)]
    return mp.pslq(basis, maxcoeff=PSLQ_MAX_COEFF)


def _interpret_pslq_relation(relation):
    """Decode the PSLQ relation as a dict {n -> q} meaning value = sum q_n pi^n.

    Convention: relation[0] is coefficient on value; relation[1..9] are
    coefficients on pi^-4, pi^-3, ..., pi^4. Then value satisfies
    relation[0] * value + sum_n relation[k] * pi^{k - 5} = 0, i.e.
    value = -(1/relation[0]) * sum_k relation[k] * pi^{k - 5}.
    """
    if relation is None:
        return None
    if relation[0] == 0:
        return None
    val_coeff = relation[0]
    powers = {}
    for idx in range(1, 10):
        n = idx - 5  # idx=1 -> n=-4, idx=9 -> n=+4
        coeff = relation[idx]
        if coeff != 0:
            powers[n] = sp.Rational(-coeff, val_coeff)
    return powers


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
class TestM1WitnessesArePiPowers:
    """PSLQ each M1 witness at 80 dps and confirm it = rational sum of pi^n."""

    @pytest.mark.parametrize("label", [
        "W1_Vol_S2_over_4",
        "W2_propinquity_4_over_pi",
        "W4_pythagorean_1_over_pi2",
    ])
    def test_witness_is_rational_pi_combination(self, label):
        witnesses = _compute_witnesses(WORKING_PRECISION_DPS)
        value, expected_weight, source = witnesses[label]

        relation = _pslq_against_pi_powers(value, WORKING_PRECISION_DPS)
        assert relation is not None, (
            f"PSLQ failed to identify {label} ({source}) as a rational "
            f"combination of pi^n, n in [-4, 4]"
        )

        powers = _interpret_pslq_relation(relation)
        assert powers is not None, (
            f"PSLQ returned a relation with zero coefficient on the value "
            f"itself for {label}"
        )

        # Verify the identification by reconstructing the value symbolically.
        mp.mp.dps = WORKING_PRECISION_DPS
        pi = mp.pi
        reconstructed = sum(mp.mpf(q.p) / mp.mpf(q.q) * pi**n
                            for n, q in powers.items())
        diff = abs(reconstructed - value)
        tolerance = mp.mpf(10)**(-TOLERANCE_DPS)
        assert diff < tolerance, (
            f"{label}: PSLQ relation reconstructs to {reconstructed}, "
            f"but witness value is {value} (diff = {diff}, tol = {tolerance})"
        )

        # Verify the Tate weight (highest |n| with nonzero coefficient).
        active_powers = list(powers.keys())
        # For pure-pi-power witnesses, the relation should be a single
        # power; the Tate weight is that power.
        nonzero_powers = [n for n in active_powers]
        assert len(nonzero_powers) >= 1, (
            f"{label}: PSLQ found no nonzero pi^n terms (relation: {relation})"
        )

        # Witnesses W1/W2/W4 are pure single-power values; the Tate weight
        # in the paper is the exponent n.
        # Check that the highest |n| matches expected_weight.
        max_n_by_abs = max(nonzero_powers, key=lambda n: abs(n))
        # Allow either single-power or matching weight; the paper's
        # Tate weight equals the unique pi exponent for pure-power witnesses.
        assert max_n_by_abs == expected_weight, (
            f"{label}: PSLQ identified pi^{max_n_by_abs} with rational "
            f"coefficient {powers[max_n_by_abs]}, but Paper 55 lists Tate "
            f"weight = {expected_weight}"
        )


class TestM1WitnessesAreInQpiInversePi:
    """Verify each witness lies in Q[pi, pi^-1] (the M1 period ring)."""

    def test_W1_vol_S2_over_4_equals_pi(self):
        """W1 = Vol(S^2)/4 = pi, an element of Q[pi, pi^-1] at Tate weight +1."""
        mp.mp.dps = WORKING_PRECISION_DPS
        pi = mp.pi
        w1 = mp.mpf(4) * pi / mp.mpf(4)
        diff = abs(w1 - pi)
        # W1 IS pi by definition; this is a tautology-level identity.
        assert diff == 0 or diff < mp.mpf(10)**(-(WORKING_PRECISION_DPS - 5)), (
            f"W1 = Vol(S^2)/4 should equal pi exactly; got diff {diff}"
        )

    def test_W2_4_over_pi_in_QpiInvPi(self):
        """W2 = 4/pi: PSLQ gives 1*W2 - 4*pi^{-1} = 0 (relation found)."""
        mp.mp.dps = WORKING_PRECISION_DPS
        w2 = mp.mpf(4) / mp.pi
        rel = _pslq_against_pi_powers(w2, WORKING_PRECISION_DPS)
        powers = _interpret_pslq_relation(rel)
        # Expected: powers == {-1: 4}, i.e. W2 = 4 * pi^-1
        assert powers == {-1: sp.Integer(4)}, (
            f"W2 = 4/pi should PSLQ-identify as 4 * pi^-1; got {powers}"
        )

    def test_W4_1_over_pi2_in_QpiInvPi(self):
        """W4 = 1/pi^2: PSLQ gives 1*W4 - 1*pi^{-2} = 0 (relation found)."""
        mp.mp.dps = WORKING_PRECISION_DPS
        w4 = mp.mpf(1) / mp.pi**2
        rel = _pslq_against_pi_powers(w4, WORKING_PRECISION_DPS)
        powers = _interpret_pslq_relation(rel)
        # Expected: powers == {-2: 1}, i.e. W4 = 1 * pi^-2
        assert powers == {-2: sp.Integer(1)}, (
            f"W4 = 1/pi^2 should PSLQ-identify as 1 * pi^-2; got {powers}"
        )


class TestM1Paper50FTheoremCrossProductHasPi2Prefactor:
    """Cross-check W5: Paper 50 Theorem 3.4 F-theorem cross-product.

    Paper 50 Thm 3.4 states (on the truncated S^3 spectral triple):
        F_s = log(2)/8 - 3*zeta(3)/(16*pi^2)
        F_D = log(2)/4 + 3*zeta(3)/(8*pi^2)
    Paper 55 line 1175 (Theorem on F-theorem dual-basis projection):
        F_D - 2*F_s = 3*zeta(3)/(4*pi^2)
    The 1/pi^2 factor in F_D - 2*F_s is the M1 component of the
    M1 x M3 cross-product (Paper 55 line 1158); zeta(3) is the M3
    component (Tate weight +3, depth 1).

    This test verifies the cross-product algebraically (at 80 dps) and
    extracts the M1 piece (= 1/pi^2) via PSLQ.
    """

    def test_F_theorem_cross_product_closed_form(self):
        """F_D - 2*F_s = 3 zeta(3)/(4 pi^2) bit-exact (Paper 50 Thm 3.4)."""
        mp.mp.dps = WORKING_PRECISION_DPS
        pi = mp.pi
        z3 = mp.zeta(3)
        l2 = mp.log(2)

        F_s = l2 / mp.mpf(8) - 3 * z3 / (16 * pi**2)
        F_D = l2 / mp.mpf(4) + 3 * z3 / (8 * pi**2)

        cross = F_D - 2 * F_s
        expected = 3 * z3 / (4 * pi**2)

        diff = abs(cross - expected)
        tol = mp.mpf(10)**(-TOLERANCE_DPS)
        assert diff < tol, (
            f"F_D - 2*F_s closed form mismatch: cross = {cross}, "
            f"expected = {expected}, diff = {diff}"
        )

    def test_F_theorem_cross_product_M1_factor_is_pi_minus_2(self):
        """The M1 factor in (F_D - 2 F_s) / (3 zeta(3)/4) is exactly pi^-2."""
        mp.mp.dps = WORKING_PRECISION_DPS
        pi = mp.pi
        z3 = mp.zeta(3)
        l2 = mp.log(2)

        F_s = l2 / mp.mpf(8) - 3 * z3 / (16 * pi**2)
        F_D = l2 / mp.mpf(4) + 3 * z3 / (8 * pi**2)
        cross = F_D - 2 * F_s

        # Strip the M3 component (3 zeta(3)/4) to isolate the M1 factor.
        m1_factor = cross / (3 * z3 / mp.mpf(4))

        # m1_factor should be exactly 1/pi^2.
        target = mp.mpf(1) / pi**2
        diff = abs(m1_factor - target)
        tol = mp.mpf(10)**(-TOLERANCE_DPS)
        assert diff < tol, (
            f"F-theorem cross-product M1 factor should be 1/pi^2; "
            f"got {m1_factor}, diff = {diff}"
        )

        # PSLQ-confirm m1_factor = 1*pi^-2 in the basis pi^-4..pi^4.
        rel = _pslq_against_pi_powers(m1_factor, WORKING_PRECISION_DPS)
        powers = _interpret_pslq_relation(rel)
        assert powers == {-2: sp.Integer(1)}, (
            f"F-theorem M1 factor should PSLQ as 1 * pi^-2; got {powers}"
        )


class TestM1AggregateClosure:
    """All five witnesses' M1 factors collected as a set of rationals-times-pi-powers."""

    def test_all_witnesses_in_QpiInvPi(self):
        """All five Paper 55 M1 witnesses are in Q[pi, pi^-1] with |Tate| <= 2."""
        mp.mp.dps = WORKING_PRECISION_DPS
        pi = mp.pi

        # Witnesses with their explicit closed forms.
        witnesses = {
            "W1_Vol(S^2)/4": (mp.mpf(4) * pi / mp.mpf(4), {1: sp.Integer(1)}),
            "W2_propinquity_4_over_pi": (mp.mpf(4) / pi, {-1: sp.Integer(4)}),
            "W3_universal_4_over_pi": (mp.mpf(4) / pi, {-1: sp.Integer(4)}),
            "W4_pythagorean_1_over_pi2": (mp.mpf(1) / pi**2, {-2: sp.Integer(1)}),
            "W5_F_theorem_M1_factor": (mp.mpf(1) / pi**2, {-2: sp.Integer(1)}),
        }

        for label, (value, expected_powers) in witnesses.items():
            rel = _pslq_against_pi_powers(value, WORKING_PRECISION_DPS)
            powers = _interpret_pslq_relation(rel)
            assert powers == expected_powers, (
                f"{label}: PSLQ expected {expected_powers}, got {powers}"
            )
            # And confirm Tate weight |n| <= 2 as Paper 55 line 393 states.
            for n in powers:
                assert abs(n) <= 2, (
                    f"{label}: Paper 55 line 393 states |n| <= 2 on empirical "
                    f"witnesses, but found pi^{n}"
                )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
