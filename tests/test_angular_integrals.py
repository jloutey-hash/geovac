"""
Tests for Wigner 3j symbols in geovac.angular_integrals.

Validates against known analytical values from standard references.
"""

import numpy as np
import pytest
from math import sqrt

from geovac.angular_integrals import wigner3j


class TestWigner3jKnownValues:
    """Test against analytically known 3j symbol values."""

    def test_110_000(self) -> None:
        """(1 1 0; 0 0 0) = -1/sqrt(3)."""
        val = wigner3j(1, 1, 0, 0, 0, 0)
        np.testing.assert_allclose(val, -1.0 / sqrt(3), atol=1e-14)

    def test_112_000(self) -> None:
        """(1 1 2; 0 0 0) = sqrt(2/15)."""
        val = wigner3j(1, 1, 2, 0, 0, 0)
        np.testing.assert_allclose(val, sqrt(2.0 / 15.0), atol=1e-14)

    def test_112_1n10(self) -> None:
        """(1 1 2; 1 -1 0) = 1/sqrt(30)."""
        val = wigner3j(1, 1, 2, 1, -1, 0)
        np.testing.assert_allclose(val, 1.0 / sqrt(30.0), atol=1e-14)

    def test_211_000(self) -> None:
        """(2 1 1; 0 0 0) = sqrt(2/15) — same as (1 1 2; 0 0 0) by symmetry."""
        val = wigner3j(2, 1, 1, 0, 0, 0)
        np.testing.assert_allclose(val, sqrt(2.0 / 15.0), atol=1e-14)

    def test_111_10n1(self) -> None:
        """(1 1 1; 1 0 -1) = -1/sqrt(6)."""
        val = wigner3j(1, 1, 1, 1, 0, -1)
        np.testing.assert_allclose(val, -1.0 / sqrt(6.0), atol=1e-14)


class TestWigner3jSelectionRules:
    """Test that selection rules return zero."""

    def test_m_sum_nonzero(self) -> None:
        """m1 + m2 + m3 != 0 gives zero."""
        assert wigner3j(1, 1, 1, 1, 1, 0) == 0.0

    def test_triangle_violated(self) -> None:
        """j3 > j1 + j2 gives zero."""
        assert wigner3j(1, 1, 3, 0, 0, 0) == 0.0

    def test_m_exceeds_j(self) -> None:
        """|m| > j gives zero."""
        assert wigner3j(1, 1, 0, 2, -2, 0) == 0.0

    def test_parity_odd(self) -> None:
        """(1 1 1; 0 0 0) = 0 because j1+j2+j3 = 3 is odd."""
        assert wigner3j(1, 1, 1, 0, 0, 0) == 0.0


class TestWigner3jSymmetries:
    """Test symmetry properties of 3j symbols."""

    def test_even_permutation(self) -> None:
        """Even permutation of columns: (j1 j2 j3) -> (j2 j3 j1)."""
        a = wigner3j(1, 2, 3, 0, 1, -1)
        b = wigner3j(2, 3, 1, 1, -1, 0)
        np.testing.assert_allclose(a, b, atol=1e-14)

    def test_odd_permutation_sign(self) -> None:
        """Odd permutation picks up (-1)^{j1+j2+j3}."""
        a = wigner3j(1, 2, 3, 0, 1, -1)
        b = wigner3j(2, 1, 3, 1, 0, -1)
        phase = (-1) ** (1 + 2 + 3)
        np.testing.assert_allclose(a, phase * b, atol=1e-14)

    def test_sign_flip_m(self) -> None:
        """Flipping all m: (j1 j2 j3; -m1 -m2 -m3) = (-1)^{j1+j2+j3} * (j1 j2 j3; m1 m2 m3)."""
        a = wigner3j(1, 1, 2, 1, -1, 0)
        b = wigner3j(1, 1, 2, -1, 1, 0)
        phase = (-1) ** (1 + 1 + 2)
        np.testing.assert_allclose(a, phase * b, atol=1e-14)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
