"""Paper 29 Corollary `int_alg` backing -- the algebraic-integer dichotomy of
the GeoVac Hopf-graph Ihara zeta.

Cor `int_alg` makes a sharp distinction between two spectra:

  * the RECIPROCAL Ihara zeros = the eigenvalues of the Hashimoto edge operator
    T.  T is a 0/1 matrix, so its characteristic polynomial is MONIC with
    integer coefficients => every reciprocal zero is an algebraic INTEGER.

  * the Ihara zeros THEMSELVES are roots of the Bass polynomial
    det(I - s A + s^2 Q),  Q = diag(d_v - 1),  whose leading coefficient is
    det(Q) = prod_v (d_v - 1).  For any graph carrying a vertex of degree >= 3
    this is != 1, so the Bass polynomial is NON-monic => the Ihara zeros are
    algebraic over Q but NOT algebraic integers.

The paper's worked example is the factor 4 s^2 + 1 (Ihara zeros +- i/2, NOT
algebraic integers) versus its reciprocal s^2 + 4 (+- 2i, algebraic integers).

This file is the backing for the corollary's keystone PRECISION -- the
monic/non-monic split -- which the algebraicity-envelope test only covers at the
weaker "algebraic over Q" level.  Everything is recomputed from the framework
(the live Hashimoto operator + Bass determinant), not hard-coded.
"""
from __future__ import annotations

import numpy as np
import sympy as sp
import pytest

from geovac.ihara_zeta import hashimoto_matrix, ihara_zeta_bass


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _s3_adjacency(max_n: int) -> np.ndarray:
    """The S^3 Coulomb Hopf graph at the given cutoff, as a 0/1 adjacency."""
    from geovac.lattice import GeometricLattice
    A = GeometricLattice(max_n=max_n).adjacency.toarray()
    return (A != 0).astype(int)


def _k4_adjacency() -> np.ndarray:
    """Complete graph K4 (3-regular): small, non-trivial Ihara zeta, fast charpoly."""
    return (np.ones((4, 4), dtype=int) - np.eye(4, dtype=int))


def _is_algebraic_integer(alpha) -> bool:
    """alpha is an algebraic integer iff its monic minimal polynomial over Q has
    integer coefficients.  sympy.minimal_polynomial returns the PRIMITIVE INTEGER
    minimal polynomial (so its coeffs are always integers); the algebraic-integer
    discriminator is therefore that its LEADING coefficient is +-1 (monic).
    E.g. i/2 -> 4x^2 + 1 (LC 4, not an algebraic integer); 2i -> x^2 + 4 (LC 1)."""
    x = sp.symbols("x")
    mp = sp.minimal_polynomial(alpha, x, polys=True)
    return abs(int(mp.LC())) == 1


# ---------------------------------------------------------------------------
# the corollary's literal worked example: 4 s^2 + 1  vs  s^2 + 4
# ---------------------------------------------------------------------------
def test_worked_example_ihara_zero_not_integer_reciprocal_is():
    """+- i/2 (Ihara zeros, roots of the NON-monic 4 s^2 + 1) are NOT algebraic
    integers; +- 2i (their reciprocals, roots of the MONIC s^2 + 4) ARE."""
    s = sp.symbols("s")
    ihara_factor = sp.Poly(4 * s**2 + 1, s)          # the Ihara-zero side
    recip_factor = sp.Poly(s**2 + 4, s)              # the reciprocal side
    assert ihara_factor.LC() == 4                    # non-monic
    assert recip_factor.LC() == 1                    # monic
    # i/2 is a root of 4s^2+1 and is NOT an algebraic integer
    assert ihara_factor.eval(sp.I / 2) == 0
    assert not _is_algebraic_integer(sp.I / 2)
    # 2i is a root of s^2+4 and IS an algebraic integer
    assert recip_factor.eval(2 * sp.I) == 0
    assert _is_algebraic_integer(2 * sp.I)
    # they are exact reciprocals: (i/2) * (-2i) = 1
    assert sp.simplify((sp.I / 2) * (-2 * sp.I)) == 1


# ---------------------------------------------------------------------------
# reciprocal zeros are algebraic integers: Hashimoto T is 0/1 => monic Z charpoly
# ---------------------------------------------------------------------------
def test_hashimoto_charpoly_monic_integer_K4():
    """On K4 (fast): the Hashimoto T is 0/1 and its characteristic polynomial is
    monic with integer coefficients => the reciprocal Ihara zeros (eigenvalues of
    T) are algebraic integers.  Full charpoly computed, not asserted by theory."""
    T = hashimoto_matrix(_k4_adjacency())
    assert set(np.unique(T)).issubset({0, 1}), "Hashimoto T must be a 0/1 matrix"
    cp = sp.Matrix([[int(v) for v in row] for row in T]).charpoly()
    coeffs = cp.all_coeffs()
    assert coeffs[0] == 1, "Hashimoto charpoly must be monic"
    assert all(int(c) == c for c in coeffs), "Hashimoto charpoly must be integer"


def test_geovac_s3_hashimoto_is_integer_0_1():
    """The GeoVac S^3 graph (max_n=3, the Ramanujan cell): its Hashimoto T is a
    0/1 integer matrix -- so its charpoly is monic over Z and the reciprocal
    Ihara zeros are algebraic integers (the property is forced by integrality of
    T; the K4 case computes the full charpoly explicitly)."""
    T = hashimoto_matrix(_s3_adjacency(3))
    assert T.shape[0] == T.shape[1] and T.shape[0] > 0
    assert set(np.unique(T)).issubset({0, 1}), "S^3 Hashimoto T must be 0/1 integer"


# ---------------------------------------------------------------------------
# Ihara zeros themselves are NOT algebraic integers: Bass poly is non-monic
# ---------------------------------------------------------------------------
def _bass_det_leading_coeff(A: np.ndarray):
    """Leading coefficient (in s) of det(I - s A + s^2 Q), Q = diag(d_v - 1)."""
    s = sp.symbols("s")
    V = A.shape[0]
    d = A.sum(axis=1).astype(int)
    M = sp.eye(V) - s * sp.Matrix(A.tolist()) + s * s * sp.diag(*[int(x - 1) for x in d])
    return sp.Poly(sp.expand(M.det()), s).LC(), int(np.prod(d - 1))


def test_bass_polynomial_is_nonmonic_K4():
    """On K4 (3-regular, no pendant vertices): det(I - sA + s^2 Q) has leading
    coefficient det(Q) = prod(d_v - 1) = 2^4 = 16 != 1 => the Bass polynomial is
    NON-monic, so the Ihara zeros are algebraic over Q but NOT algebraic integers.
    Verified both as prod(d_v - 1) and as the symbolic determinant's leading coeff.
    (A graph with degree-1 vertices has det(Q)=0 and the top term vanishes; K4 is
    the clean witness for the non-monic half.)"""
    A = _k4_adjacency()
    lc, prod_q = _bass_det_leading_coeff(A)
    assert prod_q == 16, "K4: det(Q) = prod(d_v - 1) = 2^4 = 16"
    assert lc == prod_q, "Bass det leading coeff must equal det(Q)"
    assert lc != 1, "Bass polynomial is non-monic => Ihara zeros not algebraic integers"


def test_bass_and_hashimoto_are_the_two_sides():
    """Sanity: the Bass route (Ihara zeros) and the Hashimoto route (reciprocal
    zeros) are reciprocal descriptions of the SAME graph -- the dichotomy is not
    an artifact of two unrelated objects.  On K4 the Bass zeta^{-1} is a genuine
    degree>0 polynomial while T is the 0/1 edge operator of the same graph."""
    A = _k4_adjacency()
    zeta = ihara_zeta_bass(A, symbolic=True)          # 1 / [(1-s^2)^{r-1} det(...)]
    assert zeta.free_symbols  # a non-trivial rational function of s
    T = hashimoto_matrix(A)
    # T acts on directed edges: dim = 2E; K4 has E=6 => 12
    assert T.shape[0] == 2 * (int(A.sum()) // 2) == 12
