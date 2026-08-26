"""Regression test for the Paper-60 molecular block-encoding 1-norm scaling.

Reproduces the paper's honest NEGATIVE resource result (``sec:manyelectron``): the
STANDARD block-encoding 1-norm ``lambda = sum|h| + sum|(pq|rs)|`` for two-electron H2 in
the Loewdin-orthonormalized two-center shared-scale Coulomb--Sturmian basis grows
POLYNOMIALLY, ``lambda ~ n_orb^2.2`` -- there is no sublinear behavior.  The atomic
isoenergetic-secular-matrix sublinearity (``eq:sublinear``) is a single-center property
that does not transfer to molecules.

Runs on a deliberately coarse integration grid so the O(N^4) grid-ERI sweep stays fast;
the measured exponent is stable against grid refinement (see the module defaults, which
reproduce the driver's finer grid).
"""
import numpy as np
import pytest

from geovac.sturmian_molecular_lambda import lambda_scaling

# Reduced grid: keeps the whole sweep well under ~2 min while the exponent stays ~2.2.
_GRID = dict(Lmax=14, nr=1000, nth=100, rmax=55.0)
_NMAX = (1, 2, 3, 4)  # n_orb = 2, 4, 6, 8


@pytest.mark.slow
def test_h2_molecular_lambda_is_polynomial() -> None:
    """lambda ~ n_orb^p with p ~ 2.2 (polynomial, NOT sublinear) and monotone growth."""
    n_orbs, lams, p = lambda_scaling(nmax_values=_NMAX, R=1.4, zeta=1.2, **_GRID)

    # n_orb ladder is the expected 2, 4, 6, 8.
    assert n_orbs == [2 * n for n in _NMAX]

    # Polynomial scaling: exponent lands in a band around the Paper-60 value 2.2.
    assert 1.9 < p < 2.6, f"exponent {p:.3f} outside the polynomial band [1.9, 2.6]"

    # Sanity: lambda strictly increases with basis size (polynomial, not
    # sublinear/flat/decreasing).  A sublinear-in-a-meaningful-sense object would not
    # more-than-double at each doubling of n_orb.
    lams_arr = np.asarray(lams)
    assert np.all(np.diff(lams_arr) > 0), f"lambda not monotone increasing: {lams}"
    assert p > 1.0, "exponent must exceed 1 for a genuinely super-linear (polynomial) law"
