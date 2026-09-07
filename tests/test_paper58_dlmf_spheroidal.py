"""Paper 58, Sec. "The continuous side": the eta-equation IS the spheroidal
equation, and DLMF 30.3 therefore constrains the angular separation constant.

Dictionary (paper Eq. dlmf_dictionary), from matching
    d/deta[(1-eta^2) G'] + (-A + c^2 eta^2) G = 0          (Paper 11 eta-eq, m=0)
against DLMF 30.2.1
    d/dz[(1-z^2) w'] + (lambda + gamma^2 (1-z^2) - mu^2/(1-z^2)) w = 0
term by term:
    gamma^2 = -c^2,      lambda^m_n(gamma^2) = c^2 - A.

Claims backed here:
  1. DLMF 30.3.3  lambda^m_n(0) = n(n+1)        ->  A(c^2=0) = -n(n+1).
  2. DLMF 30.3.4  -1 < dlambda/dgamma^2 < 0     ->  0 < dA/d(c^2) < 1.
  3. Independent route: A = -obl_cv(m, n, c) (scipy), asserted by the
     molecular_sturmian docstring but previously untested.

What a wrong answer looks like (the reason each assertion is here):
  - a sign error on the c^2 term flips the slope negative       -> test 2 fires
  - a wrong Legendre diagonal breaks the united-atom limit      -> test 1 fires
  - a basis/normalisation error detunes us from scipy           -> test 3 fires
Claim 2 is the load-bearing one: it is the only assertion that constrains the
RATE of l-decompactification rather than an endpoint.
"""
import numpy as np
import pytest
from scipy.special import obl_cv

from geovac.molecular_sturmian import _angular_sep_const

# (m, n_sph) pairs; n = m + n_sph is the DLMF degree
MODES = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1)]


def _A(c2: float, m: int, n_sph: int, n_basis: int = 80) -> float:
    """Angular separation constant A at a given c^2."""
    return _angular_sep_const(m, n_sph, np.sqrt(max(c2, 0.0)), b=0.0,
                              n_basis=n_basis)


@pytest.mark.parametrize("m,n_sph", MODES)
def test_united_atom_limit_is_legendre(m: int, n_sph: int) -> None:
    """DLMF 30.3.3 via the dictionary: A(0) = -n(n+1) exactly."""
    n = m + n_sph
    assert _A(0.0, m, n_sph) == pytest.approx(-n * (n + 1), abs=1e-9)


@pytest.mark.parametrize("m,n_sph", MODES)
def test_dlmf_rate_bound_on_separation_constant(m: int, n_sph: int) -> None:
    """DLMF 30.3.4 via the dictionary: 0 < dA/d(c^2) < 1, strictly.

    This is the rate of l-decompactification. It may not stall (slope > 0)
    and may not outrun c^2 (slope < 1).
    """
    h = 1e-4
    slopes = []
    for c2 in np.linspace(0.25, 60.0, 12):
        d = (_A(c2 + h, m, n_sph) - _A(c2 - h, m, n_sph)) / (2.0 * h)
        slopes.append(d)
    lo, hi = min(slopes), max(slopes)
    assert lo > 0.0, f"slope not strictly positive: min={lo!r}"
    assert hi < 1.0, f"slope not strictly below 1: max={hi!r}"


@pytest.mark.parametrize("m,n_sph", MODES[:6])
def test_independent_route_scipy_oblate(m: int, n_sph: int) -> None:
    """A(c^2) = -obl_cv(m, n, c): second algorithmic route on the solver."""
    n = m + n_sph
    for c2 in (0.5, 2.0, 8.0, 20.0, 50.0):
        c = np.sqrt(c2)
        assert _A(c2, m, n_sph) == pytest.approx(-obl_cv(m, n, c), abs=1e-9)


@pytest.mark.parametrize("m", [0, 1, 2])
def test_dlmf_ordering_relation(m: int) -> None:
    """DLMF 30.3.1: lambda^m_m < lambda^m_{m+1} < ... at fixed m, every c^2.

    Wrong answer this rejects: an eigenvalue selector that picks from the
    wrong end of the spectrum, or crosses branches as c^2 grows, would
    break strict ordering.
    """
    for c2 in (0.0, 1.0, 8.0, 30.0, 60.0):
        lams = [c2 - _A(c2, m, ns) for ns in range(5)]
        assert all(b > a for a, b in zip(lams, lams[1:])), (m, c2, lams)


def test_separated_atom_limit_is_doubly_degenerate() -> None:
    """The OTHER end of the correlation diagram: adjacent levels pair up.

    As c^2 grows the spheroidal eigenvalues collapse into near-degenerate
    doublets -- the sigma_g/sigma_u pairs of the separated-atom limit, the
    R -> infinity end that DLMF 30.3.3 (A(0) = -n(n+1)) does not reach.
    At c^2 = 140 the within-doublet gap is ~1e-7 while the between-doublet
    gap is O(10): the spectrum has visibly doubled.
    """
    c2 = 140.0
    lams = [c2 - _A(c2, 0, ns) for ns in range(4)]
    within_01 = lams[1] - lams[0]
    within_23 = lams[3] - lams[2]
    between = lams[2] - lams[1]
    assert within_01 < 1e-6, within_01
    assert within_23 < 1.0, within_23
    assert between > 10.0, between
    assert between / within_01 > 1e6


def test_doublet_gap_decays_as_exp_minus_2c() -> None:
    """Within-doublet gap ~ c^2 exp(-2c): exchange-splitting class.

    Two independent constraints, because either alone is weak:
      (a) the LOCAL slope d ln(gap)/dc must approach -2 from above (a naive
          single-exponential fit returns ~-1.67 -- the c^2 prefactor
          contaminating the slope, which is why the raw fit is not asserted);
      (b) after dividing out exp(-2c) the residual must be a clean power law.
    Wrong answer this rejects: a power-law or Gaussian decay, either of
    which would leave the local slope wandering rather than tending to -2.
    """
    c2s = np.array([40.0, 60.0, 100.0, 140.0, 200.0, 300.0])
    c = np.sqrt(c2s)
    gap = np.array([(c2 - _A(c2, 0, 1)) - (c2 - _A(c2, 0, 0)) for c2 in c2s])
    assert np.all(np.diff(gap) < 0), gap

    slopes = np.diff(np.log(gap)) / np.diff(c)
    assert np.all(slopes > -2.0), slopes           # approaches -2 from above
    assert np.all(slopes < -1.6), slopes
    assert slopes[-1] < slopes[0], slopes          # monotonically toward -2

    resid = np.log(gap) + 2.0 * c                  # divide out exp(-2c)
    fit = np.polyfit(np.log(c), resid, 1)
    r2 = np.corrcoef(np.log(c), resid)[0, 1] ** 2
    assert 1.8 < fit[0] < 2.4, fit                 # prefactor ~ c^2
    assert r2 > 0.99, r2


def test_monotone_increasing_across_front_sweep() -> None:
    """A is strictly increasing in c^2 over the Paper 58 front range.

    A corollary of the rate bound, checked on the actual sweep rather than
    on the parametrized grid: R in [1, 16] a_0 at a representative bound
    energy, i.e. the range over which R* runs from 1.565 to 17.15.
    """
    E = -1.1
    c2s = [(-R ** 2 * E / 2.0) for R in (1.0, 2.0, 4.0, 8.0, 16.0)]
    vals = [_A(c2, 0, 0) for c2 in c2s]
    assert all(b > a for a, b in zip(vals, vals[1:])), vals
