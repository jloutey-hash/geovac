"""Paper 60 -- the SCALE LOCK between the metric-free isoenergetic posing and the
variational CI over the same Goscinskian span.

Four claims, measured 2026-09-08 (drivers ``debug/p60_*.py``; the assembly they
and this file share now lives in ``geovac/sturmian_variational.py``):

  C1  the one-body ``1/r`` metric is EXACTLY diagonal,
      ``W_{mu nu} = <Phi_mu| sum_j 1/r_j |Phi_nu> = R_nu delta_{mu nu}``.
      This is the lemma that makes the posing metric-free.
  C2  SCALE LOCK: the isoenergetic energy equals the variational energy of the
      SAME span evaluated at the locked scale ``lambda = p_kappa``.
  C3  the span is NOT the limitation -- ``min_lambda E_var < E_iso`` strictly, a
      POSITIVE posing cost, exceeding 3 mHa by ``K >= 78`` (s-only).
  C4  the posing cost is STATE dependent -- at ``n_max = 14`` the ground-state
      cost exceeds 4x the ``2^1S`` cost.

Mechanism (why C2 is a theorem and not a coincidence), all three legs asserted
below rather than assumed.  With ``pk_ref = 1`` the configuration charge obeys
``Q_nu R_nu = 1``, so C1 gives ``W[:, nu] Q_nu = I`` and hence ``T = -S/2 + I``.
Then

    H(lam) C = E S C   with   E = -lam^2 / 2
      <=>  ( lam^2 I - lam Z diag(R) + lam G ) C = 0
      <=>  ( Z diag(R) - G ) C = lam C
      <=>  M C = lam C ,

i.e. the variational problem at scale ``lam`` has ``E = -lam^2/2`` exactly when
``lam`` is an eigenvalue of the isoenergetic secular matrix ``M``.  C2 is the
numerical statement of that equivalence; C1 is its load-bearing lemma.

ROUTE INDEPENDENCE (see ``test_c2_scale_lock_routes_are_independent``).  The two
sides are computed by disjoint assemblies:

  route A (E_iso)   geovac.sturmian_secular.solve -> build_M = diag(Z R_nu)
                    + T' -> eigh -> E = -p^2/2.  Never forms S or W.
  route B (E_var)   geovac.sturmian_variational.build -> (S, T, W, G) ->
                    var_energy solves the whitened generalized eigenproblem.
                    Never forms M.

They share only the primitive ``SS.repulsion_terms``.  The test below asserts the
independence operationally (route B's answer moves when ANY of its four matrix
arguments is perturbed) so that a refactor collapsing them into one call fails.

Grid: the Paper-60 box rule ``set_grid(max(80, 5 n_max^2), 24000, "grade", 2.0)``
and ``Z = 2``.  ``set_grid`` mutates ``geovac.sturmian_secular`` module globals; a
module-scoped autouse fixture restores them so the rest of the suite is unaffected.

Both routes import from ``geovac/`` only: the route-B assembly was promoted from
the transient ``debug/p60_*.py`` drivers into ``geovac.sturmian_variational`` on
2026-09-08 (those drivers are now thin shims re-exporting it), so the four guards
below do not die silently when ``debug/`` is pruned.
"""
import math
import os
import sys

import numpy as np
import pytest
from scipy.optimize import minimize_scalar

# Ensure project root on path (mirrors tests/conftest.py).
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import geovac.sturmian_secular as SS                                  # noqa: E402
import geovac.sturmian_variational as SV                              # noqa: E402
from geovac.sturmian_variational import (                             # noqa: E402
    build, var_energy, var_levels)

Z = 2.0
NPTS = 24000

# Textbook single-exponent (Kellner) helium: E = -(Z - 5/16)^2, exponent Z - 5/16.
KELLNER_EXPONENT = Z - 5.0 / 16.0
KELLNER_ENERGY = -KELLNER_EXPONENT ** 2          # -2.84765625


# --------------------------------------------------------------------------
# grid hygiene: set_grid patches module globals of geovac.sturmian_secular
# --------------------------------------------------------------------------
@pytest.fixture(scope="module", autouse=True)
def _restore_secular_grid():
    saved = dict(r=SS.r, dr=SS.dr, r2=SS.r2, R_MAX=SS.R_MAX, N_GRID=SS.N_GRID,
                 fwd=SS._ctrap_fwd, rev=SS._ctrap_rev)
    yield
    SS.r, SS.dr, SS.r2 = saved['r'], saved['dr'], saved['r2']
    SS.R_MAX, SS.N_GRID = saved['R_MAX'], saved['N_GRID']
    SS._ctrap_fwd, SS._ctrap_rev = saved['fwd'], saved['rev']
    SS._GAUNT_CACHE.clear()
    SS.reset_caches()


# --------------------------------------------------------------------------
# one build per (n_max, l_max); both routes are evaluated under the SAME grid
# --------------------------------------------------------------------------
class Case:
    """Both routes for one ``(n_max, l_max)`` span, built under one grid setting."""

    def __init__(self, nmax: int, lmax: int) -> None:
        SV.set_grid(max(80.0, 5.0 * nmax * nmax), NPTS, "grade", 2.0)
        self.nmax, self.lmax = nmax, lmax
        self.tuples = SV.family(nmax, lmax)
        # route B
        self.S, self.T, self.W, self.G, self.K, self.asym = build(nmax, lmax)
        # route A (independent assembly; same grid)
        cfgs = SS.build_configs(self.tuples)
        self.Rnu = np.array([c.Rnu for c in cfgs])
        self.Q = np.array([c.Q for c in cfgs])
        self.M = SS.build_M(cfgs, Z=Z)
        self.p = np.sort(np.linalg.eigvalsh(self.M))[::-1]   # descending roots
        self.E_iso = -self.p[0] ** 2 / 2

    def E_iso_k(self, k: int) -> float:
        return -self.p[k] ** 2 / 2

    def levels(self, lam: float, nlev: int, tol: float = 1e-10) -> np.ndarray:
        """The ``nlev`` lowest roots of ``H(lam) C = E S C`` (whitened)."""
        return var_levels(self.S, self.T, self.W, self.G, Z, lam, nlev, tol)


_CASES: dict = {}


def case(nmax: int, lmax: int) -> Case:
    key = (nmax, lmax)
    if key not in _CASES:
        _CASES[key] = Case(nmax, lmax)
    else:
        # re-install this case's grid; a later Case() may have moved it.
        SV.set_grid(max(80.0, 5.0 * nmax * nmax), NPTS, "grade", 2.0)
    return _CASES[key]


def level_energy(c: "Case", lam: float, level: int) -> float:
    """Root ``level`` of the variational problem at scale ``lam``."""
    if level == 0:
        return var_energy(c.S, c.T, c.W, c.G, Z, lam)
    return float(c.levels(lam, level + 1)[level])


def min_over_scale(c: "Case", level: int = 0, lo: float = 0.5, hi: float = 40.0,
                   npts: int = 400) -> tuple:
    """GLOBAL minimum of ``E_var(lambda)`` for root ``level``: grid scan, then a
    local refinement bracketed by the neighbouring grid points.

    A bare bounded optimizer over the whole window is NOT usable here: at
    n_max = 4, ``minimize_scalar(..., bounds=(0.3, 40), method='bounded')``
    converges to a LOCAL minimum at lambda = 4.663 and reports
    ``E = -2.8721637 > E_iso``, i.e. a NEGATIVE posing cost, which the
    variational bound forbids.  The grid finds the true optimum at lambda = 2.62
    (``E = -2.8756340``).  The refinement is bracketed by the grid-selected
    basin, so it cannot wander back out to a spurious local minimum -- and it is
    what makes the K = 1 cost vanish to 1e-12 rather than to the grid spacing.
    """
    grid = np.linspace(lo, hi, npts)
    vals = np.array([level_energy(c, L, level) for L in grid])
    i = int(vals.argmin())
    a, b = grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]
    res = minimize_scalar(lambda L: level_energy(c, L, level),
                          bounds=(a, b), method="bounded",
                          options=dict(xatol=1e-10))
    if res.fun < vals[i]:
        return float(res.fun), float(res.x)
    return float(vals[i]), float(grid[i])


# ==========================================================================
# pipeline unit test -- absolute calibration of BOTH routes at K = 1
# ==========================================================================
def test_pipeline_single_config_reproduces_kellner_helium():
    """n_max = 1, l_max = 0 (K = 1): both routes must give ``-(Z - 5/16)^2``.

    Not evidence for C2 -- at K = 1 both routes are 1x1 and the scale lock is
    trivial.  This pins the ABSOLUTE calibration of the pipeline (Slater
    integrals, normalization, Z handling, box rule): a broken radial integral or
    a mis-scaled charge moves this number in the 3rd decimal, not the 8th.
    """
    c = case(1, 0)
    assert c.K == 1
    assert SS.R_MAX >= 60.0
    assert round(c.E_iso, 6) == round(KELLNER_ENERGY, 6) == -2.847656

    e_var, lam = min_over_scale(c)
    assert round(e_var, 6) == -2.847656

    # the two postings COINCIDE at K = 1 -> the posing cost vanishes
    cost = c.E_iso - e_var
    assert abs(cost) < 1e-9, f"K=1 posing cost should vanish, got {cost:.3e} Ha"

    # and the optimal scale reproduces the Kellner exponent: lam* Q_nu = Z - 5/16
    assert abs(lam * c.Q[0] - KELLNER_EXPONENT) < 1e-6
    # locked scale is p_kappa = sqrt(-2 E), NOT the orbital exponent itself
    assert abs(lam - math.sqrt(-2 * c.E_iso)) < 1e-6


# ==========================================================================
# C1 -- the one-body 1/r metric is exactly diagonal
# ==========================================================================
@pytest.mark.parametrize("nmax,lmax,tol", [
    (6, 0, 1e-9),
    (10, 0, 1e-9),     # the measured point: off-diag 3.9e-11, |diag - R| 4.9e-11
    (3, 1, 1e-8),
    (4, 3, 1e-8),      # graded-grid error grows with l; still 7 OoM below the L2 scale
])
def test_c1_one_body_metric_is_exactly_diagonal(nmax, lmax, tol):
    """``W_{mu nu} = R_nu delta_{mu nu}``.

    Wrong answers this excludes: (a) a metric that is merely SMALL off the
    diagonal rather than zero -- the same builder's L2 overlap S is asserted to
    carry O(0.1) off-diagonals, so "everything this builder returns is tiny" is
    ruled out; (b) a diagonal equal to Z R_nu, R_nu^2, or Q_nu instead of R_nu.
    """
    c = case(nmax, lmax)
    off_W = np.abs(c.W - np.diag(np.diag(c.W))).max()
    dev = np.abs(np.diag(c.W) - c.Rnu).max()
    assert off_W < tol, f"W off-diagonal {off_W:.3e} (K={c.K})"
    assert dev < tol, f"max|diag(W) - R_nu| = {dev:.3e}"

    # ANTI-TAUTOLOGY: the same builder's L2 metric is emphatically NOT diagonal.
    off_S = np.abs(c.S - np.diag(np.diag(c.S))).max()
    assert off_S > 0.1, (
        f"L2 overlap off-diagonal is only {off_S:.3e}; the C1 assertion above "
        "would then be measuring 'the builder returns near-zero matrices', not "
        "potential-weighted orthogonality")
    # and the diagonal itself is nowhere near the tolerance: R_nu runs from
    # sqrt(2) (the 1s^2 config) down to sqrt(2)/n_max, all >> tol.
    assert np.abs(np.diag(c.W)).min() > 1e4 * tol
    assert np.abs(np.diag(c.W)).max() > 1.0

    # the wrong diagonals the assertion excludes are all far away
    assert np.abs(np.diag(c.W) - Z * c.Rnu).max() > 0.5
    assert np.abs(np.diag(c.W) - c.Rnu ** 2).max() > 0.1
    assert np.abs(np.diag(c.W) - c.Q).max() > 0.1


# ==========================================================================
# C2 -- the scale-lock identity
# ==========================================================================
@pytest.mark.parametrize("nmax,lmax", [(6, 0), (4, 3)])
def test_c2_scale_lock_identity(nmax, lmax):
    """``E_var(lambda = p_kappa) == E_iso`` to ~1e-9, at l_max = 0 and l_max = 3.

    Wrong answers this excludes: any OTHER locked scale (lam = 1, lam = Z,
    lam = the Kellner exponent), and a lambda-independent E_var.  The sharpness
    assertion below makes the second explicit: a 10% move off p_kappa costs
    > 1 mHa, so the equality at p_kappa is a lock, not a plateau.
    """
    c = case(nmax, lmax)
    pk = math.sqrt(-2 * c.E_iso)
    e_lock = var_energy(c.S, c.T, c.W, c.G, Z, pk)
    assert abs(e_lock - c.E_iso) < 1e-9, (
        f"scale lock broken: E_var(p_kappa) - E_iso = {e_lock - c.E_iso:.3e}")

    # SHARPNESS: E_var is strongly lambda-dependent, so the lock is informative.
    # NOT-THE-ENERGY-SHELL: E_var must be a genuine variational energy, not the
    # parametrization -lam^2/2 (which satisfies the identity at p_kappa for free
    # and is what a collapsed route would return).  Measured gap ~0.5-0.6 Ha.
    for f in (0.9, 1.1):
        lam = f * pk
        e_off = var_energy(c.S, c.T, c.W, c.G, Z, lam)
        assert abs(e_off - c.E_iso) > 1e-3, (
            f"E_var is flat in lambda near p_kappa (f={f}); the identity would "
            "then hold at any scale and assert nothing")
        assert abs(e_off - (-lam ** 2 / 2)) > 0.1, (
            f"E_var(lam) == -lam^2/2 at f={f}: route B is returning the energy "
            "shell rather than solving a variational problem, which makes the "
            "lock at p_kappa an identity instead of a result")

    # ... and the other candidate scales genuinely fail the identity
    for bad in (1.0, Z, KELLNER_EXPONENT):
        e_bad = var_energy(c.S, c.T, c.W, c.G, Z, bad)
        assert abs(e_bad - c.E_iso) > 1e-6, f"lambda={bad} also satisfies the lock"


@pytest.mark.parametrize("nmax,lmax", [(6, 0), (4, 3)])
def test_c2_scale_lock_mechanism(nmax, lmax):
    """The three algebraic legs that MAKE C2 a theorem, asserted separately.

    If a refactor breaks one of them the lock degenerates into a coincidence;
    this test says which leg went.
    """
    c = case(nmax, lmax)
    # leg 1: reference posing pk_ref = 1  =>  Q_nu R_nu = 1
    assert np.abs(c.Q * c.Rnu - 1.0).max() < 1e-12
    # leg 2 (needs C1): W[:, nu] Q_nu = I  =>  T = -S/2 + I
    assert np.abs(c.T - (-0.5 * c.S + np.eye(c.K))).max() < 1e-8
    # leg 3: route B's ingredients reassemble route A's secular matrix
    M_from_B = Z * np.diag(np.diag(c.W)) - c.G
    assert np.abs(c.M - M_from_B).max() < 1e-8


def test_c2_scale_lock_routes_are_independent():
    """The lock is not a tautology of a shared code path.

    Route A (``sturmian_secular.solve`` -> ``build_M`` -> ``eigh``) never forms S
    or W; route B (``sturmian_variational.build`` -> ``var_energy``) never forms
    M.  Asserted operationally: route B's answer must MOVE when any one of its
    four matrix arguments is perturbed.  A refactor that made ``var_energy``
    delegate to ``solve`` (or memoize E_iso) would leave the perturbed answers
    unchanged and fail here.
    """
    c = case(6, 0)
    pk = math.sqrt(-2 * c.E_iso)
    base = var_energy(c.S, c.T, c.W, c.G, Z, pk)
    assert abs(base - c.E_iso) < 1e-9

    eps = 1e-4
    mats = dict(S=c.S, T=c.T, W=c.W, G=c.G)
    for name in mats:
        pert = {k: v.copy() for k, v in mats.items()}
        pert[name][0, 0] += eps
        pert[name][1, 1] -= eps
        e = var_energy(pert['S'], pert['T'], pert['W'], pert['G'], Z, pk)
        assert abs(e - base) > 1e-8, (
            f"var_energy ignored its `{name}` argument -- the routes have "
            "collapsed into one code path and the scale-lock assertion is vacuous")

    # route A does not consume S/W/G at all: its matrix is diag(Z R_nu) + T'
    Tprime = SS.build_Tprime(SS.build_configs(c.tuples))
    assert np.abs(c.M - (Tprime + np.diag(Z * c.Rnu))).max() < 1e-12


# ==========================================================================
# C3 -- the span is not the limitation: a strictly POSITIVE posing cost
# ==========================================================================
@pytest.mark.parametrize("nmax", [4, 6, 8])
def test_c3_posing_cost_is_strictly_positive(nmax):
    """``min_lambda E_var < E_iso``, i.e. the variational bound is respected and
    the isoenergetic posing pays a strictly positive price for its scale lock.

    n_max = 4 is the regression case: a bare bounded optimizer reports
    E = -2.8721637 there (a local minimum at lambda = 4.663), i.e. a posing cost
    of -1.058 mHa, which VIOLATES the variational bound.  The assertion
    ``e_var <= E_iso`` is exactly what that wrong answer fails.
    """
    c = case(nmax, 0)
    e_var, lam = min_over_scale(c)
    assert e_var <= c.E_iso + 1e-12, (
        f"variational bound violated: min_lambda E_var = {e_var:.9f} > "
        f"E_iso = {c.E_iso:.9f} -- the scale scan found only a LOCAL minimum")
    cost_mha = (c.E_iso - e_var) * 1000
    assert cost_mha > 1.0, f"posing cost only {cost_mha:.4f} mHa at n_max={nmax}"
    assert 0.5 < lam < 40.0


@pytest.mark.slow
def test_c3_posing_cost_exceeds_3mha_at_K78():
    """C3 headline: the cost exceeds 3 mHa by K = 78 (s-only, n_max = 12).

    Wrong answer excluded: "the span is the limitation", which predicts a posing
    cost of ~0 -- i.e. that the isoenergetic answer already saturates its own
    span.  Measured 4.125 mHa; pinned to +-0.02 (the value is stable to
    1e-3 mHa across lambda grids of 200..2000 points and radial grids of
    12k..40k points).
    """
    c = case(12, 0)
    assert c.K == 78
    e_var, _ = min_over_scale(c)
    assert e_var <= c.E_iso + 1e-12
    cost_mha = (c.E_iso - e_var) * 1000
    assert cost_mha > 3.0, f"posing cost only {cost_mha:.4f} mHa at K=78"
    assert abs(cost_mha - 4.125) < 0.02, f"posing cost drifted to {cost_mha:.4f} mHa"


# ==========================================================================
# C4 -- the posing cost is state dependent
# ==========================================================================
@pytest.mark.slow
def test_c4_posing_cost_is_state_dependent():
    """n_max = 14, l_max = 0 (K = 105): ground-state posing cost > 4x the 2^1S cost.

    Measured 4.212 mHa (ground) vs 0.983 mHa (2^1S), ratio 4.287 -- stable to
    ~1e-3 across lambda-grid and radial-grid refinement.

    K-trend, measured 2026-09-08 (gnd / 2^1S / ratio / gap, mHa):

        K= 36   3.524  0.681   5.174   2.843
        K= 55   3.943  0.771   5.114   3.172
        K= 78   4.125  0.885   4.660   3.240
        K=105   4.212  0.983   4.287   3.230

    The RATIO narrows monotonically with K; what grows is the absolute GAP (and
    it flattens by K=105).  The 4x assertion below is therefore pinned at
    K = 105 and must NOT be assumed to hold a fortiori at larger K.

    Wrong answers excluded: (a) a state-INdependent posing cost (ratio ~ 1),
    which is what "the posing cost is just a basis-size artifact" predicts;
    (b) mis-identified roots -- E_iso(k) comes from the k-th LARGEST eigenvalue
    of M and both E_iso values are pinned below, so taking the k-th smallest
    eigenvalue, or the wrong variational root, fails here.
    """
    c = case(14, 0)
    assert c.K == 105

    # state identification is pinned, so the ratio cannot come out right by accident
    assert abs(c.E_iso_k(0) - (-2.8745946)) < 1e-5
    assert abs(c.E_iso_k(1) - (-2.1429339)) < 1e-5

    # the scale lock holds root by root: level k of H(p_k) is E_iso(k)
    for k in (0, 1):
        lk = float(c.levels(math.sqrt(-2 * c.E_iso_k(k)), k + 1)[k])
        assert abs(lk - c.E_iso_k(k)) < 1e-9

    costs = []
    for k in (0, 1):
        e_var, _ = min_over_scale(c, level=k)
        assert e_var <= c.E_iso_k(k) + 1e-12, f"bound violated at root k={k}"
        costs.append((c.E_iso_k(k) - e_var) * 1000)

    gnd, s2 = costs
    # the CLAIM first ...
    assert gnd > 4.0 * s2, (
        f"ground/2^1S posing-cost ratio is only {gnd / s2:.3f}x "
        f"({gnd:.4f} vs {s2:.4f} mHa); measured 4.287x")
    # ... then the regression pins, which are strictly tighter than the claim
    # (they imply it) and exist to catch drift rather than to state the result.
    assert abs(gnd - 4.212) < 0.02, f"ground posing cost {gnd:.4f} mHa"
    assert abs(s2 - 0.983) < 0.01, f"2^1S posing cost {s2:.4f} mHa"
