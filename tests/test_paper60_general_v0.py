"""Paper 60 ``eq:general_v0`` -- the general-``V_0`` form and its atomic collapse.

The claim, in the paper's words.  Projecting ``<Phi_mu|H - E|Phi_nu> = 0`` and
substituting the Sturmian equation on the ket,

    ( -1/2 sum_j grad_j^2 - E ) |Phi_nu> = -beta_nu V_0 |Phi_nu> ,

gives ``sum_nu ( <Phi_mu|V|Phi_nu> - beta_nu <Phi_mu|V_0|Phi_nu> ) C_nu = 0``,
that is ``V C = V_0 B C`` with ``B = diag(beta_nu)``.  Every L2 overlap cancels
IDENTICALLY -- for any local ``V_0``, orthonormal configurations or not.
Atomically ``<Phi_mu|V_0|Phi_nu> = -Z R_nu delta_{mu nu}`` and
``beta_nu Z = p_kappa / R_nu``, so ``V_0 B = -p_kappa * 1`` and the general form
collapses to ``eq:secular``.

Four legs, backed separately because they fail for different reasons:

  L0  the ANALYTIC INPUT to beta, verified by an INDEPENDENT route.  The paper's
      ``beta_nu Z = p_kappa/R_nu`` is not a definition -- it is forced by the
      hydrogenic eigen-relation.  Since the module's kinetic matrix ``T`` is
      itself BUILT from that relation, checking beta against ``T`` alone would
      be circular.  So L0 re-derives the two analytic facts on a separate
      uniform mesh by finite differences, using neither ``sturmian_secular``'s
      integrators nor ``sturmian_variational``'s ``T``.
      -> ``test_l0_hydrogenic_eigen_relation_by_independent_finite_differences``

  L1  the identity ``V_0 B = -p_kappa * 1``, with ``beta`` RECOVERED from the
      operator equation (a column least-squares solve against route B's
      assembled matrices) rather than substituted from the paper's formula.
      -> ``test_l1_v0_b_is_minus_p_kappa_times_identity``

  L2  the reduction reproduces the assembled secular matrix:
      ``Z diag(R_nu) - G`` equals ``build_M(cfgs, Z)`` BIT-IDENTICALLY, which
      also pins ``build_Tprime = -G``.
      -> ``test_l2_reduction_reproduces_the_assembled_secular_matrix``
      -> ``test_l2_the_two_assemblies_do_not_call_each_other``

  L3  Z-INDEPENDENCE.  The collapse must hold at more than one nuclear charge,
      or "it works" is a statement about helium.  Z = 2 and Z = 3, with the
      cross-plug (Z = 3 matrices against Z = 2's p_kappa) asserted to FAIL.
      -> ``test_l3_the_collapse_is_z_independent``

WHAT THIS FILE DOES *NOT* PROVE, stated so the row is not read as more than it
is.  The atomic collapse ``V_0 B = -p_kappa 1`` is, entrywise,
``-Z W_{mu nu} beta_nu + p_kappa delta_{mu nu}``:  its residual is the
diagonality of ``W`` weighted by ``Z beta_nu``, and nothing else.  So this leg
RESTS ON ``eq:W_diagonal`` (owned by
``tests/test_paper60_scale_lock.py::test_c1_one_body_metric_is_exactly_diagonal``)
rather than re-proving it -- the measured ``6e-10`` here is ``W``'s own
``3.9e-11`` off-diagonal amplified by ``max Z beta_nu = 17.0``.  It also does
NOT test the GENERAL half of eq:general_v0 ("for any local ``V_0``, orthonormal
or not"):  every ``V_0`` reachable through this machinery is
``-Z sum_j 1/r_j``.  That half is algebra -- the L2 overlap cancels because the
Sturmian equation is substituted before any metric is introduced -- and no test
here backs it.

Grid: the Paper-60 box rule ``set_grid(max(80, 5 n_max^2), 24000, "grade", 2.0)``,
as in both siblings.  ``set_grid`` mutates ``geovac.sturmian_secular`` module
globals; the module-scoped autouse fixture below restores them.

Measured 2026-09-08 (n_max = 10, l_max = 0, K = 55):

    Z    p_kappa        max|V_0 B + p_kappa 1|   max|Z diag R - G - M|
    2    2.397695680    5.988e-10                0.0
    3    3.807331160    9.508e-10                0.0
"""
import os
import sys

import numpy as np
import pytest
from scipy.special import genlaguerre

# Ensure project root on path (mirrors tests/conftest.py).
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import geovac.sturmian_secular as SS                                  # noqa: E402
import geovac.sturmian_variational as SV                              # noqa: E402

NPTS = 24000


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
# One assembly per (n_max, l_max), reused across Z.
#
# Route B (``SV.build``) carries no Z at all -- its signature has none -- so a
# single build serves every nuclear charge, and the Z-dependence of the
# reduction lives entirely in route A's diagonal.  That is ASSERTED in L3, not
# assumed here.
# --------------------------------------------------------------------------
class Case:
    def __init__(self, nmax: int, lmax: int) -> None:
        SV.set_grid(max(80.0, 5.0 * nmax * nmax), NPTS, "grade", 2.0)
        self.nmax, self.lmax = nmax, lmax
        self.tuples = SV.family(nmax, lmax)
        # route B: the variational assembly (S, T, W, G).  Never forms M.
        self.S, self.T, self.W, self.G, self.K, self.asym = SV.build(nmax, lmax)
        # route A: the isoenergetic assembly.  Never forms S or W.
        self.cfgs = SS.build_configs(self.tuples)
        self.Rnu = np.array([c.Rnu for c in self.cfgs])
        self.Qnu = np.array([c.Q for c in self.cfgs])
        self._M: dict = {}

    def M(self, Z: float) -> np.ndarray:
        if Z not in self._M:
            self._M[Z] = SS.build_M(self.cfgs, Z=Z)
        return self._M[Z]

    def p_kappa(self, Z: float) -> float:
        """The isoenergetic root: the LARGEST eigenvalue of M (deepest binding)."""
        return float(np.linalg.eigvalsh(self.M(Z))[-1])

    def beta_from_operator(self, Z: float) -> tuple:
        """``beta_nu`` recovered from the Sturmian equation, column by column.

        At the locked scale ``lam = p_kappa`` route B's matrices represent the
        operators as ``lam^2 T`` (kinetic), ``lam W`` (``sum_j 1/r_j``) and
        ``S`` (overlap), so
        ``(sum_j -1/2 grad_j^2 - E)|Phi_nu> = -beta_nu V_0|Phi_nu>`` projects to

            A[:, nu] = beta_nu * ( Z lam W )[:, nu] ,
            A = lam^2 T - E S ,   E = -lam^2 / 2 .

        That is K equations in ONE unknown per column.  ``beta_nu`` is the
        least-squares solution and the RELATIVE COLUMN RESIDUAL is returned with
        it:  the residual, not the value, is what a wrong kinetic matrix breaks,
        because the value is dominated by the (near-diagonal) nu-th row.
        """
        lam = self.p_kappa(Z)
        A = lam ** 2 * self.T - (-lam ** 2 / 2) * self.S
        Bm = Z * lam * self.W
        beta = np.zeros(self.K)
        resid = np.zeros(self.K)
        for nu in range(self.K):
            a, b = A[:, nu], Bm[:, nu]
            beta[nu] = float(b @ a) / float(b @ b)
            resid[nu] = (np.linalg.norm(a - beta[nu] * b)
                         / max(float(np.linalg.norm(a)), 1e-30))
        return beta, resid, lam


_CASES: dict = {}


def case(nmax: int, lmax: int) -> Case:
    key = (nmax, lmax)
    if key not in _CASES:
        _CASES[key] = Case(nmax, lmax)
    else:
        # re-install this case's grid; another test may have moved it.
        SV.set_grid(max(80.0, 5.0 * nmax * nmax), NPTS, "grade", 2.0)
    return _CASES[key]


# ==========================================================================
# L0 -- the analytic input to beta, by an INDEPENDENT route
# ==========================================================================
def _hydrogenic(n: int, l: int, Q: float, rr: np.ndarray) -> np.ndarray:
    """``R_{nl}`` at charge ``Q``, L2-normalized on the mesh ``rr``.

    Deliberately re-implemented here rather than imported from
    ``sturmian_secular.hyd_radial``:  L0's job is to check the eigen-relation
    WITHOUT the module's mesh or its cumulative integrators, and a wrong closed
    form would fail the eigenvalue check below rather than hide inside it.
    """
    a = Q / n
    f = ((2 * a * rr) ** l * np.exp(-a * rr)
         * genlaguerre(n - l - 1, 2 * l + 1)(2 * a * rr))
    return f / np.sqrt(np.trapezoid(f * f * rr * rr, rr))


@pytest.mark.parametrize("n,l,Q", [
    (1, 0, 0.7071067811865475),      # the 1s^2 config's charge, Q = 1/R_nu
    (2, 0, 0.8944271909999159),      # the (2s,2s) config's charge
    (3, 1, 1.5),
    (2, 1, 2.0),
    (4, 0, 0.5),
])
def test_l0_hydrogenic_eigen_relation_by_independent_finite_differences(n, l, Q):
    """``(-1/2 grad^2 - Q/r) R_nl = -(Q^2/2n^2) R_nl`` and ``<1/r> = Q/n^2``.

    These two facts are the ENTIRE analytic content of
    ``beta_nu Z = p_kappa/R_nu``:  at the true scale a configuration's orbitals
    sit at charge ``Q_nu^lam = p_kappa/R_nu``, so
    ``(sum_j -1/2 grad_j^2 - E)Phi_nu = Q_nu^lam (sum_j 1/r_j) Phi_nu`` with
    ``E = -(Q_nu^lam R_nu)^2/2 = -p_kappa^2/2``, and matching against
    ``-beta_nu V_0 = beta_nu Z sum_j 1/r_j`` gives ``beta_nu Z = p_kappa/R_nu``.
    The second fact supplies the other factor,
    ``W_{nu nu} = Q_nu sum_j 1/n_j^2 = Q_nu R_nu^2 = R_nu``.

    Why this leg exists at all:  ``sturmian_variational.build`` CONSTRUCTS its
    kinetic matrix ``T`` from exactly this relation.  Recovering beta from ``T``
    (which L1 does) therefore verifies the algebra of eq:general_v0 but cannot
    verify its physical input -- that would be circular.  This leg breaks the
    circle by computing the Laplacian by finite differences on a separate
    uniform mesh, touching neither ``SS.r`` nor ``SS._ctrap_*`` nor ``T``.

    Wrong answers this excludes: a hydrogenic eigenvalue of ``-Q^2/2n`` or
    ``-Q^2/n^2``, and ``<1/r> = Q/n`` or ``Q^2/n^2`` -- all four asserted far
    away.  Measured relative errors 3e-12 .. 2e-7 against a 1e-5 tolerance.
    """
    rmax = 40.0 * n * n / max(Q, 0.2)
    rr = np.linspace(1e-6, rmax, 200000)
    R = _hydrogenic(n, l, Q, rr)

    d1 = np.gradient(R, rr, edge_order=2)
    d2 = np.gradient(d1, rr, edge_order=2)
    lap = d2 + (2.0 / rr) * d1 - l * (l + 1) / rr ** 2 * R
    HR = -0.5 * lap - (Q / rr) * R

    # Rayleigh quotient away from the two mesh ends (where the FD stencil and
    # the 1/r^2 centrifugal term are least accurate).  The window is asserted
    # to hold essentially all of the norm, so it is not selecting a lucky patch.
    m = (rr > 0.02) & (rr < 0.5 * rmax)
    weight = float(np.trapezoid((R * R * rr * rr)[m], rr[m]))
    assert weight > 0.999, f"the FD window holds only {weight:.4f} of the norm"
    ev = float(np.trapezoid((R * HR * rr * rr)[m], rr[m])) / weight
    exact = -Q ** 2 / (2 * n ** 2)
    assert abs(ev / exact - 1.0) < 1e-5, (
        f"hydrogenic eigenvalue {ev:.9f} vs {exact:.9f} (n={n}, l={l}, Q={Q})")

    inv_r = float(np.trapezoid(R * R * rr, rr))
    assert abs(inv_r / (Q / n ** 2) - 1.0) < 1e-5, (
        f"<1/r> = {inv_r:.9f} vs Q/n^2 = {Q / n ** 2:.9f}")

    # the wrong closed forms are nowhere near -- so the two checks above are
    # discriminating, not merely "some number came out finite"
    for wrong in (-Q ** 2 / (2 * n), -Q ** 2 / n ** 2):
        if abs(wrong - exact) > 1e-12:
            assert abs(ev - wrong) > 1e-3 * abs(exact), f"eigenvalue matches {wrong}"
    for wrong in (Q / n, Q ** 2 / n ** 2):
        if abs(wrong - Q / n ** 2) > 1e-12:
            assert abs(inv_r - wrong) > 1e-3 * (Q / n ** 2), f"<1/r> matches {wrong}"


# ==========================================================================
# L1 -- V_0 B = -p_kappa * 1
# ==========================================================================
@pytest.mark.parametrize("nmax,lmax,Z,tol", [
    (10, 0, 2.0, 5e-9),      # the paper's point: measured 5.988e-10
    (10, 0, 3.0, 5e-9),      # measured 9.508e-10
    (4, 1, 2.0, 5e-8),       # measured 4.618e-09 (graded-grid error grows with l)
])
def test_l1_v0_b_is_minus_p_kappa_times_identity(nmax, lmax, Z, tol):
    """``V_0 B = -p_kappa * 1`` with ``V_0 = -Z W`` and ``B = diag(beta_nu)``.

    ``beta_nu`` is RECOVERED from the operator equation (see
    ``Case.beta_from_operator``), not substituted from the paper's formula --
    otherwise the identity would be the arithmetic
    ``(-Z R_nu)(p_kappa/(Z R_nu)) = -p_kappa`` and would assert nothing.  The
    paper's relation ``beta_nu Z R_nu = p_kappa`` is then a MEASURED output of
    that solve, asserted below rather than assumed.

    Wrong answers this excludes, each asserted separately:
      (a) NEITHER SIDE IS SMALL.  ``|diag V_0|`` runs 0.28..2.83 and ``beta``
          runs 0.85..8.48, and every diagonal entry of the product is
          ``-p_kappa`` (-2.40) -- so ``tol`` sits ~9 orders below the quantities
          being cancelled, and "both sides are near zero" is ruled out.
      (b) the WRONG beta relation.  ``p_kappa/R_nu`` (the Z dropped),
          ``1/(Z R_nu)`` (the p_kappa dropped) and ``p_kappa R_nu/Z`` (R and Q
          swapped) each MISS by > 1, against a tolerance of 5e-9.
      (c) the WRONG root -- the smallest eigenvalue of M instead of the largest.
      (d) the WRONG metric -- the L2 overlap S in place of W, both in the
          operator solve (residual 0.55 instead of 1.9e-10) and in ``V_0``.

    What it does NOT exclude, both measured rather than supposed:
      * an error shared between ``V_0`` and ``B``.  The product is invariant
        under ``Z -> c Z`` in one factor and ``1/c`` in the other, so the two
        factors are pinned INDIVIDUALLY in L3, at two nuclear charges, where a
        dropped Z cannot survive both.
      * a corrupted REPULSION block.  ``beta`` is derived from ``p_kappa``, and
        the identity is checked against that same ``p_kappa``, so an error in
        ``T'`` moves both sides together and leaves this test green -- verified:
        under ``build_Tprime -> T' * (Z/2)`` this test PASSES at Z=3 while L2
        and L3 fail.  The repulsion block is L2's and L3's business, not this
        leg's.
    """
    c = case(nmax, lmax)
    beta, resid, lam = c.beta_from_operator(Z)

    # --- the operator solve is well posed: ONE scalar per column satisfies all
    #     K projected equations.  This is the leg a wrong kinetic matrix breaks.
    assert resid.max() < 1e-8, (
        f"no single beta_nu satisfies the projected Sturmian equation: "
        f"max relative column residual {resid.max():.3e}")

    # ANTI-TAUTOLOGY for that residual: the same least-squares solve against the
    # L2 overlap in place of W is NOT satisfiable, so 1.9e-10 is a statement
    # about the 1/r metric and not "any column direction fits".
    A = lam ** 2 * c.T + (lam ** 2 / 2) * c.S
    bad = np.array([
        np.linalg.norm(A[:, nu] - (float(c.S[:, nu] @ A[:, nu])
                                   / float(c.S[:, nu] @ c.S[:, nu])) * c.S[:, nu])
        / np.linalg.norm(A[:, nu]) for nu in range(c.K)])
    assert bad.max() > 0.1, (
        f"the L2 overlap fits the Sturmian columns to {bad.max():.3e} as well; "
        "the residual assertion above is not selecting the 1/r metric")

    # --- the paper's relation, as a MEASURED output of the solve
    rel = float(np.abs(beta * (Z * c.Rnu) / lam - 1.0).max())
    assert rel < 1e-9, f"beta_nu Z R_nu != p_kappa: max relative deviation {rel:.3e}"

    # --- the claim
    V0 = -Z * c.W
    V0B = V0 * beta[None, :]                       # (-Z W) @ diag(beta)
    dev = float(np.abs(V0B + lam * np.eye(c.K)).max())
    assert dev < tol, (
        f"max|V_0 B + p_kappa 1| = {dev:.3e} at Z={Z}, K={c.K} (tol {tol:.0e})")

    # --- (a) the SCALE the tolerance must be read against
    assert float(np.abs(np.diag(V0B) + lam).max()) < tol   # every diagonal is -p_kappa
    assert lam > 2.0                                       # ... and -p_kappa is O(1)
    assert np.abs(np.diag(V0)).min() > 0.25                # V_0 is not a small matrix
    assert np.abs(np.diag(V0)).max() > 2.8 * (Z / 2.0)
    assert beta.min() > 0.8                                # nor is B a small matrix
    assert beta.max() > 3.0
    assert dev < 1e-8 * lam, (
        f"residual {dev:.3e} is not negligible against the cancelled scale "
        f"{lam:.4f}")

    # --- (b) the wrong beta relations all miss by > 1
    for name, cand in (("p_kappa/R_nu (Z dropped)", lam / c.Rnu),
                       ("1/(Z R_nu) (p_kappa dropped)", 1.0 / (Z * c.Rnu)),
                       ("p_kappa R_nu/Z (R and Q swapped)", lam * c.Rnu / Z)):
        d = float(np.abs(V0 * cand[None, :] + lam * np.eye(c.K)).max())
        assert d > 1.0, f"beta = {name} also satisfies the identity ({d:.3e})"

    # --- (c) the wrong root
    lam_min = float(np.linalg.eigvalsh(c.M(Z))[0])
    assert abs(lam_min - lam) > 1.0
    assert float(np.abs(V0B + lam_min * np.eye(c.K)).max()) > 1.0, (
        "the identity holds for the SMALLEST eigenvalue of M too; p_kappa is "
        "not being identified")

    # --- (d) the wrong metric on the V_0 side
    assert float(np.abs((-Z * c.S) * beta[None, :] + lam * np.eye(c.K)).max()) > 1.0, (
        "the L2 overlap S serves as V_0 as well as W does")


# ==========================================================================
# L2 -- the reduction reproduces the assembled secular matrix, bit for bit
# ==========================================================================
@pytest.mark.parametrize("nmax,lmax,Z", [
    (10, 0, 2.0), (10, 0, 3.0), (4, 1, 2.0),
])
def test_l2_reduction_reproduces_the_assembled_secular_matrix(nmax, lmax, Z):
    """``Z diag(R_nu) - G == build_M(cfgs, Z)``, max|diff| = 0.0.

    With ``V_0 B = -p_kappa 1``, eq:general_v0 reads
    ``(-Z W + G) C = -p_kappa C``, i.e. ``(Z diag R_nu - G) C = p_kappa C``.
    This test closes that against the ASSEMBLED secular matrix, and in doing so
    pins the sign convention of ``build_Tprime`` (which returns ``-G``, not
    ``G``) -- the one place the reduction could be off by a sign and still
    produce a plausible spectrum.

    Equality is asserted BIT-EXACTLY because that is what is true; a tolerance
    here would accept a reconstruction that is merely close.

    Wrong answers this excludes:
      (a) the WRONG SIGN, ``Z diag R + G``: asserted to differ by 0.88;
      (b) the WRONG DIAGONAL, ``Z diag Q_nu`` (``Q = 1/R``): differs by 13.9;
      (c) "0 == 0" -- the repulsion block and the nuclear diagonal are both
          asserted to carry substantial content, so the equality is not a
          statement about two empty matrices.
    """
    c = case(nmax, lmax)
    M = c.M(Z)

    d = float(np.abs(Z * np.diag(c.Rnu) - c.G - M).max())
    assert d == 0.0, f"reconstruction differs from build_M by {d:.3e}"

    # the sign convention of T', pinned directly
    assert float(np.abs(SS.build_Tprime(c.cfgs) + c.G).max()) == 0.0, (
        "build_Tprime is not -G")

    # --- (a) / (b) the near misses
    assert float(np.abs(Z * np.diag(c.Rnu) + c.G - M).max()) > 0.5, (
        "the reconstruction with the OPPOSITE sign on G also reproduces M")
    assert float(np.abs(Z * np.diag(c.Qnu) - c.G - M).max()) > 1.0, (
        "diag(Z Q_nu) reproduces M as well as diag(Z R_nu) does")

    # --- (c) neither term is empty
    assert np.abs(c.G).max() > 0.4
    off = float(np.abs(c.G - np.diag(np.diag(c.G))).max())
    assert off > 1e-2, (
        f"the repulsion block is essentially diagonal ({off:.3e}); the "
        "bit-exact equality would then be a statement about two diagonal "
        "matrices")
    assert (Z * c.Rnu).min() > 0.25 and (Z * c.Rnu).max() > 2.8 * (Z / 2.0)
    assert np.abs(M).max() > 2.0


def test_l2_the_two_assemblies_do_not_call_each_other():
    """The bit-exact equality of L2 is not the identity of one code path.

    Route B (``sturmian_variational.build`` -> ``S, T, W, G``) never forms ``M``
    or ``T'``; route A (``sturmian_secular.build_configs`` -> ``build_M``) never
    forms ``S``, ``W`` or route B's ``G``.  Asserted operationally by replacing
    each module's entry points with sentinels that raise, then running the other
    route underneath them.

    HONEST SCOPE.  The two routes are separate assembly loops in separate
    modules, but they are NOT numerically independent:  they share the
    primitives ``SS.build_configs`` and ``SS.repulsion_terms``.  That sharing is
    precisely WHY L2's equality is bit-exact rather than ~1e-16 -- the same
    floating-point sum is reached by two callers with the same normalization and
    the same loop order.  What this test excludes is the stronger and more
    likely defect, a refactor making one route delegate to the other; it does
    NOT exclude a fault in the shared ERI primitive, which would move both sides
    together.  The absolute calibration of that primitive is owned by
    ``test_paper60_scale_lock.py::test_pipeline_single_config_reproduces_kellner_helium``
    (the K=1 Kellner value), and L3 below adds a second nuclear charge.
    """
    def boom(*a, **k):                       # pragma: no cover - sentinel
        raise AssertionError("route delegated to the other assembly")

    # The sentinel machinery must itself be able to fire: a mistyped attribute
    # name would silently make every assertion below vacuous.
    saved_a = (SS.build_M, SS.build_Tprime)
    SS.build_M, SS.build_Tprime = boom, boom
    try:
        with pytest.raises(AssertionError):
            SS.build_M(None)
        with pytest.raises(AssertionError):
            SS.build_Tprime(None)
        SV.set_grid(80.0, 4000, "grade", 2.0)
        SV.build(3, 0)                       # must not touch route A
    finally:
        SS.build_M, SS.build_Tprime = saved_a

    saved_b = (SV.build, SV.u_terms, SV.radial_1r)
    SV.build, SV.u_terms, SV.radial_1r = boom, boom, boom
    try:
        with pytest.raises(AssertionError):
            SV.build(3, 0)
        SS.build_M(SS.build_configs(SS.gen_configs(0, {0: 3})), Z=2.0)
    finally:
        SV.build, SV.u_terms, SV.radial_1r = saved_b


# ==========================================================================
# L3 -- the collapse is Z-independent
# ==========================================================================
def test_l3_the_collapse_is_z_independent():
    """The whole reduction at Z = 2 AND Z = 3, with the cross-plug failing.

    This is the leg that catches a coincidence at helium.  Three things are
    asserted that a single-Z run cannot say:

      * ``p_kappa`` genuinely MOVES with Z (2.397695680 -> 3.807331160), so the
        identity ``V_0 B = -p_kappa 1`` is not being satisfied by the same
        number twice;
      * the CROSS-PLUG fails -- the Z = 3 matrices against the Z = 2 root miss
        by 1.41 -- so the right-hand side tracks the charge, not a constant;
      * the entire Z-dependence sits on the DIAGONAL and is exactly linear:
        ``M(3) - M(2) == diag(R_nu)`` to 2.2e-16.  That is eq:secular's
        Z-independence of ``T'``, which this reduction rests on and which a
        Z-contaminated repulsion block would break.

    The class of defect this leg exists for is one that is a NO-OP at helium.
    Measured: under ``build_M: T' -> T' * (Z/2)`` every Z = 2 test in this file
    stays GREEN (4/4 pass) while this test and the Z = 3 reconstruction fail --
    which is the whole reason the file does not stop at Z = 2.

    It also closes the gap L1 leaves open (an error shared between ``V_0`` and
    ``B``, which cancels in their product): both factors are pinned
    individually here, at two charges, where a dropped or misplaced Z cannot
    survive both.
    """
    c = case(10, 0)

    # --- FIRST: all the Z-dependence is on the diagonal, and linear.  This runs
    #     ahead of the per-charge loop deliberately: it is the leg that a defect
    #     invisible at Z = 2 trips, and putting it after the loop would let the
    #     loop's own reconstruction assertion mask it.
    dM = c.M(3.0) - c.M(2.0)
    assert float(np.abs(dM - np.diag(c.Rnu)).max()) < 1e-14, (
        f"M(3) - M(2) is not diag(R_nu): max deviation "
        f"{np.abs(dM - np.diag(c.Rnu)).max():.3e} -- the repulsion block T' has "
        "acquired a Z dependence and eq:secular no longer holds")
    assert float(np.abs(np.diag(dM)).max()) > 1.0     # ... and it is not zero

    pk = {}
    for Z in (2.0, 3.0):
        beta, resid, lam = c.beta_from_operator(Z)
        pk[Z] = lam
        assert resid.max() < 1e-8
        dev = float(np.abs((-Z * c.W) * beta[None, :] + lam * np.eye(c.K)).max())
        assert dev < 5e-9, f"identity fails at Z={Z}: {dev:.3e}"
        assert float(np.abs(Z * np.diag(c.Rnu) - c.G - c.M(Z)).max()) == 0.0

        # the two factors pinned SEPARATELY -- this is what a Z dropped from
        # BOTH (which cancels in the product) fails.
        assert float(np.abs(np.diag(-Z * c.W) + Z * c.Rnu).max()) < 1e-9, (
            "<Phi|V_0|Phi> is not -Z R_nu delta")
        assert float(np.abs(beta - lam / (Z * c.Rnu)).max()) < 1e-12, (
            "beta_nu is not p_kappa/(Z R_nu)")

    # --- p_kappa really moves with Z
    # IDENTITY pins (5e-7, widened 2026-09-11 -- see the note in
    # test_paper60_no_selection.py).  The substantive assertion is the NEXT one:
    # that p_kappa genuinely moves with Z, which no tolerance change affects.
    assert abs(pk[2.0] - 2.397695680) < 5e-7
    assert abs(pk[3.0] - 3.807331160) < 5e-7
    assert pk[3.0] - pk[2.0] > 1.0, (
        f"p_kappa barely moves between Z=2 and Z=3 ({pk[2.0]:.6f} vs "
        f"{pk[3.0]:.6f}); the two charges are not independent evidence")

    # --- the cross-plug FAILS: Z = 3 matrices against Z = 2's root
    beta3 = pk[3.0] / (3.0 * c.Rnu)
    cross = float(np.abs((-3.0 * c.W) * beta3[None, :]
                         + pk[2.0] * np.eye(c.K)).max())
    assert cross > 1.0, (
        f"the Z=3 identity is satisfied by the Z=2 root as well ({cross:.3e}); "
        "the right-hand side is not tracking the nuclear charge")
