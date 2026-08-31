"""
Paper 34 projection spot-checks --- batch 1 (load-bearing rows).

Companion to tests/test_paper34_projection_spot_checks.py (6 of 28
projections covered). This file adds 8 load-bearing rows from
followon_register.md A9 batch 1:

  §III.1   Fock conformal (sec:proj_fock)
  §III.5   Sturmian (sec:proj_sturmian)
  §III.11  Vector-photon promotion (sec:proj_vector_photon)
  §III.13  Drake--Swainson asymptotic subtraction (sec:proj_drake_swainson)
  §III.16  Two-body Dirac / Breit retardation (sec:proj_breit_retardation)
  §III.17  Foldy/Friar charge density (sec:proj_charge_density)
  §III.18  Zemach magnetization density (sec:proj_magnetization_density)
  §III.19  Tensor multipole (sec:proj_tensor_multipole)

After batch 1: 14 of 28 Paper 34 projections covered.

Per CLAUDE.md §13.4a, each test verifies the projection's stated
transcendental signature and/or its load-bearing identity through
analytical limit, symbolic identity, or numerical cross-check.
"""

from __future__ import annotations

import functools
import math
import numpy as np
import pytest
import sympy as sp


# ----------------------------------------------------------------------------
# Shared symbolic helpers (added by the 2026-08-28 adversarial audit).
#
# That audit found several checks in these batch files that BUILT a quantity
# from a formula and then asserted the same formula (tautologies with no
# discriminating power).  The helpers below let the repaired tests DERIVE the
# quantities they check: the hydrogenic radial function is normalized by an
# explicit integral rather than by quoting a normalization constant, so a
# wrong downstream value cannot be hidden by a matching hardcoded prefactor.
# ----------------------------------------------------------------------------

_R_SYM = sp.Symbol('r', positive=True)
_Z_SYM = sp.Symbol('Z', positive=True)


@functools.lru_cache(maxsize=None)
def _hydrogenic_radial_normalized(n: int, l: int):
    """Normalized hydrogenic radial function R_{nl}(r) at charge Z.

    DERIVED, not quoted: the associated-Laguerre shape is built and then
    divided by sqrt(int_0^inf R^2 r^2 dr), computed symbolically.
    """
    r, Z = _R_SYM, _Z_SYM
    rho = 2 * Z * r / n
    shape = rho ** l * sp.exp(-rho / 2) * sp.assoc_laguerre(n - l - 1, 2 * l + 1, rho)
    norm_sq = sp.simplify(sp.integrate(shape ** 2 * r ** 2, (r, 0, sp.oo)))
    return sp.simplify(shape / sp.sqrt(norm_sq))


@functools.lru_cache(maxsize=None)
def _hydrogenic_density_at_origin(n: int, l: int):
    """|psi_{nlm}(0)|^2 for the normalized hydrogenic orbital.

    |psi(0)|^2 = |R_{nl}(0)|^2 |Y_{lm}(0)|^2; only l = 0 survives, with
    |Y_00|^2 = 1/(4 pi).  Derived from _hydrogenic_radial_normalized.
    """
    R0 = sp.simplify(sp.limit(_hydrogenic_radial_normalized(n, l), _R_SYM, 0))
    if l != 0:
        return sp.simplify(R0 ** 2)
    return sp.simplify(R0 ** 2 / (4 * sp.pi))


# ----------------------------------------------------------------------------
# §III.1  Fock conformal: kappa = -1/16, Vol(S^3) = 2 pi^2, n^2 - 1 spectrum
# ----------------------------------------------------------------------------

def test_paper34_III1_kappa_rational_prefactor():
    """Paper 34 §III.1 (sec:proj_fock): transcendental signature 'rational
    prefactor kappa = -1/16 (Rydberg-to-graph-eigenvalue conversion).'

    Verifies the universal topological constant kappa = -1/16 used to
    map graph eigenvalues to Rydberg energies (CLAUDE.md §4).

    Per CLAUDE.md §8 ('-1/16 is the universal topological constant
    can be used directly'), this value is the framework's free-to-use
    rational prefactor. Production-side symbolic version lives at
    geovac.graph_qed_propagator.KAPPA_SCALAR = sympy Rational(-1, 16).
    """
    from geovac.graph_qed_propagator import KAPPA_SCALAR
    from sympy import Rational

    # Symbolic exactness check
    assert KAPPA_SCALAR == Rational(-1, 16), (
        f"KAPPA_SCALAR = {KAPPA_SCALAR} != -1/16"
    )

    # Numerical exactness check
    assert math.isclose(float(KAPPA_SCALAR), -1.0 / 16.0,
                        rel_tol=0.0, abs_tol=0.0), (
        f"float(KAPPA_SCALAR) = {float(KAPPA_SCALAR)} != -1/16"
    )

    # 2026-08-28 audit: tie the constant to the place it is USED, so the test
    # fails if the solver stops applying it.  Previously this test only
    # compared the symbol to a second literal copy of itself.
    from geovac.atomic_solver import AtomicSolver
    for Z in (1, 2, 3):
        solver = AtomicSolver(3, Z)
        assert math.isclose(solver.kinetic_scale, float(KAPPA_SCALAR) * Z ** 2,
                            rel_tol=0.0, abs_tol=0.0), (
            f"AtomicSolver(Z={Z}).kinetic_scale = {solver.kinetic_scale} "
            f"!= kappa * Z^2 = {float(KAPPA_SCALAR) * Z ** 2}"
        )


def test_paper34_III1_S3_volume_2pi_squared():
    """Paper 34 §III.1: 'pi enters through Vol(S^3) = 2 pi^2 when
    subsequent spectral integrals are performed.'

    Confirms the same Vol(S^3) constant referenced by §III.2 (Hopf
    bundle) and §III.6 (spectral action) -- so all three projections
    share a single Vol(S^3) measure source.
    """
    from geovac.hopf_bundle import VOL_S3

    # 2026-08-28 audit: derive the value from the general
    # Vol(S^n) = 2 pi^{(n+1)/2} / Gamma((n+1)/2) instead of comparing the
    # production constant to a second literal copy of itself.
    vol_sym = sp.simplify(2 * sp.pi ** sp.Rational(4, 2) / sp.gamma(sp.Rational(4, 2)))
    assert sp.simplify(vol_sym - 2 * sp.pi ** 2) == 0, (
        f"Gamma-function derivation gives Vol(S^3) = {vol_sym}, not 2 pi^2"
    )
    assert math.isclose(VOL_S3, float(vol_sym), rel_tol=1e-15, abs_tol=1e-15), (
        f"Production Vol(S^3) = {VOL_S3} != derived {float(vol_sym)}"
    )


def _s3_laplace_beltrami_radial(R_chi, l, chi):
    """Unit-S^3 Laplace-Beltrami operator on R(chi) Y_lm, divided by Y_lm.

    With ds^2 = dchi^2 + sin^2(chi) dOmega_2^2 on the unit S^3,

        Delta_{S^3} f = sin^{-2}(chi) d_chi( sin^2(chi) d_chi f )
                        + sin^{-2}(chi) Delta_{S^2} f,

    and Delta_{S^2} Y_lm = -l(l+1) Y_lm.
    """
    return (sp.diff(sp.sin(chi) ** 2 * sp.diff(R_chi, chi), chi) / sp.sin(chi) ** 2
            - l * (l + 1) * R_chi / sp.sin(chi) ** 2)


@pytest.mark.parametrize(
    "n,l", [(n, l) for n in range(1, 7) for l in range(n)]
)
def test_paper34_III1_laplacian_spectrum_n2_minus_1(n, l):
    """Paper 34 §III.1 / CLAUDE.md §4: 'Eigenvalues of the Laplace-Beltrami
    operator on unit S^3 are pure integers: lambda_n = -(n^2 - 1).'

    REWRITTEN 2026-08-28 (adversarial audit).  The previous body was

        lam = -(n**2 - 1); assert isinstance(lam, int); assert lam == 1 - n**2

    i.e. it asserted -(n^2-1) == 1-n^2 and never touched an operator at all;
    it could not fail for any spectrum whatsoever.  This version APPLIES the
    unit-S^3 Laplace-Beltrami operator to the actual hyperspherical harmonic
    Y_{nlm} ~ sin^l(chi) C^{l+1}_{n-l-1}(cos chi) Y_lm and checks that the
    eigenvalue is -(n^2 - 1), symbolically, for every (n, l) with n <= 6 --
    including its l-independence (the SO(4) degeneracy).

    HONEST SCOPE (consistent with docs/claim_test_matrix.md, Papers 1/7 row):
    this is a CONTINUUM property of the round-S^3 Laplace-Beltrami operator.
    The discrete graph Laplacian L = D - A is positive-semidefinite and does
    NOT carry -(n^2-1); that distinction is deliberate.
    """
    chi = sp.Symbol('chi')
    Y_radial = sp.sin(chi) ** l * sp.gegenbauer(n - l - 1, l + 1, sp.cos(chi))

    residual = sp.simplify(sp.expand_trig(sp.expand(sp.trigsimp(
        _s3_laplace_beltrami_radial(Y_radial, l, chi) + (n ** 2 - 1) * Y_radial
    ))))
    assert residual == 0, (
        f"Delta_S3 Y_(n={n},l={l}) != -(n^2-1) Y; residual = {residual}"
    )

    # Non-tautology guard: shifting the eigenvalue by 1 must NOT vanish, so
    # the assertion above is genuinely sensitive to the eigenvalue's value.
    wrong = sp.simplify(sp.expand_trig(sp.expand(sp.trigsimp(
        _s3_laplace_beltrami_radial(Y_radial, l, chi) + (n ** 2) * Y_radial
    ))))
    assert wrong != 0, (
        f"guard failed: eigenvalue -(n^2) also 'works' at (n={n}, l={l})"
    )


# ----------------------------------------------------------------------------
# §III.5  Sturmian relabeling at lambda = Z/n preserves rationality
# ----------------------------------------------------------------------------

def test_paper34_III5_sturmian_rationality_preserved():
    """Paper 34 §III.5 (sec:proj_sturmian): transcendental signature
    'preserves rationality of Layer 1. No new pi or other
    transcendentals enter at this step.'

    Tests by symbolic relabeling: the Sturmian focal length lam = Z/n
    is rational in (Z, n), and the graph eigenvalues -(n^2 - 1) remain
    integer under this relabeling. This is the categorical statement
    Paper 34 §III.5 claims.
    """
    Z, n = sp.symbols('Z n', positive=True, integer=True)
    lam_sturmian = Z / n
    # lam should be a rational function of (Z, n)
    assert lam_sturmian.is_rational_function(Z, n), (
        f"Sturmian focal length lam = Z/n = {lam_sturmian} is not rational"
    )

    # The energy formula E_n = -Z^2/(2 n^2) under Sturmian relabeling
    # E -> -lam^2 / 2 = -(Z/n)^2 / 2 stays rational
    E_sturmian = -lam_sturmian ** 2 / 2
    E_direct = -Z ** 2 / (2 * n ** 2)
    assert sp.simplify(E_sturmian - E_direct) == 0, (
        f"Sturmian energy mismatch: {E_sturmian} != {E_direct}"
    )

    # No pi or other transcendental enters at this level (Layer 1 only)
    # NOTE: `sp.pi in expr.free_symbols` is ALWAYS False (pi is a NumberSymbol,
    # not a Symbol), so the old form of this guard never fired.  Use .has().
    assert not E_sturmian.has(sp.pi), (
        f"Unexpected pi in Sturmian Layer 1: {E_sturmian}"
    )

    # 2026-08-28 audit: the guard above is only informative if it CAN fire.
    # A downstream projection that DOES inject pi (the Hopf/Fock measure
    # factor Vol(S^3) = 2 pi^2) must be caught by exactly the same predicate.
    assert (E_sturmian * 2 * sp.pi ** 2).has(sp.pi), (
        "non-tautology guard failed: .has(sp.pi) cannot detect an injected pi"
    )

    # 2026-08-28 audit: and the operational statement, on PRODUCTION code
    # rather than on a hand-built symbol.  Relabelling at lam = Z/n is an
    # exact rational rescaling, so the solved spectrum must satisfy
    # spec(Z) = Z^2 * spec(1) to machine precision at every Z.  Any
    # transcendental injected by the relabelling would be Z-dependent and
    # would break this.
    from geovac.atomic_solver import AtomicSolver

    ref = None
    for Z_val in (1, 2, 3, 5):
        solver = AtomicSolver(4, Z_val)
        eigs, _ = solver.compute_ground_state(n_states=8)
        scaled = np.sort(np.real(eigs)) / float(Z_val) ** 2
        if ref is None:
            ref = scaled
        else:
            dev = float(np.max(np.abs(scaled - ref)))
            assert dev < 1e-13, (
                f"Sturmian Z^2 relabelling not exact at Z={Z_val}: max dev = {dev}"
            )


def _velocity_form_closure(n: int, l: int):
    """Velocity-form closure I_v(n,l), DERIVED rather than quoted.

    The Bethe-logarithm denominator is the Thomas-Reiche-Kuhn-type sum

        I_v(n,l) = sum_m |<nl|p|m>|^2 (E_m - E_n)
                 = (1/2) <nl| [p,[H,p]] |nl>
                 = (1/2) <nl| grad^2 V |nl>
                 = 2 pi Z |psi_{nl}(0)|^2

    using grad^2 (-Z/r) = 4 pi Z delta^3(r).  |psi_{nl}(0)|^2 comes from the
    SYMBOLICALLY NORMALIZED hydrogenic orbital, so nothing about the value
    2 Z^4/n^3 is assumed anywhere in the chain.
    """
    return sp.simplify(2 * sp.pi * _Z_SYM * _hydrogenic_density_at_origin(n, l))


@pytest.mark.parametrize("n", [1, 2, 3])
def test_paper34_III5_sturmian_iv_closure_at_ell0(n):
    """Paper 34 §III.5 + §III.13: the velocity-form closure
    I_v(nl) = (Z^4/n^3) delta_{l,0} vanishes exactly for l > 0, and at l = 0
    is recovered by the Drake-Swainson structural denominator
    D_drake(n, 0) = 2 Z^4 / n^3.

    REWRITTEN 2026-08-28 (adversarial audit).  The previous body computed
    2*(2*0+1)*Z**4/n**3 and asserted it equals 2*Z**4/n**3 -- the same
    expression twice, a restatement of its own construction with no
    discriminating power.  Here I_v is DERIVED from the closure identity
    I_v = 2 pi Z |psi_nl(0)|^2 with |psi_nl(0)|^2 obtained from an explicitly
    normalized hydrogenic orbital, so both the l = 0 value and the l > 0
    vanishing are computed results that can disagree with the paper.
    """
    Z = _Z_SYM

    # l = 0: closure gives 2 Z^4 / n^3, matching D_drake(n, 0).
    I_v_ell0 = _velocity_form_closure(n, 0)
    D_drake_ell0 = 2 * (2 * 0 + 1) * Z ** 4 / sp.Integer(n) ** 3
    assert sp.simplify(I_v_ell0 - D_drake_ell0) == 0, (
        f"n={n}: derived I_v(nS) = {I_v_ell0} != D_drake(n,0) = {D_drake_ell0}"
    )

    # l > 0: the closure vanishes identically -- which is WHY Drake-Swainson
    # needs a structural denominator at all.
    for l in range(1, n):
        val = _velocity_form_closure(n, l)
        assert sp.simplify(val) == 0, (
            f"n={n}, l={l}: I_v should vanish for l > 0, got {val}"
        )


# ----------------------------------------------------------------------------
# §III.11 Vector-photon promotion: 1/(4 pi) per loop is S^2 Weyl factor
# ----------------------------------------------------------------------------

def test_paper34_III11_vector_photon_1_over_4pi():
    """Paper 34 §III.11 (sec:proj_vector_photon): transcendental signature
    '1/(4 pi) per loop, identified as the S^2 Weyl exchange constant
    of the Hopf base (Paper 33).'

    Verifies the structural identity 1/Vol(S^2) = 1/(4 pi), tying this
    projection to the same Hopf base measure source as §III.2.
    """
    from geovac.hopf_bundle import VOL_S2

    # Vol(S^2) = 4 pi (standard) -> 1/Vol(S^2) = 1/(4 pi)
    assert math.isclose(VOL_S2, 4.0 * math.pi, rel_tol=1e-15), (
        f"Vol(S^2) = {VOL_S2} != 4 pi"
    )

    one_over_4pi = 1.0 / VOL_S2
    expected = 1.0 / (4.0 * math.pi)
    assert math.isclose(one_over_4pi, expected, rel_tol=1e-15), (
        f"1/Vol(S^2) = {one_over_4pi} != 1/(4 pi) = {expected}"
    )

    # Symbolic cross-check
    vol_S2_sym = 4 * sp.pi
    assert sp.simplify(1 / vol_S2_sym - sp.Rational(1, 4) / sp.pi) == 0

    # 2026-08-28 audit: the checks above only relate two spellings of 4 pi.
    # The projection claim is that the 1/(4 pi) the framework carries PER
    # LOOP is this same S^2 measure factor.  Read it off the PRODUCTION
    # vertex amplitude: geovac.vector_qed.vertex_coupling returns
    #   sqrt((2l_a+1)(2q+1)(2l_b+1)/(4 pi)) * (-1)^(l_a-m_a) * 3j,
    # so squaring and dividing out the (integer) degeneracy product and the
    # (algebraic) 3j^2 must leave exactly 1/Vol(S^2) -- one factor of
    # 1/(4 pi) per closed vertex pair, i.e. per loop.
    from sympy.physics.wigner import wigner_3j
    from geovac.vector_qed import vertex_coupling

    checked = 0
    for l_a in range(0, 3):
        for l_b in range(0, 3):
            for q in range(1, 3):
                for m_a in range(-l_a, l_a + 1):
                    for m_b in range(-l_b, l_b + 1):
                        m_q = m_a - m_b
                        if abs(m_q) > q:
                            continue
                        V = vertex_coupling(1, l_a, m_a, 2, l_b, m_b, q, m_q)
                        threej = float(wigner_3j(l_a, q, l_b, -m_a, m_q, m_b))
                        if V == 0.0 or abs(threej) < 1e-14:
                            continue
                        degeneracy = (2 * l_a + 1) * (2 * q + 1) * (2 * l_b + 1)
                        residual = V ** 2 / (degeneracy * threej ** 2)
                        assert math.isclose(residual, 1.0 / VOL_S2,
                                            rel_tol=1e-12), (
                            f"vertex ({l_a},{m_a}|{l_b},{m_b}; q={q},{m_q}) "
                            f"leaves {residual}, expected 1/Vol(S^2) = "
                            f"{1.0 / VOL_S2}"
                        )
                        checked += 1
    assert checked >= 6, f"only {checked} non-vanishing vertices exercised"


# ----------------------------------------------------------------------------
# §III.13 Drake-Swainson: D_drake(n, l) = 2 (2l + 1) Z^4 / n^3
# ----------------------------------------------------------------------------

@pytest.mark.parametrize("n,l", [
    (1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2),
])
def test_paper34_III13_drake_swainson_structural_denominator(n, l):
    """Paper 34 §III.13 (sec:proj_drake_swainson): structural denominator
    D_drake(n, l) = 2 (2l + 1) Z^4 / n^3, read as
    (spin) x (angular degeneracy 2l+1) x (hydrogenic density Z^4/n^3), which
    'recovers D_drake(nS) = I_v(nS) = 2 Z^4/n^3 for l = 0'.

    REWRITTEN 2026-08-28 (adversarial audit).  The previous body built
    sp.Rational(2*(2l+1)*Z**4, n**3), asserted `.is_rational` (always True
    for a sp.Rational -- an always-true predicate), and at l=0 compared it to
    the same expression.  It could not fail for any denominator.  This
    version checks the two things the paper's factorization actually
    asserts, against an INDEPENDENTLY DERIVED hydrogenic density:

      (a) the l = 0 anchor equals the derived velocity-form closure
          I_v(nS) = 2 pi Z |psi_nS(0)|^2 (see _velocity_form_closure), and
      (b) the entire l-dependence is the angular degeneracy 2l+1, with no
          residual n or Z dependence.
    """
    Z = _Z_SYM
    D = 2 * (2 * l + 1) * Z ** 4 / sp.Integer(n) ** 3
    D0 = 2 * Z ** 4 / sp.Integer(n) ** 3

    # (a) l = 0 anchor is the derived closure, not a restated formula.
    closure = _velocity_form_closure(n, 0)
    assert sp.simplify(D0 - closure) == 0, (
        f"n={n}: D_drake(n,0) = {D0} != derived I_v(nS) = {closure}"
    )

    # (b) the whole l-dependence is the angular degeneracy.
    ratio = sp.simplify(D / D0)
    assert ratio == 2 * l + 1, (
        f"D_drake(n={n},l={l})/D_drake(n,0) = {ratio}, expected 2l+1 = {2*l+1}"
    )
    assert ratio.free_symbols == set(), (
        f"degeneracy ratio should be a pure number, got {ratio}"
    )

    # Non-tautology guard: a bare spin factor 2 (instead of the angular
    # 2l+1) must be distinguishable whenever l > 0.
    if l > 0:
        assert ratio != 2, (
            "guard failed: 2l+1 is indistinguishable from a constant 2"
        )


# ----------------------------------------------------------------------------
# §III.16 Breit retardation: R^0_BP closed forms (alpha^4 * Q[log 2, log 3])
# ----------------------------------------------------------------------------

def test_paper34_III16_breit_R0_1s1s_1s1s():
    """Paper 34 §III.16 (sec:proj_breit_retardation): closed form
    R^0_BP(1s, 1s; 1s, 1s) = -5 + 8 log 2 (Z=1, k=0).

    Validates the Z=1 closed form via the production
    geovac.breit_integrals.compute_radial routine.
    """
    from geovac.breit_integrals import compute_radial

    val = compute_radial(
        n1=1, l1=0, n3=1, l3=0,
        n2=1, l2=0, n4=1, l4=0,
        k=0, kernel_type="breit", Z=1,
    )
    expected = -5 + 8 * sp.log(2)
    diff = sp.simplify(val - expected)
    assert diff == 0, (
        f"R^0_BP(1s,1s;1s,1s) at Z=1: production = {val}, "
        f"Paper 34 stated = {expected}, diff = {diff}"
    )


def test_paper34_III16_breit_both_1s2s_orderings():
    """Paper 34 §III.16 + Appendix Breit table: the two distinct 1s/2s
    BP-retarded radial integrals at Z=1, k=0, pinned by explicit quantum
    numbers so that no label convention is needed to read the test.

    Production (geovac.breit_integrals.compute_radial, whose own docstring
    defines R^k_BP(n1 l1, n3 l3; n2 l2, n4 l4) = int int P_13(r1) P_24(r2)
    K(r1,r2), i.e. slots 1-2 of the LABEL are electron 1's density and
    slots 3-4 are electron 2's) gives:

      * mixed x mixed densities   P_13 = 1s.2s , P_24 = 1s.2s  ->  4/81
        (a pure rational -- this is the exchange-type combination)
      * pure  x pure  densities   P_13 = 1s.1s , P_24 = 2s.2s
        ->  -19/9 + log(81 sqrt(3)/16) = -4 log 2 - 19/9 + 9 log(3)/2
        (carries the log content -- the direct-type combination)

    PAPER DEFECT SURFACED 2026-08-28 (adversarial audit).  Paper 34 states
    BOTH values but attaches them to labels that contradict each other
    between two locations:

      §III.16 body:   "R^0_BP(1s,2s;1s,2s) = 4/81 ... a pure rational,
                       no log content"
      Appendix table: "(1s,1s; 2s,2s) = 4/81"  and
                      "(1s,2s; 1s,2s) = -4 log 2 - 19/9 + 9 log(3)/2"

    The same label (1s,2s;1s,2s) therefore carries two different values in
    one paper.  Each location is self-consistent under a DIFFERENT reading
    of R^k(a,b;c,d) (module convention = a,b on electron 1; Condon-Shortley
    = a,c on electron 1), and the paper never states which it uses.  Both
    NUMBERS are correct; the label convention needs to be stated once and
    applied in both places.  RESOLVED 2026-08-28: the appendix table's two
    labels were swapped back to the module convention, which is now stated
    inline beside that table.  This test pins the production values that
    settled which locus was right.

    (This supersedes the 2026-06-04 note formerly carried in this
    docstring, which read the mismatch as a wrong VALUE in the paper.  It
    is not: it is a convention collision.)
    """
    from geovac.breit_integrals import compute_radial

    mixed_x_mixed = compute_radial(
        n1=1, l1=0, n3=2, l3=0,      # electron 1 density: 1s * 2s
        n2=1, l2=0, n4=2, l4=0,      # electron 2 density: 1s * 2s
        k=0, kernel_type="breit", Z=1,
    )
    assert sp.simplify(mixed_x_mixed - sp.Rational(4, 81)) == 0, (
        f"P_13 = P_24 = 1s.2s: production = {mixed_x_mixed}, expected 4/81"
    )
    assert not mixed_x_mixed.has(sp.log), (
        f"the mixed-density integral should be log-free, got {mixed_x_mixed}"
    )

    pure_x_pure = compute_radial(
        n1=1, l1=0, n3=1, l3=0,      # electron 1 density: 1s * 1s
        n2=2, l2=0, n4=2, l4=0,      # electron 2 density: 2s * 2s
        k=0, kernel_type="breit", Z=1,
    )
    expected_log = -4 * sp.log(2) - sp.Rational(19, 9) + 9 * sp.log(3) / 2
    assert sp.simplify(pure_x_pure - expected_log) == 0, (
        f"P_13 = 1s.1s, P_24 = 2s.2s: production = {pure_x_pure}, "
        f"expected {expected_log}"
    )
    assert pure_x_pure.has(sp.log), (
        "the pure-density integral is supposed to carry the Q[log 2, log 3] "
        f"content, got {pure_x_pure}"
    )

    # The two are genuinely different objects -- this is what makes the
    # label collision above a real defect rather than a notational nit.
    assert sp.simplify(mixed_x_mixed - pure_x_pure) != 0

    # Appendix row 4, checked here because nothing else covers it.
    val_2s = compute_radial(2, 0, 2, 0, 2, 0, 2, 0, k=0,
                            kernel_type="breit", Z=1)
    assert sp.simplify(val_2s - (-sp.Rational(175, 256) + sp.log(2))) == 0, (
        f"R^0_BP(2s,2s;2s,2s) = {val_2s}, expected -175/256 + log 2"
    )


def test_paper34_III16_breit_Z3_scaling():
    """Paper 34 §III.16: Breit retardation integrals scale as Z^3
    (per geovac.breit_integrals docstring: 'Breit integrals scale as Z^3
    while Coulomb integrals scale as Z^1').

    This is the rest-mass-free Z-scaling that the alpha^4 * Q
    ring-preserving claim relies on.
    """
    from geovac.breit_integrals import compute_radial

    Z1 = compute_radial(1, 0, 1, 0, 1, 0, 1, 0, k=0, kernel_type="breit", Z=1)
    Z3 = compute_radial(1, 0, 1, 0, 1, 0, 1, 0, k=0, kernel_type="breit", Z=3)
    ratio = sp.simplify(Z3 / Z1)
    assert ratio == 27, (
        f"Z^3 scaling check failed: Z=3 / Z=1 = {ratio} (expected 27 = 3^3)"
    )


# ----------------------------------------------------------------------------
# §III.17 Foldy/Friar: (2 pi / 3) Z alpha <r^2>_E delta^3(r) contact term
# ----------------------------------------------------------------------------

def test_paper34_III17_foldy_friar_prefactor_2pi_over_3():
    """Paper 34 §III.17 (sec:proj_charge_density): Foldy/Friar contact term
    Delta V = +(2 pi / 3) Z alpha <r^2>_E delta^3(r).

    REWRITTEN 2026-08-28 (adversarial audit).  The previous body asserted
    math.isclose(2*pi/3, 2.0944) and (2/3)*pi == 4*pi/6 -- two restatements
    of arithmetic about pi that say nothing about the Foldy/Friar
    projection, and could not fail for any prefactor the paper might have
    written.  This version DERIVES the 2 pi / 3 from the two ingredients it
    is actually built from:

      (i) grad^2 (1/r) = -4 pi delta^3(r).  Verified, not quoted, via the
          divergence theorem: the flux of grad(1/r) through a sphere of any
          radius R is a radius-independent -4 pi.
      (ii) the spherical average <r'_i r'_j> = (<r^2>/3) delta_ij, so the
           second-order term of the Taylor expansion of the convolved
           Coulomb potential is (1/2)(<r^2>/3) grad^2 (1/r).

    Chaining: V = -Z int rho(r') / |r - r'| d^3r'
                = -Z/r - Z (<r^2>/6) grad^2(1/r) + ...
                = -Z/r + (2 pi/3) Z <r^2> delta^3(r) + ...
    so the sign is + and the coefficient is 2 pi / 3.
    """
    R, th = sp.symbols('R theta', positive=True)

    # (i) flux of grad(1/r) through a sphere of radius R (divergence theorem).
    #     grad(1/r) . n_hat = d/dr (1/r) |_{r=R} = -1/R^2 ; dS = R^2 sin(th) dth dphi
    flux = sp.integrate(
        sp.integrate((-1 / R ** 2) * R ** 2 * sp.sin(th), (th, 0, sp.pi)),
        (sp.Symbol('phi'), 0, 2 * sp.pi),
    )
    flux = sp.simplify(flux)
    assert flux == -4 * sp.pi, (
        f"flux of grad(1/r) = {flux}, expected -4 pi (so grad^2(1/r) = -4 pi delta^3)"
    )
    assert sp.diff(flux, R) == 0, (
        "flux must be radius-independent for the delta identification to hold"
    )
    lap_inv_r_delta_coeff = flux          # = -4 pi

    # (ii) spherical average of r'_i r'_j over the unit sphere, weighted by
    #      the radial second moment: coefficient of delta_ij is 1/3.  Derive
    #      the 1/3 from <n_z^2> = (1/2) int_0^pi cos^2(th) sin(th) dth.
    iso = sp.simplify(sp.Rational(1, 2)
                      * sp.integrate(sp.cos(th) ** 2 * sp.sin(th), (th, 0, sp.pi)))
    assert iso == sp.Rational(1, 3), f"<n_z^2> = {iso}, expected 1/3"

    # Chain: second-order Taylor term = (1/2) * iso * <r^2> * grad^2(1/r)
    # so Delta V = -Z * (1/2) * iso * <r^2> * (-4 pi) delta^3(r).
    r_sq = sp.Symbol('r_E_sq', positive=True)
    Z = sp.Symbol('Z', positive=True)
    delta_V_coeff = sp.simplify(
        -Z * sp.Rational(1, 2) * iso * r_sq * lap_inv_r_delta_coeff
    )
    expected = sp.Rational(2, 3) * sp.pi * Z * r_sq
    assert sp.simplify(delta_V_coeff - expected) == 0, (
        f"derived Foldy/Friar contact coefficient = {delta_V_coeff}, "
        f"expected +(2 pi/3) Z <r^2> = {expected}"
    )

    # Sign check: the correction is REPULSIVE (positive) -- the paper writes
    # a leading '+'.  A finite-size nucleus binds an s electron less.
    assert sp.simplify(delta_V_coeff / (Z * r_sq)) > 0


def test_paper34_III17_hydrogenic_1s_contact_density():
    """Paper 34 §III.17 cross-check: hydrogenic 1s contact density at
    origin is |psi_1s(0)|^2 = Z^3 / pi (standard result; rational over
    pi, ring-preserving over Q(alpha) when combined with the Foldy/Friar
    pi prefactor).

    The Foldy/Friar contact term Delta E = (2 pi / 3) Z alpha
    <r^2>_E |psi(0)|^2 then evaluates to
    (2 pi / 3) Z alpha <r^2>_E * Z^3 / pi = (2/3) Z^4 alpha <r^2>_E,
    which is the canonical Lamb-shift r_p contribution -- the pi
    cancels, leaving a rational coefficient (the 'ring-preserving over
    Q(alpha)' claim).
    """
    # 2026-08-28 audit: |psi_1s(0)|^2 is now DERIVED from the symbolically
    # normalized hydrogenic orbital (see _hydrogenic_density_at_origin), not
    # asserted.  The previous body wrote `psi_1s_sq_origin = Z**3/sp.pi` and
    # then "verified" that (2 pi/3)*Z*alpha*<r^2>*(Z^3/pi) = (2/3)Z^4 alpha
    # <r^2> -- true by construction for ANY value of |psi(0)|^2 of the form
    # (rational)/pi, so it tested nothing about the hydrogenic density.
    Z = _Z_SYM
    psi_1s_sq_origin = _hydrogenic_density_at_origin(1, 0)
    assert sp.simplify(psi_1s_sq_origin - Z ** 3 / sp.pi) == 0, (
        f"derived |psi_1s(0)|^2 = {psi_1s_sq_origin}, expected Z^3/pi"
    )
    # And the higher-n densities the paper's n-scaling relies on.
    for n in (2, 3):
        assert sp.simplify(_hydrogenic_density_at_origin(n, 0)
                           - Z ** 3 / (sp.pi * n ** 3)) == 0

    # Foldy/Friar evaluated at 1s
    r_sq_E = sp.symbols('r_E_sq', positive=True)
    alpha = sp.symbols('alpha', positive=True)
    contact_shift = (sp.Rational(2, 3) * sp.pi * Z * alpha * r_sq_E) * psi_1s_sq_origin

    # pi must cancel
    contact_simpl = sp.simplify(contact_shift)
    # see the note above: free_symbols never contains sp.pi; use .has()
    assert not contact_simpl.has(sp.pi), (
        f"pi did not cancel in Foldy/Friar 1s contact: {contact_simpl}"
    )

    # Expected coefficient: (2/3) Z^4 alpha <r^2>_E
    expected = sp.Rational(2, 3) * Z ** 4 * alpha * r_sq_E
    diff = sp.simplify(contact_simpl - expected)
    assert diff == 0, (
        f"Foldy/Friar 1s evaluation: {contact_simpl} != (2/3) Z^4 alpha r^2_E "
        f"= {expected}, diff = {diff}"
    )


# ----------------------------------------------------------------------------
# §III.18 Zemach: A_hf (1 - 2 Z alpha m_e r_Z + O(r_Z^2)) leading-order
# ----------------------------------------------------------------------------

def test_paper34_III18_zemach_leading_order_linear_in_rZ():
    """Paper 34 §III.18 (sec:proj_magnetization_density): leading-order
    Zemach correction is A_hf_contact (1 - 2 Z alpha m_e r_Z + O(r_Z^2)).

    Test: the production module's regression scales linearly with r_Z
    at leading order. The delta_LO_ppm at r_Z=1.045 fm should be ~2x
    the delta_LO_ppm at r_Z=0.5225 fm (factor of 2 in r_Z -> factor of 2
    in correction at LO).
    """
    from geovac.magnetization_density import hydrogen_zemach_eides_leading_order

    key = "delta_LO_ppm"
    rz_panel = [5.225e-6, 1.045e-5, 2.09e-5, 3.0e-5]
    vals = {}
    for rz in rz_panel:
        res = hydrogen_zemach_eides_leading_order(r_Z_bohr=rz, profile="gaussian")
        assert key in res, f"Expected key {key} in result, got {list(res.keys())}"
        vals[rz] = res[key]

    # 2026-08-28 audit: the old guard was |ratio - 2| < 0.05 at a single
    # pair of r_Z values.  Production is EXACTLY linear (delta_LO = -2 Z m_e
    # M_1 with M_1 = r_Z), so a 5% window is ~14 orders of magnitude looser
    # than the truth and would not catch a genuine quadratic leak.  Tighten
    # to the exact statement: delta_LO / r_Z is constant to double precision
    # across a 6x span of r_Z.
    slopes = [vals[rz] / rz for rz in rz_panel]
    spread = max(slopes) - min(slopes)
    assert abs(spread / slopes[0]) < 1e-14, (
        f"Zemach LO not exactly linear in r_Z: delta_LO/r_Z = {slopes} "
        f"(relative spread {abs(spread / slopes[0]):.3e})"
    )

    # The slope is the paper's -2 Z m_e coefficient, read off production
    # rather than hardcoded: doubling the lepton mass must double it.
    res_m1 = hydrogen_zemach_eides_leading_order(r_Z_bohr=1.045e-5, lepton_mass=1.0)
    res_m2 = hydrogen_zemach_eides_leading_order(r_Z_bohr=1.045e-5, lepton_mass=2.0)
    assert math.isclose(res_m2[key] / res_m1[key], 2.0, rel_tol=1e-12), (
        f"delta_LO not linear in lepton mass: ratio = {res_m2[key] / res_m1[key]}"
    )

    # Sign: the Zemach correction reduces the hyperfine splitting.
    assert all(v < 0 for v in vals.values()), (
        f"Zemach LO should be negative (binding reduction), got {vals}"
    )


def test_paper34_III18_zemach_profile_independence_leading_order():
    """Paper 34 §III.18: 'profile-dependence at sub-leading order admits
    Gaussian, exponential, and dipole forms with structurally identical
    leading -2 Z alpha m_e r_Z behaviour (profile independence at
    leading order verified in Sprint HF Track 4 / Sprint MH Track C).'

    Verifies Gaussian vs exponential profile agreement at leading order.
    """
    from geovac.magnetization_density import hydrogen_zemach_eides_leading_order

    rz = 1.045e-5  # bohr
    res_gauss = hydrogen_zemach_eides_leading_order(r_Z_bohr=rz, profile="gaussian")
    res_exp = hydrogen_zemach_eides_leading_order(r_Z_bohr=rz, profile="exponential")

    key = "delta_LO_ppm"
    delta_gauss = res_gauss[key]
    delta_exp = res_exp[key]

    # 2026-08-28 audit -- UPGRADE.  The old guard allowed a 5% profile leak.
    # Production is BIT-EXACT: delta_LO depends on rho_M only through its
    # FIRST moment M_1 = r_Z, which the two profiles share by calibration.
    # Assert the bit-exact statement, which is stronger than the paper's
    # "structurally identical leading behaviour".
    assert delta_gauss == delta_exp, (
        f"Profile independence at LO is not bit-exact: gaussian = "
        f"{delta_gauss!r}, exponential = {delta_exp!r}"
    )

    # NON-TAUTOLOGY GUARD (the reason this is not a false positive): the
    # `profile` argument must actually reach the density.  The two profiles
    # agree on M_1 by construction but DISAGREE on M_2 by ~13%, so the
    # bit-exactness above is a real leading-order cancellation and not the
    # signature of an ignored keyword.
    m2_g = res_gauss['rho_M_moments']['M_2']
    m2_e = res_exp['rho_M_moments']['M_2']
    assert m2_g != m2_e, (
        "guard failed: the two profiles produce identical second moments, so "
        "the `profile` argument may be ignored and the LO agreement vacuous"
    )
    assert abs(m2_e / m2_g - 1.0) > 0.05, (
        f"guard failed: M_2 differs by only {abs(m2_e/m2_g - 1.0):.2%}; the "
        "profiles are too close to certify a genuine LO cancellation"
    )
    # M_1 is the shared calibration -- this is WHAT makes LO profile-blind.
    assert res_gauss['rho_M_moments']['M_1'] == res_exp['rho_M_moments']['M_1'] == rz

    # And the complementary control: a profile with a DIFFERENT M_1 must give
    # a different LO value (the 'delta' profile is the point limit, M_1 = 0).
    res_delta = hydrogen_zemach_eides_leading_order(r_Z_bohr=rz, profile="delta")
    assert res_delta['rho_M_moments']['M_1'] == 0.0
    assert abs(res_delta[key]) < 1e-12 < abs(delta_gauss), (
        f"delta-profile LO = {res_delta[key]}, expected ~0 (M_1 = 0)"
    )


# ----------------------------------------------------------------------------
# §III.19 Nuclear tensor multipole: rank-2 Wigner 3j selection rules
# ----------------------------------------------------------------------------

def test_paper34_III19_rank2_multipole_triangle_inequality():
    r"""Paper 34 §III.19 (sec:proj_tensor_multipole): rank-2 quadrupole
    coupling H_Q = -e Q_N T^{(2)}_{ij} (\partial_i E_j) / 6 satisfies
    Wigner 3j triangle inequality at every step ('ring-preserving over
    Q at the angular level').

    Test: rank-2 spherical tensor selection rule -- non-zero matrix
    element between |l_1, m_1> and |l_2, m_2> requires |l_1 - l_2| <= 2
    <= l_1 + l_2 and m_1 - m_2 = M with -2 <= M <= 2.
    """
    from sympy.physics.wigner import wigner_3j

    # Rank-2 multipole: triangle inequality requires |l1 - l2| <= 2 <= l1 + l2

    # Outside upper bound -> vanishes
    # s -> g (l=0 to l=4) under rank-2: 0 + 4 = 4 >> 2, vanishes
    assert wigner_3j(0, 2, 4, 0, 0, 0) == 0
    # f -> s (l=3 to l=0) under rank-2: |3-0| = 3 > 2, vanishes
    assert wigner_3j(3, 2, 0, 0, 0, 0) == 0

    # Inside triangle, non-zero
    # s -> d (l=0 to l=2): |0-2| = 2 = 2, on boundary
    assert wigner_3j(0, 2, 2, 0, 0, 0) != 0
    # p -> p (l=1 to l=1) under rank-2: 0 <= 2 <= 2, allowed
    assert wigner_3j(1, 2, 1, 0, 0, 0) != 0
    # d -> d under rank-2
    assert wigner_3j(2, 2, 2, 0, 0, 0) != 0


def test_paper34_III19_rank2_M_quantum_number_conservation():
    """Paper 34 §III.19: rank-2 spherical tensor M_q in {-2, -1, 0, 1, 2}
    couples m_1 - m_2 = M (the magnetic-quantum-number conservation
    that gives the quadrupole hyperfine M selection rule).

    Test: wigner_3j(l1, 2, l2, -m1, M, m2) is zero unless m1 - m2 = M.
    """
    from sympy.physics.wigner import wigner_3j

    # Take l1 = l2 = 2 (d-d coupling under rank-2)
    # Non-zero requires m1 + M + m2 = 0 -> m2 = -m1 - M -> m1 - m2 = ... wait
    # wigner_3j(l1, j, l2, m1, M, m2) is non-zero only if m1 + M + m2 = 0
    # So m1 - (-m2) = m1 + m2 = -M -> equivalent statement.

    # Conservation: m1 + M + m2 = 0
    assert wigner_3j(2, 2, 2, 1, -1, 0) != 0   # 1 + (-1) + 0 = 0 OK
    assert wigner_3j(2, 2, 2, 1, 0, -1) != 0   # 1 + 0 + (-1) = 0 OK
    assert wigner_3j(2, 2, 2, 2, 0, -2) != 0   # 2 + 0 + (-2) = 0 OK
    # Violation
    assert wigner_3j(2, 2, 2, 1, 1, 1) == 0    # 1 + 1 + 1 = 3 != 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
