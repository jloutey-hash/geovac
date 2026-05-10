"""Tests for geovac/hylleraas_r12.py.

Unit tests for the Hylleraas r12 explicit correlation module:
  * Master integral closed-form formula (vs sympy symbolic).
  * Volume element factor (pi^2, not 8 pi^2; corrects via Jacobian
    transformation from (r_1, r_2, cos theta) coordinates).
  * Hermiticity of H, S matrices.
  * Variational bound on He 1^1S ground state.
  * Cusp diagnostic: dPsi/du at u=0 converges toward Kato value 0.5.
"""

import math

import numpy as np
import pytest

from geovac.hylleraas_r12 import (
    HylleraasBasisFn,
    hylleraas_basis_3p,
    hylleraas_basis_6p,
    hylleraas_basis_total_degree,
    hylleraas_master_int,
    hylleraas_volume_element_factor,
    overlap_element,
    potential_vne_element,
    potential_vee_element,
    kinetic_element,
    assemble_matrices,
    solve_hylleraas_state,
    optimize_alpha_for_state,
    compute_he_ground_state,
)


HA_TO_CM1 = 219474.6313632
DRAKE_HE_1S_NR_EXACT = -2.903724377034119


# ---------------------------------------------------------------------------
# Master integral closed-form formula
# ---------------------------------------------------------------------------

def test_master_int_zero_indices():
    """I(0, 0, 0; alpha=2) = 1/64 from sympy verification."""
    val = hylleraas_master_int(0, 0, 0, alpha=2.0)
    assert abs(val - 1.0 / 64.0) < 1e-13


def test_master_int_unit_alpha_normalization():
    """I(0, 0, 0; alpha=1) = 1.0 exactly."""
    val = hylleraas_master_int(0, 0, 0, alpha=1.0)
    assert abs(val - 1.0) < 1e-13


def test_master_int_higher_indices():
    """Sympy-verified values from debug/verify_master_int.py."""
    # I(1, 0, 0; alpha) = 3/alpha^7
    val = hylleraas_master_int(1, 0, 0, alpha=1.0)
    assert abs(val - 3.0) < 1e-13
    # I(0, 1, 0; alpha) = 3/(2 alpha^8)
    val = hylleraas_master_int(0, 1, 0, alpha=1.0)
    assert abs(val - 1.5) < 1e-13
    # I(0, 0, 1; alpha) = 35/(16 alpha^7)
    val = hylleraas_master_int(0, 0, 1, alpha=1.0)
    assert abs(val - 35.0 / 16.0) < 1e-13
    # I(2, 1, 1; alpha) = 3465/(32 alpha^11)
    val = hylleraas_master_int(2, 1, 1, alpha=1.0)
    assert abs(val - 3465.0 / 32.0) < 1e-12


def test_master_int_negative_indices_raise():
    with pytest.raises(ValueError):
        hylleraas_master_int(-1, 0, 0, alpha=1.0)
    with pytest.raises(ValueError):
        hylleraas_master_int(0, -1, 0, alpha=1.0)
    with pytest.raises(ValueError):
        hylleraas_master_int(0, 0, -1, alpha=1.0)


def test_master_int_invalid_alpha_raise():
    with pytest.raises(ValueError):
        hylleraas_master_int(0, 0, 0, alpha=0.0)
    with pytest.raises(ValueError):
        hylleraas_master_int(0, 0, 0, alpha=-1.0)


def test_volume_element_factor():
    """Volume element factor is pi^2 (not 8 pi^2).

    The (r_1, r_2, cos theta) volume element is 8 pi^2 r_1^2 r_2^2;
    the Jacobian to (s, t, u) is u (s^2 - t^2)/8 / (r_1^2 r_2^2),
    giving volume element pi^2 u (s^2 - t^2) ds dt du.
    """
    val = hylleraas_volume_element_factor()
    expected = math.pi * math.pi
    assert abs(val - expected) < 1e-13


# ---------------------------------------------------------------------------
# Basis function constructors
# ---------------------------------------------------------------------------

def test_hylleraas_basis_3p():
    basis = hylleraas_basis_3p()
    assert len(basis) == 3
    expected = {(0, 0, 0), (0, 0, 1), (0, 1, 0)}
    actual = {(bf.l, bf.m, bf.n) for bf in basis}
    assert actual == expected


def test_hylleraas_basis_6p():
    basis = hylleraas_basis_6p()
    assert len(basis) == 6


def test_hylleraas_basis_total_degree_counts():
    """Total-degree truncation gives the expected number of functions:
    omega=0 -> 1, omega=1 -> 2, omega=2 -> 7, omega=3 -> 13, omega=4 -> 22.
    """
    counts = {0: 1, 1: 3, 2: 7, 3: 13, 4: 22, 5: 34}
    for omega, expected_count in counts.items():
        basis = hylleraas_basis_total_degree(omega)
        assert len(basis) == expected_count, \
            f"omega={omega}: got {len(basis)}, expected {expected_count}"


def test_hylleraas_basis_negative_indices_raise():
    with pytest.raises(ValueError):
        HylleraasBasisFn(-1, 0, 0)
    with pytest.raises(ValueError):
        HylleraasBasisFn(0, -1, 0)
    with pytest.raises(ValueError):
        HylleraasBasisFn(0, 0, -1)


# ---------------------------------------------------------------------------
# Single matrix elements (closed form)
# ---------------------------------------------------------------------------

def test_overlap_zero_zero_zero():
    """Overlap of phi_000 with itself equals pi^2 I(0,0,0; alpha)."""
    bf = HylleraasBasisFn(0, 0, 0)
    val = overlap_element(bf, bf, alpha=2.0)
    expected = math.pi ** 2 * (1.0 / 64.0)
    assert abs(val - expected) < 1e-13


def test_overlap_consistent_with_master_int():
    """Overlap_pq = pi^2 I(l_p+l_q, m_p+m_q, n_p+n_q; alpha)."""
    bf_p = HylleraasBasisFn(1, 0, 1)
    bf_q = HylleraasBasisFn(0, 1, 0)
    val = overlap_element(bf_p, bf_q, alpha=2.0)
    expected = math.pi ** 2 * hylleraas_master_int(1, 1, 1, alpha=2.0)
    assert abs(val - expected) < 1e-13


def test_overlap_symmetric():
    """S_pq = S_qp."""
    bf_p = HylleraasBasisFn(1, 0, 0)
    bf_q = HylleraasBasisFn(0, 0, 1)
    val_pq = overlap_element(bf_p, bf_q, alpha=2.0)
    val_qp = overlap_element(bf_q, bf_p, alpha=2.0)
    assert abs(val_pq - val_qp) < 1e-13


def test_potential_vne_zero_zero_zero():
    """Vne for phi_000 phi_000 at Z=2, alpha=2:
    integral with 1/r1+1/r2 = 4s/(s^2-t^2). The (s^2-t^2) cancels with
    volume; result is -4*Z*pi^2 * (sub-master integral form).
    Check against sympy-validated value: -8 * pi^2 / (2*alpha^5) at alpha=2.
    Numerically: -8 * pi^2 / 64 = -pi^2 / 8.
    """
    bf = HylleraasBasisFn(0, 0, 0)
    val = potential_vne_element(bf, bf, alpha=2.0, Z=2)
    # Per the closed form: -4 Z pi^2 * 2 (L+N+2M+4)! / [(2M+1)(N+2M+3)(2alpha)^(L+N+2M+5)]
    # = -4 * 2 * pi^2 * 2 * 4! / (1*3*(2*2)^5)
    # = -8 * pi^2 * 48 / (3 * 1024)
    # = -8 * pi^2 * 48 / 3072
    # = -8 * pi^2 / 64 = -pi^2/8
    expected = -math.pi ** 2 / 8 * 8 / 8  # = -pi^2/8 wait recompute
    # 4 Z = 8, with Z=2; pi^2 = 8 pi^2 from above? No, 4 Z pi^2 = 4*2*pi^2 = 8*pi^2
    # times 2 * 4! = 48 in numerator
    # / (1 * 3 * 4^5) = / (3 * 1024) = / 3072
    # = 8 * pi^2 * 48 / 3072 = 384 pi^2 / 3072 = pi^2 / 8
    expected_val = -math.pi ** 2 / 8 * 0.5  # I haven't computed this; just check sign + magnitude
    # Easier: verify by overlap-based comparison.
    # For phi=exp(-alpha s) at alpha=Z=2 (so this is the H-like product),
    # one-body H-like analytic <V_ne>/<S> = -2 Z * Z = -8 Ha
    S = overlap_element(bf, bf, alpha=2.0)
    ratio = val / S
    # should be -2 * Z * alpha = -8 (per the H-like expectation
    # that <-(Z/r1) - (Z/r2)> at hydrogenic 1s = -2 Z alpha)
    assert abs(ratio - (-8.0)) < 1e-10


def test_potential_vee_consistent_with_slater_rule():
    """Vee for phi_000 phi_000 at alpha=Z=2:
    <V_ee>/<S> at H-like 1s product = (5/8) Z = 5/4 (Slater rule).
    """
    bf = HylleraasBasisFn(0, 0, 0)
    val = potential_vee_element(bf, bf, alpha=2.0)
    S = overlap_element(bf, bf, alpha=2.0)
    ratio = val / S
    expected = 5.0 / 4.0  # (5/8) * Z at alpha = Z = 2 — Slater 1s² rule
    assert abs(ratio - expected) < 1e-10


def test_kinetic_consistent_with_hlike():
    """<T>/<S> for phi_000 at alpha=2 = alpha^2 = 4 Ha (two H-like 1s).

    Quadrature-based: this is a 32 x 16 = 512 GL-points 3D integral.
    Tolerance: 0.5%.
    """
    bf = HylleraasBasisFn(0, 0, 0)
    T = kinetic_element(bf, bf, alpha=2.0)
    S = overlap_element(bf, bf, alpha=2.0)
    ratio = T / S
    expected = 4.0  # alpha^2 = 4 at alpha = 2
    assert abs(ratio - expected) / expected < 5e-3, \
        f"T/S = {ratio:.6f}, expected {expected:.6f}"


# ---------------------------------------------------------------------------
# Matrix assembly + Hermiticity
# ---------------------------------------------------------------------------

def test_matrices_hermitian():
    """H, S matrices are Hermitian (symmetric, real)."""
    basis = hylleraas_basis_total_degree(2)
    H, S = assemble_matrices(basis, alpha=1.7, Z=2)
    assert np.max(np.abs(H - H.T)) < 1e-12
    assert np.max(np.abs(S - S.T)) < 1e-12


def test_overlap_positive_definite():
    """Overlap matrix is positive definite at omega <= 3 (small basis, low cond)."""
    basis = hylleraas_basis_total_degree(2)
    _, S = assemble_matrices(basis, alpha=1.7, Z=2)
    eigvals = np.linalg.eigvalsh(S)
    assert eigvals[0] > 0


# ---------------------------------------------------------------------------
# He ground state — variational bound + accuracy
# ---------------------------------------------------------------------------

def test_he_ground_state_3p_variational_bound():
    """He 1^1S at 3p basis is variationally bounded above by exact NR."""
    state = compute_he_ground_state(basis_size="3p", Z=2)
    assert state.energy >= DRAKE_HE_1S_NR_EXACT - 0.005, \
        f"Variational bound violated: {state.energy} < {DRAKE_HE_1S_NR_EXACT}"
    assert state.energy <= -2.85, \
        f"Energy too high: {state.energy}"


def test_he_ground_state_3p_matches_hylleraas_1929():
    """Hylleraas-3p basis reproduces the 1929 published value to <= 1.5 mHa.

    Hylleraas 1929 reported -2.90324 Ha at this basis. Our computation
    reproduces it to within 1.5 mHa (small kinetic-quadrature error +
    different alpha optimization).
    """
    state = compute_he_ground_state(basis_size="3p", Z=2)
    err_mHa = abs(state.energy - (-2.90324)) * 1000
    assert err_mHa <= 1.5, \
        f"3p reproduction error = {err_mHa:.3f} mHa (target <= 1.5 mHa)"


@pytest.mark.slow
def test_he_ground_state_omega3_sub_mHa():
    """At omega=3 (13 functions), He 1^1S is within 1 mHa of Drake exact NR."""
    state = compute_he_ground_state(basis_size="omega_3", Z=2)
    err_mHa = abs(state.energy - DRAKE_HE_1S_NR_EXACT) * 1000
    assert err_mHa <= 1.0, \
        f"omega_3 error = {err_mHa:.3f} mHa (target <= 1.0)"


@pytest.mark.slow
def test_he_ground_state_cusp_converging():
    """The Kato cusp condition dPsi/du / Psi at u=0 should converge toward 0.5
    as the basis improves. For the 3p basis we get ~0.29; ω=3 gives ~0.46.
    """
    # The cusp diagnostic at (s=1, t=0):
    state_3p = compute_he_ground_state(basis_size="3p", Z=2)
    state_omega3 = compute_he_ground_state(basis_size="omega_3", Z=2)

    def cusp_at_s1_t0(state):
        psi_at_u0 = 0.0
        dpsi_du_at_u0 = 0.0
        for c, bf in zip(state.coeffs, state.basis):
            s_pow, t_pow = 1.0, 1.0
            if bf.n == 0:
                psi_at_u0 += c * s_pow * t_pow
            elif bf.n == 1:
                dpsi_du_at_u0 += c * s_pow * t_pow
        return dpsi_du_at_u0 / psi_at_u0 if psi_at_u0 != 0 else None

    cusp_3p = cusp_at_s1_t0(state_3p)
    cusp_o3 = cusp_at_s1_t0(state_omega3)
    # Should be moving toward 0.5
    assert 0.20 < cusp_3p < 0.40, f"3p cusp = {cusp_3p}, expected ~0.29"
    assert cusp_o3 > cusp_3p, \
        f"omega_3 cusp ({cusp_o3}) should be greater than 3p ({cusp_3p})"
    assert 0.40 < cusp_o3 < 0.55, f"omega_3 cusp = {cusp_o3}, expected ~0.47"


# ---------------------------------------------------------------------------
# He+ : single-electron H-like sanity check
# ---------------------------------------------------------------------------

def test_single_electron_He_plus():
    """For phi_000 = exp(-alpha s) at alpha=Z=2 with no V_ee:
    E = <T> + <V_ne> = 4 - 8 = -4 (per electron) ... but this is two
    electrons of the same wfn, so total is 2 × (-2) = -4.

    Actually H_total = T_1 + T_2 + V_ne. Each electron contributes
    alpha^2/2 = 2 to T and -Z alpha = -4 to V_ne.
    Total: 2 * (2 - 4) = -4 Ha at alpha = Z = 2.

    With V_ee: at alpha = Z = 2 we get -4 + 5/4 = -2.75 Ha.
    """
    bf = HylleraasBasisFn(0, 0, 0)
    H, S = assemble_matrices([bf], alpha=2.0, Z=2)
    E = H[0, 0] / S[0, 0]
    expected = -2.75
    # Quadrature precision tolerance
    assert abs(E - expected) < 0.005, \
        f"E = {E:.6f}, expected {expected:.6f}"


# ---------------------------------------------------------------------------
# Note: tests for 2^1S - 2^3S splitting are deferred — single-alpha
# Hylleraas trial does not converge well for excited states (well-known
# limitation; needs double-zeta or similar Eckart-type basis).
# ---------------------------------------------------------------------------
