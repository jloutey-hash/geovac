"""Tests for the Sprint L2-D Connes axiom audit at signature (3, 1).

Tests cover:
  - BBB Table 1 sign determination at (m, n) = (4, 6).
  - 4-spinor (bare Cl(3,1)) verification of BBB signs via U_4 = i*gamma^2.
  - J_L^2 = +I bit-exact at finite n_max (LOAD-BEARING).
  - {J_L, gamma^5} = 0 bit-exact (BBB eps'' = -1).
  - {J_L, gamma^0} = 0 bit-exact (BBB eps*kappa = -1).
  - J_L D_L = +D_L J_L bit-exact (BBB universal Sec 5(v)).
  - chi D + D chi residual reported (NOT zero on truthful D_GV --
    documented structural finding).
  - Order-zero / order-one sampled residuals (expected ~5-10%
    finite-resolution per Paper 32 §IV).
  - M3 trivialization both readings (n_fock-parity vs chirality-pairing).
"""

from __future__ import annotations

import numpy as np
import pytest

from geovac.connes_axiom_audit_31 import (
    apply_J_L,
    apply_J_L_to_operator,
    audit_at_4_6,
    bbb_signs_at_4_6,
    J_L_squared_matrix,
    lorentzian_J_spatial_matrix,
    lorentzian_real_structure_matrix,
    m3_vertex_parity_chirality_pairing,
    m3_vertex_parity_sum_lorentzian_chirality_symmetric,
    m3_vertex_parity_sum_riemannian,
    verify_bbb_signs_at_4_spinor_level,
    verify_chi_D_anticommutes,
    verify_J_L_anticommutes_chi,
    verify_J_L_anticommutes_eta,
    verify_J_L_D_relation,
    verify_J_L_squared,
    verify_m3_trivialization,
    verify_order_one,
    verify_order_zero,
)
from geovac.full_dirac_operator_system import (
    FullDiracTruncatedOperatorSystem,
    full_dirac_basis,
)
from geovac.krein_space_construction import KreinSpace
from geovac.lorentzian_dirac import (
    lorentzian_dirac_matrix,
    spacetime_chirality_lifted,
)


# ---------------------------------------------------------------------------
# BBB Table 1 sign determination
# ---------------------------------------------------------------------------


def test_bbb_signs_at_4_6_primitive_values():
    """BBB Table 1 primitive signs at (m, n) = (4, 6)."""
    s = bbb_signs_at_4_6()
    assert s['eps'] == +1
    assert s['eps_pp'] == -1
    assert s['kappa'] == -1
    assert s['kappa_pp'] == +1


def test_bbb_signs_at_4_6_derived_relations():
    """Derived BBB relation signs at (m, n) = (4, 6)."""
    s = bbb_signs_at_4_6()
    assert s['J_squared_sign'] == +1                  # J^2 = +I
    assert s['J_chi_sign'] == -1                       # Jχ = -χJ
    assert s['J_eta_sign'] == s['eps'] * s['kappa']   # = -1
    assert s['J_eta_sign'] == -1
    assert s['eta_chi_sign'] == s['eps_pp'] * s['kappa_pp']  # = -1
    assert s['eta_chi_sign'] == -1


def test_bbb_signs_at_4_6_universal_relations():
    """BBB Sec 5 item (v) universal (signature-independent) relations."""
    s = bbb_signs_at_4_6()
    assert s['J_D_sign'] == +1     # J D = +D J  universal
    assert s['chi_D_sign'] == -1   # chi D = -D chi  universal


def test_bbb_signs_at_4_6_signature_label():
    """Signature and (m, n) labels."""
    s = bbb_signs_at_4_6()
    assert '(3, 1)' in s['signature']
    assert '(4, 6)' in s['mn']
    assert '1611.07062' in s['source']


# ---------------------------------------------------------------------------
# 4-spinor level: bit-exact verification via i*gamma^2
# ---------------------------------------------------------------------------


def test_bbb_signs_at_4_spinor_level_all_pass():
    """U_4 = i*gamma^2 satisfies all four BBB (4, 6) signs bit-exact."""
    ok, residuals = verify_bbb_signs_at_4_spinor_level(tol=1e-14)
    assert ok
    # All four residuals should be bit-exact zero
    for key, res in residuals.items():
        assert res == 0.0, f"{key} residual {res} != 0"


def test_bbb_signs_at_4_spinor_J_squared_bit_exact():
    """(i*gamma^2) conj(i*gamma^2) = +I_4 bit-exact."""
    _, residuals = verify_bbb_signs_at_4_spinor_level()
    assert residuals['J_squared_residual'] == 0.0


def test_bbb_signs_at_4_spinor_J_chi_anticommutes_bit_exact():
    """(i*gamma^2) gamma^5 + gamma^5 (i*gamma^2) = 0 bit-exact."""
    _, residuals = verify_bbb_signs_at_4_spinor_level()
    assert residuals['J_chi_anticomm_residual'] == 0.0


def test_bbb_signs_at_4_spinor_J_eta_anticommutes_bit_exact():
    """(i*gamma^2) gamma^0 + gamma^0 (i*gamma^2) = 0 bit-exact."""
    _, residuals = verify_bbb_signs_at_4_spinor_level()
    assert residuals['J_eta_anticomm_residual'] == 0.0


# ---------------------------------------------------------------------------
# J_L construction on H_GV (spatial U_L is bit-exact unitary involution-like)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_lorentzian_J_spatial_unitary(n_max):
    """U_L_spatial U_L_spatial^dagger = I bit-exact at each n_max."""
    basis = full_dirac_basis(n_max)
    U = lorentzian_J_spatial_matrix(basis)
    dim = U.shape[0]
    UUH = U @ U.conj().T
    assert np.linalg.norm(UUH - np.eye(dim)) < 1e-12


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_lorentzian_J_spatial_real_entries(n_max):
    """U_L_spatial entries are real (chiral basis property)."""
    basis = full_dirac_basis(n_max)
    U = lorentzian_J_spatial_matrix(basis)
    assert np.allclose(U.imag, 0.0, atol=1e-14)


# ---------------------------------------------------------------------------
# Full Krein-space J_L on K = H_GV (x) C^{N_t}
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max,N_t", [(1, 1), (1, 11), (2, 1), (2, 11), (3, 1)])
def test_lorentzian_real_structure_dim(n_max, N_t):
    """U_L has the correct Krein-space dimension."""
    krein = KreinSpace(n_max=n_max, N_t=N_t)
    U_L = lorentzian_real_structure_matrix(krein)
    assert U_L.shape == (krein.dim, krein.dim)


@pytest.mark.parametrize("n_max,N_t", [(1, 1), (1, 11), (2, 1), (2, 11), (3, 1)])
def test_lorentzian_real_structure_unitary(n_max, N_t):
    """U_L U_L^dagger = I bit-exact."""
    krein = KreinSpace(n_max=n_max, N_t=N_t)
    U_L = lorentzian_real_structure_matrix(krein)
    UUH = U_L @ U_L.conj().T
    assert np.linalg.norm(UUH - np.eye(krein.dim)) < 1e-12


# ---------------------------------------------------------------------------
# L2D-FALS-1 LOAD-BEARING: J_L^2 = +I bit-exact at every panel cell
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max,N_t", [
    (1, 1), (1, 11), (1, 21),
    (2, 1), (2, 11), (2, 21),
    (3, 1), (3, 11), (3, 21),
])
def test_FALS_1_J_L_squared_plus_I_bit_exact(n_max, N_t):
    """L2D-FALS-1 (LOAD-BEARING): J_L^2 = +I at (m, n) = (4, 6)."""
    krein = KreinSpace(n_max=n_max, N_t=N_t)
    U_L = lorentzian_real_structure_matrix(krein)
    ok, residual = verify_J_L_squared(U_L, expected_sign=+1, tol=1e-12)
    assert ok, f"J_L^2 != +I at (n_max={n_max}, N_t={N_t}): residual={residual}"
    assert residual == 0.0, (
        f"J_L^2 = +I should be bit-exact, got residual={residual}"
    )


def test_FALS_1_J_L_squared_not_minus_I():
    """Sanity: J_L^2 != -I (distinguishes (4, 6) from KO-dim 3 J_GV)."""
    krein = KreinSpace(n_max=2, N_t=1)
    U_L = lorentzian_real_structure_matrix(krein)
    J2 = J_L_squared_matrix(U_L)
    # J^2 should be +I, NOT -I
    target_minus = -np.eye(krein.dim, dtype=np.complex128)
    residual_to_minus = float(np.linalg.norm(J2 - target_minus))
    # Far from -I (should be ~ sqrt(2 * dim) ~ 5.66 at dim=16)
    assert residual_to_minus > 1.0


# ---------------------------------------------------------------------------
# Axiom (ii): {J_L, gamma^5} = 0  (BBB eps'' = -1)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max,N_t", [
    (1, 1), (1, 11), (2, 1), (2, 11), (3, 1), (3, 11),
])
def test_axiom_ii_J_L_anticommutes_chi(n_max, N_t):
    """{J_L, gamma^5} = 0 bit-exact (BBB eps'' = -1 at (4, 6))."""
    krein = KreinSpace(n_max=n_max, N_t=N_t)
    U_L = lorentzian_real_structure_matrix(krein)
    gamma5 = spacetime_chirality_lifted(krein)
    ok, residual = verify_J_L_anticommutes_chi(U_L, gamma5, tol=1e-12)
    assert ok, f"{{J, chi}} != 0 at (n_max={n_max}, N_t={N_t}): residual={residual}"


# ---------------------------------------------------------------------------
# Axiom (iii): {J_L, gamma^0} = 0  (BBB eps*kappa = -1)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max,N_t", [
    (1, 1), (1, 11), (2, 1), (2, 11), (3, 1), (3, 11),
])
def test_axiom_iii_J_L_anticommutes_eta(n_max, N_t):
    """{J_L, gamma^0} = 0 bit-exact (BBB eps*kappa = -1 at (4, 6))."""
    krein = KreinSpace(n_max=n_max, N_t=N_t)
    U_L = lorentzian_real_structure_matrix(krein)
    eta = krein.J  # gamma^0 on the Krein space
    ok, residual = verify_J_L_anticommutes_eta(U_L, eta, tol=1e-12)
    assert ok, f"{{J, eta}} != 0 at (n_max={n_max}, N_t={N_t}): residual={residual}"


# ---------------------------------------------------------------------------
# Axiom (iv): J_L D_L = +D_L J_L  (BBB universal)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max,N_t", [
    (1, 1), (1, 11), (2, 1), (2, 11), (3, 1), (3, 11),
])
def test_axiom_iv_J_L_D_commutation(n_max, N_t):
    """J_L D_L = +D_L J_L bit-exact (BBB universal Sec 5(v) at (4, 6))."""
    krein = KreinSpace(n_max=n_max, N_t=N_t)
    U_L = lorentzian_real_structure_matrix(krein)
    D_L = lorentzian_dirac_matrix(krein)
    ok, residual = verify_J_L_D_relation(U_L, D_L, expected_sign=+1, tol=1e-12)
    assert ok, f"J D != +D J at (n_max={n_max}, N_t={N_t}): residual={residual}"


# ---------------------------------------------------------------------------
# Axiom (v): {gamma^5, D_L} = 0  (BBB universal)
# STRUCTURAL FINDING: fails on truthful D_GV (documented).
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max,N_t", [(1, 1), (2, 1), (3, 1)])
def test_axiom_v_chi_D_anticomm_fails_at_Nt1_on_truthful_DGV(n_max, N_t):
    """{gamma^5, D_L} != 0 at N_t = 1 on truthful D_GV (structural finding).

    Mechanism: at N_t = 1 the temporal part vanishes (d/dt = 0), so
    D_L = i * D_GV.  Since chirality = gamma^5 eigenvalue and D_GV is
    chirality-diagonal, {chi, D_GV} = 2 chi D_GV which is non-zero
    (norm = 2 |D_GV|_F).

    The non-zero residual is the documented STRUCTURAL FINDING for the
    L2-D memo: the BBB axiom chi D = -D chi does NOT hold with truthful
    chirality-diagonal D_GV.  Use offdiag CH for an axiom-compatible
    Dirac (R3 in the memo).
    """
    krein = KreinSpace(n_max=n_max, N_t=N_t)
    D_L = lorentzian_dirac_matrix(krein)
    gamma5 = spacetime_chirality_lifted(krein)
    ok, residual = verify_chi_D_anticommutes(gamma5, D_L, tol=1e-12)
    # We EXPECT this to fail (structural finding)
    assert not ok, "Unexpectedly passed: BBB axiom holds, structural finding inverted"
    assert residual > 1.0, (
        f"Expected large residual for structural finding, got {residual}"
    )


# ---------------------------------------------------------------------------
# Axioms (vi)-(vii): order-zero and order-one
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max,N_t", [(1, 1), (2, 1), (3, 1)])
def test_axiom_vi_order_zero_small_at_finite_resolution(n_max, N_t):
    """[a, J_L b J_L^{-1}] residual is bounded (Paper 32 §IV scope).

    Paper 32 §IV reports 5-20% residual on Riemannian side as finite-
    resolution artifact analogous to multiplicative-closure failure of
    the Connes-vS truncated operator system.  Lorentzian side should
    preserve this scope (residual <= 0.5 max as a generous bound).
    """
    krein = KreinSpace(n_max=n_max, N_t=N_t)
    U_L = lorentzian_real_structure_matrix(krein)
    I_t = np.eye(N_t, dtype=np.complex128)
    op_sys = FullDiracTruncatedOperatorSystem(n_max)
    A_basis = [
        np.kron(m, I_t)
        for m in op_sys.multiplier_matrices[:3]
    ]
    _, n_failures, max_res = verify_order_zero(U_L, A_basis, tol=1e-10)
    # Allow finite-resolution residual up to 0.5 (paper 32 §IV scope says ~5-20%)
    assert max_res < 0.5, (
        f"Order-zero residual {max_res} > 0.5 at (n_max={n_max}, N_t={N_t})"
    )


@pytest.mark.parametrize("n_max,N_t", [(1, 1), (2, 1), (3, 1)])
def test_axiom_vii_order_one_small_at_finite_resolution(n_max, N_t):
    """[[D_L, a], J_L b J_L^{-1}] residual is bounded (Paper 32 §IV scope)."""
    krein = KreinSpace(n_max=n_max, N_t=N_t)
    U_L = lorentzian_real_structure_matrix(krein)
    D_L = lorentzian_dirac_matrix(krein)
    I_t = np.eye(N_t, dtype=np.complex128)
    op_sys = FullDiracTruncatedOperatorSystem(n_max)
    A_basis = [
        np.kron(m, I_t)
        for m in op_sys.multiplier_matrices[:3]
    ]
    _, n_failures, max_res = verify_order_one(U_L, A_basis, D_L, tol=1e-10)
    assert max_res < 0.5, (
        f"Order-one residual {max_res} > 0.5 at (n_max={n_max}, N_t={N_t})"
    )


# ---------------------------------------------------------------------------
# Antilinear operations sanity
# ---------------------------------------------------------------------------


def test_apply_J_L_to_real_vector():
    """J_L on a real psi: J(psi) = U conj(psi) = U psi (for real psi)."""
    krein = KreinSpace(n_max=1, N_t=1)
    U_L = lorentzian_real_structure_matrix(krein)
    psi = np.zeros(krein.dim, dtype=np.complex128)
    psi[0] = 1.0
    Jpsi = apply_J_L(U_L, psi)
    # Should equal U_L[:, 0]
    assert np.allclose(Jpsi, U_L[:, 0])


def test_apply_J_L_squared_is_identity():
    """J_L^2 (psi) = psi for any psi at finite n_max (J^2 = +I)."""
    krein = KreinSpace(n_max=2, N_t=1)
    U_L = lorentzian_real_structure_matrix(krein)
    rng = np.random.default_rng(0)
    psi = (rng.standard_normal(krein.dim)
           + 1j * rng.standard_normal(krein.dim))
    J_psi = apply_J_L(U_L, psi)
    JJ_psi = apply_J_L(U_L, J_psi)
    assert np.allclose(JJ_psi, psi, atol=1e-12)


def test_apply_J_L_antilinear():
    """J_L(alpha * psi) = conj(alpha) * J_L(psi) for any complex alpha."""
    krein = KreinSpace(n_max=1, N_t=1)
    U_L = lorentzian_real_structure_matrix(krein)
    rng = np.random.default_rng(0)
    psi = (rng.standard_normal(krein.dim)
           + 1j * rng.standard_normal(krein.dim))
    alpha = 2.0 + 3.0j
    LHS = apply_J_L(U_L, alpha * psi)
    RHS = np.conj(alpha) * apply_J_L(U_L, psi)
    assert np.allclose(LHS, RHS, atol=1e-12)


def test_apply_J_L_to_operator_inversion():
    """J_L (J_L op J_L^-1) J_L^-1 = op for any op."""
    krein = KreinSpace(n_max=2, N_t=1)
    U_L = lorentzian_real_structure_matrix(krein)
    rng = np.random.default_rng(0)
    op = (rng.standard_normal((krein.dim, krein.dim))
          + 1j * rng.standard_normal((krein.dim, krein.dim)))
    J_op_Jinv = apply_J_L_to_operator(U_L, op)
    J_J_op_Jinv_Jinv = apply_J_L_to_operator(U_L, J_op_Jinv)
    # For J^2 = +I, J (J op J^-1) J^-1 = op
    assert np.allclose(J_J_op_Jinv_Jinv, op, atol=1e-12)


# ---------------------------------------------------------------------------
# M3 trivialization (Sprint L0 prediction)
# ---------------------------------------------------------------------------


def test_m3_riemannian_d_even_d_odd_nontrivial_at_nmax_2():
    """Riemannian-side D_even(4) - D_odd(4) != 0 at n_max = 2."""
    DeR, DoR, diffR = m3_vertex_parity_sum_riemannian(2)
    assert abs(diffR) > 0.1


def test_m3_lorentzian_chirality_symmetric_equals_riemannian():
    """Reading A (n_fock-parity) gives identical Riemannian and Lorentzian.

    On chirality-symmetric truncation, both chirality blocks contribute
    equally to the same parity bucket (Paper 28 §QED-vertex reading).
    The Lorentzian D_even - D_odd is the SAME as the Riemannian one
    because parity is chirality-independent.
    """
    for n_max in [1, 2, 3, 5]:
        DeR, DoR, diffR = m3_vertex_parity_sum_riemannian(n_max)
        DeL, DoL, diffL = m3_vertex_parity_sum_lorentzian_chirality_symmetric(n_max)
        assert abs(diffR - diffL) < 1e-12, (
            f"Reading A: R diff {diffR} != L diff {diffL} at n_max={n_max}"
        )


def test_m3_lorentzian_chirality_pairing_trivializes_bit_exact():
    """Reading B (chirality-pairing): D_+ = D_- bit-exact."""
    for n_max in [1, 2, 3, 5, 10]:
        Dp, Dm, diffB = m3_vertex_parity_chirality_pairing(n_max)
        assert diffB == 0.0, (
            f"Reading B: D_+ - D_- != 0 at n_max={n_max}: {diffB}"
        )


def test_m3_trivialization_reading_A_falsified():
    """L0 prediction FALSIFIED under Reading A at every tested n_max."""
    for n_max in [1, 2, 3, 4, 5]:
        r = verify_m3_trivialization(n_max)
        assert not r['reading_A_trivializes'], (
            f"L0 prediction (Reading A) unexpectedly passed at n_max={n_max}"
        )


def test_m3_trivialization_reading_B_confirmed():
    """L0 prediction CONFIRMED under Reading B (chirality-pairing)."""
    for n_max in [1, 2, 3, 4, 5]:
        r = verify_m3_trivialization(n_max)
        assert r['reading_B_trivializes']


def test_m3_verdict_string_includes_convention_dependent():
    """The M3 verdict text mentions the convention-dependent finding."""
    r = verify_m3_trivialization(3)
    assert 'convention-dependent' in r['verdict'].lower() \
        or 'reading b' in r['verdict'].lower()


# ---------------------------------------------------------------------------
# Comprehensive audit_at_4_6 driver
# ---------------------------------------------------------------------------


def test_audit_at_4_6_n_max_1_load_bearing_pass():
    """audit_at_4_6 load-bearing falsifier (J_L^2 = +I) passes at n_max=1."""
    r = audit_at_4_6(n_max=1, N_t=1)
    assert r['all_load_bearing_pass']
    assert r['(i)_J_L_squared_plus_I']['pass']


def test_audit_at_4_6_four_bbb_axioms_pass_at_all_panel_cells():
    """The four BBB-predicted-sign axioms pass at every (n_max, N_t)."""
    for n_max in [1, 2, 3]:
        for N_t in [1, 11]:
            r = audit_at_4_6(n_max=n_max, N_t=N_t, sample_size=2)
            for ax in ['(i)_J_L_squared_plus_I',
                       '(ii)_J_L_anticomm_chi',
                       '(iii)_J_L_anticomm_eta',
                       '(iv)_J_L_D_commutation']:
                assert r[ax]['pass'], (
                    f"Axiom {ax} failed at (n_max={n_max}, N_t={N_t}): "
                    f"residual={r[ax]['residual']}"
                )


def test_audit_at_4_6_structural_finding_chi_D_documented():
    """Axiom (v) fails on truthful D_GV; STRUCTURAL_FINDING key present."""
    r = audit_at_4_6(n_max=2, N_t=1)
    assert not r['(v)_chi_D_anticomm']['pass']  # expected fail
    assert 'STRUCTURAL_FINDING' in r['(v)_chi_D_anticomm']


def test_audit_at_4_6_bbb_signs_captured():
    """audit result captures the BBB signs dict."""
    r = audit_at_4_6(n_max=1, N_t=1)
    assert r['BBB_signs']['eps'] == +1
    assert r['BBB_signs']['eps_pp'] == -1
