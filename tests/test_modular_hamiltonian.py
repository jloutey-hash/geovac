"""Tests for geovac.modular_hamiltonian: Sprint L1 modular Hamiltonian closure.

Verifies:
  - HemisphericWedge: idempotent (P^2 = P), Hermitian, dim = dim_H/2
  - TomitaConjugation: J_TT^2 = +I (NOT -I; signature distinct from J_GV)
  - KMSState: trace cyclicity, KMS condition at beta
  - ModularHamiltonian:
      * BW canonical sigma_{2*pi}(O) = O at machine precision (n_max=2,3,4)
      * Witness collapse: HH, Sew, Unruh all match BW residual
      * Propinquity rate cross-check: residual << gamma_{n_max} from L2
      * KMS condition at expectation level

Sprint L1 (2026-05-16): principal falsifier named in Paper 32 SVIII rem,
Paper 34 SIII.27, Paper 38.
"""

from __future__ import annotations

import numpy as np
import pytest

from geovac.modular_hamiltonian import (
    HemisphericWedge,
    KMSState,
    ModularHamiltonian,
    TomitaConjugation,
    for_bisognano_wichmann,
    for_hartle_hawking,
    for_sewell,
    for_unruh,
    propinquity_rate_check,
    verify_cross_witness_collapse,
)
from geovac.full_dirac_operator_system import (
    FullDiracTruncatedOperatorSystem,
    camporesi_higuchi_full_dirac_matrix,
    full_dirac_basis,
)


# ---------------------------------------------------------------------------
# HemisphericWedge
# ---------------------------------------------------------------------------


def test_hemispheric_wedge_idempotent_n_max_2():
    basis = full_dirac_basis(n_max=2)
    wedge = HemisphericWedge(axis="hopf")
    P = wedge.projection_matrix(basis)
    residual = np.linalg.norm(P @ P - P)
    assert residual < 1e-12, f"P^2 = P failed: ||P^2-P|| = {residual:.2e}"


def test_hemispheric_wedge_idempotent_n_max_3():
    basis = full_dirac_basis(n_max=3)
    wedge = HemisphericWedge(axis="hopf")
    P = wedge.projection_matrix(basis)
    residual = np.linalg.norm(P @ P - P)
    assert residual < 1e-12, f"P^2 = P failed: {residual:.2e}"


def test_hemispheric_wedge_idempotent_n_max_4():
    basis = full_dirac_basis(n_max=4)
    wedge = HemisphericWedge(axis="hopf")
    P = wedge.projection_matrix(basis)
    residual = np.linalg.norm(P @ P - P)
    assert residual < 1e-12, f"P^2 = P failed: {residual:.2e}"


def test_hemispheric_wedge_hermitian():
    """P_W must be Hermitian (it is a real orthogonal projector)."""
    for n_max in (2, 3, 4):
        basis = full_dirac_basis(n_max=n_max)
        wedge = HemisphericWedge(axis="hopf")
        P = wedge.projection_matrix(basis)
        residual = np.linalg.norm(P - P.conj().T)
        assert residual < 1e-13, (
            f"n_max={n_max}: P_W not Hermitian, ||P - P^*|| = {residual:.2e}"
        )


def test_hemispheric_wedge_reflection_is_involution():
    """R_axis^2 = I for any axis."""
    for n_max in (2, 3, 4):
        basis = full_dirac_basis(n_max=n_max)
        wedge = HemisphericWedge(axis="hopf")
        R = wedge.reflection_matrix(basis)
        I = np.eye(R.shape[0], dtype=np.complex128)
        residual = np.linalg.norm(R @ R - I)
        assert residual < 1e-13, (
            f"n_max={n_max}: R^2 = I failed, ||R^2-I|| = {residual:.2e}"
        )


def test_hemispheric_wedge_dim_half_of_total():
    """The hemisphere has dimension dim_H / 2 (m_j -> -m_j has no fixed points
    on the spinor basis since two_m_j is always odd)."""
    for n_max in (2, 3, 4):
        basis = full_dirac_basis(n_max=n_max)
        wedge = HemisphericWedge(axis="hopf")
        dim_W = wedge.wedge_dim(basis)
        expected = len(basis) // 2
        assert dim_W == expected, (
            f"n_max={n_max}: wedge_dim={dim_W}, expected={expected}"
        )


def test_hemispheric_wedge_axis_validation():
    """Only axis='hopf' is supported."""
    with pytest.raises(ValueError, match="axis"):
        HemisphericWedge(axis="invalid")


def test_hemispheric_wedge_restrict_to_wedge_kills_complement():
    """P_W O P_W applied to a state in the -1 eigenspace gives zero."""
    basis = full_dirac_basis(n_max=2)
    wedge = HemisphericWedge(axis="hopf")
    P = wedge.projection_matrix(basis)
    R = wedge.reflection_matrix(basis)
    # The -1 eigenspace projector is Q = (1/2)(I - R)
    Q = 0.5 * (np.eye(R.shape[0], dtype=np.complex128) - R)
    # P_W Q P_W should be zero (different eigenspaces)
    PQP = P @ Q @ P
    assert np.linalg.norm(PQP) < 1e-12, (
        f"||PQP|| = {np.linalg.norm(PQP):.2e}, expected zero"
    )


# ---------------------------------------------------------------------------
# TomitaConjugation — the J_TT != J_GV signature
# ---------------------------------------------------------------------------


def test_tomita_J_squared_is_plus_identity():
    """The CRITICAL test: J_TT^2 = +I (NOT -I).

    This distinguishes the Tomita-Takesaki conjugation from the
    Connes KO-dim 3 charge conjugation J_GV.
    """
    for n_max in (2, 3, 4):
        bw = for_bisognano_wichmann(n_max=n_max)
        J_TT = bw.construct_tomita_J()
        ok, residual = J_TT.verify_J_squared_positive_identity()
        assert ok, (
            f"n_max={n_max}: J_TT^2 != +I, ||J^2 - I|| = {residual:.2e}"
        )


def test_tomita_J_NOT_equal_to_J_GV_signature():
    """Verify J_TT is categorically different from J_GV.

    J_GV has J^2 = -I (KO-dim 3); J_TT has J^2 = +I (Tomita).
    For our canonical Tomita conjugation U=I, J^2 = +I, so ||J^2-(-I)||
    should be 2*sqrt(dim).
    """
    for n_max in (2, 3, 4):
        bw = for_bisognano_wichmann(n_max=n_max)
        J_TT = bw.construct_tomita_J()
        not_J_GV, distance = J_TT.verify_J_squared_negative_identity_DOES_NOT_HOLD()
        assert not_J_GV, (
            f"n_max={n_max}: J_TT^2 looks like -I (would conflate with J_GV)"
        )
        # ||I - (-I)|| = ||2I|| = 2*sqrt(dim)
        expected_distance = 2.0 * np.sqrt(bw.dim)
        assert abs(distance - expected_distance) < 1e-12, (
            f"n_max={n_max}: distance={distance:.4f}, expected"
            f" {expected_distance:.4f}"
        )


def test_tomita_apply_antilinear():
    """J_TT(alpha psi) = conj(alpha) J_TT(psi)."""
    J = TomitaConjugation(U=np.eye(4, dtype=np.complex128))
    psi = np.array([1.0, 2.0 + 1.0j, 3.0, -1.0j])
    alpha = 2.0 + 3.0j
    lhs = J.apply(alpha * psi)
    rhs = np.conj(alpha) * J.apply(psi)
    assert np.allclose(lhs, rhs, atol=1e-14), "J_TT antilinearity failed"


def test_tomita_J_apply_to_operator_consistency():
    """J O J^{-1} = U conj(O) U^T."""
    U = np.eye(3, dtype=np.complex128)
    J = TomitaConjugation(U=U)
    O = np.random.RandomState(0).randn(3, 3) + 1j * np.random.RandomState(1).randn(3, 3)
    JOJinv = J.apply_to_operator(O)
    expected = np.conj(O)  # since U = I, U conj(O) U^T = conj(O)
    assert np.allclose(JOJinv, expected, atol=1e-14)


# ---------------------------------------------------------------------------
# KMSState
# ---------------------------------------------------------------------------


def test_kms_state_density_matrix_normalized():
    """Tr(rho_beta) = 1."""
    H = np.diag([1.0, 2.0, 3.0, 4.0]).astype(np.complex128)
    kms = KMSState(beta=2.0, H=H)
    rho = kms.density_matrix()
    tr = np.trace(rho)
    assert abs(tr - 1.0) < 1e-12, f"Tr(rho) = {tr}, expected 1.0"


def test_kms_state_partition_function_positive():
    H = np.diag([0.0, 1.0, 2.0]).astype(np.complex128)
    kms = KMSState(beta=1.0, H=H)
    Z = kms.partition_function()
    assert np.real(Z) > 0, f"Z = {Z}, should be positive"


def test_kms_state_expectation_self_adjoint():
    """omega(O) is real if O is Hermitian."""
    H = np.diag([0.0, 1.0, 2.0, 3.0]).astype(np.complex128)
    kms = KMSState(beta=1.0, H=H)
    O = np.diag([1.0, 2.0, 3.0, 4.0]).astype(np.complex128)
    val = kms.expectation(O)
    assert abs(np.imag(val)) < 1e-12, f"Im(omega(O)) = {np.imag(val)}"


# ---------------------------------------------------------------------------
# ModularHamiltonian — THE primary L1 falsifier
# ---------------------------------------------------------------------------


def test_modular_periodicity_n_max_2_BW():
    """L1 PRIMARY FALSIFIER: sigma_{2*pi}(O) = O at n_max=2.

    Expected: STRONG_IDENTIFICATION (residual ~ machine eps).
    """
    bw = for_bisognano_wichmann(n_max=2)
    results = bw.verify_witness()
    assert results["verdict"] == "STRONG_IDENTIFICATION", (
        f"n_max=2 verdict: {results['verdict']}, "
        f"max_residual={results['max_periodicity_residual']:.2e}"
    )


def test_modular_periodicity_n_max_3_BW():
    """L1 falsifier at n_max=3."""
    bw = for_bisognano_wichmann(n_max=3)
    results = bw.verify_witness()
    assert results["verdict"] == "STRONG_IDENTIFICATION", (
        f"n_max=3 verdict: {results['verdict']}, "
        f"max_residual={results['max_periodicity_residual']:.2e}"
    )


def test_modular_periodicity_n_max_4_BW():
    """L1 falsifier at n_max=4."""
    bw = for_bisognano_wichmann(n_max=4)
    results = bw.verify_witness()
    assert results["verdict"] == "STRONG_IDENTIFICATION", (
        f"n_max=4 verdict: {results['verdict']}, "
        f"max_residual={results['max_periodicity_residual']:.2e}"
    )


def test_K_geometric_integer_spectrum():
    """K_boost has integer spectrum (two_m_j is odd integer).

    This is the structural reason sigma_{2*pi}(O) = O closes exactly:
    e^{i*2*pi*n} = 1 for any integer n.
    """
    for n_max in (2, 3, 4):
        bw = for_bisognano_wichmann(n_max=n_max)
        K = bw.K_geometric
        # K should be diagonal
        K_diag = np.diag(K)
        K_off = K - np.diag(K_diag)
        assert np.linalg.norm(K_off) < 1e-13, (
            f"K not diagonal at n_max={n_max}"
        )
        # K eigenvalues are odd integers
        for k_val in K_diag:
            assert abs(np.imag(k_val)) < 1e-14, f"K not real: {k_val}"
            k_real = np.real(k_val)
            assert abs(k_real - np.round(k_real)) < 1e-14, (
                f"K eigenvalue {k_real} not integer"
            )
            assert int(np.round(k_real)) % 2 != 0, (
                f"K eigenvalue {int(np.round(k_real))} not odd"
            )


def test_BW_beta_equals_2pi():
    bw = for_bisognano_wichmann(n_max=2)
    assert abs(bw.beta - 2.0 * np.pi) < 1e-14


def test_BW_kappa_g_equals_1():
    bw = for_bisognano_wichmann(n_max=2)
    assert bw.kappa_g == 1.0


# ---------------------------------------------------------------------------
# Witness factories (HH, Sew, Unruh)
# ---------------------------------------------------------------------------


def test_HH_beta_equals_8piM():
    M = 1.0
    hh = for_hartle_hawking(n_max=2, M=M)
    assert abs(hh.beta - 8.0 * np.pi * M) < 1e-14
    assert abs(hh.kappa_g - 1.0 / (4.0 * M)) < 1e-14


def test_HH_M_2_doubles_beta():
    M = 2.0
    hh = for_hartle_hawking(n_max=2, M=M)
    assert abs(hh.beta - 8.0 * np.pi * M) < 1e-14


def test_Sewell_equals_HH():
    """Sewell IS Hartle-Hawking at the operator-system level."""
    sw = for_sewell(n_max=2, M=1.0)
    hh = for_hartle_hawking(n_max=2, M=1.0)
    assert sw.kappa_g == hh.kappa_g
    assert sw.beta == hh.beta


def test_Unruh_beta_equals_2pi_over_a():
    a = 2.0
    un = for_unruh(n_max=2, a=a)
    assert abs(un.beta - 2.0 * np.pi / a) < 1e-14
    assert un.kappa_g == a


def test_HH_modular_periodicity_strong():
    """Hartle-Hawking at M=1 also closes at machine precision."""
    hh = for_hartle_hawking(n_max=2, M=1.0)
    results = hh.verify_witness()
    assert results["verdict"] == "STRONG_IDENTIFICATION", (
        f"HH n_max=2 verdict: {results['verdict']}, "
        f"max_res={results['max_periodicity_residual']:.2e}"
    )


def test_Unruh_modular_periodicity_strong_at_a_1():
    un = for_unruh(n_max=2, a=1.0)
    results = un.verify_witness()
    assert results["verdict"] == "STRONG_IDENTIFICATION"


def test_Unruh_modular_periodicity_strong_at_a_2():
    """Test Unruh at a=2 (different beta), still period closure at 2*pi."""
    un = for_unruh(n_max=2, a=2.0)
    results = un.verify_witness()
    assert results["verdict"] == "STRONG_IDENTIFICATION", (
        f"Unruh a=2 verdict: {results['verdict']}, "
        f"max_res={results['max_periodicity_residual']:.2e}"
    )


# ---------------------------------------------------------------------------
# Cross-witness collapse
#
# NOTE: the `verify_cross_witness_collapse` tests below check that the
# period-closure RESIDUALS agree across witnesses. Those residuals are
# kappa_g/beta-independent by construction (the flow uses K_geometric only), so
# these are weak consistency checks, NOT the load-bearing collapse backing. The
# genuine cross-witness collapse (bit-identical modular OPERATORS rho/Delta/K_TT
# across genuinely-distinct-beta witnesses, with a beta-cancellation negative
# control) is proved by `test_cross_witness_modular_operator_bit_identical` and
# `test_cross_witness_collapse_requires_beta_cancellation` above.
# ---------------------------------------------------------------------------


def test_cross_witness_collapse_n_max_2():
    """All four witnesses give bit-identical period-closure residuals (weak
    consistency check; see the section note — the operator-level test is the
    load-bearing collapse backing)."""
    results = verify_cross_witness_collapse(n_max=2)
    assert results["max_consistency_residual"] < 1e-14, (
        f"Cross-witness inconsistency at n_max=2: "
        f"{results['max_consistency_residual']:.2e}"
    )


def test_cross_witness_collapse_n_max_3():
    results = verify_cross_witness_collapse(n_max=3)
    assert results["max_consistency_residual"] < 1e-14


def test_cross_witness_collapse_n_max_4():
    results = verify_cross_witness_collapse(n_max=4)
    assert results["max_consistency_residual"] < 1e-14


def test_cross_witness_kappa_g_values_distinct():
    """Make sure the witness factories actually use distinct kappa_g."""
    results = verify_cross_witness_collapse(n_max=2)
    # BW kappa_g = 1, HH kappa_g = 0.25, Sew = 0.25, Unruh_a1 = 1, Unruh_a2 = 2
    assert results["BW"]["kappa_g"] == 1.0
    assert results["HH"]["kappa_g"] == 0.25
    assert results["Unruh_a2"]["kappa_g"] == 2.0


def test_cross_witness_modular_operator_bit_identical():
    """The cross-witness collapse is a NON-TRIVIAL consequence of the BW choice
    H_local = K_alpha^W / beta, not a code-identity.

    Each witness has a GENUINELY DIFFERENT beta (= 2*pi/kappa_g), and beta enters
    the Tomita construction through rho_W = e^{-beta H_local}/Z.  Because
    beta H_local = K_alpha^W is beta-independent, the density matrix rho_W (and
    hence the full modular operator Delta and modular Hamiltonian K_TT) is
    bit-identical across the witnesses.  We compare the OPERATORS (not just the
    period-closure residuals, which are all ~0 and would match trivially)."""
    n_max = 3
    witnesses = {
        "BW": for_bisognano_wichmann(n_max=n_max),       # beta = 2*pi
        "HH_M1": for_hartle_hawking(n_max=n_max, M=1.0),  # beta = 8*pi
        "Unruh_a2": for_unruh(n_max=n_max, a=2.0),        # beta = pi
    }
    betas = {name: w.beta for name, w in witnesses.items()}
    # the betas must be genuinely distinct (else the collapse is vacuous)
    bvals = sorted(betas.values())
    assert all(b2 - b1 > 0.5 for b1, b2 in zip(bvals, bvals[1:])), (
        f"witness betas not distinct enough: {betas}"
    )

    tms = {name: w.tomita_structure() for name, w in witnesses.items()}
    ref = tms["BW"]
    for name, t in tms.items():
        assert np.linalg.norm(t.rho - ref.rho) < 1e-13, (
            f"{name}: rho_W not bit-identical to BW (beta={betas[name]:.3f} vs "
            f"{betas['BW']:.3f}); ||drho||={np.linalg.norm(t.rho - ref.rho):.2e}"
        )
        assert np.linalg.norm(t.Delta_matrix - ref.Delta_matrix) < 1e-13, (
            f"{name}: Delta not bit-identical to BW"
        )
        assert np.linalg.norm(t.K_TT - ref.K_TT) < 1e-13, (
            f"{name}: K_TT not bit-identical to BW"
        )


def test_cross_witness_collapse_requires_beta_cancellation():
    """Negative control: WITHOUT the 1/beta in H_local the collapse FAILS.

    If one (wrongly) set H_local = K_alpha^W (no division by beta), then
    rho_W = e^{-beta K_alpha^W}/Z is genuinely beta-dependent, so witnesses with
    different beta give DIFFERENT density matrices.  This confirms the collapse
    tested above is a real consequence of the H_local = K_alpha^W/beta
    construction, not an artifact of comparing identical code paths."""
    n_max = 3
    un = for_unruh(n_max=n_max, a=2.0)               # beta = pi   (widest spread)
    hh = for_hartle_hawking(n_max=n_max, M=1.0)       # beta = 8*pi
    assert abs(un.beta - hh.beta) > 0.5

    # beta-independent H_local (the WRONG choice) -> beta-dependent rho_W
    tms_un = un.tomita_structure(H_local=un.restrict_K_alpha_to_wedge())
    tms_hh = hh.tomita_structure(H_local=hh.restrict_K_alpha_to_wedge())
    drho = np.linalg.norm(tms_un.rho - tms_hh.rho)
    # ~4.7e-4 here vs < 1e-13 for the bit-identical (correct-/beta) collapse: a
    # ~9-order separation confirms the collapse is non-vacuous.
    assert drho > 1e-5, (
        f"control FAILED to discriminate: rho_W should be beta-dependent without "
        f"the 1/beta, but ||drho||={drho:.2e} (collapse test would be vacuous)"
    )


def test_operator_level_period_lift_DW_vs_Kalpha():
    """Paper 42 §5.5(II): the operator-level period lift distinguishes the intrinsic
    wedge Dirac D_W (half-integer spectrum -> e^{i2pi D_W} = -I, the spinor double
    cover) from the boost-class K_alpha^W (odd-integer spectrum -> e^{i2pi K_alpha} =
    +I, the identity).  The CONJUGATION modular flow sigma_2pi(O) = U O U^{-1} closes
    for BOTH (a global scalar cancels), so the distinction is at the operator level,
    NOT the flow-closure level.  Backs the 'derived structural finding' of §5.5(II)
    that was prose-only before (the v4.43.x cert coverage gap)."""
    import scipy.linalg
    for n_max in (2, 3, 4):
        mh = for_bisognano_wichmann(n_max=n_max, axis="hopf")
        K_alpha = mh.restrict_K_alpha_to_wedge()
        D_W = mh.restrict_to_wedge_block(mh.D)
        dim = D_W.shape[0]
        Iden = np.eye(dim, dtype=np.complex128)

        lift_DW = scipy.linalg.expm(1j * 2.0 * np.pi * D_W)
        lift_Ka = scipy.linalg.expm(1j * 2.0 * np.pi * K_alpha)
        # operator-level: D_W -> -I (double cover), K_alpha -> +I (identity)
        assert np.linalg.norm(lift_DW + Iden) < 1e-12, (
            f"n_max={n_max}: e^(i2pi D_W) should be -I; "
            f"||lift+I|| = {np.linalg.norm(lift_DW + Iden):.2e}"
        )
        assert np.linalg.norm(lift_Ka - Iden) < 1e-12, (
            f"n_max={n_max}: e^(i2pi K_alpha) should be +I; "
            f"||lift-I|| = {np.linalg.norm(lift_Ka - Iden):.2e}"
        )
        # the lifts are genuinely DIFFERENT operators (not vacuously equal)
        assert np.linalg.norm(lift_DW - lift_Ka) > 1.0, (
            f"n_max={n_max}: the +I and -I lifts must differ (||diff||~2*sqrt(dim))"
        )
        # but the conjugation flow sigma_2pi(O) = U O U^{-1} closes for D_W too:
        # the global -I cancels, (-I) O (-I) = O, for an arbitrary deterministic O
        O = (np.arange(dim * dim, dtype=np.float64).reshape(dim, dim)
             + 1j * np.arange(dim * dim, 0, -1, dtype=np.float64).reshape(dim, dim))
        O = O / np.linalg.norm(O)   # scale-free so the residual reflects the lift
        sigma_DW = lift_DW @ O @ lift_DW.conj().T
        assert np.linalg.norm(sigma_DW - O) < 1e-12, (
            f"n_max={n_max}: conjugation flow with D_W must close (scalar cancels); "
            f"||sigma_2pi(O)-O|| = {np.linalg.norm(sigma_DW - O):.2e}"
        )


# ---------------------------------------------------------------------------
# KMS condition (expectation level)
# ---------------------------------------------------------------------------


def test_kms_condition_at_beta_n_max_2():
    """omega(A * sigma_{i*beta}(B)) = omega(B * A) at beta = 2*pi."""
    bw = for_bisognano_wichmann(n_max=2)
    op_sys = FullDiracTruncatedOperatorSystem(n_max=2)
    # Pick first two non-identity multipliers
    gens = [M for (lbl, M) in op_sys.basis_matrices if lbl != (1, 0, 0)][:2]
    A_W = bw.wedge.restrict_to_wedge(gens[0], bw.basis)
    B_W = bw.wedge.restrict_to_wedge(gens[1], bw.basis)
    ok, residual = bw.verify_kms_condition_at_beta(A_W, B_W, tol=1e-10)
    assert ok, f"KMS condition failed: residual = {residual:.2e}"


def test_kms_condition_at_beta_n_max_3():
    bw = for_bisognano_wichmann(n_max=3)
    op_sys = FullDiracTruncatedOperatorSystem(n_max=3)
    gens = [M for (lbl, M) in op_sys.basis_matrices if lbl != (1, 0, 0)][:2]
    A_W = bw.wedge.restrict_to_wedge(gens[0], bw.basis)
    B_W = bw.wedge.restrict_to_wedge(gens[1], bw.basis)
    ok, residual = bw.verify_kms_condition_at_beta(A_W, B_W, tol=1e-10)
    assert ok, f"KMS condition failed: residual = {residual:.2e}"


# ---------------------------------------------------------------------------
# Propinquity rate consistency with L2 (cross-check to Paper 38)
# ---------------------------------------------------------------------------


def test_propinquity_rate_residual_dominated_by_machine_eps():
    """L1 modular residual is at machine precision; gamma_n is O(1).

    This confirms the STRONG_IDENTIFICATION verdict — the residual is
    bit-exact, not merely scaling as gamma_n -> 0. The ratio
    residual/gamma_n should be << 1 (much smaller than the qualitative-
    rate expectation), which is the "stronger than soft" closure.
    """
    for n_max in (2, 3, 4):
        pr = propinquity_rate_check(n_max=n_max)
        # The residual / gamma_n ratio should be at machine precision
        assert pr["max_residual"] < 1e-12, (
            f"n_max={n_max}: max_residual={pr['max_residual']:.2e}"
        )
        # gamma_n is O(1) at small n_max
        assert pr["gamma_n_L2"] > 0.1, f"gamma_n unexpectedly small: {pr['gamma_n_L2']}"
        # Ratio should be machine-eps over O(1)
        assert pr["residual_over_gamma"] < 1e-12, (
            f"n_max={n_max}: ratio = {pr['residual_over_gamma']:.2e}"
        )


# ---------------------------------------------------------------------------
# Leakage quantification (per L1-A §10a obstruction)
# ---------------------------------------------------------------------------


def test_operator_system_leakage_per_n_max():
    """L1-A §10a obstruction: operator-system non-multiplicative-closure
    creates 'leakage' when modular flow is applied. Quantify the leakage
    at n_max=2,3,4 by measuring how much of sigma_{2*pi}(O) sits outside
    the operator system.

    Since we get STRONG_IDENTIFICATION (sigma_{2*pi}(O) = O bit-exact),
    the leakage should also be zero by construction (O is in O_system,
    sigma_{2*pi}(O) = O, so sigma_{2*pi}(O) is also in O_system).
    """
    for n_max in (2, 3, 4):
        bw = for_bisognano_wichmann(n_max=n_max)
        op_sys = FullDiracTruncatedOperatorSystem(n_max=n_max)
        # Take first non-identity multiplier
        M = None
        for label, mat in op_sys.basis_matrices:
            if label != (1, 0, 0):
                M = mat
                break
        # Verify M is in the operator system (it should be, by construction)
        is_in_O, residual_in_O = op_sys.contains(M)
        assert is_in_O, f"Generator not in O at n_max={n_max}"
        # Apply sigma_{2*pi}
        sigma_M = bw.modular_flow_real_time(2.0 * np.pi, M)
        # Verify sigma_M is also in O (leakage = 0)
        is_in_O_after, leak = op_sys.contains(sigma_M, tol=1e-8)
        assert is_in_O_after, (
            f"sigma_{{2pi}}(M) NOT in O at n_max={n_max}: leak={leak:.2e}"
        )


# ---------------------------------------------------------------------------
# Modular flow basic properties
# ---------------------------------------------------------------------------


def test_modular_flow_at_t_0_is_identity():
    """sigma_{t=0}(O) = O for any O."""
    bw = for_bisognano_wichmann(n_max=2)
    O = np.random.RandomState(42).randn(bw.dim, bw.dim) + 1j * np.random.RandomState(43).randn(bw.dim, bw.dim)
    sigma_O = bw.modular_flow_real_time(0.0, O)
    assert np.allclose(sigma_O, O, atol=1e-13), "sigma_0 != id"


def test_modular_flow_unitary_real_time():
    """sigma_t is unitary, so it preserves Frobenius norm of Hermitian O."""
    bw = for_bisognano_wichmann(n_max=2)
    O = np.random.RandomState(7).randn(bw.dim, bw.dim) + 1j * np.random.RandomState(8).randn(bw.dim, bw.dim)
    O = 0.5 * (O + O.conj().T)  # Hermitian
    for t in (0.5, 1.0, np.pi):
        sigma_O = bw.modular_flow_real_time(t, O)
        # Frobenius norm should be preserved
        assert abs(np.linalg.norm(sigma_O) - np.linalg.norm(O)) < 1e-12, (
            f"||sigma_t(O)|| != ||O|| at t={t}"
        )


def test_modular_flow_group_law():
    """sigma_t * sigma_s = sigma_{t+s}."""
    bw = for_bisognano_wichmann(n_max=2)
    O = np.random.RandomState(11).randn(bw.dim, bw.dim) + 1j * np.random.RandomState(12).randn(bw.dim, bw.dim)
    t = 0.7
    s = 1.3
    lhs = bw.modular_flow_real_time(t, bw.modular_flow_real_time(s, O))
    rhs = bw.modular_flow_real_time(t + s, O)
    assert np.allclose(lhs, rhs, atol=1e-12), (
        "Modular flow group law failed"
    )


# ---------------------------------------------------------------------------
# Slow tests (n_max = 5 if budget permits)
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_modular_periodicity_n_max_5_slow():
    """n_max=5 verification (dim_H = 110, slow but tractable)."""
    bw = for_bisognano_wichmann(n_max=5)
    results = bw.verify_witness()
    assert results["verdict"] == "STRONG_IDENTIFICATION", (
        f"n_max=5 verdict: {results['verdict']}, "
        f"max_res={results['max_periodicity_residual']:.2e}"
    )


@pytest.mark.slow
def test_cross_witness_collapse_n_max_5_slow():
    results = verify_cross_witness_collapse(n_max=5)
    assert results["max_consistency_residual"] < 1e-12


# ===========================================================================
# Sprint L1-tighten: Load-bearing BW-gamma Tomita-Takesaki construction
# ===========================================================================
#
# These tests verify the Tomita-Takesaki polar decomposition is itself
# load-bearing (NOT just an expectation-level cross-check) and that the
# BW-gamma modular flow closes sigma_{2*pi}^TT(O) = O bit-exactly.
#
# Outcome 1 prediction (PI's prior): with H_local = K_alpha^W / beta,
# K_TT = ad(K_alpha^W) and sigma_t^TT(a) = e^{-it K_alpha^W} a e^{+it K_alpha^W}
# has the same period closure as BW-alpha (up to flow sign, which
# vanishes at the period). UNIFIED_STRONG verdict.

from geovac.modular_hamiltonian import TomitaModularStructure


# ----- TomitaModularStructure core construction -----


def test_tomita_structure_construct_n_max_2():
    """TomitaModularStructure constructs cleanly at n_max=2."""
    bw = for_bisognano_wichmann(n_max=2)
    tms = bw.tomita_structure()
    assert tms.dim_W == bw.dim // 2
    assert tms.dim_GNS == tms.dim_W ** 2
    assert tms.rho is not None
    assert tms.K_TT is not None


def test_tomita_structure_rho_is_density_matrix_n_max_2():
    """rho is positive Hermitian with unit trace."""
    bw = for_bisognano_wichmann(n_max=2)
    tms = bw.tomita_structure()
    # Hermitian
    assert np.allclose(tms.rho, tms.rho.conj().T, atol=1e-14)
    # Positive
    eigs = np.linalg.eigvalsh(tms.rho)
    assert (eigs >= -1e-14).all()
    # Trace = 1
    assert abs(np.trace(tms.rho) - 1.0) < 1e-12


def test_tomita_rho_sqrt_squares_to_rho():
    """rho_sqrt @ rho_sqrt = rho (rho_sqrt is the principal square root)."""
    for n_max in (2, 3, 4):
        bw = for_bisognano_wichmann(n_max=n_max)
        tms = bw.tomita_structure()
        recon = tms.rho_sqrt @ tms.rho_sqrt
        assert np.allclose(recon, tms.rho, atol=1e-12), (
            f"n_max={n_max}: ||rho_sqrt^2 - rho|| = "
            f"{np.linalg.norm(recon - tms.rho):.2e}"
        )


def test_tomita_rho_invsqrt_inverts_rho_sqrt():
    """rho_invsqrt @ rho_sqrt = I."""
    for n_max in (2, 3, 4):
        bw = for_bisognano_wichmann(n_max=n_max)
        tms = bw.tomita_structure()
        prod = tms.rho_invsqrt @ tms.rho_sqrt
        I = np.eye(tms.dim_W, dtype=np.complex128)
        assert np.allclose(prod, I, atol=1e-10), (
            f"n_max={n_max}: ||rho_invsqrt rho_sqrt - I|| = "
            f"{np.linalg.norm(prod - I):.2e}"
        )


def test_tomita_delta_positive_n_max_2():
    """Delta is positive Hermitian on H_GNS."""
    bw = for_bisognano_wichmann(n_max=2)
    tms = bw.tomita_structure()
    # Hermitian
    assert np.allclose(
        tms.Delta_matrix, tms.Delta_matrix.conj().T, atol=1e-12,
    )
    # Positive
    eigs = np.linalg.eigvalsh(tms.Delta_matrix)
    assert (eigs > -1e-12).all()


def test_tomita_K_TT_hermitian_n_max_2():
    """K_TT = -log Delta is Hermitian."""
    bw = for_bisognano_wichmann(n_max=2)
    tms = bw.tomita_structure()
    K = tms.K_TT
    assert np.allclose(K, K.conj().T, atol=1e-12), (
        f"K_TT not Hermitian: ||K - K^*|| = "
        f"{np.linalg.norm(K - K.conj().T):.2e}"
    )


# ----- J_TT^2 = +I at the FULL Tomita-Takesaki level (CRITICAL) -----


def test_tomita_J_TT_squared_is_plus_identity_full_construction_n_max_2():
    """J_TT^2 = +I at the FULL Tomita-Takesaki level on H_GNS at n_max=2.

    This is the critical test distinguishing J_TT (state-dependent) from
    J_GV (KO-dim 3 charge conjugation, J_GV^2 = -I). The L1 version of
    this test (test_tomita_J_squared_is_plus_identity) only checks the
    canonical sanity-check U = I; L1-tighten checks the actual Tomita
    J_TT acting on the GNS Hilbert-Schmidt space via apply_J_TT.
    """
    bw = for_bisognano_wichmann(n_max=2)
    tms = bw.tomita_structure()
    residual = tms.J_TT_squared_residual()
    assert residual < 1e-10, (
        f"J_TT^2 != +I at full TT level: residual = {residual:.2e}"
    )


def test_tomita_J_TT_squared_is_plus_identity_full_construction_n_max_3():
    bw = for_bisognano_wichmann(n_max=3)
    tms = bw.tomita_structure()
    residual = tms.J_TT_squared_residual()
    assert residual < 1e-10, (
        f"J_TT^2 != +I at full TT level n_max=3: residual = {residual:.2e}"
    )


# ----- Modular flow on the algebra (BW-gamma) -----


def test_tomita_modular_flow_at_t_0_identity():
    """sigma_t^TT(a) = a at t = 0."""
    bw = for_bisognano_wichmann(n_max=2)
    tms = bw.tomita_structure()
    rng = np.random.RandomState(0)
    a = rng.randn(tms.dim_W, tms.dim_W) + 1j * rng.randn(tms.dim_W, tms.dim_W)
    sigma_a = tms.modular_flow_on_algebra(0.0, a)
    assert np.allclose(sigma_a, a, atol=1e-13)


def test_tomita_modular_flow_unitary():
    """sigma_t^TT(a) preserves the Frobenius norm of Hermitian a."""
    bw = for_bisognano_wichmann(n_max=2)
    tms = bw.tomita_structure()
    rng = np.random.RandomState(1)
    a = rng.randn(tms.dim_W, tms.dim_W) + 1j * rng.randn(tms.dim_W, tms.dim_W)
    a = 0.5 * (a + a.conj().T)
    for t in (0.5, 1.0, np.pi):
        sigma_a = tms.modular_flow_on_algebra(t, a)
        norm_a = np.linalg.norm(a)
        norm_sigma_a = np.linalg.norm(sigma_a)
        assert abs(norm_a - norm_sigma_a) < 1e-12, (
            f"t={t}: ||sigma_t(a)|| = {norm_sigma_a}, ||a|| = {norm_a}"
        )


# ----- Modular periodicity at the FULL Tomita-Takesaki level -----


def test_tomita_modular_periodicity_n_max_2():
    """sigma_{2*pi}^TT(a) = a at the full Tomita-Takesaki level at n_max=2.

    THE LOAD-BEARING L1-TIGHTEN TEST: this is the closure of the
    principal Sprint L1 falsifier via the Tomita-Takesaki construction
    (not just BW-alpha geometric).
    """
    bw = for_bisognano_wichmann(n_max=2)
    tms = bw.tomita_structure()
    # Test on first 5 wedge-restricted multipliers
    op_sys = FullDiracTruncatedOperatorSystem(n_max=2)
    gens = [M for (lbl, M) in op_sys.basis_matrices if lbl != (1, 0, 0)][:5]
    max_res = 0.0
    for M in gens:
        a_W = bw.restrict_to_wedge_block(M)
        ok, res = tms.verify_modular_periodicity_tomita(a_W, tol=1e-10)
        assert ok, (
            f"sigma_{{2pi}}^TT(O) != O at n_max=2: residual = {res:.2e}"
        )
        max_res = max(max_res, res)
    assert max_res < 1e-10


def test_tomita_modular_periodicity_n_max_3():
    bw = for_bisognano_wichmann(n_max=3)
    tms = bw.tomita_structure()
    op_sys = FullDiracTruncatedOperatorSystem(n_max=3)
    gens = [M for (lbl, M) in op_sys.basis_matrices if lbl != (1, 0, 0)][:5]
    for M in gens:
        a_W = bw.restrict_to_wedge_block(M)
        ok, res = tms.verify_modular_periodicity_tomita(a_W, tol=1e-10)
        assert ok, (
            f"sigma_{{2pi}}^TT(O) != O at n_max=3: residual = {res:.2e}"
        )


def test_tomita_modular_periodicity_n_max_4():
    bw = for_bisognano_wichmann(n_max=4)
    tms = bw.tomita_structure()
    op_sys = FullDiracTruncatedOperatorSystem(n_max=4)
    gens = [M for (lbl, M) in op_sys.basis_matrices if lbl != (1, 0, 0)][:5]
    for M in gens:
        a_W = bw.restrict_to_wedge_block(M)
        ok, res = tms.verify_modular_periodicity_tomita(a_W, tol=1e-10)
        assert ok, (
            f"sigma_{{2pi}}^TT(O) != O at n_max=4: residual = {res:.2e}"
        )


# ----- Tomita vs alpha comparison (Outcome 1: UNIFIED_STRONG) -----


def test_tomita_vs_alpha_comparison_n_max_2():
    """Verify UNIFIED_STRONG verdict: BW-alpha and BW-gamma both close,
    and the flows agree up to sign (sigma_t^TT = sigma_{-t}^alpha).

    This is the PI's prior "Outcome 1" that L1-tighten is verifying.
    """
    bw = for_bisognano_wichmann(n_max=2)
    comp = bw.compare_alpha_vs_tomita(n_test_operators=5)
    assert comp["verdict"] == "UNIFIED_STRONG", (
        f"n_max=2 verdict: {comp['verdict']}, "
        f"max_res_alpha={comp['max_residual_alpha_2pi']:.2e}, "
        f"max_res_TT={comp['max_residual_TT_2pi']:.2e}, "
        f"avg_diff_at_t1={comp['avg_TT_vs_alpha_t1']:.2e}, "
        f"avg_diff_at_neg_t1={comp['avg_TT_vs_alpha_neg_t1']:.2e}"
    )


def test_tomita_vs_alpha_comparison_n_max_3():
    bw = for_bisognano_wichmann(n_max=3)
    comp = bw.compare_alpha_vs_tomita(n_test_operators=5)
    assert comp["verdict"] == "UNIFIED_STRONG"


def test_tomita_vs_alpha_comparison_n_max_4():
    bw = for_bisognano_wichmann(n_max=4)
    comp = bw.compare_alpha_vs_tomita(n_test_operators=5)
    assert comp["verdict"] == "UNIFIED_STRONG"


def test_tomita_flow_equals_alpha_flow_at_neg_t():
    """sigma_t^TT(a) = sigma_{-t}^alpha(a) bit-exact (the conjugate
    relation between Tomita and geometric flows)."""
    bw = for_bisognano_wichmann(n_max=2)
    tms = bw.tomita_structure()
    K_alpha_W = bw.restrict_K_alpha_to_wedge()
    op_sys = FullDiracTruncatedOperatorSystem(n_max=2)
    gens = [M for (lbl, M) in op_sys.basis_matrices if lbl != (1, 0, 0)][:3]
    for M in gens:
        a_W = bw.restrict_to_wedge_block(M)
        for t in (0.5, 1.0, np.pi):
            cmp = tms.compare_K_TT_action_to_K_alpha(
                K_alpha_W, bw.beta, a_W, t=t,
            )
            assert cmp["TT_vs_alpha_neg_t_diff"] < 1e-10, (
                f"t={t}: sigma_t^TT != sigma_{{-t}}^alpha, "
                f"diff = {cmp['TT_vs_alpha_neg_t_diff']:.2e}"
            )


def test_tomita_flow_distinct_from_alpha_flow_at_t():
    """sigma_t^TT(a) != sigma_t^alpha(a) at general t (the flows are
    conjugate, not equal, except at the period)."""
    bw = for_bisognano_wichmann(n_max=2)
    tms = bw.tomita_structure()
    K_alpha_W = bw.restrict_K_alpha_to_wedge()
    op_sys = FullDiracTruncatedOperatorSystem(n_max=2)
    gens = [M for (lbl, M) in op_sys.basis_matrices if lbl != (1, 0, 0)][:3]
    # Pick non-trivial t (and a generator that's not commuting with K_alpha)
    found_distinct = False
    for M in gens:
        a_W = bw.restrict_to_wedge_block(M)
        cmp = tms.compare_K_TT_action_to_K_alpha(
            K_alpha_W, bw.beta, a_W, t=1.0,
        )
        if cmp["TT_vs_alpha_diff"] > 1e-8:
            found_distinct = True
            break
    assert found_distinct, (
        "Expected at least one generator where sigma_t^TT != sigma_t^alpha"
        " at t=1"
    )


# ----- Cross-witness Tomita-level collapse -----


def test_tomita_cross_witness_collapse_n_max_2():
    """All four witnesses give bit-identical Tomita-level period closure.

    For the unit-normalization U-1 (kappa_g varies, beta = 2*pi/kappa_g
    varies), the Tomita modular flow at any beta has period 2*pi at the
    level of the BW canonical setup, because K_TT = ad(K_alpha) is
    beta-independent (K_alpha is the integer-spectrum generator).
    """
    n_max = 2
    bws = {
        "BW": for_bisognano_wichmann(n_max=n_max),
        "HH_M1": for_hartle_hawking(n_max=n_max, M=1.0),
        "Sew_M1": for_sewell(n_max=n_max, M=1.0),
        "Unruh_a1": for_unruh(n_max=n_max, a=1.0),
        "Unruh_a2": for_unruh(n_max=n_max, a=2.0),
    }
    op_sys = FullDiracTruncatedOperatorSystem(n_max=n_max)
    M = [m for (lbl, m) in op_sys.basis_matrices if lbl != (1, 0, 0)][0]
    residuals = {}
    for name, bw in bws.items():
        tms = bw.tomita_structure()
        a_W = bw.restrict_to_wedge_block(M)
        ok, res = tms.verify_modular_periodicity_tomita(a_W, tol=1e-10)
        residuals[name] = res
        assert ok, (
            f"Witness {name} Tomita periodicity failed: res = {res:.2e}"
        )
    # All witnesses should have bit-similar residuals
    res_values = list(residuals.values())
    assert max(res_values) < 1e-10
    # Cross-consistency (residuals should match to machine precision)
    bw_res = residuals["BW"]
    for name, r in residuals.items():
        assert abs(r - bw_res) < 1e-10, (
            f"Cross-witness Tomita mismatch: {name} = {r:.2e} vs BW = "
            f"{bw_res:.2e}"
        )


def test_tomita_cross_witness_collapse_n_max_3():
    n_max = 3
    bws = {
        "BW": for_bisognano_wichmann(n_max=n_max),
        "HH_M1": for_hartle_hawking(n_max=n_max, M=1.0),
        "Unruh_a2": for_unruh(n_max=n_max, a=2.0),
    }
    op_sys = FullDiracTruncatedOperatorSystem(n_max=n_max)
    M = [m for (lbl, m) in op_sys.basis_matrices if lbl != (1, 0, 0)][0]
    for name, bw in bws.items():
        tms = bw.tomita_structure()
        a_W = bw.restrict_to_wedge_block(M)
        ok, res = tms.verify_modular_periodicity_tomita(a_W, tol=1e-10)
        assert ok, (
            f"n_max=3 witness {name} Tomita: res = {res:.2e}"
        )


# ----- K_TT spectrum verification -----


def test_K_TT_spectrum_is_log_lambda_ratios():
    """K_TT eigenvalues are log(lambda_i / lambda_j) for rho eigenvalues
    lambda_k. Diagonal pairs (i = j) give zero; cross-pairs give
    log(lambda_i / lambda_j)."""
    bw = for_bisognano_wichmann(n_max=2)
    tms = bw.tomita_structure()
    rho_eigs = np.linalg.eigvalsh(tms.rho)
    # K_TT eigenvalues should be {log(rho_eigs[i]/rho_eigs[j]) : i, j}
    K_eigs = np.linalg.eigvalsh(tms.K_TT)
    K_eigs_sorted = np.sort(K_eigs)
    # Expected: all pairwise log-ratios
    expected_eigs = []
    for i in range(len(rho_eigs)):
        for j in range(len(rho_eigs)):
            expected_eigs.append(np.log(rho_eigs[i] / rho_eigs[j]))
    expected_eigs = np.sort(np.array(expected_eigs))
    assert np.allclose(K_eigs_sorted, expected_eigs, atol=1e-9), (
        "K_TT spectrum does not match log(lambda_i/lambda_j) pattern"
    )


def test_K_TT_has_dim_W_zero_eigenvalues():
    """K_TT has at least dim_W zero eigenvalues (from the diagonal i=j
    pairs in the log(lambda_i/lambda_j) spectrum)."""
    for n_max in (2, 3, 4):
        bw = for_bisognano_wichmann(n_max=n_max)
        tms = bw.tomita_structure()
        K_eigs = np.linalg.eigvalsh(tms.K_TT)
        n_zeros = int(np.sum(np.abs(K_eigs) < 1e-10))
        assert n_zeros >= tms.dim_W, (
            f"n_max={n_max}: K_TT has {n_zeros} zero eigs, expected at "
            f"least dim_W = {tms.dim_W}"
        )


# ----- ModularHamiltonian.k_alpha_geometric and k_tomita accessors -----


def test_modular_hamiltonian_k_alpha_geometric_accessor():
    """k_alpha_geometric() returns the BW-alpha (L1) generator."""
    bw = for_bisognano_wichmann(n_max=2)
    K_alpha = bw.k_alpha_geometric()
    # Should equal K_geometric
    assert np.allclose(K_alpha, bw.K_geometric)


def test_modular_hamiltonian_k_tomita_accessor():
    """k_tomita() returns the BW-gamma load-bearing modular Hamiltonian
    on H_GNS."""
    bw = for_bisognano_wichmann(n_max=2)
    K_TT = bw.k_tomita()
    # Should have shape (dim_W^2, dim_W^2) = (64, 64) at n_max=2
    assert K_TT.shape == (64, 64)
    # Should be Hermitian
    assert np.allclose(K_TT, K_TT.conj().T, atol=1e-12)


# ----- Slow tests at n_max = 5 -----


@pytest.mark.slow
def test_tomita_modular_periodicity_n_max_5_slow():
    """Tomita-Takesaki period closure at n_max=5."""
    bw = for_bisognano_wichmann(n_max=5)
    tms = bw.tomita_structure()
    op_sys = FullDiracTruncatedOperatorSystem(n_max=5)
    gens = [M for (lbl, M) in op_sys.basis_matrices if lbl != (1, 0, 0)][:3]
    for M in gens:
        a_W = bw.restrict_to_wedge_block(M)
        ok, res = tms.verify_modular_periodicity_tomita(a_W, tol=1e-10)
        assert ok, f"n_max=5: res = {res:.2e}"


@pytest.mark.slow
def test_tomita_vs_alpha_n_max_5_slow():
    """UNIFIED_STRONG verdict at n_max=5."""
    bw = for_bisognano_wichmann(n_max=5)
    comp = bw.compare_alpha_vs_tomita(n_test_operators=3)
    assert comp["verdict"] == "UNIFIED_STRONG"
