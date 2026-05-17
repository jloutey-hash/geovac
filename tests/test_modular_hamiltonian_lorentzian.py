"""Tests for geovac.modular_hamiltonian_lorentzian: Sprint L2-E Krein-level
Paper 42 redo.

Verifies the four LOAD-BEARING falsifiers + standard structural checks:

  L2E-FALS-1 (LOAD-BEARING): sigma_{2*pi}^{L,alpha}(O) = O bit-exact
  L2E-FALS-2: sigma_{2*pi}^{L,TT}(O) = O bit-exact
  L2E-FALS-3: flow conjugacy sigma_t^{L,TT}(O) = sigma_{-t}^{L,alpha}(O)
  RIE-LIMIT (LOAD-BEARING #4): at N_t = 1, K_L_alpha = K_alpha and
                                K_L_TT = K_TT bit-identically.

Plus the Paper 42 §7.2 / O3 H_local verdict at (3, 1) and six-witness
collapse at the Krein level (Paper 42 §8 lifted).

Sprint L2-E (2026-05-16): four-witness Wick-rotation theorem at the
Krein level on a hemispheric wedge of S^3 x R at signature (3, 1).
"""

from __future__ import annotations

import numpy as np
import pytest

from geovac.full_dirac_operator_system import (
    FullDiracTruncatedOperatorSystem,
    full_dirac_basis,
)
from geovac.krein_space_construction import KreinSpace
from geovac.modular_hamiltonian_lorentzian import (
    LorentzianModularHamiltonian,
    LorentzianTomitaStructure,
    LorentzianWedge,
    build_K_L_alpha_full,
    for_bisognano_wichmann_lorentzian,
    for_hartle_hawking_lorentzian,
    for_sewell_lorentzian,
    for_unruh_lorentzian,
    restrict_K_L_alpha_to_wedge,
    restrict_operator_to_wedge_block,
    temporal_positive_half_projector,
    verify_cross_witness_collapse_lorentzian,
)


# ---------------------------------------------------------------------------
# Temporal projector
# ---------------------------------------------------------------------------


def test_temporal_positive_half_projector_N_t_1_is_identity():
    """At N_t = 1, the temporal grid is {t=0} and P_t = I_1."""
    P_t = temporal_positive_half_projector(1)
    assert P_t.shape == (1, 1)
    assert np.allclose(P_t, np.eye(1))


def test_temporal_positive_half_projector_N_t_11_dim():
    """At N_t = 11 (symmetric grid t = linspace(-1, 1, 11)), 6 points have
    t >= 0 (indices 5..10), so dim(P_t) image is 6."""
    P_t = temporal_positive_half_projector(11)
    assert P_t.shape == (11, 11)
    assert abs(np.trace(P_t).real - 6.0) < 1e-12


def test_temporal_positive_half_projector_N_t_21():
    """N_t = 21 symmetric grid: 11 points have t >= 0."""
    P_t = temporal_positive_half_projector(21)
    assert abs(np.trace(P_t).real - 11.0) < 1e-12


def test_temporal_projector_idempotent():
    """P_t^2 = P_t for any N_t."""
    for N_t in (1, 5, 11, 21):
        P_t = temporal_positive_half_projector(N_t)
        residual = np.linalg.norm(P_t @ P_t - P_t)
        assert residual < 1e-14, (
            f"N_t={N_t}: P_t^2 - P_t residual = {residual:.2e}"
        )


def test_temporal_projector_invalid_input():
    with pytest.raises(ValueError, match="N_t"):
        temporal_positive_half_projector(0)


# ---------------------------------------------------------------------------
# LorentzianWedge — wedge construction
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max,N_t", [(1, 1), (1, 11), (2, 1), (2, 11), (3, 1), (3, 11)])
def test_wedge_idempotent(n_max, N_t):
    krein = KreinSpace(n_max=n_max, N_t=N_t, T_max=1.0)
    wedge = LorentzianWedge(krein=krein)
    ok, residual = wedge.verify_idempotent()
    assert ok, f"n_max={n_max}, N_t={N_t}: P_W_L^2 != P_W_L, res={residual:.2e}"


@pytest.mark.parametrize("n_max,N_t", [(1, 1), (1, 11), (2, 1), (2, 11), (3, 1), (3, 11)])
def test_wedge_hermitian(n_max, N_t):
    krein = KreinSpace(n_max=n_max, N_t=N_t, T_max=1.0)
    wedge = LorentzianWedge(krein=krein)
    ok, residual = wedge.verify_hermitian()
    assert ok, f"n_max={n_max}, N_t={N_t}: P_W_L not Hermitian, res={residual:.2e}"


def test_wedge_dim_factorizes():
    """dim_W_L = dim_W_spatial * N_t_positive."""
    for n_max in (1, 2, 3):
        for N_t in (1, 11, 21):
            krein = KreinSpace(n_max=n_max, N_t=N_t, T_max=1.0)
            wedge = LorentzianWedge(krein=krein)
            expected = wedge.dim_W_spatial * wedge.N_t_positive
            assert wedge.dim_W_L == expected, (
                f"n_max={n_max}, N_t={N_t}: dim_W_L={wedge.dim_W_L}, "
                f"expected={expected}"
            )


def test_wedge_at_N_t_1_is_spatial_wedge():
    """At N_t = 1, P_W_L collapses to spatial P_W * 1."""
    krein = KreinSpace(n_max=2, N_t=1, T_max=1.0)
    wedge = LorentzianWedge(krein=krein)
    assert wedge.dim_W_L == wedge.dim_W_spatial
    # Trace check
    assert abs(np.trace(wedge.P_W_L).real - wedge.dim_W_spatial) < 1e-12


# ---------------------------------------------------------------------------
# K_L_alpha geometric boost generator
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max,N_t", [(1, 1), (1, 11), (2, 1), (2, 11), (3, 1)])
def test_K_L_alpha_full_Hermitian_diagonal(n_max, N_t):
    """K_L_alpha is Hermitian and diagonal."""
    krein = KreinSpace(n_max=n_max, N_t=N_t, T_max=1.0)
    K_L = build_K_L_alpha_full(krein)
    # Hermitian
    assert np.linalg.norm(K_L - K_L.conj().T) < 1e-14
    # Diagonal
    off_diag = K_L - np.diag(np.diag(K_L))
    assert np.linalg.norm(off_diag) < 1e-14


def test_K_L_alpha_full_integer_spectrum():
    """K_L_alpha has odd-integer eigenvalues (inherited from two_m_j)."""
    krein = KreinSpace(n_max=2, N_t=11, T_max=1.0)
    K_L = build_K_L_alpha_full(krein)
    eigs = np.diag(K_L).real
    for e in eigs:
        assert abs(e - round(e)) < 1e-14, f"non-integer eigenvalue {e}"
        # Odd integer
        assert int(round(e)) % 2 != 0, f"non-odd eigenvalue {e}"


@pytest.mark.parametrize("n_max,N_t", [(2, 1), (2, 11), (3, 1)])
def test_K_L_alpha_wedge_diag_integer(n_max, N_t):
    """K_L_alpha^W has positive odd-integer eigenvalues on the wedge."""
    krein = KreinSpace(n_max=n_max, N_t=N_t, T_max=1.0)
    wedge = LorentzianWedge(krein=krein)
    K_L_W = restrict_K_L_alpha_to_wedge(krein, wedge)
    eigs = np.diag(K_L_W).real
    for e in eigs:
        assert e > 0, f"non-positive eigenvalue {e} on wedge"
        assert int(round(e)) % 2 != 0, f"non-odd eigenvalue {e}"


def test_K_L_alpha_at_N_t_1_equals_spatial_K_alpha():
    """At N_t = 1, K_L_alpha bit-identically = K_alpha (Paper 42 §5)."""
    n_max = 2
    krein = KreinSpace(n_max=n_max, N_t=1, T_max=1.0)
    K_L = build_K_L_alpha_full(krein)
    # Reference: spatial K_alpha = diag(two_m_j) on the full Dirac basis
    basis = full_dirac_basis(n_max)
    K_spatial = np.diag(
        np.array([float(b.two_m_j) for b in basis], dtype=np.complex128)
    )
    residual = np.linalg.norm(K_L - K_spatial)
    assert residual < 1e-14, (
        f"K_L != K_alpha at N_t=1, residual={residual:.2e}"
    )


# ---------------------------------------------------------------------------
# L2E-FALS-1 (LOAD-BEARING): BW-alpha period closure
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max,N_t", [(2, 1), (2, 11), (3, 1)])
def test_L2E_FALS_1_BW_alpha_period_closure(n_max, N_t):
    """L2E-FALS-1 (LOAD-BEARING): sigma_{2*pi}^{L,alpha}(O) = O bit-exact.

    Direct Lorentzian-side analog of Paper 42 Theorem 5.4.
    """
    lmh = for_bisognano_wichmann_lorentzian(n_max=n_max, N_t=N_t, T_max=1.0)
    results = lmh.verify_witness_lorentzian()
    assert results["verdict"] == "STRONG_IDENTIFICATION_LORENTZIAN", (
        f"n_max={n_max}, N_t={N_t}: verdict={results['verdict']}, "
        f"max_alpha={results['max_alpha_residual']:.2e}"
    )
    # Bit-exact threshold: residual ~ O(sqrt(dim) * eps_machine)
    assert results["max_alpha_residual"] < 1e-12, (
        f"n_max={n_max}, N_t={N_t}: max_alpha={results['max_alpha_residual']:.2e}"
    )


@pytest.mark.slow
def test_L2E_FALS_1_BW_alpha_period_closure_n_max_3_N_t_11():
    """L2E-FALS-1 at n_max=3, N_t=11 (dim_W_L=120, GNS=14400 — slow)."""
    lmh = for_bisognano_wichmann_lorentzian(n_max=3, N_t=11, T_max=1.0)
    results = lmh.verify_witness_lorentzian()
    assert results["max_alpha_residual"] < 1e-12


# ---------------------------------------------------------------------------
# L2E-FALS-2: BW-gamma Tomita-Takesaki period closure
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max,N_t", [(2, 1), (2, 11), (3, 1)])
def test_L2E_FALS_2_BW_gamma_period_closure(n_max, N_t):
    """L2E-FALS-2: sigma_{2*pi}^{L,TT}(O) = O bit-exact via Tomita polar
    decomposition on Krein-GNS Hilbert-Schmidt space.

    Lorentzian-side analog of Paper 42 Theorem 6.3.
    """
    lmh = for_bisognano_wichmann_lorentzian(n_max=n_max, N_t=N_t, T_max=1.0)
    results = lmh.verify_witness_lorentzian()
    assert results["max_TT_residual"] < 1e-12, (
        f"n_max={n_max}, N_t={N_t}: max_TT={results['max_TT_residual']:.2e}"
    )


# ---------------------------------------------------------------------------
# L2E-FALS-3: Flow conjugacy
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max,N_t", [(2, 1), (2, 11), (3, 1)])
def test_L2E_FALS_3_flow_conjugacy(n_max, N_t):
    """L2E-FALS-3: sigma_t^{L,TT}(O) = sigma_{-t}^{L,alpha}(O) bit-exact at general t.

    Lorentzian-side analog of Paper 42 Theorem 7.1.
    """
    lmh = for_bisognano_wichmann_lorentzian(n_max=n_max, N_t=N_t, T_max=1.0)
    results = lmh.verify_witness_lorentzian()
    assert results["max_conjugacy_residual"] < 1e-12, (
        f"n_max={n_max}, N_t={N_t}: max_conj={results['max_conjugacy_residual']:.2e}"
    )


def test_L2E_FALS_3_conjugacy_at_multiple_t():
    """Flow conjugacy at t in {1, pi, 2*pi}.

    Per Paper 42 §7.1: conjugacy is bit-exact at all real t (not just at
    the period 2*pi).
    """
    lmh = for_bisognano_wichmann_lorentzian(n_max=2, N_t=11, T_max=1.0)
    gens = lmh._default_generators(n_test=3)
    for t in (1.0, np.pi, 2.0 * np.pi):
        for O in gens:
            O_W = restrict_operator_to_wedge_block(O, lmh.krein, lmh.wedge)
            ok, residual = lmh.verify_flow_conjugacy_lorentzian(O_W, t=t)
            assert ok, (
                f"conjugacy at t={t}: residual={residual:.2e}"
            )


# ---------------------------------------------------------------------------
# Riemannian-limit recovery (LOAD-BEARING falsifier #4)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_riemannian_limit_recovery_LOAD_BEARING(n_max):
    """At N_t = 1, K_L_alpha and K_L_TT reduce bit-identically to K_alpha
    and K_TT of geovac.modular_hamiltonian.

    LOAD-BEARING falsifier #4: if this fails, the Lorentzian extension
    is structurally different from the Riemannian construction even at
    the reduction limit. Would re-open Paper 42 closure at the
    operator-system level.
    """
    lmh = for_bisognano_wichmann_lorentzian(n_max=n_max, N_t=11, T_max=1.0)
    ok, details = lmh.riemannian_limit_recovery()
    assert ok, (
        f"n_max={n_max}: Riemannian-limit recovery failed. "
        f"k_alpha_residual={details['k_alpha_W_residual']:.2e}, "
        f"k_tomita_residual={details['k_tomita_residual']:.2e}"
    )
    # Bit-identical (zero) is the load-bearing target
    assert details["k_alpha_W_residual"] == 0.0
    assert details["k_tomita_residual"] == 0.0


@pytest.mark.parametrize("n_max,N_t", [(2, 11), (2, 21), (3, 11)])
def test_riemannian_limit_recovery_at_N_t_gt_1(n_max, N_t):
    """Even at N_t > 1, building the N_t = 1 reduction recovers the
    Riemannian construction bit-exactly (the reduction takes a fresh
    KreinSpace at N_t = 1)."""
    lmh = for_bisognano_wichmann_lorentzian(n_max=n_max, N_t=N_t, T_max=1.0)
    ok, details = lmh.riemannian_limit_recovery()
    assert ok
    assert details["N_t_test"] == 1


# ---------------------------------------------------------------------------
# H_local verdict (Paper 42 §7.2 / O3 at (3, 1))
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_h_local_verdict_at_N_t_1_matches_riemannian(n_max):
    """At N_t = 1, the Lorentzian H_local-vs-D_L_W residual EQUALS the
    Riemannian H_local-vs-D_W residual bit-exactly.

    THE HEADLINE STRUCTURAL FINDING (§9 of memo): Paper 42 §7.2 O3
    holds at (3, 1) as a signature-INDEPENDENT baseline at the Riemannian
    reduction. The framework's intrinsic Dirac is NOT the right local
    Hamiltonian for the BW vacuum at beta = 2*pi at either signature.
    """
    lmh = for_bisognano_wichmann_lorentzian(n_max=n_max, N_t=1, T_max=1.0)
    verdict = lmh.h_local_verdict_at_3_1()
    assert verdict["verdict"] == "ii", (
        f"n_max={n_max}: expected verdict 'ii', got {verdict['verdict']}; "
        f"residual_full={verdict['residual_full']:.4e}"
    )
    assert verdict["lorentzian_eq_riemannian_baseline"], (
        f"n_max={n_max}: Lorentzian residual at N_t=1 should match "
        f"Riemannian baseline bit-exactly"
    )


@pytest.mark.parametrize("n_max,N_t", [(1, 11), (2, 11), (3, 11)])
def test_h_local_verdict_at_N_t_gt_1_refined_by_temporal(n_max, N_t):
    """At N_t > 1, the Lorentzian residual EXCEEDS the Riemannian baseline
    because D_L has the temporal-derivative piece i * gamma^0 (x) d/dt
    which H_local does not capture.

    This is verdict (iii) under the refined classification: Paper 42 §7.2
    O3 is PRESERVED at the Riemannian limit and REFINED at N_t > 1 by the
    temporal-derivative contribution.
    """
    lmh = for_bisognano_wichmann_lorentzian(n_max=n_max, N_t=N_t, T_max=1.0)
    verdict = lmh.h_local_verdict_at_3_1()
    assert verdict["verdict"] == "iii", (
        f"n_max={n_max}, N_t={N_t}: expected verdict 'iii' "
        f"(refined by temporal-derivative), got {verdict['verdict']}"
    )
    # Lorentzian residual normalized > Riemannian residual normalized
    assert (
        verdict["residual_full_normalized"]
        > verdict["residual_RIE_normalized"]
    ), (
        f"n_max={n_max}, N_t={N_t}: residual_full_normalized "
        f"{verdict['residual_full_normalized']:.4e} should exceed "
        f"Riemannian baseline {verdict['residual_RIE_normalized']:.4e}"
    )


@pytest.mark.parametrize("n_max,N_t", [(1, 1), (2, 1), (2, 11), (3, 1)])
def test_h_local_NOT_eq_D_L_W(n_max, N_t):
    """H_local := K_L_alpha^W / beta is NOT D_L_W bit-exactly at any
    panel cell. This is the structural content of Paper 42 §7.2 lifted
    to (3, 1).
    """
    lmh = for_bisognano_wichmann_lorentzian(
        n_max=n_max, N_t=N_t, T_max=1.0,
    )
    verdict = lmh.h_local_verdict_at_3_1()
    assert verdict["verdict"] != "i", (
        f"n_max={n_max}, N_t={N_t}: verdict 'i' (H_local = D_L_W) "
        f"would refine Paper 42 §7.2 to signature-DEPENDENT; "
        f"this would be a surprise."
    )


# ---------------------------------------------------------------------------
# Six-witness collapse (Paper 42 §8 Krein-level)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max,N_t", [(1, 1), (2, 1), (2, 11), (3, 1)])
def test_six_witness_collapse_lorentzian(n_max, N_t):
    """Six witnesses (BW, HH_M1, HH_M2, Sew_M1, Unruh_a1, Unruh_a2) collapse
    to a single Krein-level construction.

    Paper 42 §8 Corollary 8.1 lifted to (3, 1): rho_W^L is beta-independent
    under the BW choice H_local = K_L_alpha^W / beta, so the modular
    operator Delta_L and K_L_TT are bit-identical across the six witness
    instantiations.
    """
    results = verify_cross_witness_collapse_lorentzian(
        n_max=n_max, N_t=N_t, T_max=1.0,
    )
    assert results["six_witness_collapse_ok"], (
        f"n_max={n_max}, N_t={N_t}: collapse failed; "
        f"max_alpha={results['max_alpha_cross_consistency']:.2e}, "
        f"max_TT={results['max_TT_cross_consistency']:.2e}"
    )
    # Each witness should give bit-identical period residual to BW canonical
    assert results["max_alpha_cross_consistency"] == 0.0
    assert results["max_TT_cross_consistency"] == 0.0


# ---------------------------------------------------------------------------
# LorentzianTomitaStructure (standalone tests)
# ---------------------------------------------------------------------------


def test_tomita_structure_construct():
    """LorentzianTomitaStructure.construct() builds Delta_L and K_L_TT."""
    lmh = for_bisognano_wichmann_lorentzian(n_max=2, N_t=1, T_max=1.0)
    ts = lmh.tomita_structure()
    assert ts.rho is not None
    assert ts.Delta_matrix is not None
    assert ts.K_TT is not None
    # Delta is Hermitian positive on H_GNS
    Delta_residual = np.linalg.norm(ts.Delta_matrix - ts.Delta_matrix.conj().T)
    assert Delta_residual < 1e-12
    eigvals = np.linalg.eigvalsh(ts.Delta_matrix)
    assert (eigvals > 0).all()


def test_tomita_J_TT_squared_is_plus_identity():
    """The Lorentzian Tomita J_L_TT satisfies J^2 = +I (Tomita signature).

    Paper 42 §6.3 Prop 6.2: J_TT^2 = +I is the categorical distinction
    from Connes J_GV (J^2 = -I, KO-dim 3). Same signature holds at (3, 1).
    """
    lmh = for_bisognano_wichmann_lorentzian(n_max=2, N_t=1, T_max=1.0)
    ts = lmh.tomita_structure()
    residual = ts.J_L_TT_squared_residual()
    assert residual < 1e-10, (
        f"J_L_TT^2 != +I, residual={residual:.2e}"
    )


def test_tomita_modular_periodicity():
    """sigma_{2*pi}^{L,TT}(a) = a bit-exact for the BW choice H_local =
    K_L_alpha^W / beta."""
    lmh = for_bisognano_wichmann_lorentzian(n_max=2, N_t=1, T_max=1.0)
    ts = lmh.tomita_structure()
    gens = lmh._default_generators(n_test=2)
    for O in gens:
        O_W = restrict_operator_to_wedge_block(O, lmh.krein, lmh.wedge)
        ok, residual = ts.verify_modular_periodicity_tomita_lorentzian(O_W)
        assert ok, f"Tomita period closure: residual={residual:.2e}"


def test_tomita_K_TT_integer_spectrum():
    """K_L_TT on Krein-GNS has integer eigenvalues (= integer differences
    of K_L_alpha^W eigenvalues).

    Paper 42 §6.5 Prop 6.4 lifted to (3, 1): Spec(K_TT) is contained in
    integer differences of two_m_j values, which are themselves integers.
    """
    lmh = for_bisognano_wichmann_lorentzian(n_max=2, N_t=1, T_max=1.0)
    ts = lmh.tomita_structure()
    eigs = np.linalg.eigvalsh(ts.K_TT).real
    for e in eigs:
        # Integer eigenvalues — bit-exactly integer (up to log-precision)
        assert abs(e - round(e)) < 1e-10, f"K_TT eigenvalue {e} not integer"


# ---------------------------------------------------------------------------
# Witness factories
# ---------------------------------------------------------------------------


def test_bw_witness_factory():
    lmh = for_bisognano_wichmann_lorentzian(n_max=2, N_t=11, T_max=1.0)
    assert lmh.kappa_g == 1.0
    assert abs(lmh.beta - 2.0 * np.pi) < 1e-14


def test_hartle_hawking_witness_factory():
    lmh = for_hartle_hawking_lorentzian(n_max=2, N_t=11, T_max=1.0, M=1.0)
    assert abs(lmh.kappa_g - 0.25) < 1e-14
    assert abs(lmh.beta - 8.0 * np.pi) < 1e-14


def test_sewell_witness_factory_equals_HH():
    """Sewell at M = m gives same modular structure as HH at M = m."""
    lmh_hh = for_hartle_hawking_lorentzian(n_max=2, N_t=1, T_max=1.0, M=2.0)
    lmh_sew = for_sewell_lorentzian(n_max=2, N_t=1, T_max=1.0, M=2.0)
    assert lmh_sew.kappa_g == lmh_hh.kappa_g
    assert lmh_sew.beta == lmh_hh.beta


def test_unruh_witness_factory():
    lmh = for_unruh_lorentzian(n_max=2, N_t=11, T_max=1.0, a=2.0)
    assert lmh.kappa_g == 2.0
    assert abs(lmh.beta - np.pi) < 1e-14


# ---------------------------------------------------------------------------
# Full witness verification battery (end-to-end)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max,N_t", [(2, 1), (2, 11), (3, 1)])
def test_full_witness_battery_BW(n_max, N_t):
    """End-to-end: all four falsifiers pass."""
    lmh = for_bisognano_wichmann_lorentzian(n_max=n_max, N_t=N_t, T_max=1.0)
    results = lmh.verify_witness_lorentzian()
    assert results["verdict"] == "STRONG_IDENTIFICATION_LORENTZIAN"
    # Wedge tests
    assert results["wedge_idempotent"]["ok"]
    assert results["wedge_hermitian"]["ok"]
    # All three FALS
    assert results["max_alpha_residual"] < 1e-12
    assert results["max_TT_residual"] < 1e-12
    assert results["max_conjugacy_residual"] < 1e-12


def test_full_witness_battery_HH_at_M_1():
    """Hartle-Hawking at M = 1 also gives strong identification."""
    lmh = for_hartle_hawking_lorentzian(n_max=2, N_t=1, T_max=1.0, M=1.0)
    results = lmh.verify_witness_lorentzian()
    assert results["verdict"] == "STRONG_IDENTIFICATION_LORENTZIAN"


def test_full_witness_battery_Unruh_at_a_2():
    """Unruh at a = 2 also gives strong identification (kappa_g-independent)."""
    lmh = for_unruh_lorentzian(n_max=2, N_t=1, T_max=1.0, a=2.0)
    results = lmh.verify_witness_lorentzian()
    assert results["verdict"] == "STRONG_IDENTIFICATION_LORENTZIAN"


# ---------------------------------------------------------------------------
# Operator wedge-block restriction
# ---------------------------------------------------------------------------


def test_restrict_operator_to_wedge_block_dim():
    """restrict_operator_to_wedge_block returns a (dim_W_L, dim_W_L) matrix."""
    krein = KreinSpace(n_max=2, N_t=11, T_max=1.0)
    wedge = LorentzianWedge(krein=krein)
    O_full = np.eye(krein.dim, dtype=np.complex128)
    O_W = restrict_operator_to_wedge_block(O_full, krein, wedge)
    assert O_W.shape == (wedge.dim_W_L, wedge.dim_W_L)


def test_restrict_operator_to_wedge_identity_preserved():
    """The identity on K restricts to the identity on the wedge."""
    krein = KreinSpace(n_max=2, N_t=11, T_max=1.0)
    wedge = LorentzianWedge(krein=krein)
    I_K = np.eye(krein.dim, dtype=np.complex128)
    I_W = restrict_operator_to_wedge_block(I_K, krein, wedge)
    I_target = np.eye(wedge.dim_W_L, dtype=np.complex128)
    residual = np.linalg.norm(I_W - I_target)
    assert residual < 1e-12, f"Identity restriction residual = {residual:.2e}"


# ---------------------------------------------------------------------------
# Cross-witness consistency at multiple cutoffs (slow check)
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_full_panel_n_max_3_N_t_21_BW_strong():
    """Full panel: n_max=3, N_t=21 BW canonical strong identification."""
    lmh = for_bisognano_wichmann_lorentzian(n_max=3, N_t=21, T_max=1.0)
    results = lmh.verify_witness_lorentzian()
    assert results["verdict"] == "STRONG_IDENTIFICATION_LORENTZIAN", (
        f"verdict={results['verdict']}, "
        f"max_alpha={results['max_alpha_residual']:.2e}, "
        f"max_TT={results['max_TT_residual']:.2e}"
    )


@pytest.mark.slow
def test_six_witness_collapse_n_max_3_N_t_11():
    """Six-witness collapse at (n_max=3, N_t=11)."""
    results = verify_cross_witness_collapse_lorentzian(
        n_max=3, N_t=11, T_max=1.0,
    )
    assert results["six_witness_collapse_ok"]
