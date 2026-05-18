"""Test suite for the joint Lichnerowicz bound on the compact-temporal
Lorentzian Dirac (L3b-2 Sub-Sprint A).

Tests cover:

  1. Bit-exact vanishing of the time-chirality cross commutator
     [gamma^0 (x) d_t, a_s (x) a_t] = 0 for every generator.
  2. Bit-exact structural identity [D_L, a_s (x) a_t] = i [D_GV, a_s] (x) a_t.
  3. Scalar-only multipliers (p = 0) reproduce Paper 38 §L3 bit-exact.
  4. Temporal-only multipliers (Y3_{1,0,0} (x) g_p) have ||[D_L, a]|| = 0.
  5. Joint Lichnerowicz bound holds at every test panel cell under L^1.
  6. Joint Lichnerowicz bound holds at every test panel cell under L^2.
  7. Asymptotic-tightness verification of the closed-form bound.
  8. C_3 ratio is independent of N_t (panel structure check).
  9. Random linear-combination multipliers respect the bound.
 10. Riemannian-limit reduction at N_t = 1 (Term A vanishes trivially).
 11. The temporal commutator [d_t, a_t] = 0 bit-exact for every p.
 12. The spatial commutator [gamma^0, a_s] = 0 bit-exact for every spatial mult.
 13. Operator-norm tensor-product identity ||A (x) B|| = ||A|| * ||B||.
 14. Riemannian limit at N_t = 1: joint commutator reduces to spatial.
"""

from __future__ import annotations

import numpy as np
import pytest

from geovac.full_dirac_operator_system import (
    camporesi_higuchi_full_dirac_matrix,
)
from geovac.krein_space_compact_temporal import CompactTemporalKreinSpace
from geovac.lorentzian_dirac_compact import (
    fourier_d_dt_matrix,
    lorentzian_dirac_compact_matrix,
)
from geovac.operator_system_compact_temporal import (
    CompactTemporalTruncatedOperatorSystem,
    compact_temporal_multiplier_matrices,
)


# ===========================================================================
# Helpers
# ===========================================================================


def _opnorm(M):
    return float(np.linalg.norm(M, ord=2))


def _commutator(A, B):
    return A @ B - B @ A


# ===========================================================================
# Test 1: Time-chirality cross commutator vanishes bit-exact
# ===========================================================================


@pytest.mark.parametrize("n_max,N_t", [(2, 3), (3, 3), (3, 5), (2, 5)])
def test_term_A_bit_exact_vanish(n_max, N_t):
    """[gamma^0 (x) d_t, a_s (x) a_t] = 0 bit-exact for every generator.

    Structural reason: a_s commutes with gamma^0 (block-diagonal multiplier
    vs chirality-swap acts trivially on identical blocks), and a_t commutes
    with d_t (both momentum-diagonal).
    """
    T = 2.0 * np.pi
    krein = CompactTemporalKreinSpace(n_max=n_max, N_t=N_t, T=T)
    op_sys = CompactTemporalTruncatedOperatorSystem(n_max=n_max, N_t=N_t, T=T)
    D_t = fourier_d_dt_matrix(N_t, T)
    time_kernel = np.kron(krein.J_spatial, D_t)  # gamma^0 (x) d_t

    max_res = 0.0
    for M in op_sys.multiplier_matrices:
        c = _commutator(time_kernel, M)
        max_res = max(max_res, _opnorm(c))

    assert max_res < 1e-12, (
        f"Term A residual {max_res:.3e} > 1e-12 at (n_max={n_max}, N_t={N_t})"
    )


# ===========================================================================
# Test 2: Structural identity (2.1) bit-exact
# ===========================================================================


@pytest.mark.parametrize("n_max,N_t", [(2, 3), (3, 3), (3, 5)])
def test_structural_identity_2_1(n_max, N_t):
    """[D_L, a_s (x) a_t] = i [D_GV, a_s] (x) a_t bit-exactly.

    Memo Eq. (2.1).
    """
    T = 2.0 * np.pi
    krein = CompactTemporalKreinSpace(n_max=n_max, N_t=N_t, T=T)
    op_sys = CompactTemporalTruncatedOperatorSystem(n_max=n_max, N_t=N_t, T=T)
    D_L = lorentzian_dirac_compact_matrix(krein)
    D_GV = camporesi_higuchi_full_dirac_matrix(krein.basis_spatial)
    temp_mats = compact_temporal_multiplier_matrices(N_t, T)
    spat_to_idx = {lbl: i for i, lbl in enumerate(op_sys.spat_labels)}

    max_res = 0.0
    for (N, L, M_lab, p), full_mat in zip(
        op_sys.multiplier_labels, op_sys.multiplier_matrices
    ):
        sidx = spat_to_idx[(N, L, M_lab)]
        a_s = op_sys._spat_matrices[sidx]
        a_t = temp_mats[p]
        lhs = _commutator(D_L, full_mat)
        rhs = 1j * np.kron(_commutator(D_GV, a_s), a_t)
        res = float(np.linalg.norm(lhs - rhs))
        max_res = max(max_res, res)

    assert max_res < 1e-10, (
        f"Structural identity residual {max_res:.3e} > 1e-10"
    )


# ===========================================================================
# Test 3: Scalar-only multipliers reproduce Paper 38 §L3
# ===========================================================================


@pytest.mark.parametrize("n_max,N_t", [(2, 3), (3, 5)])
def test_scalar_only_p0_multipliers_match_spatial(n_max, N_t):
    """For p = 0 multipliers (a_t = I_{N_t}), ||[D_L, a]||_op = ||[D_GV, a_s]||_op."""
    T = 2.0 * np.pi
    krein = CompactTemporalKreinSpace(n_max=n_max, N_t=N_t, T=T)
    op_sys = CompactTemporalTruncatedOperatorSystem(n_max=n_max, N_t=N_t, T=T)
    D_L = lorentzian_dirac_compact_matrix(krein)
    D_GV = camporesi_higuchi_full_dirac_matrix(krein.basis_spatial)
    spat_to_idx = {lbl: i for i, lbl in enumerate(op_sys.spat_labels)}

    for (N, L, M_lab, p), full_mat in zip(
        op_sys.multiplier_labels, op_sys.multiplier_matrices
    ):
        if p != 0:
            continue
        sidx = spat_to_idx[(N, L, M_lab)]
        a_s = op_sys._spat_matrices[sidx]
        joint_norm = _opnorm(_commutator(D_L, full_mat))
        spatial_norm = _opnorm(_commutator(D_GV, a_s))
        # ||a_t||_op = 1 (a_t = I_{N_t} at p=0), so joint = spatial
        assert abs(joint_norm - spatial_norm) < 1e-10, (
            f"p=0 mismatch at (N,L,M)={(N,L,M_lab)}: {joint_norm} vs {spatial_norm}"
        )


# ===========================================================================
# Test 4: Temporal-only (Y3_{1,0,0} (x) g_p) gives zero commutator
# ===========================================================================


@pytest.mark.parametrize("n_max,N_t", [(2, 3), (3, 5)])
def test_temporal_only_multiplier_zero_commutator(n_max, N_t):
    """For Y3_{1,0,0} = constant on S^3, the multiplier a_s = I on H_GV.

    Hence [D_GV, I] = 0 and the joint commutator vanishes.
    """
    T = 2.0 * np.pi
    krein = CompactTemporalKreinSpace(n_max=n_max, N_t=N_t, T=T)
    op_sys = CompactTemporalTruncatedOperatorSystem(n_max=n_max, N_t=N_t, T=T)
    D_L = lorentzian_dirac_compact_matrix(krein)

    n_temporal_only_tested = 0
    for (N, L, M_lab, p), full_mat in zip(
        op_sys.multiplier_labels, op_sys.multiplier_matrices
    ):
        if N != 1:
            continue
        c = _commutator(D_L, full_mat)
        norm_c = _opnorm(c)
        assert norm_c < 1e-12, (
            f"Temporal-only multiplier ({N},{L},{M_lab},p={p}) "
            f"has [D_L, a] = {norm_c}"
        )
        n_temporal_only_tested += 1
    # At least one temporal-only multiplier should exist (p=0..N_t-1 for N=1)
    assert n_temporal_only_tested >= 1, "No temporal-only multipliers found"


# ===========================================================================
# Test 5: Joint Lichnerowicz under L^1 (single-tensor cases)
# ===========================================================================


@pytest.mark.parametrize("n_max,N_t", [(2, 3), (3, 5)])
def test_joint_lichnerowicz_L1_single_tensor(n_max, N_t):
    """For pure-tensor a = a_s (x) a_t: ||[D_L, a]||_op <= bound * RHS_L1.

    Where bound = sup_{N<=n_max} (N-1)/sqrt(N^2-1) and
    RHS_L1 = ||grad_x f_s||_inf * ||f_t||_inf + ||f_s||_inf * ||d_t f_t||_inf.

    Since the LHS = spatial Paper 38 bound x ||a_t||_op, and ||a_t||_op =
    ||f_t||_inf for momentum-diagonal a_t, the L^1 RHS upper bounds the LHS.
    """
    from geovac.r25_l3_lipschitz_bound import lipschitz_norm_inf_y3

    T = 2.0 * np.pi
    krein = CompactTemporalKreinSpace(n_max=n_max, N_t=N_t, T=T)
    op_sys = CompactTemporalTruncatedOperatorSystem(n_max=n_max, N_t=N_t, T=T)
    D_L = lorentzian_dirac_compact_matrix(krein)

    paper_38_bound = max(
        (N - 1) / float(np.sqrt(N**2 - 1)) for N in range(2, n_max + 1)
    )

    # Test a few low-N tensor products
    for label, full_mat in zip(
        op_sys.multiplier_labels[:15], op_sys.multiplier_matrices[:15]
    ):
        N, L, M_lab, p = label
        if N == 1:
            continue  # constant spatial multiplier, ratio is 0/0
        # LHS
        lhs = _opnorm(_commutator(D_L, full_mat))
        if lhs < 1e-14:
            continue  # vacuous
        # Spatial Lipschitz of Y3
        spat_lip = float(lipschitz_norm_inf_y3(N, L, abs(M_lab), prec=20))
        # f_t L^infinity (= ||a_t||_op for momentum-diagonal)
        if p == 0:
            f_t_inf = 1.0
            f_t_deriv = 0.0
        else:
            omegas = np.array([
                2.0 * np.pi * k / T for k in krein.momentum_grid
            ])
            f_t_inf = float(np.max(np.abs(omegas) ** p))
            f_t_deriv = float(np.max(
                p * (np.abs(omegas) ** (p - 1)) * (2.0 * np.pi / T)
            ))
        # ||f_s||_inf upper-bounded by 1 (Y3 unit harmonic, rough)
        f_s_inf = 1.0
        rhs_L1 = spat_lip * f_t_inf + f_s_inf * f_t_deriv
        if rhs_L1 < 1e-14:
            continue
        ratio = lhs / rhs_L1
        # Bound holds with margin (the rough ||f_s||_inf ≤ 1 makes the bound
        # only slightly tighter; the L^1 bound is by construction looser
        # than the spatial-only bound)
        assert ratio <= paper_38_bound + 1e-8, (
            f"L^1 bound violated at {label}: {ratio:.4f} > {paper_38_bound:.4f}"
        )


# ===========================================================================
# Test 6: Joint Lichnerowicz under L^2 (single-tensor cases)
# ===========================================================================


@pytest.mark.parametrize("n_max,N_t", [(2, 3), (3, 5)])
def test_joint_lichnerowicz_L2_single_tensor(n_max, N_t):
    """For pure-tensor a = a_s (x) a_t under the L^2 Pythagorean norm."""
    from geovac.r25_l3_lipschitz_bound import lipschitz_norm_inf_y3

    T = 2.0 * np.pi
    krein = CompactTemporalKreinSpace(n_max=n_max, N_t=N_t, T=T)
    op_sys = CompactTemporalTruncatedOperatorSystem(n_max=n_max, N_t=N_t, T=T)
    D_L = lorentzian_dirac_compact_matrix(krein)

    paper_38_bound = max(
        (N - 1) / float(np.sqrt(N**2 - 1)) for N in range(2, n_max + 1)
    )

    for label, full_mat in zip(
        op_sys.multiplier_labels[:15], op_sys.multiplier_matrices[:15]
    ):
        N, L, M_lab, p = label
        if N == 1:
            continue
        lhs = _opnorm(_commutator(D_L, full_mat))
        if lhs < 1e-14:
            continue
        spat_lip = float(lipschitz_norm_inf_y3(N, L, abs(M_lab), prec=20))
        if p == 0:
            f_t_inf = 1.0
            f_t_deriv = 0.0
        else:
            omegas = np.array([
                2.0 * np.pi * k / T for k in krein.momentum_grid
            ])
            f_t_inf = float(np.max(np.abs(omegas) ** p))
            f_t_deriv = float(np.max(
                p * (np.abs(omegas) ** (p - 1)) * (2.0 * np.pi / T)
            ))
        f_s_inf = 1.0
        rhs_L2 = float(np.sqrt(
            (spat_lip * f_t_inf) ** 2 + (f_s_inf * f_t_deriv) ** 2
        ))
        if rhs_L2 < 1e-14:
            continue
        ratio = lhs / rhs_L2
        assert ratio <= paper_38_bound + 1e-8, (
            f"L^2 bound violated at {label}: {ratio:.4f} > {paper_38_bound:.4f}"
        )


# ===========================================================================
# Test 7: Asymptotic tightness of the closed-form bound
# ===========================================================================


def test_asymptotic_tightness_closed_form_bound():
    """The closed-form bound (N-1)/sqrt(N^2-1) is monotone increasing in N
    and approaches 1 from below as N -> infinity.
    """
    bounds = [
        (N - 1) / float(np.sqrt(N**2 - 1)) for N in range(2, 200)
    ]
    # Strictly increasing
    for i in range(1, len(bounds)):
        assert bounds[i] > bounds[i - 1], (
            f"Bound not monotone at N={i + 2}"
        )
    # All < 1
    assert all(b < 1.0 for b in bounds)
    # Approaches 1
    assert bounds[-1] > 0.99
    assert bounds[-1] < 1.0
    # Numerical values
    assert abs(bounds[0] - 1.0 / np.sqrt(3)) < 1e-10  # N=2: 1/sqrt(3)
    assert abs(bounds[1] - 1.0 / np.sqrt(2)) < 1e-10  # N=3: 1/sqrt(2)
    assert abs(bounds[2] - 3.0 / np.sqrt(15)) < 1e-10  # N=4: 3/sqrt(15)


# ===========================================================================
# Test 8: C_3 bound is INDEPENDENT of N_t (joint constant is dictated by SU(2))
# ===========================================================================


@pytest.mark.slow
@pytest.mark.parametrize("n_max", [2, 3])
def test_C3_independent_of_Nt(n_max):
    """For fixed n_max, the empirical spatial-only C_3 must be the same
    across different N_t values, because the joint commutator factors as
    [D_L, a] = i [D_GV, a_s] (x) a_t (memo Eq. 2.1).

    This test verifies the structural identity at the *ratio* level:
    the spatial-only ratio ||[D_GV, a_s]|| / ||grad f_s||_inf is the
    same regardless of N_t (since the spatial factor is N_t-independent).
    """
    from geovac.r25_l3_lipschitz_bound import lipschitz_norm_inf_y3

    T = 2.0 * np.pi
    paper_38_bound = max(
        (N - 1) / float(np.sqrt(N**2 - 1)) for N in range(2, n_max + 1)
    )
    # Cache Lipschitz norms
    lip_cache = {}
    max_ratios = []
    for N_t in [1, 3, 5]:
        krein = CompactTemporalKreinSpace(n_max=n_max, N_t=N_t, T=T)
        op_sys = CompactTemporalTruncatedOperatorSystem(n_max=n_max, N_t=N_t, T=T)
        D_GV = camporesi_higuchi_full_dirac_matrix(krein.basis_spatial)
        spat_to_idx = {lbl: i for i, lbl in enumerate(op_sys.spat_labels)}

        # Only test scalar (p=0) generators to get the spatial-only ratio
        max_ratio = 0.0
        for label, full_mat in zip(
            op_sys.multiplier_labels, op_sys.multiplier_matrices
        ):
            N, L, M_lab, p = label
            if N == 1 or p != 0:
                continue
            sidx = spat_to_idx[(N, L, M_lab)]
            a_s = op_sys._spat_matrices[sidx]
            spatial_lhs = _opnorm(_commutator(D_GV, a_s))
            key = (N, L, abs(M_lab))
            if key not in lip_cache:
                lip_cache[key] = float(
                    lipschitz_norm_inf_y3(N, L, abs(M_lab), prec=20)
                )
            spat_lip = lip_cache[key]
            if spat_lip < 1e-14:
                continue
            ratio = spatial_lhs / spat_lip
            max_ratio = max(max_ratio, ratio)
        max_ratios.append(max_ratio)

    # The N_t-independent fact: max_ratios should be identical at every N_t
    for i in range(1, len(max_ratios)):
        assert abs(max_ratios[i] - max_ratios[0]) < 1e-10, (
            f"Spatial ratio at N_t={[1,3,5][i]} differs from N_t=1: "
            f"{max_ratios[i]} vs {max_ratios[0]}"
        )
    # And bounded by asymptotic 1
    for ratio in max_ratios:
        assert ratio <= 1.0 + 1e-8, (
            f"Spatial-only ratio {ratio} exceeds asymptotic 1 at n_max={n_max}"
        )


# ===========================================================================
# Test 9: Random linear-combination multipliers respect the bound
# ===========================================================================


@pytest.mark.parametrize("n_max,N_t,seed", [(2, 3, 7), (3, 5, 11)])
def test_random_combinations_bound(n_max, N_t, seed):
    """Random linear combinations of generators respect the bound via the
    triangle-inequality form.
    """
    from geovac.r25_l3_lipschitz_bound import lipschitz_norm_inf_y3

    T = 2.0 * np.pi
    rng = np.random.default_rng(seed)
    krein = CompactTemporalKreinSpace(n_max=n_max, N_t=N_t, T=T)
    op_sys = CompactTemporalTruncatedOperatorSystem(n_max=n_max, N_t=N_t, T=T)
    D_L = lorentzian_dirac_compact_matrix(krein)

    paper_38_bound = max(
        (N - 1) / float(np.sqrt(N**2 - 1)) for N in range(2, n_max + 1)
    )

    # Cache Lipschitz norms
    spat_lip_cache = {}
    for label in op_sys.multiplier_labels:
        N, L, M_lab, p = label
        key = (N, L, abs(M_lab))
        if key not in spat_lip_cache:
            if N >= 2:
                spat_lip_cache[key] = float(
                    lipschitz_norm_inf_y3(N, L, abs(M_lab), prec=20)
                )
            else:
                spat_lip_cache[key] = 0.0

    n_gens = len(op_sys.multiplier_matrices)
    omegas = np.array([
        2.0 * np.pi * k / T for k in krein.momentum_grid
    ])

    for trial in range(5):
        # Random 3-5 sparsity linear combination
        n_nonzero = int(rng.integers(3, 6))
        indices = rng.choice(n_gens, size=n_nonzero, replace=False)
        coeffs = rng.normal(0, 1, n_nonzero) + 1j * rng.normal(0, 1, n_nonzero)

        # Build matrix
        M = np.zeros_like(op_sys.multiplier_matrices[0])
        labels_used = []
        for j, idx in enumerate(indices):
            M = M + coeffs[j] * op_sys.multiplier_matrices[idx]
            labels_used.append((op_sys.multiplier_labels[idx], coeffs[j]))

        lhs = _opnorm(_commutator(D_L, M))
        if lhs < 1e-14:
            continue

        # RHS_L1 via triangle inequality
        rhs_L1 = 0.0
        for (label, c) in labels_used:
            N, L, M_lab, p = label
            spat_lip = spat_lip_cache.get((N, L, abs(M_lab)), 0.0)
            f_s_inf = 1.0
            if p == 0:
                f_t_inf = 1.0
                f_t_deriv = 0.0
            else:
                f_t_inf = float(np.max(np.abs(omegas) ** p))
                f_t_deriv = float(np.max(
                    p * (np.abs(omegas) ** (p - 1)) * (2.0 * np.pi / T)
                ))
            rhs_L1 += abs(c) * (spat_lip * f_t_inf + f_s_inf * f_t_deriv)

        if rhs_L1 < 1e-14:
            continue
        ratio = lhs / rhs_L1
        assert ratio <= paper_38_bound + 1e-6, (
            f"Random combo ratio {ratio:.4f} > bound {paper_38_bound:.4f}"
        )


# ===========================================================================
# Test 10: Riemannian-limit reduction at N_t = 1
# ===========================================================================


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_riemannian_limit_term_A_vanish(n_max):
    """At N_t = 1: gamma^0 (x) d_t = gamma^0 (x) 0 = 0, so Term A is
    trivially zero. The joint commutator reduces to i [D_GV, a_s] (x) 1.
    """
    T = 2.0 * np.pi
    N_t = 1
    krein = CompactTemporalKreinSpace(n_max=n_max, N_t=N_t, T=T)
    op_sys = CompactTemporalTruncatedOperatorSystem(n_max=n_max, N_t=N_t, T=T)
    D_t = fourier_d_dt_matrix(N_t, T)
    # At N_t=1, D_t = 0 trivially
    assert D_t.shape == (1, 1)
    assert abs(D_t[0, 0]) < 1e-15, "D_t at N_t=1 should be 0"
    # Term A residual is zero at every generator
    time_kernel = np.kron(krein.J_spatial, D_t)
    for M in op_sys.multiplier_matrices:
        assert _opnorm(_commutator(time_kernel, M)) < 1e-15


# ===========================================================================
# Test 11: Temporal commutator [d_t, a_t] = 0 bit-exact for every p
# ===========================================================================


@pytest.mark.parametrize("N_t", [3, 5, 7])
def test_temporal_commutator_vanishes(N_t):
    """Both d_t and g_p are diagonal in the momentum basis, so their
    commutator vanishes identically.
    """
    T = 2.0 * np.pi
    D_t = fourier_d_dt_matrix(N_t, T)
    temp_mats = compact_temporal_multiplier_matrices(N_t, T)
    for p in range(N_t):
        c = D_t @ temp_mats[p] - temp_mats[p] @ D_t
        norm = _opnorm(c)
        assert norm < 1e-14, (
            f"[d_t, g_{p}] = {norm} at N_t={N_t} (expected 0)"
        )


# ===========================================================================
# Test 12: Spatial commutator [gamma^0, a_s] = 0 bit-exact for every a_s
# ===========================================================================


@pytest.mark.parametrize("n_max", [2, 3])
def test_spatial_chirality_commutator_vanishes(n_max):
    """a_s = blkdiag(W, W) commutes with gamma^0 = [[0,I],[I,0]].

    Verifies the structural claim in §2.1 Term A analysis.
    """
    T = 2.0 * np.pi
    N_t = 1  # only spatial structure matters
    krein = CompactTemporalKreinSpace(n_max=n_max, N_t=N_t, T=T)
    op_sys = CompactTemporalTruncatedOperatorSystem(n_max=n_max, N_t=N_t, T=T)
    gamma0 = krein.J_spatial
    for a_s in op_sys._spat_matrices:
        c = gamma0 @ a_s - a_s @ gamma0
        assert _opnorm(c) < 1e-14, (
            f"[gamma^0, a_s] = {_opnorm(c)} (expected 0)"
        )


# ===========================================================================
# Test 13: Operator-norm tensor-product identity
# ===========================================================================


def test_opnorm_tensor_product_identity():
    """For any matrices A, B: ||A (x) B||_op = ||A||_op * ||B||_op."""
    rng = np.random.default_rng(42)
    A = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    B = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))

    n_A = _opnorm(A)
    n_B = _opnorm(B)
    n_AB = _opnorm(np.kron(A, B))
    assert abs(n_AB - n_A * n_B) < 1e-10, (
        f"||A (x) B|| = {n_AB} != {n_A} * {n_B} = {n_A * n_B}"
    )


# ===========================================================================
# Test 14: Riemannian limit at N_t=1 — joint commutator IS spatial
# ===========================================================================


@pytest.mark.slow
@pytest.mark.parametrize("n_max,N_t", [(4, 7)])
def test_joint_lichnerowicz_asymptotic_bound(n_max, N_t):
    """Asymptotic bound: every panel ratio is <= 1 + eps.

    This is the rigorous claim of Theorem 3.1 / 4.1 at the strong
    asymptotic level. At finite cutoff (n_max <= 4 tested here), the
    empirical C_3 may exceed the closed-form per-harmonic SUP bound
    sup_{N <= n_max} (N-1)/sqrt(N^2-1) on multi-shell combinations
    while staying strictly below 1.
    """
    from geovac.r25_l3_lipschitz_bound import lipschitz_norm_inf_y3

    T = 2.0 * np.pi
    krein = CompactTemporalKreinSpace(n_max=n_max, N_t=N_t, T=T)
    op_sys = CompactTemporalTruncatedOperatorSystem(n_max=n_max, N_t=N_t, T=T)
    D_L = lorentzian_dirac_compact_matrix(krein)

    omegas = np.array([2.0 * np.pi * k / T for k in krein.momentum_grid])

    for label, full_mat in zip(
        op_sys.multiplier_labels[:25], op_sys.multiplier_matrices[:25]
    ):
        N, L, M_lab, p = label
        if N == 1:
            continue
        lhs = _opnorm(_commutator(D_L, full_mat))
        if lhs < 1e-14:
            continue
        spat_lip = float(lipschitz_norm_inf_y3(N, L, abs(M_lab), prec=20))
        if p == 0:
            f_t_inf = 1.0
            f_t_deriv = 0.0
        else:
            f_t_inf = float(np.max(np.abs(omegas) ** p))
            f_t_deriv = float(np.max(
                p * (np.abs(omegas) ** (p - 1)) * (2.0 * np.pi / T)
            ))
        f_s_inf = 1.0
        rhs_L1 = spat_lip * f_t_inf + f_s_inf * f_t_deriv
        if rhs_L1 < 1e-14:
            continue
        ratio = lhs / rhs_L1
        # Asymptotic claim: ratio < 1 + eps
        assert ratio <= 1.0 + 1e-8, (
            f"Asymptotic bound violated at {label}: {ratio:.4f} > 1"
        )


@pytest.mark.parametrize("n_max", [2, 3])
def test_riemannian_limit_joint_equals_spatial(n_max):
    """At N_t = 1: [D_L, a_s (x) 1] = i [D_GV, a_s] (x) 1.

    So ||[D_L, a]||_op = ||[D_GV, a_s]||_op.

    This reduces the joint Lichnerowicz to the spatial Paper 38 §L3 EXACTLY
    at the Riemannian boundary.
    """
    T = 2.0 * np.pi
    N_t = 1
    krein = CompactTemporalKreinSpace(n_max=n_max, N_t=N_t, T=T)
    op_sys = CompactTemporalTruncatedOperatorSystem(n_max=n_max, N_t=N_t, T=T)
    D_L = lorentzian_dirac_compact_matrix(krein)
    D_GV = camporesi_higuchi_full_dirac_matrix(krein.basis_spatial)
    spat_to_idx = {lbl: i for i, lbl in enumerate(op_sys.spat_labels)}

    for label, full_mat in zip(
        op_sys.multiplier_labels, op_sys.multiplier_matrices
    ):
        N, L, M_lab, p = label
        if p != 0:
            continue  # at N_t=1 only p=0 is non-degenerate
        sidx = spat_to_idx[(N, L, M_lab)]
        a_s = op_sys._spat_matrices[sidx]
        joint_norm = _opnorm(_commutator(D_L, full_mat))
        spatial_norm = _opnorm(_commutator(D_GV, a_s))
        assert abs(joint_norm - spatial_norm) < 1e-12, (
            f"At N_t=1, joint != spatial at {label}: {joint_norm} vs {spatial_norm}"
        )
