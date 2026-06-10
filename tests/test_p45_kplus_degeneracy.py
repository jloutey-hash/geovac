"""Frozen falsifier for the Paper 45 v2 annihilation theorem (2026-06-09).

Asserts the structural facts established in the P45 hardening sprint
(debug/sprint_p45_hardening_phase1_memo.md, Theorem thm:kplus_annihilation):

  T1. {J, D_GV (x) I} = 0 exactly (forced by Krein-self-adjointness of
      D_L = i*(J_s (x) D_t + D_GV (x) I)).
  T2. P+ (i D_GV (x) I) P+ = 0 exactly (compression annihilates the
      spatial Dirac) while the uncompressed norm is > 1.
  T3. The K+-restricted Lipschitz seminorm ||[P+ D_L P+, P+ a P+]||
      vanishes on EVERY operator-system multiplier (kernel = 100%).
  T4. The full-space seminorm ||[D_L, a]|| is nonzero for some a
      (the unrestricted system does see the Dirac), but its kernel
      exceeds the scalars (kernel-condition failure pre-restriction).
  T5. Spatial kernel condition at N_t = 1: truthful chirality-diagonal
      CH Dirac FAILS (kernel 10/14 at n_max=2); the engineered offdiag
      Dirac PASSES (kernel = 1/14, identity only).

If any of these ever flips, the v2 erratum analysis is wrong and
Paper 45 v2 / the option-B repair design must be revisited.
"""

from __future__ import annotations

import numpy as np
import pytest

from geovac.krein_space_compact_temporal import CompactTemporalKreinSpace
from geovac.lorentzian_dirac_compact import (
    fourier_d_dt_matrix,
    lorentzian_dirac_compact_matrix,
)
from geovac.operator_system_compact_temporal import (
    CompactTemporalTruncatedOperatorSystem,
)
from geovac.full_dirac_operator_system import (
    FullDiracTruncatedOperatorSystem,
    camporesi_higuchi_full_dirac_matrix,
    camporesi_higuchi_offdiag_dirac_matrix,
)

TOL = 1e-12


def _build_cell(n_max: int, N_t: int):
    op_sys = CompactTemporalTruncatedOperatorSystem(n_max=n_max, N_t=N_t)
    krein = op_sys.krein
    D_L = lorentzian_dirac_compact_matrix(krein)
    D_t = fourier_d_dt_matrix(N_t, krein.T)
    time_part = np.kron(krein.J_spatial, D_t)
    space_part = D_L / 1j - time_part  # = D_GV (x) I_t
    I_K = np.eye(krein.dim, dtype=np.complex128)
    P_plus = 0.5 * (I_K + krein.J)
    return op_sys, krein, D_L, space_part, P_plus


@pytest.mark.parametrize("n_max,N_t", [(2, 3)])
def test_anticommutator_and_compression_annihilation(n_max, N_t):
    _, krein, _, space_part, P_plus = _build_cell(n_max, N_t)
    anti = krein.J @ space_part + space_part @ krein.J
    assert np.linalg.norm(anti) < TOL  # T1
    compressed = P_plus @ (1j * space_part) @ P_plus
    assert np.linalg.norm(compressed) < TOL  # T2 (zero)
    assert np.linalg.norm(1j * space_part) > 1.0  # T2 (nontrivial input)


@pytest.mark.parametrize("n_max,N_t", [(2, 3)])
def test_restricted_seminorm_vanishes_identically(n_max, N_t):
    op_sys, _, D_L, _, P_plus = _build_cell(n_max, N_t)
    D_restr = P_plus @ D_L @ P_plus
    restr = []
    full = []
    for a in op_sys.multiplier_matrices:
        a_r = P_plus @ a @ P_plus
        restr.append(np.linalg.norm(D_restr @ a_r - a_r @ D_restr, 2))
        full.append(np.linalg.norm(D_L @ a - a @ D_L, 2))
    restr = np.array(restr)
    full = np.array(full)
    # T3: restricted seminorm identically zero on the whole system.
    assert restr.max() < TOL
    # T4: unrestricted seminorm nonzero somewhere, kernel > 1 element.
    assert full.max() > 1e-6
    n_kernel_full = int(np.sum(full < TOL))
    assert n_kernel_full > 1  # kernel strictly exceeds the scalars
    assert n_kernel_full < len(full)  # but is not everything


def test_spatial_kernel_condition_truthful_vs_offdiag():
    osys = FullDiracTruncatedOperatorSystem(2)
    D_truth = camporesi_higuchi_full_dirac_matrix(osys.basis)
    D_off = camporesi_higuchi_offdiag_dirac_matrix(osys.basis)

    def kernel_count(D):
        vals = [
            np.linalg.norm(D @ a - a @ D, 2)
            for a in osys.multiplier_matrices
        ]
        return sum(v < TOL for v in vals), len(vals)

    k_truth, n = kernel_count(D_truth)
    k_off, _ = kernel_count(D_off)
    # T5: truthful CH fails the kernel condition (10/14 at n_max=2)...
    assert k_truth == 10 and n == 14
    # ...the engineered offdiag Dirac passes (identity only).
    assert k_off == 1
