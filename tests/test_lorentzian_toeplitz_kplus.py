"""Frozen falsifier: the Lorentzian-propinquity rescue probe (Toeplitz under K+).

Sprint 2026-06-18 (PI-directed, during /qa group1 Bite B-3). Sharpens Paper 45's
K+ annihilation theorem and confirms the WH7 "weakens-to-convention" lean with an
exact space/time sector decomposition.

Two facts, both bit-exact on the real P45 Krein machinery:

  v1 (J = gamma0_spatial (x) I_t, the Paper-45 involution):
    - SPATIAL multipliers: K+-compressed Lipschitz seminorm = 0 (P45 annihilation
      reproduced -- it is spatial-multiplier-only).
    - TOEPLITZ TEMPORAL multipliers a = I_spatial (x) S_q (the WH7 algebra):
      K+-compressed seminorm = omega_q = 2*pi*q/T EXACTLY (survives, nonzero).
      => the annihilation is NOT system-wide; a nondegenerate seminorm exists.

  v2 (J = gamma0_spatial (x) sign(D_t), a genuine temporal Wick involution):
    - same Toeplitz temporal seminorm = omega_q EXACTLY.
      => the surviving metric carries NO indefinite-signature content even when
      the Krein involution genuinely acts on time. The Lorentzian signature is
      convention relative to the propinquity Lipschitz seminorm.

If either invariant breaks (spatial seminorm survives compression; or the Toeplitz
temporal seminorm departs from omega_q under v1 or v2), the recorded conclusion in
Paper 45 / CLAUDE.md WH7 must be revisited.

Memo: debug/sprint_lorentzian_toeplitz_kplus_probe_memo.md.
"""

from __future__ import annotations

import numpy as np
import pytest

from geovac.krein_space_compact_temporal import CompactTemporalKreinSpace
from geovac.lorentzian_dirac_compact import (
    lorentzian_dirac_compact_matrix,
    fourier_momentum_grid,
    fourier_d_dt_matrix,
)
from geovac.operator_system_compact_temporal import (
    CompactTemporalTruncatedOperatorSystem,
)

TWO_PI = 2.0 * np.pi


def _opn(X: np.ndarray) -> float:
    return float(np.linalg.norm(X, 2)) if X.size else 0.0


def _temporal_shift(N_t: int, q: int) -> np.ndarray:
    """S_q = P_K M_{e_q} P_K in the fourier_momentum_grid basis (mode k -> k+q)."""
    grid = list(fourier_momentum_grid(N_t))
    pos = {k: i for i, k in enumerate(grid)}
    S = np.zeros((N_t, N_t), dtype=complex)
    for i, k in enumerate(grid):
        if (k + q) in pos:
            S[pos[k + q], i] = 1.0
    return S


CELLS = [(2, 3), (3, 5)]


@pytest.mark.parametrize("n_max,N_t", CELLS)
def test_kplus_annihilates_spatial_multipliers(n_max: int, N_t: int) -> None:
    """P45 control: K+ compression kills the SPATIAL multiplier seminorm."""
    T = TWO_PI
    krein = CompactTemporalKreinSpace(n_max=n_max, N_t=N_t, T=T)
    D_L = lorentzian_dirac_compact_matrix(krein)
    I = np.eye(krein.dim, dtype=complex)
    Pp = 0.5 * (I + krein.J)
    D_restr = Pp @ D_L @ Pp
    op = CompactTemporalTruncatedOperatorSystem(n_max=n_max, N_t=N_t, T=T)
    max_restr = 0.0
    max_full = 0.0
    for a in op.multiplier_matrices:
        ar = Pp @ a @ Pp
        max_full = max(max_full, _opn(D_L @ a - a @ D_L))
        max_restr = max(max_restr, _opn(D_restr @ ar - ar @ D_restr))
    assert max_full > 1e-6, "spatial multipliers must be visible BEFORE compression"
    assert max_restr < 1e-9, (
        f"K+ should annihilate spatial seminorm; got {max_restr:.2e}"
    )


@pytest.mark.parametrize("n_max,N_t", CELLS)
def test_kplus_spares_toeplitz_temporal_v1(n_max: int, N_t: int) -> None:
    """v1 rescue: Toeplitz temporal seminorm survives K+ = omega_q exactly."""
    T = TWO_PI
    krein = CompactTemporalKreinSpace(n_max=n_max, N_t=N_t, T=T)
    D_L = lorentzian_dirac_compact_matrix(krein)
    I = np.eye(krein.dim, dtype=complex)
    Pp = 0.5 * (I + krein.J)
    D_restr = Pp @ D_L @ Pp
    dim_s = krein.dim // N_t
    Is = np.eye(dim_s, dtype=complex)
    for q in range(1, (N_t - 1) // 2 + 1):
        a = np.kron(Is, _temporal_shift(N_t, q))
        ar = Pp @ a @ Pp
        s_restr = _opn(D_restr @ ar - ar @ D_restr)
        omega_q = TWO_PI * q / T
        assert s_restr == pytest.approx(omega_q, abs=1e-9), (
            f"v1 Toeplitz seminorm should survive = omega_q={omega_q}; got {s_restr}"
        )


@pytest.mark.parametrize("n_max", [2, 3, 4, 5])
def test_truncated_modular_flow_is_compact(n_max: int) -> None:
    """STRUCTURAL CLOSURE (2026-06-19): the truncated Bisognano-Wichmann boost
    is COMPACT — integer spectrum, exp(2*pi*i*K) = I — at every finite cutoff.

    K_boost = diag(two_m_j) (odd integers; the modular generator of
    geovac/modular_hamiltonian.py _build_geometric_K_boost). Finite truncation
    => finite Hermitian K => discrete (integer) spectrum => modular flow
    e^{itK} is 2*pi-periodic = a circle = the KMS beta=2*pi loop. A compact
    circle carries only a Euclidean (forward-triangle) metric, which is WHY
    the chain inequality is the only universal direction (Paper 45 Q1) and why
    every Lipschitz/boost seminorm is signature-blind. The Lorentzian
    (reverse-triangle) signature is the Wick rotation t->it of this circle,
    needing the non-compact continuum boost (continuous spectrum), unreachable
    at any finite cutoff. See debug/sprint_lorentzian_toeplitz_kplus_probe_memo.md.
    """
    from geovac.full_dirac_operator_system import full_dirac_basis
    basis = full_dirac_basis(n_max)
    two_mj = np.array([b.two_m_j for b in basis], dtype=float)
    # integer (in fact odd-integer) spectrum
    assert np.allclose(two_mj, np.round(two_mj))
    # compact: e^{2 pi i K} = I (period-2pi closure) bit-exact
    U = np.exp(1j * 2.0 * np.pi * two_mj)
    assert np.allclose(U, 1.0 + 0.0j, atol=1e-12), (
        f"modular flow must close at 2*pi (compact) at n_max={n_max}"
    )


def test_continuum_decompactification_trend() -> None:
    """The de-compactification (toward the non-compact Lorentzian boost) is
    STRICTLY an n->inf limit: under the continuum range-normalization
    K/max|spec|, the spectral spacing -> 0 and the flow period -> infinity
    monotonically, but the flow still CLOSES (compact) at every finite cutoff.
    """
    from geovac.full_dirac_operator_system import full_dirac_basis
    spacings, periods = [], []
    for n_max in [2, 3, 4, 5]:
        spec = np.array(sorted({b.two_m_j for b in full_dirac_basis(n_max)}), dtype=float)
        M = np.abs(spec).max()
        g = int(np.gcd.reduce([int(abs(x)) for x in spec if x != 0]))
        spacings.append(float(np.min(np.diff(spec))) / M)   # normalized spacing
        periods.append(2.0 * np.pi * M / g)                 # normalized period
        # still closes (compact) at finite cutoff even after rescaling
        assert np.allclose(np.exp(1j * periods[-1] * spec / M), 1.0, atol=1e-9)
    # spacing strictly decreasing toward 0; period strictly increasing toward inf
    assert all(spacings[i] > spacings[i + 1] for i in range(len(spacings) - 1))
    assert all(periods[i] < periods[i + 1] for i in range(len(periods) - 1))


@pytest.mark.parametrize("n_max,N_t", CELLS)
def test_kplus_temporal_signature_blind_v2(n_max: int, N_t: int) -> None:
    """v2: genuine temporal Wick involution J_t=sign(D_t) leaves seminorm = omega_q.

    The indefinite signature contributes NOTHING to the surviving Lipschitz
    seminorm even when the Krein involution genuinely acts on time.
    """
    T = TWO_PI
    krein = CompactTemporalKreinSpace(n_max=n_max, N_t=N_t, T=T)
    D_L = lorentzian_dirac_compact_matrix(krein)
    Js = krein.J_spatial
    dim_s = Js.shape[0]
    Is = np.eye(dim_s, dtype=complex)
    grid = fourier_momentum_grid(N_t)
    # v2 temporal Krein involution: J_t = sign(omega_k), sign(0)=+1
    Jt = np.diag([1.0 if k >= 0 else -1.0 for k in grid]).astype(complex)
    Jv2 = np.kron(Js, Jt)
    I = np.eye(krein.dim, dtype=complex)
    # legitimacy: genuine involution, preserves spatial annihilation, acts on time
    assert _opn(Jv2 @ Jv2 - I) < 1e-12
    D_t = fourier_d_dt_matrix(N_t, T)
    D_GV = D_L / 1j - np.kron(Js, D_t)  # = D_GV (x) I_t
    assert _opn(Jv2 @ D_GV + D_GV @ Jv2) < 1e-9, "spatial annihilation must persist"
    assert _opn(Jv2 @ D_L - D_L @ Jv2) > 1e-6, "v2 J must genuinely act (not commute)"
    Pp = 0.5 * (I + Jv2)
    D_restr = Pp @ D_L @ Pp
    for q in range(1, (N_t - 1) // 2 + 1):
        a = np.kron(Is, _temporal_shift(N_t, q))
        ar = Pp @ a @ Pp
        s_restr = _opn(D_restr @ ar - ar @ D_restr)
        omega_q = TWO_PI * q / T
        assert s_restr == pytest.approx(omega_q, abs=1e-9), (
            f"v2 seminorm should stay Euclidean = omega_q={omega_q}; got {s_restr} "
            "(signature would be visible if this departed from omega_q)"
        )
