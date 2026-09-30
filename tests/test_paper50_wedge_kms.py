"""Genuine backing for Paper 50 Sec. 3 (wedge KMS entropy) and Sec. 4 /
Thm 4.2 (Casini-Huerta-Myers identification).

Two load-bearing headlines that were previously NO-TEST (debug-memo-only):

  Gap 4 (CHM, thm:chm_identification): the framework wedge modular Hamiltonian
    K_alpha = J_polar on the truncated Camporesi-Higuchi triple has spectrum
    {1, 3, 5, ..., 2 n_max - 1} (the odd integers 2 m_j) with multiplicity
    g(2k+1) = (n_max - k)(n_max - k + 1), and generates a 2*pi-periodic
    modular flow.

  Gap 3 (sec:wedge_kms): the von Neumann entropy of the unit-beta Gibbs
    state rho_W = exp(-K_alpha)/Z (beta = 1 in the units where spec K_alpha
    is the odd integers -- NOT the Bisognano-Wichmann normalisation, which is
    beta = pi for K_alpha = 2 J_z; see rem:bw_normalisation and the b-family
    test below) scales as
        S(rho_W) = 2 log n_max + (coth 1 - log(2 sinh 1)) + O(1/n_max),
        coth 1 - log(2 sinh 1) = 0.458448...
    The slope 2 is the growth exponent of the lowest-shell degeneracy
    g(1) = n_max(n_max+1) ~ n_max^2 (the Gibbs weight concentrates on the
    lowest shells, so S ~ log g(1)); it is a counting fact about the
    truncation, not a dimension of a wedge boundary.  The area-law fit
    S ~ n_max^2 is rejected.

Everything is RECOMPUTED from the live modular_hamiltonian machinery
(for_bisognano_wichmann -> restrict_K_alpha_to_wedge), not read from the
debug/bh_phase0 memo.  The full-Dirac construction is O(dim^2) and cannot
reach the n_max needed for the asymptotic constant, so the large-n tests use
the analytical shell-degeneracy formula -- which test_entropy_formula_matches_machinery
validates bit-exactly against the live machinery at n_max = 2..5.
"""
from __future__ import annotations

import math
from collections import Counter

import numpy as np
import pytest

from geovac.modular_hamiltonian import for_bisognano_wichmann


# ---------------------------------------------------------------------------
# live-machinery helpers (n_max <= ~10; full Dirac is O(dim^2))
# ---------------------------------------------------------------------------

def _wedge_K_diag(n_max: int) -> np.ndarray:
    """Diagonal of the wedge-restricted modular Hamiltonian K_alpha^W, from
    the live BW construction on the truncated Camporesi-Higuchi triple."""
    mh = for_bisognano_wichmann(n_max)
    return np.real(np.diag(mh.restrict_K_alpha_to_wedge()))


def _entropy_from_diag(k_diag: np.ndarray, s: float = 1.0) -> float:
    """S(rho) = -Tr[rho log rho], rho = exp(-s K)/Z, from a diagonal K."""
    w = np.exp(-s * k_diag)
    Z = w.sum()
    p = w / Z
    return float(-(p * np.log(p)).sum())


# ---------------------------------------------------------------------------
# analytical shell-degeneracy formula (validated vs machinery below)
# ---------------------------------------------------------------------------

def _degeneracy(k: int, n_max: int) -> int:
    """Multiplicity of eigenvalue two_m_j = k (positive odd) on the wedge."""
    assert k > 0 and k % 2 == 1
    m = (k - 1) // 2
    return (n_max - m) * (n_max - m + 1) if n_max > m else 0


def _analytic_entropy(n_max: int, s: float = 1.0) -> float:
    k_vals = range(1, 2 * n_max, 2)
    degs = np.array([_degeneracy(k, n_max) for k in k_vals], dtype=float)
    w = degs * np.exp(-s * np.array(list(k_vals), dtype=float))
    Z = w.sum()
    # per-state probability p(k) = exp(-s k)/Z; entropy sums over all states
    S = 0.0
    for k, d in zip(k_vals, degs):
        if d == 0:
            continue
        p = math.exp(-s * k) / Z
        if p <= 0.0:          # high-k shells underflow to 0; contribute 0 log 0 = 0
            continue
        S -= d * p * math.log(p)
    return S


CHM_CONSTANT = 1.0 / math.tanh(1.0) - math.log(2.0 * math.sinh(1.0))  # 0.458448...


# ===========================================================================
# Gap 4 -- CHM identification (thm:chm_identification)
# ===========================================================================

@pytest.mark.parametrize("n_max", [2, 3, 4, 5])
def test_chm_spectrum_is_odd_integers(n_max):
    """spec(K_alpha^W) = {1, 3, ..., 2 n_max - 1} (the odd integers 2 m_j)."""
    k = _wedge_K_diag(n_max)
    distinct = sorted({int(round(x)) for x in k})
    assert distinct == list(range(1, 2 * n_max, 2))
    # genuinely integer-valued (canonical (2pi)^-1 units), not rounded noise
    assert np.allclose(k, np.round(k), atol=1e-9)


@pytest.mark.parametrize("n_max", [2, 3, 4, 5])
def test_chm_multiplicities(n_max):
    """Eigenvalue 2k+1 has multiplicity (n_max - k)(n_max - k + 1)."""
    k = _wedge_K_diag(n_max)
    counts = Counter(int(round(x)) for x in k)
    for m in range(n_max):
        assert counts[2 * m + 1] == (n_max - m) * (n_max - m + 1)
    # wedge dimension = sum_{n=1}^{n_max} n(n+1)
    assert len(k) == sum(n * (n + 1) for n in range(1, n_max + 1))


@pytest.mark.parametrize("n_max", [2, 3, 4])
def test_modular_flow_2pi_closure(n_max):
    """sigma_{2pi}(O) = O bit-exactly (integer K spectrum -> exact closure)."""
    mh = for_bisognano_wichmann(n_max)
    dim = mh.dim
    rng = np.random.default_rng(0)
    O = rng.standard_normal((dim, dim)) + 1j * rng.standard_normal((dim, dim))
    O = O + O.conj().T
    ok, residual = mh.verify_modular_periodicity(O)
    assert ok and residual < 1e-10, f"2pi closure residual {residual:.2e}"


# ===========================================================================
# Gap 3 -- wedge KMS entropy (sec:wedge_kms)
# ===========================================================================

@pytest.mark.parametrize("n_max,S_memo", [
    (2, 1.922212), (3, 2.691039), (4, 3.250103), (5, 3.688641),
])
def test_entropy_matches_memo_and_thermo_identity(n_max, S_memo):
    """S(rho_W) reproduces the BH-Phase0 values AND equals the thermodynamic
    closed form log Z + <K_alpha> (the paper's 'reproduces direct
    diagonalization to <= 4e-16')."""
    k = _wedge_K_diag(n_max)
    S_direct = _entropy_from_diag(k)              # -Tr rho log rho
    w = np.exp(-k); Z = w.sum(); p = w / Z
    S_thermo = float(np.log(Z) + (p * k).sum())   # log Z + <K>
    assert abs(S_direct - S_memo) < 1e-5
    assert abs(S_direct - S_thermo) < 1e-13


@pytest.mark.parametrize("n_max", [2, 3, 4, 5])
def test_entropy_formula_matches_machinery(n_max):
    """The analytical shell-degeneracy entropy equals the live-machinery
    entropy bit-exactly -- this VALIDATES the formula used for the large-n
    asymptotic tests below."""
    assert abs(_analytic_entropy(n_max) - _entropy_from_diag(_wedge_K_diag(n_max))) < 1e-12


def test_chm_constant_closed_form():
    """coth 1 - log(2 sinh 1) = 0.45844874... (the eq:S_W_fit additive
    constant; the paper quotes it rounded as 0.458448)."""
    assert abs(CHM_CONSTANT - 0.45844874336819064) < 1e-12
    assert abs(CHM_CONSTANT - 0.458448) < 1e-5


@pytest.mark.slow
def test_asymptotic_constant_converges_to_coth():
    """S(rho_W) - 2 log n_max -> coth 1 - log(2 sinh 1) as n_max -> infinity
    (analytical formula, validated against machinery at small n)."""
    deviations = []
    for n in (50, 200, 1000, 4000):
        deviations.append(_analytic_entropy(n) - 2.0 * math.log(n))
    # monotone approach to the closed-form constant
    errs = [abs(d - CHM_CONSTANT) for d in deviations]
    assert errs[-1] < 5e-3, f"S-2ln n -> {deviations[-1]:.5f}, want {CHM_CONSTANT:.5f}"
    assert errs[0] > errs[-1], "should converge (error decreasing in n_max)"


@pytest.mark.slow
def test_slope_converges_to_2():
    """Windowed least-squares slope of S vs log n_max climbs to ~2 at large
    n_max: the growth exponent of the lowest-shell degeneracy
    g(1) = n_max(n_max+1) ~ n_max^2 on which the unit-beta Gibbs state
    concentrates (a counting fact, not a boundary dimension)."""
    ns = np.arange(350, 2001, 50, dtype=float)
    S = np.array([_analytic_entropy(int(n)) for n in ns])
    A = np.column_stack([np.log(ns), np.ones_like(ns)])
    slope, _ = np.linalg.lstsq(A, S, rcond=None)[0]
    assert abs(slope - 2.0) < 1e-2, f"slope {slope:.5f} not ~2"


def test_area_law_rejected_log_preferred():
    """S ~ n_max^2 (area law) is rejected (R^2 < 0.9) while S ~ a log n + b
    is strong (R^2 > 0.999), over n_max = 2..12 (BH-Phase0 finding)."""
    ns = np.arange(2, 13, dtype=float)
    S = np.array([_analytic_entropy(int(n)) for n in ns])
    SS_tot = np.sum((S - S.mean()) ** 2)

    def _r2(cols):
        A = np.column_stack(cols + [np.ones_like(ns)])
        c = np.linalg.lstsq(A, S, rcond=None)[0]
        return 1.0 - np.sum((S - A @ c) ** 2) / SS_tot

    r2_log = _r2([np.log(ns)])
    r2_area = _r2([ns ** 2])
    assert r2_log > 0.999, f"log fit R^2={r2_log:.5f}"
    assert r2_area < 0.9, f"area-law fit R^2={r2_area:.5f} not rejected"


@pytest.mark.parametrize("n_max", [6])
def test_s_deformation_limits(n_max):
    """rho(s)=exp(-s K)/Z: s->0 gives S->log dim_W (maximally mixed on the
    wedge); s->infty gives S->log n_eq (maximally mixed on the equator shell,
    n_eq = n_max(n_max+1))."""
    k = _wedge_K_diag(n_max)
    dim_W = len(k)
    n_eq = n_max * (n_max + 1)
    S_s0 = _entropy_from_diag(k, s=1e-6)
    S_sinf = _entropy_from_diag(k, s=40.0)
    assert abs(S_s0 - math.log(dim_W)) < 1e-3
    assert abs(S_sinf - math.log(n_eq)) < 1e-6

# ===========================================================================
# rem:bw_normalisation -- K_alpha = 2 J_z; the b-family of Gibbs states
# (guard pass, 2026-09-29, /qa trunk DELTA #16)
# ===========================================================================

def _b_family_constant(b: float) -> float:
    """Closed-form n_max -> infinity constant of S(e^{-bK}/Z) - 2 log n_max."""
    return b / math.tanh(b) - math.log(2.0 * math.sinh(b))


@pytest.mark.parametrize("b", [0.5, 1.0, math.pi])
def test_b_family_entropy_constant(b):
    """S(rho_b) = 2 log n_max + b coth b - log(2 sinh b) + O(1/n_max) for
    rho_b = e^{-bK_alpha}/Z.  Rejects: 'the constant coth 1 - log(2 sinh 1)
    is b-independent / the Bisognano-Wichmann-normalised value' (at the BW
    normalisation b = pi the constant is 0.0136..., not 0.4584...)."""
    # signed residual e(n) = S - 2 log n - C(b) = c/n + O(1/n^2); Richardson
    # (4 e(4000) - e(1000))/3 cancels the 1/n term (measured |.| <= 6.3e-7
    # at b = 0.5, 1, pi), so a +4e-4 error in C(b) is caught at every b.
    e = {n: _analytic_entropy(n, b) - 2.0 * math.log(n) - _b_family_constant(b)
         for n in (1000, 4000)}
    rich = (4.0 * e[4000] - e[1000]) / 3.0
    assert abs(rich) < 5e-6, f"b={b}: Richardson residual {rich:.2e}"
    assert abs(e[4000]) < abs(e[1000]), "residual must shrink with n_max"
    if b == math.pi:
        assert abs(_b_family_constant(b) - 0.0136246) < 1e-6
        assert abs(_b_family_constant(b) - CHM_CONSTANT) > 0.4


@pytest.mark.parametrize("n_max", [2, 3, 4])
def test_adjoint_flow_closes_at_pi(n_max):
    """spec K_alpha = odd integers (K_alpha = 2 J_z): Ad e^{i pi K} = id while
    e^{i pi K} = -1 (the spinor sign), so the adjoint flow's minimal period
    is pi, half the 2pi the paper's closure statement uses.  Rejects:
    'the minimal period of the modular-type flow is 2pi'."""
    k = _wedge_K_diag(n_max)
    ki = np.rint(k).astype(int)
    assert np.allclose(k, ki, atol=1e-10) and np.all(ki % 2 == 1)
    u = np.exp(1j * math.pi * k)
    assert np.allclose(u, -1.0, atol=1e-10)
    phase = np.exp(1j * math.pi * (k[:, None] - k[None, :]))
    assert np.allclose(phase, 1.0, atol=1e-10)
    half = np.exp(1j * (math.pi / 2) * (k[:, None] - k[None, :]))
    assert not np.allclose(half, 1.0, atol=1e-6)
