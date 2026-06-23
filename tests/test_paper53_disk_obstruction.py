"""Paper 53 — finite-disk OBSTRUCTION + prop:prop2 (permanent backing).

Backfill of the two Paper-53 coverage gaps surfaced by the v4.43.0 /qa code
dimension (the disk-obstruction NEGATIVES and the prop:prop2 Proposition had
backing only under debug/archive/, which is pruning-scheduled).  These tests
RECOMPUTE the paper's own constructions from scratch and assert the *current*
(corrected, 2026-06) claims:

  (1) The finite disk-with-cone is OBSTRUCTED by its Dirichlet boundary, but the
      obstruction is the NON-DECAYING interior approximation rate (weight-robust,
      ~Λ^{+0.07}, plateaus ~0.3), NOT positivity.  The Cesàro weight (1-λ/Λ)_+^s
      at s>=2 (Def 4.1) PRESERVES positivity on the disk (min g > 0); only the
      truncated HEAT weight e^{-tλ} fails positivity (Gibbs, min g ≈ −0.13).  The
      non-decaying rate is what pivots the paper to the boundaryless plane (where
      it decays Λ^{-0.6..-0.9}).  (Corrected 2026-06-23: the positivity leg was
      previously mis-attributed to the Cesàro construction; the −0.13 dip is the
      heat weight.  Was test_disk_markov_berezin_positivity_fails.)

  (2) prop:prop2 — the assembled disk operator system O_n = P_n C(disk) P_n has
      propagation number 2: O² already fills the full M_N(C) envelope at every
      tested (K,J).  Computed via the live geovac.operator_system machinery.
      (Porting debug/archive/gravity_arc/l6_disk_prop2.py.)

The plane-convergence POSITIVE results live in test_paper53_plane_bochner_riesz.py.
"""

from __future__ import annotations

import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss
from scipy.special import jn_zeros, jv

from geovac.operator_system import operator_system_dim, _extract_matrix_basis

R = 1.0


# --------------------------------------------------------------------------
# (1) Finite-disk Berezin — the OBSTRUCTION is the NON-DECAYING interior rate
#     (weight-robust), NOT positivity.  Cesàro (Def 4.1) preserves positivity;
#     only the heat weight fails it.  k=0 (radial) sector: f=(ρ/R)² is radial.
# --------------------------------------------------------------------------
def _k0_modes(Jmax: int, Lam: float):
    """k=0 Dirichlet radial modes R_j(ρ)=√2/(R|J_1(a_j)|)·J_0(a_j ρ/R),
    normalized ∫_0^R R_j² ρ dρ = 1; eigenvalue λ_j=(a_j/R)²; keep λ_j ≤ Λ."""
    modes = []
    for a in jn_zeros(0, Jmax):
        lam = (a / R) ** 2
        if lam <= Lam:
            modes.append((a, lam, np.sqrt(2.0) / (R * abs(jv(1, a)))))
    return modes


def _Rj(rho, a, norm):
    return norm * jv(0, a * rho / R)


def _smooth_ratio(f_fn, modes, weight, rho):
    """Markov smoothed ratio g = (Σ w_a<R_a,f>R_a)/(Σ w_a<R_a,1>R_a) on the
    disk, for an arbitrary spectral weight w_a = weight(λ_a)."""
    xg, wg = leggauss(600)
    s = 0.5 * R * (xg + 1.0)
    ws = 0.5 * R * wg
    num = np.zeros_like(rho)
    den = np.zeros_like(rho)
    for (a, lam, norm) in modes:
        Rj_s = _Rj(s, a, norm)
        cf = np.sum(Rj_s * f_fn(s) * s * ws)   # <R_j, f>
        c1 = np.sum(Rj_s * 1.0 * s * ws)       # <R_j, 1>
        w = weight(lam)
        Rj_rho = _Rj(rho, a, norm)
        num += w * cf * Rj_rho
        den += w * c1 * Rj_rho
    return num / den


def _heat_weight(Lam):
    """Truncated heat weight w_a = e^{-t λ_a}, t = 1/Λ (Gibbs-non-positive)."""
    return lambda lam: np.exp(-(1.0 / Lam) * lam)


def _cesaro_weight(Lam, s=2):
    """Cesàro weight w_a = (1 - λ_a/Λ)_+^s of Def 4.1 (positivity-preserving)."""
    return lambda lam: max(0.0, (1.0 - lam / Lam)) ** s


def test_disk_cesaro_positivity_holds_heat_fails():
    """Positivity leg (corrected 2026-06-23): the Cesàro weight (Def 4.1) at
    s>=2 PRESERVES positivity on the finite disk (min g > 0); the truncated
    HEAT weight e^{-tλ} fails it (min g ≈ −0.13, Gibbs).  So positivity is NOT
    the disk obstruction — the −0.13 dip belongs to the heat weight, not the
    Cesàro construction."""
    f1 = lambda r: (r / R) ** 2
    rho = np.linspace(1e-4, R, 2000)
    Jmax = 80
    Lambdas = [100.0, 200.0, 400.0, 800.0, 1600.0]

    cesaro_mins, heat_mins = [], []
    for Lam in Lambdas:
        modes = _k0_modes(Jmax, Lam)
        g_c = _smooth_ratio(f1, modes, _cesaro_weight(Lam, s=2), rho)
        g_h = _smooth_ratio(f1, modes, _heat_weight(Lam), rho)
        cesaro_mins.append(float(g_c.min()))
        heat_mins.append(float(g_h.min()))

    # Cesàro (Def 4.1) PRESERVES positivity on the disk: min g > 0 at every Λ.
    assert min(cesaro_mins) > 0.0, (
        f"Cesàro s=2 must preserve positivity (min g > 0); got per-Lambda mins "
        f"{[round(m,4) for m in cesaro_mins]}"
    )
    # the HEAT weight does NOT: it dips to ≈ −0.13 (the Gibbs failure).
    assert min(heat_mins) < -0.05, (
        f"heat weight should fail positivity (min g < -0.05); got "
        f"{[round(m,4) for m in heat_mins]}"
    )
    # specifically the Λ=200 heat cell reproduces the paper's min g ≈ −0.13.
    assert heat_mins[1] < -0.10, (
        f"heat Lambda=200 min g = {heat_mins[1]:.4f}, expected < -0.10"
    )


def test_disk_interior_rate_does_not_decay():
    """The genuine finite-disk OBSTRUCTION (weight-robust): the
    approximate-identity coefficient does NOT decay toward 0 — it plateaus
    ~0.3 across the Λ sweep, for BOTH the Cesàro and heat weights.  This is the
    Dirichlet-boundary obstruction that pivots the paper to the plane (where the
    rate decays Λ^{-0.6..-0.9})."""
    f1 = lambda r: (r / R) ** 2
    gradf1 = 2.0 / R
    rho = np.linspace(1e-4, R, 2000)
    Jmax = 80
    Lambdas = [100.0, 200.0, 400.0, 800.0, 1600.0]

    for name in ("cesaro", "heat"):
        e1s = []
        for Lam in Lambdas:
            modes = _k0_modes(Jmax, Lam)
            w = _cesaro_weight(Lam, s=2) if name == "cesaro" else _heat_weight(Lam)
            g = _smooth_ratio(f1, modes, w, rho)
            e1s.append(float(np.max(np.abs(g - f1(rho))) / gradf1))
        # no decay: the coefficient stays bounded away from 0 (plateaus ~0.3),
        # i.e. it does NOT approach 0 like the plane's Λ^{-0.6..-0.9}.
        assert min(e1s) > 0.15, (
            f"{name}: interior rate must NOT decay (e1 plateaus, min > 0.15); "
            f"got e1={[round(e,4) for e in e1s]}"
        )
        # and the last value is not an order of magnitude below the first.
        assert e1s[-1] > 0.5 * e1s[0], (
            f"{name}: e1 must not converge; first={e1s[0]:.4f}, last={e1s[-1]:.4f}"
        )


# --------------------------------------------------------------------------
# (2) prop:prop2 — assembled disk operator system O_n = P_n C(disk) P_n has
#     propagation number 2 (O² fills the full M_N(C) envelope).
# --------------------------------------------------------------------------
def _radial_modes(m: float, R_: float, N_rho: int, J: int):
    from scipy.linalg import eigh_tridiagonal
    a = R_ / N_rho
    k = np.arange(1, N_rho + 1)
    rho = k * a
    diag = 2.0 / a**2 + (m**2 - 0.25) / rho**2
    off = -np.ones(N_rho - 1) / a**2
    _, v = eigh_tridiagonal(diag, off, select="i", select_range=(0, J - 1))
    return v, rho, a


def _build_disk_operator_system(K: int, J: int, R_: float, N_rho: int, max_deg: int):
    ks = list(range(-K, K))           # 2K azimuthal modes
    n_modes = len(ks) * J
    mode_vecs, rho, a = {}, None, None
    for k in ks:
        v, rho, a = _radial_modes(abs(k + 0.5), R_, N_rho, J)
        mode_vecs[k] = v

    def idx(k, j):
        return (k + K) * J + j

    def mult_matrix(p: int, q: int):
        M = np.zeros((n_modes, n_modes), dtype=complex)
        for k in ks:
            kp = k + q
            if kp not in ks:
                continue
            wgt = (rho ** p) * a
            Ovl = (mode_vecs[kp] * wgt[:, None]).T @ mode_vecs[k]
            for jp in range(J):
                for j in range(J):
                    M[idx(kp, jp), idx(k, j)] = Ovl[jp, j]
        return M

    gens = [np.eye(n_modes, dtype=complex)]
    seen = set()
    for a_ in range(max_deg + 1):
        for b_ in range(max_deg + 1):
            if a_ + b_ == 0 or a_ + b_ > max_deg:
                continue
            q, p = a_ - b_, a_ + b_
            if abs(q) >= 2 * K or (p, q) in seen:
                continue
            seen.add((p, q))
            gens.append(mult_matrix(p, q))
    return gens, n_modes


def _prop_from_generators(gens, max_k=4, tol=1e-9):
    N = gens[0].shape[0]
    target = N * N
    cur = list(gens)
    dims = [operator_system_dim(cur, tol=tol)]
    if dims[0] == target:
        return 1, dims
    for k in range(2, max_k + 1):
        basis = _extract_matrix_basis(cur, tol=tol)
        newg = [A @ B for A in basis for B in gens]
        dk = operator_system_dim(newg, tol=tol)
        dims.append(dk)
        if dk == target:
            return k, dims
        if dk == dims[-2]:
            return -1, dims
        cur = newg
    return -1, dims


@pytest.mark.parametrize("K,J,N", [(2, 2, 8), (2, 3, 12), (3, 2, 12), (3, 3, 18)])
def test_disk_operator_system_prop2(K, J, N):
    """prop:prop2 — O² already spans the full M_N(C) envelope (propagation 2)."""
    gens, n_modes = _build_disk_operator_system(K, J, R_=8.0, N_rho=120, max_deg=8)
    assert n_modes == N, f"(K={K},J={J}) expected N={N}, got {n_modes}"
    prop, dims = _prop_from_generators(gens)
    assert prop == 2, f"(K={K},J={J}) expected prop=2; got {prop} (dims={dims})"
    # the 2-step generating signature: O² == full M_N(C) envelope (N²)
    assert dims[1] == N * N, (
        f"(K={K},J={J}) O² should fill M_{N}(C) envelope {N*N}; got {dims[1]}"
    )
