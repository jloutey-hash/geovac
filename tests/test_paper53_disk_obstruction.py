"""Paper 53 — finite-disk OBSTRUCTION + prop:prop2 (permanent backing).

Backfill of the two Paper-53 coverage gaps surfaced by the v4.43.0 /qa code
dimension (the disk-obstruction NEGATIVES and the prop:prop2 Proposition had
backing only under debug/archive/, which is pruning-scheduled).  These tests
RECOMPUTE the paper's own constructions from scratch and assert the *current*
(corrected, 2026-06) claims:

  (1) The finite disk-with-cone is OBSTRUCTED by its Dirichlet boundary.  The
      Markov–Cesàro Berezin g_t = (e^{-tΔ}f)/(e^{-tΔ}1) of a *positive* radial
      function f=(ρ/R)² is NOT positivity-preserving (min g ≈ −0.13 < 0) and its
      interior approximation rate does NOT decay (Λ^{+0.07}, p>0).  This is the
      negative result that pivots the paper to the boundaryless plane.
      (Porting debug/archive/spectral_triple_arc/b3_markov_berezin.py.)

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
# (1) Finite-disk Markov–Cesàro Berezin — the OBSTRUCTION (positivity fails,
#     interior rate does not decay).  k=0 (radial) sector: f=(ρ/R)² is radial.
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


def _smooth_ratio(f_fn, modes, t, rho):
    """Markov–Cesàro smoothed ratio g_t = (e^{-tΔ}f)/(e^{-tΔ}1) on the disk."""
    xg, wg = leggauss(600)
    s = 0.5 * R * (xg + 1.0)
    ws = 0.5 * R * wg
    num = np.zeros_like(rho)
    den = np.zeros_like(rho)
    for (a, lam, norm) in modes:
        Rj_s = _Rj(s, a, norm)
        cf = np.sum(Rj_s * f_fn(s) * s * ws)   # <R_j, f>
        c1 = np.sum(Rj_s * 1.0 * s * ws)       # <R_j, 1>
        e = np.exp(-t * lam)
        Rj_rho = _Rj(rho, a, norm)
        num += e * cf * Rj_rho
        den += e * c1 * Rj_rho
    return num / den


def test_disk_markov_berezin_positivity_fails():
    """Disk obstruction: f=(ρ/R)² ≥ 0 but g_t goes NEGATIVE (min ≈ −0.13) and the
    interior rate does NOT decay (Λ^{+p}, p>0).  The opposite of the plane, where
    positivity holds (min g ≈ −2e-4) and the rate converges (Λ^{-0.6..-0.9})."""
    f1 = lambda r: (r / R) ** 2
    gradf1 = 2.0 / R
    rho = np.linspace(1e-4, R, 2000)
    Jmax = 80
    Lambdas = [100.0, 200.0, 400.0, 800.0, 1600.0]

    mins, e1s = [], []
    for Lam in Lambdas:
        modes = _k0_modes(Jmax, Lam)
        g1 = _smooth_ratio(f1, modes, 1.0 / Lam, rho)
        mins.append(float(g1.min()))
        e1s.append(float(np.max(np.abs(g1 - f1(rho))) / gradf1))

    # (a) positivity FAILS: a positive f produces a Berezin image dipping well
    #     below zero (the Dirichlet-boundary 0/0 artifact at ρ=R).
    assert min(mins) < -0.05, (
        f"expected disk positivity failure (min g < -0.05); got min over sweep "
        f"= {min(mins):.4f} (per-Lambda mins={[round(m,4) for m in mins]})"
    )
    # specifically the Λ=200 cell reproduces the paper's min g ≈ −0.13
    assert mins[1] < -0.10, f"Lambda=200 min g = {mins[1]:.4f}, expected < -0.10"

    # (b) interior rate does NOT decay: log-log slope p = +slope > 0
    #     (Λ^{+0.07} no-decay), unlike the plane's p < 0 convergence.
    p = float(np.polyfit(np.log(Lambdas), np.log(e1s), 1)[0])
    assert p > 0.0, (
        f"disk interior rate must NOT decay (slope > 0); got Lambda^{{{p:.3f}}} "
        f"(e1={[round(e,4) for e in e1s]})"
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
