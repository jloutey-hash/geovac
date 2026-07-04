"""Sturmian Bethe-logarithm machinery for the Paper 36 Lamb-chain tests.

Faithful port of the shared machinery of the archived sprint drivers
  debug/archive/qed_arc/ls3_bethe_log_regularized.py       (LS-3)
  debug/archive/precision_arc/ls4_bethe_log_drake.py       (LS-4)
into the durable tests/ tree (see README.md for the provenance policy).

Two Bethe-log forms are exposed:

* ``bethe_log_acceleration(n, l, N, lam)`` -- LS-3's acceleration form
  (matrix elements of nabla V = Z r-hat / r^2), which converges for
  s-states (ln k0(2S) = 2.786 at N=40, -0.92% vs Drake--Swainson).

* ``bethe_log_drake(n, l, N, lam)`` -- LS-4's velocity form with the
  Drake--Swainson structural normalization D_drake = 2(2l+1) Z^4 / n^3,
  which closes the l>0 structural floor (ln k0(2P) = -0.0307 and
  ln k0(3D) = -0.005236 at N=40).

Basis: chi_{k,l}(r) = r^{l+1} e^{-lam r} L_k^{(2l+1)}(2 lam r), k=0..N-1,
with all matrix elements exact polynomial-times-exponential integrals
evaluated in mpmath at high precision, then cast to float for the
generalized eigenproblem.

Reference values (Drake & Swainson, Phys. Rev. A 41, 1243 (1990)).
"""

from __future__ import annotations

import math
from typing import Dict, List, Tuple

import mpmath as mp
import numpy as np
import scipy.linalg as sla

__all__ = [
    "DRAKE_REF",
    "sturmian_basis_matrices",
    "sturmian_dipole_r_matrix",
    "sturmian_acceleration_matrix",
    "diagonalize_h0_sturmian",
    "build_velocity_spectrum",
    "bethe_log_drake",
    "bethe_log_acceleration",
]

# Drake & Swainson 1990 tabulated ln k_0(n, l) reference values
DRAKE_REF: Dict[Tuple[int, int], float] = {
    (1, 0): 2.9841285558,
    (2, 0): 2.8117698931,
    (2, 1): -0.0300167089,
    (3, 0): 2.7676636450,
    (3, 1): -0.0381902694,
    (3, 2): -0.0052321481,
    (4, 0): 2.7498118405,
}


# ---------------------------------------------------------------------------
# Sturmian basis primitives (LS-2/LS-3/LS-4 shared machinery, ported verbatim)
# ---------------------------------------------------------------------------

def _laguerre_poly_coeffs_mp(n: int, alpha: int) -> List[mp.mpf]:
    """Coefficients of L_n^{(alpha)}(x) as an mpmath list."""
    coeffs = []
    for k in range(n + 1):
        c = math.comb(n + alpha, n - k)
        coeffs.append(mp.mpf((-1) ** k * c) / mp.mpf(math.factorial(k)))
    return coeffs


def _multiply_polys_mp(p: List[mp.mpf], q: List[mp.mpf]) -> List[mp.mpf]:
    out = [mp.mpf(0)] * (len(p) + len(q) - 1)
    for i, pi in enumerate(p):
        for j, qj in enumerate(q):
            out[i + j] += pi * qj
    return out


def _sturmian_radial_poly_mp(k: int, l: int, lam: float) -> List[mp.mpf]:
    """Polynomial coefficients of chi_{k,l}(r) / exp(-lam r):
    chi_{k,l}(r) = r^{l+1} L_k^{(2l+1)}(2 lam r) exp(-lam r)."""
    lam_mp = mp.mpf(lam)
    L_coeffs = _laguerre_poly_coeffs_mp(k, 2 * l + 1)
    L_in_r = [L_coeffs[m] * (mp.mpf(2) * lam_mp) ** m for m in range(k + 1)]
    return [mp.mpf(0)] * (l + 1) + L_in_r


def _laguerre_integral_polyexp(coeffs: List[mp.mpf], alpha) -> mp.mpf:
    """int_0^inf (sum_k c_k r^k) exp(-alpha r) dr = sum_k c_k k!/alpha^{k+1}."""
    val = mp.mpf(0)
    a_pow = mp.mpf(1) / mp.mpf(alpha)
    for k, c in enumerate(coeffs):
        val += mp.mpf(c) * mp.mpf(math.factorial(k)) * a_pow
        a_pow /= mp.mpf(alpha)
    return val


def sturmian_basis_matrices(N: int, l: int, lam: float, Z: int = 1,
                            mp_dps: int = 50) -> dict:
    """Build T, V, S in the Sturmian basis chi_{k,l}, k = 0..N-1."""
    mp.mp.dps = max(mp_dps, 30)
    lam_mp = mp.mpf(lam)
    decay = 2 * lam_mp

    polys = [_sturmian_radial_poly_mp(k, l, lam) for k in range(N)]
    poly_derivs = []
    for poly in polys:
        d_poly = [mp.mpf(0)] * len(poly)
        for m in range(1, len(poly)):
            d_poly[m - 1] = poly[m] * m
        chi_deriv = [d_poly[m] - lam_mp * poly[m] for m in range(len(poly))]
        poly_derivs.append(chi_deriv)

    S_mp = [[mp.mpf(0)] * N for _ in range(N)]
    T_mp = [[mp.mpf(0)] * N for _ in range(N)]
    V_mp = [[mp.mpf(0)] * N for _ in range(N)]

    for i in range(N):
        for j in range(i, N):
            P_prod = _multiply_polys_mp(polys[i], polys[j])
            S_ij = _laguerre_integral_polyexp(P_prod, decay)

            D_prod = _multiply_polys_mp(poly_derivs[i], poly_derivs[j])
            T1 = mp.mpf("0.5") * _laguerre_integral_polyexp(D_prod, decay)
            if len(P_prod) > 2:
                P_div_r2 = P_prod[2:]
                T2_int = mp.mpf(0)
                for k, c in enumerate(P_div_r2):
                    T2_int += c * mp.mpf(math.factorial(k)) / decay ** (k + 1)
            else:
                T2_int = mp.mpf(0)
            T2 = mp.mpf(l * (l + 1)) / 2 * T2_int
            T_ij = T1 + T2

            if len(P_prod) > 1:
                P_div_r = P_prod[1:]
                V_int = mp.mpf(0)
                for k, c in enumerate(P_div_r):
                    V_int += c * mp.mpf(math.factorial(k)) / decay ** (k + 1)
            else:
                V_int = mp.mpf(0)
            V_ij = -mp.mpf(Z) * V_int

            S_mp[i][j] = S_ij
            S_mp[j][i] = S_ij
            T_mp[i][j] = T_ij
            T_mp[j][i] = T_ij
            V_mp[i][j] = V_ij
            V_mp[j][i] = V_ij

    S = np.array([[float(S_mp[i][j]) for j in range(N)] for i in range(N)])
    T = np.array([[float(T_mp[i][j]) for j in range(N)] for i in range(N)])
    V = np.array([[float(V_mp[i][j]) for j in range(N)] for i in range(N)])

    return {"T": T, "V": V, "S": S, "basis_l": l, "lam": lam, "N": N}


def sturmian_dipole_r_matrix(N1: int, l1: int, lam1: float,
                             N2: int, l2: int, lam2: float,
                             mp_dps: int = 50) -> np.ndarray:
    """Dipole matrix R_{ij} = <chi_i^{l1,lam1} | r | chi_j^{l2,lam2}>."""
    mp.mp.dps = max(mp_dps, 30)
    decay = mp.mpf(lam1) + mp.mpf(lam2)
    polys1 = [_sturmian_radial_poly_mp(k, l1, lam1) for k in range(N1)]
    polys2 = [_sturmian_radial_poly_mp(k, l2, lam2) for k in range(N2)]
    R = np.zeros((N1, N2))
    for i in range(N1):
        for j in range(N2):
            P_prod = _multiply_polys_mp(polys1[i], polys2[j])
            val = mp.mpf(0)
            for k, c in enumerate(P_prod):
                val += c * mp.mpf(math.factorial(k + 1)) / decay ** (k + 2)
            R[i, j] = float(val)
    return R


def sturmian_acceleration_matrix(N1: int, l1: int, lam1: float,
                                 N2: int, l2: int, lam2: float,
                                 Z: int = 1, mp_dps: int = 50) -> np.ndarray:
    """Acceleration matrix W_{ij} = <chi_i^{l1,lam1} | 1/r^2 | chi_j^{l2,lam2}>."""
    mp.mp.dps = max(mp_dps, 30)
    decay = mp.mpf(lam1) + mp.mpf(lam2)
    polys1 = [_sturmian_radial_poly_mp(k, l1, lam1) for k in range(N1)]
    polys2 = [_sturmian_radial_poly_mp(k, l2, lam2) for k in range(N2)]
    W = np.zeros((N1, N2))
    for i in range(N1):
        for j in range(N2):
            P_prod = _multiply_polys_mp(polys1[i], polys2[j])
            if len(P_prod) > 2:
                P_div_r2 = P_prod[2:]
                val = mp.mpf(0)
                for k, c in enumerate(P_div_r2):
                    val += c * mp.mpf(math.factorial(k)) / decay ** (k + 1)
            else:
                val = mp.mpf(0)
            W[i, j] = float(val)
    return W


def diagonalize_h0_sturmian(N: int, l: int, lam: float, Z: int = 1,
                            sv_cutoff: float = 1e-10):
    """Diagonalize H_0 = T - Z/r in the Sturmian basis (canonical
    orthogonalization with singular-value cutoff)."""
    mats = sturmian_basis_matrices(N, l, lam, Z)
    H = mats["T"] + mats["V"]
    S = mats["S"]
    s_eigs, U_S = sla.eigh(S)
    keep = s_eigs > sv_cutoff * s_eigs.max()
    s_kept = s_eigs[keep]
    U_kept = U_S[:, keep]
    X = U_kept * (1.0 / np.sqrt(s_kept))
    H_orth = X.T @ H @ X
    H_orth = 0.5 * (H_orth + H_orth.T)
    eigvals, V_orth = sla.eigh(H_orth)
    eigvecs = X @ V_orth
    return eigvals, eigvecs, mats


# ---------------------------------------------------------------------------
# LS-4: velocity-form spectrum + Drake--Swainson structural normalization
# ---------------------------------------------------------------------------

def build_velocity_spectrum(n_target: int, l_target: int, N: int,
                            lam: float, Z: int = 1) -> dict:
    """Velocity-form pseudostate spectrum {(Delta_m, |<nl|p|m>|^2_msum)}.

    Velocity form: |<nl|p|m>|^2_msum = (E_m - E_n)^2 * l_> * |<r>|^2.
    """
    Ry = 0.5
    E_target = -0.5 * Z ** 2 / n_target ** 2
    eigvals_t, eigvecs_t, _ = diagonalize_h0_sturmian(N, l_target, lam, Z)
    idx = int(np.argmin(np.abs(eigvals_t - E_target)))
    E_t = float(eigvals_t[idx])
    psi = eigvecs_t[:, idx]

    channels = []
    if l_target >= 1:
        channels.append(l_target - 1)
    channels.append(l_target + 1)

    DEs_list: List[float] = []
    psqs_list: List[float] = []
    for lp in channels:
        eigvals_p, eigvecs_p, _ = diagonalize_h0_sturmian(N, lp, lam, Z)
        R_basis = sturmian_dipole_r_matrix(N, l_target, lam, N, lp, lam)
        R_pseudo = psi @ R_basis @ eigvecs_p
        l_gt = max(l_target, lp)
        for alpha in range(len(eigvals_p)):
            dE = eigvals_p[alpha] - E_t
            if abs(dE) < 1e-12:
                continue
            R_sq = float(R_pseudo[alpha] ** 2)
            p_sq = (dE ** 2) * l_gt * R_sq
            DEs_list.append(float(dE))
            psqs_list.append(p_sq)

    DEs = np.array(DEs_list)
    psqs = np.array(psqs_list)

    I_v = float((psqs * DEs).sum())
    log_args = np.log(np.abs(2.0 * DEs / Ry))
    J_v = float((psqs * DEs * log_args).sum())

    return {
        "n_target": n_target, "l_target": l_target, "N": N, "lam": lam,
        "E_target": E_t, "E_target_exact": E_target,
        "DEs": DEs, "psqs": psqs,
        "I_v": I_v, "J_v": J_v,
        "n_states": int(len(DEs)),
    }


def bethe_log_drake(n_target: int, l_target: int, N: int,
                    lam: float = None, Z: int = 1) -> dict:
    """LS-4 Drake--Swainson regularized Bethe log:
    ln k_0(n, l) = J_v(nl) / D_drake(n, l),
    D_drake(n, l) = 2 (2l+1) Z^4 / n^3.

    For l > 0 the spectral sum is cutoff-independent (I_v = 0 by closure),
    so the direct J_v assembly IS the regularized value; no K sweep needed
    (verified K-independent over 3.6 orders of magnitude in LS-4).
    """
    if lam is None:
        lam = Z / n_target
    spec = build_velocity_spectrum(n_target, l_target, N, lam, Z)
    D_drake_nl = 2.0 * (2 * l_target + 1) * Z ** 4 / n_target ** 3
    ln_k0 = spec["J_v"] / D_drake_nl
    ln_k0_natural = (spec["J_v"] / spec["I_v"]
                     if abs(spec["I_v"]) > 1e-12 else float("nan"))
    return {
        "n_target": n_target, "l_target": l_target, "N": N, "lam": lam,
        "I_v": spec["I_v"], "J_v": spec["J_v"],
        "D_drake_nl": D_drake_nl,
        "ln_k0_drake": ln_k0,
        "ln_k0_natural_J_over_I_v": ln_k0_natural,
    }


# ---------------------------------------------------------------------------
# LS-3: acceleration-form Bethe log (converges for s-states)
# ---------------------------------------------------------------------------

def bethe_log_acceleration(n_target: int, l_target: int, N: int = 20,
                           lam: float = None, Z: int = 1) -> dict:
    """LS-3 acceleration-form Bethe log:
    I_a = sum_m |<nl|nabla V|m>|^2 / (E_m - E_n)
    J_a = sum_m |<nl|nabla V|m>|^2 ln|2(E_m - E_n)/Ry| / (E_m - E_n)
    ln k_0 = J_a / I_a
    with |<n,l|nabla V|m>|^2_msum_radial = l_> Z^2 |W_radial|^2.
    """
    if lam is None:
        lam = Z / n_target
    Ry = 0.5
    E_target = -0.5 * Z ** 2 / n_target ** 2

    eigvals_t, eigvecs_t, _ = diagonalize_h0_sturmian(N, l_target, lam, Z)
    idx_target = int(np.argmin(np.abs(eigvals_t - E_target)))
    E_t_basis = eigvals_t[idx_target]
    psi_target = eigvecs_t[:, idx_target]

    channels = []
    if l_target >= 1:
        channels.append(l_target - 1)
    channels.append(l_target + 1)

    I_total = 0.0
    J_total = 0.0
    n_used = 0

    for lp in channels:
        eigvals_p, eigvecs_p, _ = diagonalize_h0_sturmian(N, lp, lam, Z)
        W_basis = sturmian_acceleration_matrix(N, l_target, lam, N, lp, lam, Z)
        W_pseudo = psi_target @ W_basis @ eigvecs_p
        l_gt = max(l_target, lp)
        for alpha in range(len(eigvals_p)):
            E_alpha = eigvals_p[alpha]
            dE = E_alpha - E_t_basis
            if abs(dE) < 1e-12:
                continue
            d2_msum = l_gt * (Z ** 2) * (W_pseudo[alpha] ** 2)
            ln_term = math.log(abs(2.0 * dE / Ry))
            I_total += d2_msum / dE
            J_total += d2_msum * ln_term / dE
            n_used += 1

    ln_k0 = J_total / I_total if abs(I_total) > 1e-30 else float("nan")

    return {
        "n_target": n_target, "l_target": l_target,
        "N_basis": N, "lam": lam,
        "E_target_exact": E_target,
        "E_target_basis": float(E_t_basis),
        "I": float(I_total), "J": float(J_total), "ln_k0": float(ln_k0),
        "n_pseudostates_used": n_used,
        "channels": channels,
        "form": "acceleration",
    }
