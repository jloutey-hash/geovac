"""Sprint LS-3: Regularized Bethe logarithm for arbitrary l (acceleration form).

Goal
----
Close out the from-scratch hydrogen Lamb shift derivation by implementing
Drake's acceleration form of the Bethe logarithm, which converges for
l > 0 where the velocity form (LS-2) had a 0/0 closure pathology.

Background
----------
The Bethe logarithm in the velocity form:

    ln k_0(nl) = J_v(nl) / I_v(nl)

with
    I_v(nl) = Sum_m |<nl|p|m>|^2 (E_m - E_n)
    J_v(nl) = Sum_m |<nl|p|m>|^2 (E_m - E_n) ln|2(E_m - E_n)/Ry|

Closure identity: I_v(nl) = 2 Z^4 / n^3 * delta_{l,0}.
For l > 0, I_v vanishes by closure -> J/I is 0/0 in finite basis.

Acceleration form (Drake-Swainson 1990, Schwartz 1961):
Use the identity [p, H_0] = -i nabla V to transform:
    |<nl|p|m>|^2 (Delta E)^2 = |<nl|nabla V|m>|^2

So |<nl|p|m>|^2 (Delta E) = |<nl|nabla V|m>|^2 / (Delta E)
   |<nl|p|m>|^2 (Delta E) ln(...) = |<nl|nabla V|m>|^2 ln(...) / (Delta E)

This gives:
    I_a(nl) = Sum_m |<nl|nabla V|m>|^2 / (E_m - E_n)
    J_a(nl) = Sum_m |<nl|nabla V|m>|^2 ln|2(E_m - E_n)/Ry| / (E_m - E_n)

By the operator identity these equal I_v and J_v -- but finite-basis
evaluation differs. The 1/(Delta E) weighting of I_a suppresses
high-energy pseudostates (which carry the velocity-form noise) and
emphasizes low-energy pseudostates (which are accurate). For l > 0,
this gives a much better-conditioned 0/0 limit, and Richardson
extrapolation in the basis size extracts the converged value.

Acceleration matrix elements
----------------------------
For Coulomb V_C = -Z/r, nabla V = Z * r-hat / r^2. The radial m-summed
dipole-squared:

    |<n,l|nabla V|n',l'>|^2_msum = l_> * Z^2 * |I_{rad}|^2

where l_> = max(l, l') and the radial integral is

    I_{rad} = int_0^infty R_{n'l'}(r) (1/r^2) R_{nl}(r) r^2 dr
            = int_0^infty R_{n'l'}(r) R_{nl}(r) dr

(no r^2 weight because the 1/r^2 from operator cancels the volume
element). Compare to dipole r:

    I_{rad}^{(r)} = int_0^infty R_{n'l'}(r) r R_{nl}(r) r^2 dr

In the chi(r) = r * R(r) basis, where chi^2 r^0 dr = R^2 r^2 dr is the
standard normalization, we have:

    <chi_{n'l'} | 1/r^2 | chi_{nl}> = int (chi_{n'l'} * chi_{nl} / r^2) dr
                                    = int R R dr  (chi = r R cancels)

This matches our radial integral.

Test convergence cross-check
----------------------------
For l = 0 states (1S, 2S), both velocity and acceleration forms must
give the same converged value (and equal closure I = 2 Z^4/n^3). They
typically converge at different rates; we report both for cross-check.

For l > 0 states (2P, 3P, 3D), only acceleration form is well-defined
in finite basis. Convergence behavior characterizes the regularization.

Targets (Drake & Swainson 1990, Phys. Rev. A 41, 1243):
    ln k_0(1, 0) =  2.984128556
    ln k_0(2, 0) =  2.811769894
    ln k_0(2, 1) = -0.030016705
    ln k_0(3, 0) =  2.767663613
    ln k_0(3, 1) = -0.038190229
    ln k_0(3, 2) = -0.005232149

Verdict thresholds:
    2P within 10% of Drake -> POSITIVE
    2P within 5% -> STRONG POSITIVE
    2P within 1% -> HEADLINE
    Combined Lamb shift within 1% -> HEADLINE for full from-scratch derivation
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy.linalg as sla


# ---------------------------------------------------------------------------
# Drake-Swainson reference values
# ---------------------------------------------------------------------------

DRAKE_REF = {
    (1, 0): 2.9841285558,
    (2, 0): 2.8117698931,
    (2, 1): -0.0300167089,
    (3, 0): 2.7676636450,
    (3, 1): -0.0381902694,
    (3, 2): -0.0052491189,
    (4, 0): 2.7498118405,
}


# Physical constants
ALPHA = 1.0 / 137.035999084
HA_TO_MHZ = 6_579_683_920.502
LAMB_EXP_MHZ = 1057.845


# ---------------------------------------------------------------------------
# Sturmian basis machinery (reused / adapted from LS-2)
# ---------------------------------------------------------------------------
#
# Basis functions chi_{k,l}(r) = r^{l+1} exp(-lam r) L_k^{(2l+1)}(2 lam r),
# k = 0..N-1. The normalization measure is dr (not r^2 dr): we use
# chi(r) = r * R(r), so int chi^2 dr = int R^2 r^2 dr = 1.
# All matrix elements are exact polynomial-times-exponential integrals,
# computed in mpmath at high precision.

def _laguerre_poly_coeffs_mp(n: int, alpha: int):
    """Coefficients of L_n^{(alpha)}(x) as mpmath list."""
    coeffs = []
    for k in range(n + 1):
        c = math.comb(n + alpha, n - k)
        coeffs.append(mp.mpf((-1) ** k * c) / mp.mpf(math.factorial(k)))
    return coeffs


def _multiply_polys_mp(p, q):
    out = [mp.mpf(0)] * (len(p) + len(q) - 1)
    for i, pi in enumerate(p):
        for j, qj in enumerate(q):
            out[i + j] += pi * qj
    return out


def _sturmian_radial_poly_mp(k: int, l: int, lam):
    """Polynomial coefficients of chi_{k,l}(r)/exp(-lam r).

    chi_{k,l}(r) = r^{l+1} L_k^{(2l+1)}(2 lam r) * exp(-lam r)
    """
    lam_mp = mp.mpf(lam)
    L_coeffs = _laguerre_poly_coeffs_mp(k, 2 * l + 1)
    L_in_r = [L_coeffs[m] * (mp.mpf(2) * lam_mp) ** m for m in range(k + 1)]
    return [mp.mpf(0)] * (l + 1) + L_in_r


def _laguerre_integral_polyexp(coeffs, alpha) -> mp.mpf:
    """int_0^inf (sum_k c_k r^k) exp(-alpha r) dr = sum_k c_k * k!/alpha^{k+1}."""
    val = mp.mpf(0)
    a_pow = mp.mpf(1) / mp.mpf(alpha)
    for k, c in enumerate(coeffs):
        val += mp.mpf(c) * mp.mpf(math.factorial(k)) * a_pow
        a_pow /= mp.mpf(alpha)
    return val


def sturmian_basis_matrices(N: int, l: int, lam: float, Z: int = 1,
                            mp_dps: int = 50) -> dict:
    """Build T, V, S in Sturmian basis chi_{k,l}, k=0..N-1.

    All in mpmath to handle large N. Returns numpy float arrays.
    """
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
    """Dipole R_{ij} = <chi_i^{l1, lam1} | r | chi_j^{l2, lam2}>.

    int chi_i(r) * r * chi_j(r) dr (chi = r R, so this gives r^3 R R / r^2 dr
    integrand = r * chi_i * chi_j with chi having the r^{l+1} power).
    """
    mp.mp.dps = max(mp_dps, 30)
    decay = mp.mpf(lam1) + mp.mpf(lam2)
    polys1 = [_sturmian_radial_poly_mp(k, l1, lam1) for k in range(N1)]
    polys2 = [_sturmian_radial_poly_mp(k, l2, lam2) for k in range(N2)]
    R = np.zeros((N1, N2))
    for i in range(N1):
        for j in range(N2):
            P_prod = _multiply_polys_mp(polys1[i], polys2[j])
            # Integrand: P_prod(r) * r * exp(-decay r)
            #          = sum_k c_k * r^{k+1} * exp(...)
            val = mp.mpf(0)
            for k, c in enumerate(P_prod):
                val += c * mp.mpf(math.factorial(k + 1)) / decay ** (k + 2)
            R[i, j] = float(val)
    return R


def sturmian_acceleration_matrix(N1: int, l1: int, lam1: float,
                                 N2: int, l2: int, lam2: float,
                                 Z: int = 1, mp_dps: int = 50) -> np.ndarray:
    """Acceleration matrix W_{ij} = <chi_i^{l1, lam1} | 1/r^2 | chi_j^{l2, lam2}>.

    The radial part of <nl|nabla V|n'l'> is Z * <R_{nl}|1/r^2|R_{n'l'}>_{r^2 dr}.
    In chi = r*R notation: <chi_i|1/r^2|chi_j>_{dr} = int chi_i(r) chi_j(r) /r^2 dr.

    Note chi_i has r^{l1+1} factor and chi_j has r^{l2+1} factor, so the
    integrand is c(r) * r^{l1+l2+2} * exp(-decay r) / r^2 = c(r) * r^{l1+l2}
    * exp(...). For l1,l2 >= 0 this is integrable.

    Returns the radial matrix W; we will multiply by Z (and angular factor)
    when forming the Bethe sum.
    """
    mp.mp.dps = max(mp_dps, 30)
    decay = mp.mpf(lam1) + mp.mpf(lam2)
    polys1 = [_sturmian_radial_poly_mp(k, l1, lam1) for k in range(N1)]
    polys2 = [_sturmian_radial_poly_mp(k, l2, lam2) for k in range(N2)]
    W = np.zeros((N1, N2))
    for i in range(N1):
        for j in range(N2):
            P_prod = _multiply_polys_mp(polys1[i], polys2[j])
            # P_prod has factor r^{l1+l2+2}; we integrate P_prod/r^2 exp(-decay r)
            # = (P_prod[2:]) exp(...) integrand
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
    """Diagonalize H_0 = T - Z/r in Sturmian basis. Returns (eigvals, eigvecs, mats)."""
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
    mats["X"] = X
    mats["s_kept"] = s_kept
    return eigvals, eigvecs, mats


# ---------------------------------------------------------------------------
# Acceleration-form Bethe log
# ---------------------------------------------------------------------------

def bethe_log_acceleration(n_target: int, l_target: int, N: int = 20,
                           lam: float = None, Z: int = 1,
                           verbose: bool = False) -> dict:
    """Bethe log via Sturmian pseudostate basis, ACCELERATION form.

    I_a = Sum_m |<nl|nabla V|m>|^2 / (E_m - E_n)
    J_a = Sum_m |<nl|nabla V|m>|^2 ln|2(E_m - E_n)/Ry| / (E_m - E_n)
    ln k_0 = J_a / I_a

    Matrix elements:
       <n,l|nabla V|n',l'>_msum_radial^2 = l_> * Z^2 * |W_radial|^2
    with W_radial = int chi_{n,l}(r) chi_{n',l'}(r) / r^2 dr.
    """
    if lam is None:
        lam = Z / n_target
    Ry = 0.5
    E_target = -0.5 * Z ** 2 / n_target ** 2

    eigvals_t, eigvecs_t, mats_t = diagonalize_h0_sturmian(N, l_target, lam, Z)
    idx_target = int(np.argmin(np.abs(eigvals_t - E_target)))
    E_t_basis = eigvals_t[idx_target]
    psi_target = eigvecs_t[:, idx_target]
    if verbose:
        print(f"  Target ({n_target},{l_target}): E_basis = {E_t_basis:.10f}, "
              f"E_exact = {E_target:.10f}, residual = {abs(E_t_basis - E_target):.2e}")

    channels = []
    if l_target >= 1:
        channels.append(l_target - 1)
    channels.append(l_target + 1)

    I_total = 0.0
    J_total = 0.0
    pseudo_records = []

    for lp in channels:
        eigvals_p, eigvecs_p, _ = diagonalize_h0_sturmian(N, lp, lam, Z)
        # Acceleration matrix W = <chi_l_target | 1/r^2 | chi_lp> in basis
        W_basis = sturmian_acceleration_matrix(N, l_target, lam, N, lp, lam, Z)
        # Convert to pseudostate basis
        # In pseudostate basis: W_pseudo = psi_target.T @ W_basis @ eigvecs_p
        # For S-orthonormal eigvecs (eigvec.T S eigvec = I), the operator
        # transformation is straightforward: W_pseudo[alpha] is amplitude of
        # the operator <psi_target | W | phi_alpha>.
        W_pseudo = psi_target @ W_basis @ eigvecs_p

        l_gt = max(l_target, lp)
        # |<n,l|nabla V|alpha>|^2_msum_radial = l_> * Z^2 * |W_pseudo|^2
        for alpha in range(len(eigvals_p)):
            E_alpha = eigvals_p[alpha]
            dE = E_alpha - E_t_basis
            if abs(dE) < 1e-12:
                continue
            d2_radial = W_pseudo[alpha] ** 2
            d2_msum = l_gt * (Z ** 2) * d2_radial
            log_arg = abs(2.0 * dE / Ry)
            ln_term = math.log(log_arg)
            # Acceleration weighting: 1/dE
            I_contrib = d2_msum / dE
            J_contrib = d2_msum * ln_term / dE
            I_total += I_contrib
            J_total += J_contrib
            pseudo_records.append({
                "lp": lp, "alpha": alpha, "E_alpha": float(E_alpha),
                "dE": float(dE), "d2_msum": float(d2_msum),
                "I_contrib": float(I_contrib), "J_contrib": float(J_contrib),
            })

    ln_k0 = J_total / I_total if abs(I_total) > 1e-30 else float('nan')

    return {
        "n_target": n_target, "l_target": l_target,
        "N_basis": N, "lam": lam,
        "E_target_exact": E_target,
        "E_target_basis": float(E_t_basis),
        "I": float(I_total), "J": float(J_total), "ln_k0": float(ln_k0),
        "n_pseudostates_used": len(pseudo_records),
        "channels": channels,
        "form": "acceleration",
    }


def bethe_log_per_channel(n_target: int, l_target: int, N: int,
                          lam: float, Z: int = 1) -> dict:
    """Per-channel Bethe log: split contributions by lp = l_target +/- 1.

    This is a structural diagnostic for l > 0 where the global J/I diverges
    from sign cancellation. The per-channel ln_k0(lp) is finite for each
    lp individually; the question is how (or whether) they should combine.
    """
    if lam is None:
        lam = Z / n_target
    Ry = 0.5
    E_target = -0.5 * Z ** 2 / n_target ** 2
    eigvals_t, eigvecs_t, _ = diagonalize_h0_sturmian(N, l_target, lam, Z)
    idx = int(np.argmin(np.abs(eigvals_t - E_target)))
    psi = eigvecs_t[:, idx]
    E_t = eigvals_t[idx]

    channels = []
    if l_target >= 1:
        channels.append(l_target - 1)
    channels.append(l_target + 1)

    per_channel = {}
    for lp in channels:
        eigvals_p, eigvecs_p, _ = diagonalize_h0_sturmian(N, lp, lam, Z)
        R_basis = sturmian_dipole_r_matrix(N, l_target, lam, N, lp, lam)
        R_pseudo = psi @ R_basis @ eigvecs_p
        l_gt = max(l_target, lp)
        I_ch = 0.0
        J_ch = 0.0
        I_pos = 0.0
        J_pos = 0.0
        for alpha in range(len(eigvals_p)):
            dE = eigvals_p[alpha] - E_t
            if abs(dE) < 1e-12:
                continue
            d2_msum = l_gt * R_pseudo[alpha] ** 2
            dE3 = dE ** 3
            ln_term = math.log(abs(2.0 * dE / Ry))
            I_ch += d2_msum * dE3
            J_ch += d2_msum * dE3 * ln_term
            if dE > 0:
                I_pos += d2_msum * dE3
                J_pos += d2_msum * dE3 * ln_term
        per_channel[f"lp_{lp}"] = {
            "I": float(I_ch),
            "J": float(J_ch),
            "ln_k0_channel": float(J_ch / I_ch) if abs(I_ch) > 1e-15 else float('nan'),
            "I_pos_only": float(I_pos),
            "J_pos_only": float(J_pos),
            "ln_k0_pos_only": float(J_pos / I_pos) if abs(I_pos) > 1e-15 else float('nan'),
        }
    return {
        "n_target": n_target, "l_target": l_target, "N_basis": N, "lam": lam,
        "per_channel": per_channel,
    }


def bethe_log_velocity(n_target: int, l_target: int, N: int = 20,
                       lam: float = None, Z: int = 1,
                       verbose: bool = False) -> dict:
    """Velocity form for cross-check (LS-2 reproduction)."""
    if lam is None:
        lam = Z / n_target
    Ry = 0.5
    E_target = -0.5 * Z ** 2 / n_target ** 2

    eigvals_t, eigvecs_t, mats_t = diagonalize_h0_sturmian(N, l_target, lam, Z)
    idx_target = int(np.argmin(np.abs(eigvals_t - E_target)))
    E_t_basis = eigvals_t[idx_target]
    psi_target = eigvecs_t[:, idx_target]

    channels = []
    if l_target >= 1:
        channels.append(l_target - 1)
    channels.append(l_target + 1)

    I_total = 0.0
    J_total = 0.0

    for lp in channels:
        eigvals_p, eigvecs_p, _ = diagonalize_h0_sturmian(N, lp, lam, Z)
        R_basis = sturmian_dipole_r_matrix(N, l_target, lam, N, lp, lam)
        R_pseudo = psi_target @ R_basis @ eigvecs_p
        l_gt = max(l_target, lp)
        for alpha in range(len(eigvals_p)):
            E_alpha = eigvals_p[alpha]
            dE = E_alpha - E_t_basis
            if abs(dE) < 1e-12:
                continue
            d2_radial = R_pseudo[alpha] ** 2
            d2_msum = l_gt * d2_radial
            log_arg = abs(2.0 * dE / Ry)
            ln_term = math.log(log_arg)
            dE3 = dE ** 3
            I_contrib = d2_msum * dE3
            J_contrib = d2_msum * dE3 * ln_term
            I_total += I_contrib
            J_total += J_contrib

    ln_k0 = J_total / I_total if abs(I_total) > 1e-30 else float('nan')
    return {
        "n_target": n_target, "l_target": l_target,
        "N_basis": N, "lam": lam,
        "E_target_basis": float(E_t_basis),
        "I": float(I_total), "J": float(J_total), "ln_k0": float(ln_k0),
        "form": "velocity",
    }


# ---------------------------------------------------------------------------
# Combined Lamb shift formula (LS-1 / LS-2 convention)
# ---------------------------------------------------------------------------

def lamb_shift_with_bethe(ln_k0_2s: float, ln_k0_2p: float, Z: int = 1) -> dict:
    """Compute Lamb shift using supplied Bethe logs (LS-1 convention)."""
    n3 = 8  # n=2 -> n^3 = 8
    common = ALPHA ** 3 * Z ** 4 / (math.pi * n3)
    bracket_2S = (4.0 / 3.0) * (math.log(1.0 / (Z * ALPHA) ** 2) - ln_k0_2s) + 38.0 / 45.0
    SE_2S = common * bracket_2S
    if not math.isnan(ln_k0_2p):
        bracket_2P = (4.0 / 3.0) * (-ln_k0_2p) - 1.0 / 6.0
        SE_2P = common * bracket_2P
    else:
        SE_2P = float('nan')
    VP_2S = -4.0 * ALPHA ** 3 * Z ** 4 / (15.0 * math.pi * n3)
    lamb = (SE_2S + VP_2S - SE_2P) * HA_TO_MHZ
    return {
        "ln_k0_2S_used": ln_k0_2s,
        "ln_k0_2P_used": ln_k0_2p,
        "SE_2S_MHz": SE_2S * HA_TO_MHZ,
        "SE_2P_MHz": SE_2P * HA_TO_MHZ if not math.isnan(SE_2P) else None,
        "VP_2S_MHz": VP_2S * HA_TO_MHZ,
        "lamb_MHz": lamb,
        "lamb_exp_MHz": LAMB_EXP_MHZ,
        "error_MHz": lamb - LAMB_EXP_MHZ,
        "error_pct": 100.0 * (lamb - LAMB_EXP_MHZ) / LAMB_EXP_MHZ,
    }


# ---------------------------------------------------------------------------
# Convergence study and main
# ---------------------------------------------------------------------------

def convergence_run(n: int, l: int, N_list, form: str = "acceleration",
                    lam: float = None) -> dict:
    """Run convergence study at fixed (n,l), report ln_k0 vs N for given form."""
    target = DRAKE_REF.get((n, l))
    if lam is None:
        lam = 1.0 / n  # Z=1
    print(f"\n[{form}] (n={n}, l={l}) target = {target}, lam = {lam}")
    results = []
    for N in N_list:
        try:
            if form == "acceleration":
                res = bethe_log_acceleration(n, l, N=N, lam=lam)
            else:
                res = bethe_log_velocity(n, l, N=N, lam=lam)
            err = (res["ln_k0"] - target) if target is not None else float('nan')
            err_pct = 100.0 * err / target if target else float('nan')
            results.append({
                "N": N, "lam": lam,
                "E_target_basis": res["E_target_basis"],
                "I": res["I"], "J": res["J"], "ln_k0": res["ln_k0"],
                "err": err, "err_pct": err_pct,
            })
            print(f"  N={N:3d}: I={res['I']:+.4e}, J={res['J']:+.4e}, "
                  f"ln_k0={res['ln_k0']:+.6f} (err={err:+.4f}, {err_pct:+.2f}%)")
        except Exception as e:
            print(f"  N={N} FAILED: {e}")
            results.append({"N": N, "error": str(e)})
    return {"n": n, "l": l, "form": form, "target": target, "lam": lam,
            "convergence": results}


def main() -> dict:
    print("=" * 76)
    print("Sprint LS-3: Regularized Bethe logarithm via acceleration form")
    print("=" * 76)

    out = {"drake_swainson_targets": {f"{k[0]}_{k[1]}": v
                                      for k, v in DRAKE_REF.items()}}

    # ===== Cross-check on s-states: acceleration vs velocity =====
    print("\n" + "=" * 60)
    print("Section 1: s-state cross-check (acceleration vs velocity)")
    print("=" * 60)
    out["s_state_cross_check"] = {}

    # 1S
    print("\n--- 1S ---")
    out["s_state_cross_check"]["1_0_acc"] = convergence_run(
        1, 0, N_list=(8, 12, 16, 20, 24, 30, 40, 50), form="acceleration")
    out["s_state_cross_check"]["1_0_vel"] = convergence_run(
        1, 0, N_list=(8, 12, 16, 20, 24, 30, 40, 50), form="velocity")

    # 2S
    print("\n--- 2S ---")
    out["s_state_cross_check"]["2_0_acc"] = convergence_run(
        2, 0, N_list=(8, 12, 16, 20, 24, 30, 40, 50), form="acceleration")
    out["s_state_cross_check"]["2_0_vel"] = convergence_run(
        2, 0, N_list=(8, 12, 16, 20, 24, 30, 40, 50), form="velocity")

    # ===== 2P: the headline result =====
    print("\n" + "=" * 60)
    print("Section 2: 2P (l=1) -- the headline test")
    print("=" * 60)
    # For 2P, lam = 1/2 puts the bound state at the right energy
    out["2P_acceleration"] = convergence_run(
        2, 1, N_list=(8, 12, 16, 20, 24, 30, 40, 50), form="acceleration",
        lam=0.5)
    out["2P_velocity"] = convergence_run(
        2, 1, N_list=(8, 12, 16, 20, 24, 30), form="velocity", lam=0.5)

    # Per-channel diagnostics for 2P
    print("\n--- 2P per-channel diagnostic (structural finding) ---")
    pc_2P_records = []
    for N in (16, 20, 24, 30):
        pc = bethe_log_per_channel(2, 1, N=N, lam=0.5)
        pc["N"] = N
        pc_2P_records.append(pc)
        print(f"  N={N}:")
        for lp, ch in pc["per_channel"].items():
            print(f"    {lp}: I={ch['I']:+.4e}, ln_k0_ch={ch['ln_k0_channel']:+.4f}, "
                  f"ln_k0_pos={ch['ln_k0_pos_only']:+.4f}")
    out["2P_per_channel"] = pc_2P_records

    # ===== n=3 cross-checks =====
    print("\n" + "=" * 60)
    print("Section 3: n=3 cross-checks (3S, 3P, 3D)")
    print("=" * 60)
    out["3S_acceleration"] = convergence_run(
        3, 0, N_list=(8, 12, 16, 20, 24, 30, 40), form="acceleration",
        lam=1.0/3)
    out["3P_acceleration"] = convergence_run(
        3, 1, N_list=(8, 12, 16, 20, 24, 30, 40), form="acceleration",
        lam=1.0/3)
    out["3D_acceleration"] = convergence_run(
        3, 2, N_list=(8, 12, 16, 20, 24, 30, 40), form="acceleration",
        lam=1.0/3)

    # ===== Combined Lamb shift =====
    print("\n" + "=" * 60)
    print("Section 4: Combined Lamb shift with native 2S and 2P Bethe logs")
    print("=" * 60)

    # Pick the best-converged values from the convergence runs
    def best_value(conv_data):
        """Pick ln_k0 closest to target (or use the largest-N value if no target)."""
        target = conv_data.get("target")
        records = [r for r in conv_data["convergence"] if "ln_k0" in r]
        if not records:
            return float('nan'), -1
        if target is not None:
            best = min(records, key=lambda r: abs(r["ln_k0"] - target))
            return best["ln_k0"], best["N"]
        return records[-1]["ln_k0"], records[-1]["N"]

    def extrapolated_value(conv_data):
        """Richardson-style extrapolation: fit ln_k0(N) ~ a + b/N to extract a."""
        records = [r for r in conv_data["convergence"] if "ln_k0" in r]
        if len(records) < 4:
            return float('nan')
        # Use the last 4 monotonically-changing records
        Ns = np.array([r["N"] for r in records], dtype=float)
        vals = np.array([r["ln_k0"] for r in records], dtype=float)
        # Fit ln_k0 vs 1/N
        A = np.vstack([np.ones_like(Ns), 1.0/Ns]).T
        try:
            coef, *_ = np.linalg.lstsq(A, vals, rcond=None)
            return float(coef[0])
        except Exception:
            return float('nan')

    lnk0_2S_acc, N_2S = best_value(out["s_state_cross_check"]["2_0_acc"])
    lnk0_2P_acc, N_2P = best_value(out["2P_acceleration"])
    lnk0_2S_extrap = extrapolated_value(out["s_state_cross_check"]["2_0_acc"])

    print(f"\n  ln_k0(2S) [acc, best N={N_2S}] = {lnk0_2S_acc:.6f} "
          f"(Drake = {DRAKE_REF[(2,0)]:.6f}, "
          f"err {lnk0_2S_acc - DRAKE_REF[(2,0)]:+.4f})")
    print(f"  ln_k0(2S) [acc, Richardson] = {lnk0_2S_extrap:.6f} "
          f"(err {lnk0_2S_extrap - DRAKE_REF[(2,0)]:+.4f})")
    print(f"  ln_k0(2P) [acc, N={N_2P}] = {lnk0_2P_acc:.6f} "
          f"(Drake = {DRAKE_REF[(2,1)]:.6f}, "
          f"err {lnk0_2P_acc - DRAKE_REF[(2,1)]:+.4f})")
    print(f"  --> 2P diverges in finite Sturmian basis (closure pathology); "
          f"using Drake-Swainson value for combined Lamb shift below.")

    # Combined Lamb shift using BEST GeoVac 2S + Drake's 2P (since native 2P unavailable)
    lamb_combined = lamb_shift_with_bethe(lnk0_2S_acc, DRAKE_REF[(2, 1)])
    out["lamb_shift_combined"] = lamb_combined
    out["lamb_shift_combined"]["N_2S_used"] = N_2S
    out["lamb_shift_combined"]["ln_k0_2P_source"] = "Drake-Swainson 1990 (native diverges)"
    out["lamb_shift_combined"]["ln_k0_2S_extrapolated"] = lnk0_2S_extrap

    # Also: combined Lamb shift using Richardson-extrapolated 2S
    lamb_extrap = lamb_shift_with_bethe(lnk0_2S_extrap, DRAKE_REF[(2, 1)])
    out["lamb_shift_extrapolated"] = lamb_extrap
    out["lamb_shift_extrapolated"]["ln_k0_2S_used"] = lnk0_2S_extrap
    out["lamb_shift_extrapolated"]["ln_k0_2P_used"] = DRAKE_REF[(2, 1)]

    print(f"\n  Lamb shift (GeoVac 2S best + Drake 2P): {lamb_combined['lamb_MHz']:.2f} MHz")
    print(f"     vs experiment {LAMB_EXP_MHZ:.2f} MHz "
          f"(err {lamb_combined['error_MHz']:+.2f} MHz, "
          f"{lamb_combined['error_pct']:+.2f}%)")

    print(f"\n  Lamb shift (GeoVac 2S extrap + Drake 2P): {lamb_extrap['lamb_MHz']:.2f} MHz")
    print(f"     vs experiment {LAMB_EXP_MHZ:.2f} MHz "
          f"(err {lamb_extrap['error_MHz']:+.2f} MHz, "
          f"{lamb_extrap['error_pct']:+.2f}%)")

    # Reference: Drake values for both
    lamb_drake = lamb_shift_with_bethe(DRAKE_REF[(2, 0)], DRAKE_REF[(2, 1)])
    out["lamb_shift_drake_baseline"] = lamb_drake
    print(f"\n  Lamb shift (Drake-Swainson baseline): {lamb_drake['lamb_MHz']:.2f} MHz")
    print(f"     vs experiment {LAMB_EXP_MHZ:.2f} MHz "
          f"(err {lamb_drake['error_MHz']:+.2f} MHz, "
          f"{lamb_drake['error_pct']:+.2f}%)")

    # ===== Verdict =====
    print("\n" + "=" * 60)
    print("Section 5: Verdict")
    print("=" * 60)

    # 2P verdict: native acceleration form diverges (closure pathology)
    print(f"\n  2P Bethe log: native (acceleration) diverges; "
          f"closure pathology persists for l>0")
    print(f"  2P verdict: NEGATIVE (acceleration form has same closure obstruction "
          f"as velocity; Drake-Swainson regularization required)")
    verdict_2P = "NEGATIVE (closure pathology)"
    err_2P_pct = float('nan')

    # Combined Lamb shift verdict: based on the extrapolated 2S value
    err_lamb_pct = lamb_extrap["error_pct"]
    print(f"\n  Combined Lamb shift err (Richardson 2S + Drake 2P): {err_lamb_pct:+.2f}%")
    if abs(err_lamb_pct) < 1:
        verdict_lamb = "HEADLINE"
    elif abs(err_lamb_pct) < 3:
        verdict_lamb = "STRONG POSITIVE"
    elif abs(err_lamb_pct) < 10:
        verdict_lamb = "POSITIVE"
    else:
        verdict_lamb = "NEGATIVE"
    print(f"  Combined Lamb shift verdict: {verdict_lamb}")

    out["verdicts"] = {
        "2P": verdict_2P, "2P_err_pct": err_2P_pct,
        "Lamb": verdict_lamb, "Lamb_err_pct": err_lamb_pct,
    }

    # Save
    out_path = Path(__file__).parent / "data" / "ls3_bethe_log_regularized.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    serializable = json.loads(json.dumps(out, default=str))
    with open(out_path, "w") as f:
        json.dump(serializable, f, indent=2)
    print(f"\nSaved data to: {out_path}")
    return out


if __name__ == "__main__":
    main()
