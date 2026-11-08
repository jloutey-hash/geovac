"""Sprint LS-2: Bethe logarithm from GeoVac spectral machinery.

Goal
----
Compute the hydrogen Bethe logarithm ln k_0(n, l) from a GeoVac-native
spectral mode sum, replacing the tabulated Drake-Swainson 1990 values
used in LS-1.

Definition
----------

    ln k_0(n, l) = J(n,l) / I(n,l)

with

    I(n,l) = Σ_m |⟨nl| p |m⟩|² (E_m - E_n)
    J(n,l) = Σ_m |⟨nl| p |m⟩|² (E_m - E_n) ln |2 (E_m - E_n) / Ry|

Sum runs over the COMPLETE Coulomb spectrum (bound + continuum). Selection
rule Δl = ±1 from p = -i∇ on a central-field eigenstate.

Equivalent (acceleration / length) form
---------------------------------------
With m=1 in atomic units, ⟨n'l'|p|nl⟩ = i (E_{n'} - E_n) ⟨n'l'|r|nl⟩, so
|⟨n'l'|p|nl⟩|² = (ΔE)² |⟨n'l'|r|nl⟩|² and the Bethe sums become

    I(n,l) = Σ_m |⟨nl|r|m⟩|² (E_m - E_n)³
    J(n,l) = Σ_m |⟨nl|r|m⟩|² (E_m - E_n)³ ln|2(E_m - E_n)/Ry|

Closure identity (velocity form) at Z=1:

    I(n, l) = 2 Z⁴ / n³ · δ_{l, 0}

(from [p, [H_0, p]] = -∇²V and |ψ_{n0}(0)|² = Z³/(π n³).)

For l = 0 this is finite; for l > 0 the velocity-form closure VANISHES,
and ln k_0 must be obtained by L'Hôpital or, equivalently, by computing
both I and J at finite cutoff and observing that their ratio is finite
(both → 0 as the basis becomes complete, but the ratio J/I tends to a
finite limit).

Routes
------
ROUTE A — bound-state-only finite sum (m runs over physical bound states
n' = max(l'+1,1)..n_max). Closes only a small fraction of the closure
sum-rule; serves as a sanity check.

ROUTE B — finite SLATER-LAGUERRE pseudostate basis. Diagonalize H_0 =
T - Z/r in a finite radial Laguerre basis at fixed exponent λ; the
eigenvalues split into a few accurate bound states plus a sequence of
positive-energy 'pseudostates' that discretize the continuum. As the
basis size N grows, ln k_0 converges to the exact value. This is the
GeoVac-native route — it's algebraic-first (matrix elements are exact
rationals from Laguerre recurrences), and it makes ln k_0 a finite
matrix expression whose limit as N → ∞ is the physical value.

ROUTE C — pure Coulomb-Sturmian basis (= Laguerre at exponent λ = Z/n
specifically). Special case of B with a fixed-exponent basis.

Both B and C are direct implementations of the GeoVac philosophy:
'continuum is approximated by discrete pseudostates in a graph-spectral
basis'.

Targets (Drake & Swainson 1990, Phys. Rev. A 41, 1243):
    ln k_0(2, 0) =  2.811769893
    ln k_0(2, 1) = -0.030016709

References
----------
- W. Gordon, Ann. Phys. 394 (1929) 1031 — closed form for bound dipole
- M. Lieber, *Phys. Rev.* 174 (1968) 2037 — hydrogenic Bethe sums
- C. Schwartz, *Phys. Rev.* 123 (1961) 1700 — finite-basis Bethe log
- S. P. Goldman & G. W. F. Drake, *Phys. Rev. A* 28 (1983) 1228 —
    pseudostate method for Bethe logs
- G. W. F. Drake, R. A. Swainson, *Phys. Rev. A* 41 (1990) 1243 — table
- F. Salvat & R. Mayol, *Comput. Phys. Commun.* 62 (1991) 65 — radial
    Coulomb wavefunctions
"""

from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy.linalg as sla
import sympy as sp

# Drake & Swainson 1990 reference values for ln k_0(n, l) at Z=1
DRAKE_REF = {
    (1, 0): 2.9841285558,
    (2, 0): 2.8117698931,
    (2, 1): -0.0300167089,
    (3, 0): 2.7676636450,
    (3, 1): -0.0381902694,
    (3, 2): -0.0052491189,
    (4, 0): 2.7498118405,
}


# ---------------------------------------------------------------------------
# Hydrogenic radial wavefunction (sympy / closed form)
# ---------------------------------------------------------------------------

def hydrogenic_radial_R(n: int, l: int, r):
    """Sympy expression for the normalized hydrogenic radial function R_{nl}(r).

    Z=1 atomic units. Convention: ∫₀^∞ R_{nl}(r)² r² dr = 1.
    """
    Z = 1
    rho = 2 * Z * r / n
    norm_sq = ((2 * Z / n) ** 3
               * sp.factorial(n - l - 1)
               / (2 * n * sp.factorial(n + l)))
    norm = sp.sqrt(norm_sq)
    L = sp.assoc_laguerre(n - l - 1, 2 * l + 1, rho)
    return norm * sp.exp(-rho / 2) * rho ** l * L


def radial_dipole_sq_msummed(n: int, l: int, np_: int, lp: int) -> float:
    """|⟨n l | r | n' l'⟩|², summed over m_l final, AVERAGED over m_l initial.

    This is the m-summed dipole that appears in standard oscillator-strength
    and Bethe-log sums for an unpolarized initial state. The radial integral
    is

        R_{n l → n' l'} = ∫₀^∞ R_{n'l'}(r) R_{nl}(r) r³ dr

    and the angular average gives

        |⟨nl|r|n'l'⟩|²_m-summed = (l_> / (2l+1)) · |R|² · (2l+1)
                                = l_> · |R|²

    where l_> = max(l, l'). Verified against the f-sum rule:
        f(1s → 2p) = (2/3)(E_2 - E_1)|⟨1s|r|2p⟩|²_msum = 0.4162 ✓.

    For an *initial* l-shell of degeneracy 2l+1, summing over m_initial
    instead of averaging would give an extra factor (2l+1). The Bethe log
    is m-AVERAGED (single state |nlm⟩ for each m_initial), so we use
    l_> · |R|².
    """
    if abs(lp - l) != 1:
        return 0.0
    r = sp.Symbol('r', positive=True)
    R1 = hydrogenic_radial_R(n, l, r)
    R2 = hydrogenic_radial_R(np_, lp, r)
    integrand = R1 * R2 * r ** 3
    R_radial = float(sp.integrate(integrand, (r, 0, sp.oo)))
    l_gt = max(l, lp)
    return l_gt * R_radial ** 2


# ---------------------------------------------------------------------------
# Energies and basic Bethe sum machinery
# ---------------------------------------------------------------------------

def hydrogen_bound_energy(n: int) -> float:
    """E_n in Hartree at Z=1: -1/(2 n²)."""
    return -0.5 / (n * n)


def closure_I_velocity(n: int, l: int, Z: int = 1) -> float:
    """Exact closure identity I = 2 Z⁴/n³ · δ_{l,0} from [p,[H_0,p]]."""
    return 2.0 * Z ** 4 / n ** 3 if l == 0 else 0.0


# ---------------------------------------------------------------------------
# ROUTE A: bound-state-only Bethe sum
# ---------------------------------------------------------------------------

def bethe_sums_bound(n: int, l: int, n_max: int) -> dict:
    """Compute I, J Bethe sums over BOUND hydrogen states up to n' = n_max.

    Selection rule Δl = ±1 enforced. The state itself excluded.
    Returns I, J, ln_k0, list of contributions.
    """
    E_n = hydrogen_bound_energy(n)
    Ry = 0.5

    channels = []
    if l >= 1:
        channels.append(l - 1)
    channels.append(l + 1)

    contributions = []
    I_total = 0.0
    J_total = 0.0

    for lp in channels:
        n_start = lp + 1
        for np_ in range(n_start, n_max + 1):
            if np_ == n and lp == l:
                continue
            E_m = hydrogen_bound_energy(np_)
            dE = E_m - E_n
            if abs(dE) < 1e-30:
                continue
            d2 = radial_dipole_sq_msummed(n, l, np_, lp)
            if d2 == 0:
                continue
            log_arg = abs(2.0 * dE / Ry)
            ln_term = math.log(log_arg)
            dE3 = dE ** 3
            I_contrib = d2 * dE3
            J_contrib = d2 * dE3 * ln_term
            I_total += I_contrib
            J_total += J_contrib
            contributions.append({
                "np": np_, "lp": lp, "dE": dE, "d2_msum": d2,
                "I_contrib": I_contrib, "J_contrib": J_contrib,
            })

    ln_k0 = J_total / I_total if abs(I_total) > 1e-30 else float('nan')

    return {
        "n": n, "l": l, "n_max": n_max,
        "I": I_total, "J": J_total, "ln_k0": ln_k0,
        "n_states_summed": len(contributions),
        "contributions": contributions,
    }


# ---------------------------------------------------------------------------
# ROUTE B: Slater-Laguerre pseudostate basis
# ---------------------------------------------------------------------------
#
# We use the basis of Slater-type orbitals (STOs) at fixed l, with a
# common Laguerre exponent λ:
#
#   φ_{k, l}(r) = N_{k l λ} · r^l · exp(-λ r) · L_k^{(2l+2)}(2 λ r)    k = 0, 1, 2, ...
#
# This is the standard 'Coulomb-Sturmian' basis when λ = Z/n_target, and
# more generally a complete orthogonal basis on L²(r²dr, [0,∞)) for any
# fixed λ > 0. (Reference: Avery & Avery 2006.)
#
# Matrix elements:
#   T_kl, V_kl = -Z<1/r>_kl, S_kl  -> T - Z/r  diagonalization gives
#   eigenvalues that include accurate bound states for n satisfying
#   n*λ = Z (bound) and pseudostates discretizing the continuum.
#
# We compute everything via direct sympy symbolic integration to keep
# all matrix elements exact rationals in λ and Z. Then float-diagonalize.

def _laguerre_integral_polyexp(coeffs, alpha) -> float:
    """Closed-form integral of (sum_k coeffs[k] r^k) * exp(-alpha r) over [0, inf).

    Uses ∫₀^∞ r^k e^{-α r} dr = k! / α^{k+1}. Coeffs may be float or mpmath.
    """
    val = mp.mpf(0)
    a_pow = mp.mpf(1) / mp.mpf(alpha)
    for k, c in enumerate(coeffs):
        val += mp.mpf(c) * mp.mpf(math.factorial(k)) * a_pow
        a_pow /= mp.mpf(alpha)
    return val


def _laguerre_poly_coeffs_mp(n: int, alpha: int):
    """Coefficients of L_n^{(alpha)}(x) as mpmath list.

    L_n^{(α)}(x) = Σ_{k=0}^n (-1)^k C(n+α, n-k) x^k / k!. Integer alpha.
    """
    coeffs = []
    for k in range(n + 1):
        c = math.comb(n + alpha, n - k)
        coeffs.append(mp.mpf((-1) ** k * c) / mp.mpf(math.factorial(k)))
    return coeffs


def _multiply_polys_mp(p, q):
    """Multiply two mpmath polynomial coefficient lists."""
    out = [mp.mpf(0)] * (len(p) + len(q) - 1)
    for i, pi in enumerate(p):
        for j, qj in enumerate(q):
            out[i + j] += pi * qj
    return out


def _sturmian_radial_poly_mp(k: int, l: int, lam):
    """Polynomial coefs of chi_k^{l, lam}(r)/exp(-lam r) = r^{l+1} L_k^{(2l+1)}(2 lam r)."""
    lam_mp = mp.mpf(lam)
    L_coeffs = _laguerre_poly_coeffs_mp(k, 2 * l + 1)
    # L in r: coeff[m] * (2 lam)^m
    L_in_r = [L_coeffs[m] * (mp.mpf(2) * lam_mp) ** m for m in range(k + 1)]
    # Shift by r^{l+1}
    return [mp.mpf(0)] * (l + 1) + L_in_r


def sturmian_basis_matrices(N: int, l: int, lam: float, Z: int = 1,
                            mp_dps: int = 50) -> dict:
    """Build T, V, S in Sturmian basis of N functions at l with exponent lam.

    Basis functions (unnormalized, dr measure):
        chi_{k, l}(r) = r^{l+1} exp(-lam r) L_k^{(2l+1)}(2 lam r),  k = 0..N-1

    Built in mpmath at high precision then converted to numpy float; this
    avoids the catastrophic loss of precision in Laguerre polynomial
    evaluation for k >~ 20 with double arithmetic.
    """
    mp.mp.dps = max(mp_dps, 30)
    lam_mp = mp.mpf(lam)
    decay = 2 * lam_mp

    polys = [_sturmian_radial_poly_mp(k, l, lam) for k in range(N)]
    # Derivatives: d/dr (P(r) exp(-lam r)) = (P'(r) - lam P(r)) exp(-lam r)
    poly_derivs = []
    for poly in polys:
        d_poly = [mp.mpf(0)] * len(poly)
        for m in range(1, len(poly)):
            d_poly[m - 1] = poly[m] * m
        chi_deriv = [d_poly[m] - lam_mp * poly[m] for m in range(len(poly))]
        poly_derivs.append(chi_deriv)

    # Build matrices
    S_mp = [[mp.mpf(0)] * N for _ in range(N)]
    T_mp = [[mp.mpf(0)] * N for _ in range(N)]
    V_mp = [[mp.mpf(0)] * N for _ in range(N)]

    for i in range(N):
        for j in range(i, N):
            P_prod = _multiply_polys_mp(polys[i], polys[j])
            S_ij = _laguerre_integral_polyexp(P_prod, decay)

            D_prod = _multiply_polys_mp(poly_derivs[i], poly_derivs[j])
            T1 = mp.mpf("0.5") * _laguerre_integral_polyexp(D_prod, decay)
            # P_prod has factor r^{2l+2}, so /r^2 → P_prod[2:]
            if len(P_prod) > 2:
                P_div_r2 = P_prod[2:]
                T2_int = mp.mpf(0)
                for k, c in enumerate(P_div_r2):
                    T2_int += c * mp.mpf(math.factorial(k)) / decay ** (k + 1)
            else:
                T2_int = mp.mpf(0)
            T2 = mp.mpf(l * (l + 1)) / 2 * T2_int
            T_ij = T1 + T2

            # V = -Z ∫ chi_i chi_j / r dr = -Z * ∫ P_prod * r^{-1} exp dr
            # P_prod has r^{2l+2} factor, /r → P_prod[1:]
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

    # Convert to numpy
    S = np.array([[float(S_mp[i][j]) for j in range(N)] for i in range(N)])
    T = np.array([[float(T_mp[i][j]) for j in range(N)] for i in range(N)])
    V = np.array([[float(V_mp[i][j]) for j in range(N)] for i in range(N)])

    return {"T": T, "V": V, "S": S, "basis_l": l, "lam": lam, "N": N,
            "S_mp": S_mp, "T_mp": T_mp, "V_mp": V_mp}


def sturmian_dipole_matrix(N1: int, l1: int, lam1: float,
                           N2: int, l2: int, lam2: float,
                           mp_dps: int = 50) -> np.ndarray:
    """Build dipole matrix R_{ij} = <chi_i^{l1, lam1} | r | chi_j^{l2, lam2}>.
    Built in mpmath then converted to float.
    """
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


def diagonalize_h0_sturmian(N: int, l: int, lam: float, Z: int = 1,
                            sv_cutoff: float = 1e-10):
    """Diagonalize H_0 = T - Z/r in the N-dim Sturmian basis at (l, lam).

    Uses canonical orthogonalization (SVD of S) to handle near-singular
    overlap matrix at large N.

    Returns:
        eigvals: array of M eigenvalues (M <= N after dropping null modes)
        eigvecs: array shape (N, M) of eigenvectors in the Sturmian basis
        mats: dict with T, V, S, transformation X (with eigvecs^T S eigvecs = I)
    """
    mats = sturmian_basis_matrices(N, l, lam, Z)
    H = mats["T"] + mats["V"]
    S = mats["S"]
    # Canonical orthogonalization: S = U_S diag(s) U_S^T
    s_eigs, U_S = sla.eigh(S)
    # Drop modes with eigenvalue below cutoff (relative to max)
    keep = s_eigs > sv_cutoff * s_eigs.max()
    s_kept = s_eigs[keep]
    U_kept = U_S[:, keep]
    X = U_kept * (1.0 / np.sqrt(s_kept))  # transformation: X^T S X = I
    H_orth = X.T @ H @ X
    # Symmetrize for numerical safety
    H_orth = 0.5 * (H_orth + H_orth.T)
    eigvals, V_orth = sla.eigh(H_orth)
    eigvecs = X @ V_orth  # back to Sturmian basis; eigvecs^T S eigvecs = I
    mats["X"] = X
    mats["s_kept"] = s_kept
    return eigvals, eigvecs, mats


def bethe_log_sturmian(n_target: int, l_target: int, N: int = 20, lam=None,
                       Z: int = 1, verbose: bool = False) -> dict:
    """Compute ln k_0(n_target, l_target) via finite Sturmian pseudostate basis.

    Strategy:
    1. Diagonalize H_0 in the (l = l_target ± 1) channels with Sturmian
       basis at exponent lam (default Z/n_target → Sturmian basis exact at
       n_target).
    2. The target state |n_target, l_target⟩ is also diagonalized in the
       l_target channel; identify it by eigenvalue ≈ -Z²/(2 n_target²).
    3. Compute dipole between target eigenvector and all pseudostates in
       l_target ± 1 channels.
    4. Form Bethe sum I, J using pseudo-eigenvalues.

    Returns dict with I, J, ln_k0, basis info, eigenvalue spectrum.
    """
    if lam is None:
        lam = Z / n_target
    Ry = 0.5
    E_target = -0.5 * Z ** 2 / n_target ** 2

    # Channel for target state
    eigvals_t, eigvecs_t, mats_t = diagonalize_h0_sturmian(N, l_target, lam, Z)
    # Find target eigenvector: eigenvalue closest to E_target
    idx_target = int(np.argmin(np.abs(eigvals_t - E_target)))
    E_t_basis = eigvals_t[idx_target]
    psi_target = eigvecs_t[:, idx_target]
    if verbose:
        print(f"  Target ({n_target}, {l_target}): E_basis = {E_t_basis:.6f}, "
              f"E_exact = {E_target:.6f}, residual = {abs(E_t_basis - E_target):.2e}")

    # Get pseudostates in l_target ± 1 channels
    channels = []
    if l_target >= 1:
        channels.append(l_target - 1)
    channels.append(l_target + 1)

    I_total = 0.0
    J_total = 0.0
    pseudo_records = []

    for lp in channels:
        eigvals_p, eigvecs_p, mats_p = diagonalize_h0_sturmian(N, lp, lam, Z)
        # Dipole matrix between Sturmian basis at l_target and Sturmian basis at lp
        R_basis = sturmian_dipole_matrix(N, l_target, lam, N, lp, lam)
        # Convert to pseudostate basis: R_pseudo[α, β] = ψ_α^T R_basis ψ_β
        # For generalized eigenvectors normalized as v^T S v = I, the operator
        # in pseudostate basis is V^T R_basis W where V, W are eigenvector
        # matrices for each channel.
        # Care: in sla.eigh(H, S) output, eigvecs satisfies eigvecs^T S eigvecs = I,
        # but the OPERATOR transformation is just eigvecs_t^T R_basis eigvecs_p
        # IF R_basis was computed in the same dr measure (it was).
        R_pseudo = psi_target @ R_basis @ eigvecs_p  # length M_p

        # Angular factor for m-summed dipole squared
        l_gt = max(l_target, lp)
        # Note: in pseudostate basis, |R_pseudo[α]|² is the m-AVERAGED radial
        # integral squared. To get |⟨nl|r|n'l'⟩|²_msum we multiply by l_>:
        M_p = len(eigvals_p)
        for alpha in range(M_p):
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
            pseudo_records.append({
                "lp": lp, "alpha": alpha, "E_alpha": E_alpha,
                "dE": dE, "d2_msum": d2_msum,
                "I_contrib": I_contrib, "J_contrib": J_contrib,
            })

    ln_k0 = J_total / I_total if abs(I_total) > 1e-30 else float('nan')

    return {
        "n_target": n_target, "l_target": l_target,
        "N_basis": N, "lam": lam,
        "E_target_exact": E_target,
        "E_target_basis": E_t_basis,
        "I": I_total, "J": J_total, "ln_k0": ln_k0,
        "n_pseudostates_used": len(pseudo_records),
        "channels": channels,
    }


# ---------------------------------------------------------------------------
# Main convergence study
# ---------------------------------------------------------------------------

def convergence_route_A(n: int, l: int,
                        n_max_list=(5, 10, 20, 30)) -> dict:
    """Route A: bound-only Bethe log convergence."""
    target = DRAKE_REF.get((n, l))
    print(f"\n[Route A] ({n},{l}) target = {target}")
    results = []
    for n_max in n_max_list:
        if n_max < n + 1:
            continue
        res = bethe_sums_bound(n, l, n_max)
        results.append({
            "n_max": n_max,
            "I": res["I"], "J": res["J"], "ln_k0": res["ln_k0"],
            "n_states": res["n_states_summed"],
        })
        err = (res["ln_k0"] - target) if target else float('nan')
        print(f"  n_max={n_max}: I={res['I']:.4e}, ln_k0={res['ln_k0']:.4f} "
              f"(err={err:+.4f})")
    return {"n": n, "l": l, "target": target, "convergence": results}


def convergence_route_B(n: int, l: int, N_list=(8, 12, 16, 20),
                        lam=None) -> dict:
    """Route B: Sturmian pseudostate Bethe log convergence."""
    target = DRAKE_REF.get((n, l))
    if lam is None:
        lam = 1.0 / n  # Z=1, target n
    print(f"\n[Route B] ({n},{l}) target = {target}, lam = {lam}")
    results = []
    for N in N_list:
        try:
            res = bethe_log_sturmian(n, l, N=N, lam=lam, verbose=True)
            err = (res["ln_k0"] - target) if target else float('nan')
            results.append({
                "N": N, "lam": lam,
                "E_target_basis": res["E_target_basis"],
                "I": res["I"], "J": res["J"], "ln_k0": res["ln_k0"],
                "n_pseudostates": res["n_pseudostates_used"],
            })
            print(f"  N={N}: ln_k0={res['ln_k0']:.4f} (err={err:+.4f})")
        except Exception as e:
            print(f"  N={N} FAILED: {e}")
    return {"n": n, "l": l, "target": target, "lam": lam,
            "convergence": results}


def lamb_shift_with_geovac_bethe(ln_k0_2s: float, ln_k0_2p: float) -> dict:
    """Compute Lamb shift using GeoVac-derived Bethe logs."""
    ALPHA = 1.0 / 137.035999084
    HA_TO_MHZ = 6_579_683_920.502
    Z = 1
    common = ALPHA ** 3 * Z ** 4 / (math.pi * 8)
    bracket_2S = (4.0 / 3.0) * (math.log(1.0 / (Z * ALPHA) ** 2)
                                - ln_k0_2s) + 38.0 / 45.0
    SE_2S = common * bracket_2S
    if not math.isnan(ln_k0_2p):
        bracket_2P = (4.0 / 3.0) * (-ln_k0_2p) - 1.0 / 6.0
        SE_2P = common * bracket_2P
    else:
        SE_2P = float('nan')
    VP_2S = -4.0 * ALPHA ** 3 * Z ** 4 / (15.0 * math.pi * 8)
    lamb = (SE_2S + VP_2S - SE_2P) * HA_TO_MHZ
    return {
        "ln_k0_2S_used": ln_k0_2s,
        "ln_k0_2P_used": ln_k0_2p,
        "SE_2S_MHz": SE_2S * HA_TO_MHZ,
        "SE_2P_MHz": SE_2P * HA_TO_MHZ if not math.isnan(SE_2P) else None,
        "VP_2S_MHz": VP_2S * HA_TO_MHZ,
        "lamb_MHz": lamb,
        "lamb_exp_MHz": 1057.845,
        "error_MHz": lamb - 1057.845,
    }


def main() -> dict:
    print("=" * 72)
    print("Sprint LS-2: Bethe logarithm from GeoVac spectral machinery")
    print("=" * 72)

    out = {"drake_swainson_targets": {f"{k[0]}_{k[1]}": v
                                      for k, v in DRAKE_REF.items()}}

    # Route A — bound only (sympy radial integrals, slow but exact)
    print("\n========= ROUTE A: bound-only sums =========")
    out["route_A"] = {}
    for nl in [(2, 0), (2, 1), (1, 0)]:
        out["route_A"][f"{nl[0]}_{nl[1]}"] = convergence_route_A(
            *nl, n_max_list=(5, 10, 15, 20))

    # Closure sanity
    print("\n--- Closure sanity check (velocity form) ---")
    out["closure_check"] = []
    for n, l in [(1, 0), (2, 0), (2, 1), (3, 0)]:
        clos = closure_I_velocity(n, l)
        bound = bethe_sums_bound(n, l, 30)
        ratio = bound["I"] / clos if clos > 0 else float('nan')
        print(f"  ({n},{l}): closure I = {clos:.6e}, "
              f"bound n_max=30 I = {bound['I']:.6e}, "
              f"frac of closure = {ratio:.4f}")
        out["closure_check"].append({
            "n": n, "l": l, "closure_I": clos, "bound_I": bound["I"],
            "fraction": ratio,
        })

    # Route B — Sturmian pseudostates
    print("\n========= ROUTE B: Sturmian pseudostate basis =========")
    out["route_B"] = {}
    # 2S and 1S — l=0 cases, well-conditioned
    for nl in [(1, 0), (2, 0)]:
        out["route_B"][f"{nl[0]}_{nl[1]}"] = convergence_route_B(
            *nl, N_list=(8, 12, 16, 20, 24, 30, 40, 50))
    # 2P — l>0, document the I→0 closure-cancellation pathology
    out["route_B"]["2_1"] = convergence_route_B(
        2, 1, N_list=(8, 12, 16, 20, 24))

    # Lamb shift downstream
    print("\n========= Lamb shift with GeoVac-derived Bethe logs =========")
    # Use the largest-N Route B values for 2S; for 2P use Drake-Swainson
    # since Route B exhibits I→0 closure pathology.
    lnk0_2s = out["route_B"]["2_0"]["convergence"][-1]["ln_k0"]
    lnk0_2p_drake = DRAKE_REF[(2, 1)]
    lamb_geovac = lamb_shift_with_geovac_bethe(lnk0_2s, lnk0_2p_drake)
    out["lamb_shift"] = lamb_geovac
    print(f"  Using ln_k0(2S) = {lnk0_2s:.4f} [GeoVac route B], "
          f"ln_k0(2P) = {lnk0_2p_drake:.4f} [Drake-Swainson]")
    print(f"  Lamb shift = {lamb_geovac['lamb_MHz']:.2f} MHz")
    print(f"  vs. experiment {lamb_geovac['lamb_exp_MHz']:.2f} MHz "
          f"(err {lamb_geovac['error_MHz']:+.2f} MHz)")

    # Save
    out_path = Path(__file__).parent / "data" / "ls2_bethe_log.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    serializable = json.loads(json.dumps(out, default=str))
    with open(out_path, "w") as f:
        json.dump(serializable, f, indent=2)
    print(f"\nSaved data to: {out_path}")
    return out


if __name__ == "__main__":
    main()
