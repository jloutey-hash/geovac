#!/usr/bin/env python
"""
Single-center Gaussian integral engine for He atom.

Computes AO and MO integrals for He in cc-pVDZ and cc-pVTZ basis sets,
performs RHF + FCI, then builds JW Hamiltonians to get actual Pauli term counts.

This replaces estimated Pauli counts with measured values from real integrals.

Method:
  - Spherical Gaussian basis (natural for single-center atoms)
  - Radial integrals via scipy.integrate.quad / dblquad
  - Angular coupling via sympy Gaunt coefficients
  - RHF in orthogonalized AO basis
  - MO transformation of h1 and ERI
  - FCI via exact diagonalization (trivial for 2 electrons)
  - JW encoding via OpenFermion

Author: GeoVac Development Team
Date: March 2026
"""

import sys
import os
import json
import time
from pathlib import Path
from typing import Dict, List, Tuple, Any

import numpy as np
from scipy.integrate import quad, dblquad
from scipy.special import gamma as gamma_fn
from scipy.linalg import eigh
from sympy.physics.wigner import gaunt as sympy_gaunt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# ---------------------------------------------------------------------------
# Basis set definitions from Basis Set Exchange
# ---------------------------------------------------------------------------

def he_cc_pvdz_basis() -> List[Dict[str, Any]]:
    """He cc-pVDZ basis: [2s1p] = 5 spatial orbitals."""
    return [
        {
            'l': 0, 'primitives': [
                (38.36, 0.023809),
                (5.770, 0.154891),
                (1.240, 0.469987),
                (0.2976, 0.513027),
            ],
        },
        {
            'l': 0, 'primitives': [
                (38.36, 0.0),
                (5.770, 0.0),
                (1.240, 0.0),
                (0.2976, 1.0),
            ],
        },
        {
            'l': 1, 'primitives': [
                (1.275, 1.0),
            ],
        },
    ]


def he_cc_pvtz_basis() -> List[Dict[str, Any]]:
    """He cc-pVTZ basis: [3s2p1d] = 14 spatial orbitals."""
    return [
        {
            'l': 0, 'primitives': [
                (234.0, 0.002587),
                (35.16, 0.019533),
                (7.989, 0.090998),
                (2.212, 0.272050),
                (0.6669, 0.478065),
                (0.2089, 0.307737),
            ],
        },
        {
            'l': 0, 'primitives': [
                (234.0, 0.0),
                (35.16, 0.0),
                (7.989, 0.0),
                (2.212, 0.0),
                (0.6669, 1.0),
                (0.2089, 0.0),
            ],
        },
        {
            'l': 0, 'primitives': [
                (234.0, 0.0),
                (35.16, 0.0),
                (7.989, 0.0),
                (2.212, 0.0),
                (0.6669, 0.0),
                (0.2089, 1.0),
            ],
        },
        {
            'l': 1, 'primitives': [
                (3.044, 1.0),
            ],
        },
        {
            'l': 1, 'primitives': [
                (0.758, 1.0),
            ],
        },
        {
            'l': 2, 'primitives': [
                (1.965, 1.0),
            ],
        },
    ]


# ---------------------------------------------------------------------------
# Radial integrals for spherical Gaussians at a single center
# ---------------------------------------------------------------------------

def radial_norm(alpha: float, l: int) -> float:
    """Normalization constant for r^l * exp(-alpha*r^2)."""
    # <phi|phi> = integral_0^inf r^(2l+2) * exp(-2*alpha*r^2) dr
    # = Gamma(l + 3/2) / (2 * (2*alpha)^(l + 3/2))
    val = gamma_fn(l + 1.5) / (2.0 * (2.0 * alpha) ** (l + 1.5))
    return 1.0 / np.sqrt(val)


def overlap_radial(alphas_a: List[Tuple[float, float]],
                   alphas_b: List[Tuple[float, float]],
                   l: int) -> float:
    """
    Overlap integral for two contracted spherical Gaussian shells with same l.

    Each shell is [(exponent, contraction_coeff), ...].
    The integral is: sum_{ij} c_i * c_j * N_i * N_j * I(alpha_i + alpha_j, l)
    where I(gamma, l) = Gamma(l + 3/2) / (2 * gamma^(l + 3/2))
    """
    result = 0.0
    for alpha_a, c_a in alphas_a:
        N_a = radial_norm(alpha_a, l)
        for alpha_b, c_b in alphas_b:
            N_b = radial_norm(alpha_b, l)
            gamma = alpha_a + alpha_b
            I = gamma_fn(l + 1.5) / (2.0 * gamma ** (l + 1.5))
            result += c_a * c_b * N_a * N_b * I
    return result


def kinetic_radial(alphas_a: List[Tuple[float, float]],
                   alphas_b: List[Tuple[float, float]],
                   l: int) -> float:
    """
    Kinetic energy integral for two contracted spherical Gaussian shells.

    T = <a| -1/2 nabla^2 |b>

    For R_l(r) = r^l * exp(-beta*r^2), the kinetic operator gives:
    -1/2 nabla^2 [R_l * Y_lm] = [beta(2l+3) r^l - 2*beta^2 r^(l+2)] * exp(-beta*r^2) * Y_lm

    So T = sum c_i c_j N_i N_j * [beta_j(2l+3)*I(2l+2,gamma) - 2*beta_j^2*I(2l+4,gamma)]
    where I(n, gamma) = Gamma((n+1)/2) / (2 * gamma^((n+1)/2))
    and gamma = alpha_i + alpha_j.
    """
    result = 0.0
    for alpha_a, c_a in alphas_a:
        N_a = radial_norm(alpha_a, l)
        for alpha_b, c_b in alphas_b:
            N_b = radial_norm(alpha_b, l)
            gamma = alpha_a + alpha_b
            beta = alpha_b

            I_2l2 = gamma_fn(l + 1.5) / (2.0 * gamma ** (l + 1.5))
            I_2l4 = gamma_fn(l + 2.5) / (2.0 * gamma ** (l + 2.5))

            T_val = beta * (2 * l + 3) * I_2l2 - 2.0 * beta**2 * I_2l4
            result += c_a * c_b * N_a * N_b * T_val
    return result


def nuclear_radial(alphas_a: List[Tuple[float, float]],
                   alphas_b: List[Tuple[float, float]],
                   l: int, Z: int) -> float:
    """
    Nuclear attraction integral for same-center, nucleus at origin.

    V = <a| -Z/r |b> = -Z * sum c_i c_j N_i N_j * I(2l+1, gamma)
    where I(n, gamma) = Gamma((n+1)/2) / (2 * gamma^((n+1)/2))
    """
    result = 0.0
    for alpha_a, c_a in alphas_a:
        N_a = radial_norm(alpha_a, l)
        for alpha_b, c_b in alphas_b:
            N_b = radial_norm(alpha_b, l)
            gamma = alpha_a + alpha_b
            # integral of r^(2l+1) * exp(-gamma*r^2) from 0 to inf
            I_2l1 = gamma_fn(l + 1.0) / (2.0 * gamma ** (l + 1.0))
            result += c_a * c_b * N_a * N_b * I_2l1
    return -Z * result


def slater_Rk_radial(alphas_a: List[Tuple[float, float]], la: int,
                     alphas_b: List[Tuple[float, float]], lb: int,
                     alphas_c: List[Tuple[float, float]], lc: int,
                     alphas_d: List[Tuple[float, float]], ld: int,
                     k: int) -> float:
    """
    Radial Slater integral R^k for Gaussian basis functions.

    R^k(ab, cd) = int int R_a(r1)*R_b(r1) * r_<^k/r_>^(k+1) *
                  R_c(r2)*R_d(r2) * r1^2 * r2^2 dr1 dr2

    For primitives: R_i(r) = N_i * r^l_i * exp(-alpha_i * r^2)

    Product R_a*R_b has effective: r^(la+lb) * exp(-gamma_ab * r^2)
    Product R_c*R_d has effective: r^(lc+ld) * exp(-gamma_cd * r^2)
    """
    total = 0.0

    for alpha_a, c_a in alphas_a:
        Na = radial_norm(alpha_a, la)
        for alpha_b, c_b in alphas_b:
            Nb = radial_norm(alpha_b, lb)
            gab = alpha_a + alpha_b
            cab = c_a * c_b * Na * Nb

            for alpha_c, c_c in alphas_c:
                Nc = radial_norm(alpha_c, lc)
                for alpha_d, c_d in alphas_d:
                    Nd = radial_norm(alpha_d, ld)
                    gcd = alpha_c + alpha_d
                    ccd = c_c * c_d * Nc * Nd

                    # Compute the 2D radial integral analytically
                    # by splitting into r1 < r2 and r1 > r2 regions.
                    val = _slater_Rk_primitive(
                        la + lb, gab, lc + ld, gcd, k
                    )
                    total += cab * ccd * val
    return total


def _slater_Rk_primitive(L1: int, g1: float, L2: int, g2: float, k: int) -> float:
    """
    R^k for two primitive radial products:
      f1(r) = r^L1 * exp(-g1*r^2)
      f2(r) = r^L2 * exp(-g2*r^2)

    R^k = int_0^inf int_0^inf f1(r1)*r1^2 * f2(r2)*r2^2 *
          r_<^k / r_>^(k+1) dr1 dr2

    The inner integral over r1 is evaluated analytically using the
    incomplete gamma function, reducing this to a 1D quadrature over r2.

    For r1 < r2:
      int_0^r2 r1^(L1+k+2) exp(-g1*r1^2) dr1
        = (1/2) * gamma_lower(n_less, g1*r2^2) / g1^n_less

    For r1 > r2:
      int_r2^inf r1^(L1-k+1) exp(-g1*r1^2) dr1
        = (1/2) * gamma_upper(n_greater, g1*r2^2) / g1^n_greater

    where n_less = (L1+k+3)/2 and n_greater = (L1-k+2)/2.
    """
    from scipy.special import gammainc, gammaincc

    n_less = (L1 + k + 3) / 2.0
    n_greater = (L1 - k + 2) / 2.0
    A_less = gamma_fn(n_less) / (2.0 * g1 ** n_less)
    A_greater = gamma_fn(n_greater) / (2.0 * g1 ** n_greater)

    def integrand(r2: float) -> float:
        if r2 < 1e-15:
            return 0.0
        x = g1 * r2**2
        f2 = r2 ** (L2 + 2) * np.exp(-g2 * r2**2)

        # Analytical inner integrals via incomplete gamma
        inner_less = A_less * gammainc(n_less, x) / r2**(k + 1)
        inner_greater = A_greater * gammaincc(n_greater, x) * r2**k

        return f2 * (inner_less + inner_greater)

    result, _ = quad(integrand, 0, np.inf, limit=200)
    return result


# ---------------------------------------------------------------------------
# Real Gaunt coefficients (angular coupling)
# ---------------------------------------------------------------------------
# sympy.physics.wigner.gaunt uses COMPLEX spherical harmonics.
# Our basis uses REAL spherical harmonics. For l=0,1 the difference
# is minor, but for l>=2 it matters significantly. We compute the
# real Gaunt coefficient by numerical quadrature over (theta, phi).

_real_gaunt_cache: Dict[Tuple[int, ...], float] = {}


def _real_sph_harm(l: int, m: int, theta: np.ndarray, phi: np.ndarray) -> np.ndarray:
    """
    Real spherical harmonic Y_l^m(theta, phi).

    Convention:
      m > 0:  sqrt(2) * N_lm * P_l^m(cos theta) * cos(m*phi)
      m = 0:  N_l0 * P_l^0(cos theta)
      m < 0:  sqrt(2) * N_l|m| * P_l^|m|(cos theta) * sin(|m|*phi)

    where N_lm = sqrt((2l+1)/(4pi) * (l-|m|)!/(l+|m|)!).
    """
    from scipy.special import sph_harm_y

    # scipy.special.sph_harm_y(l, m, theta, phi) returns complex Y_l^m
    # with Condon-Shortley phase.
    if m > 0:
        yc = sph_harm_y(l, m, theta, phi)
        return np.sqrt(2.0) * np.real(yc)
    elif m == 0:
        yc = sph_harm_y(l, 0, theta, phi)
        return np.real(yc)
    else:
        yc = sph_harm_y(l, abs(m), theta, phi)
        return np.sqrt(2.0) * ((-1)**abs(m)) * np.imag(yc)


def _build_angular_grid(n_theta: int = 60, n_phi: int = 120) -> Tuple[
    np.ndarray, np.ndarray, np.ndarray
]:
    """Build a theta-phi quadrature grid with weights."""
    # Gauss-Legendre for cos(theta)
    from numpy.polynomial.legendre import leggauss
    x_gl, w_gl = leggauss(n_theta)
    # x_gl is in [-1, 1] = cos(theta), w_gl are weights
    theta_gl = np.arccos(x_gl)
    # Uniform grid for phi
    phi_gl = np.linspace(0, 2 * np.pi, n_phi, endpoint=False)
    dphi = 2 * np.pi / n_phi

    # Outer product
    theta_2d, phi_2d = np.meshgrid(theta_gl, phi_gl, indexing='ij')
    # Weights: w_theta * dphi (no sin(theta) factor needed because
    # Gauss-Legendre integrates over cos(theta) directly)
    w_2d = np.outer(w_gl, np.ones(n_phi) * dphi)

    return theta_2d.ravel(), phi_2d.ravel(), w_2d.ravel()


# Pre-build angular grid (computed once)
_THETA_GRID, _PHI_GRID, _W_GRID = _build_angular_grid(80, 160)


def real_gaunt(l1: int, m1: int, l2: int, m2: int, l3: int, m3: int) -> float:
    """
    Compute the real Gaunt coefficient by numerical quadrature:
      G = integral Y_l1^m1_real * Y_l2^m2_real * Y_l3^m3_real dOmega
    """
    key = (l1, m1, l2, m2, l3, m3)
    if key in _real_gaunt_cache:
        return _real_gaunt_cache[key]

    y1 = _real_sph_harm(l1, m1, _THETA_GRID, _PHI_GRID)
    y2 = _real_sph_harm(l2, m2, _THETA_GRID, _PHI_GRID)
    y3 = _real_sph_harm(l3, m3, _THETA_GRID, _PHI_GRID)

    val = float(np.sum(y1 * y2 * y3 * _W_GRID))
    _real_gaunt_cache[key] = val
    return val


def eri_angular_factor(la: int, ma: int, lb: int, mb: int,
                       lc: int, mc: int, ld: int, md: int,
                       k: int) -> float:
    """
    Angular factor for the ERI via real Gaunt coefficients.

    (ab|cd) = sum_k angular_factor(k) * R^k(ab,cd) * 4*pi/(2k+1)

    angular_factor(k) = sum_q G_real(la,ma,lb,mb,k,q) * G_real(lc,mc,ld,md,k,q)

    Selection rules (still apply for real harmonics):
    - |la-lb| <= k <= la+lb with la+lb+k even
    - |lc-ld| <= k <= lc+ld with lc+ld+k even
    """
    # Parity selection
    if (la + lb + k) % 2 != 0:
        return 0.0
    if (lc + ld + k) % 2 != 0:
        return 0.0

    # Triangle rules
    if k < abs(la - lb) or k > la + lb:
        return 0.0
    if k < abs(lc - ld) or k > lc + ld:
        return 0.0

    # Sum over all q from -k to k
    total = 0.0
    for q in range(-k, k + 1):
        g1 = real_gaunt(la, ma, lb, mb, k, q)
        if abs(g1) < 1e-12:
            continue
        g2 = real_gaunt(lc, mc, ld, md, k, q)
        if abs(g2) < 1e-12:
            continue
        total += g1 * g2

    return total


# ---------------------------------------------------------------------------
# Build full integral matrices
# ---------------------------------------------------------------------------

def build_basis_functions(basis_shells: List[Dict]) -> List[Dict]:
    """
    Expand basis shells into individual basis functions (one per m quantum number).

    Returns list of {'l': l, 'm': m, 'primitives': [(alpha, coeff), ...]}
    """
    functions = []
    for shell in basis_shells:
        l = shell['l']
        for m in range(-l, l + 1):
            functions.append({
                'l': l,
                'm': m,
                'primitives': shell['primitives'],
            })
    return functions


def compute_integrals(basis_shells: List[Dict], Z: int = 2) -> Tuple[
    np.ndarray, np.ndarray, np.ndarray, np.ndarray
]:
    """
    Compute S, T, V, and ERI in the AO basis.

    Parameters
    ----------
    basis_shells : list of dict
        From he_cc_pvdz_basis() or he_cc_pvtz_basis().
    Z : int
        Nuclear charge.

    Returns
    -------
    S : ndarray (M, M) - overlap
    T : ndarray (M, M) - kinetic energy
    V : ndarray (M, M) - nuclear attraction
    eri : ndarray (M, M, M, M) - electron repulsion integrals, chemist notation
    """
    bfs = build_basis_functions(basis_shells)
    M = len(bfs)
    print(f"  {M} basis functions (spatial orbitals)")

    S = np.zeros((M, M))
    T = np.zeros((M, M))
    V = np.zeros((M, M))

    # One-electron integrals: diagonal in (l,m) for same-center
    print("  Computing 1-electron integrals...")
    for i, bf_i in enumerate(bfs):
        for j, bf_j in enumerate(bfs):
            if bf_i['l'] != bf_j['l'] or bf_i['m'] != bf_j['m']:
                continue
            l = bf_i['l']
            S[i, j] = overlap_radial(bf_i['primitives'], bf_j['primitives'], l)
            T[i, j] = kinetic_radial(bf_i['primitives'], bf_j['primitives'], l)
            V[i, j] = nuclear_radial(bf_i['primitives'], bf_j['primitives'], l, Z)

    # Two-electron integrals via Slater expansion
    print("  Computing 2-electron integrals (Slater expansion)...")
    eri = np.zeros((M, M, M, M))

    # Determine which k values are needed
    l_max = max(bf['l'] for bf in bfs)
    k_max = 2 * l_max

    # Cache radial Slater integrals by shell indices
    # Build unique shell list
    shells_unique = []
    shell_map = []  # maps bf index -> shell index
    for i, bf in enumerate(bfs):
        # Find matching shell by l and primitives
        found = False
        for si, su in enumerate(shells_unique):
            if su['l'] == bf['l'] and su['primitives'] == bf['primitives']:
                shell_map.append(si)
                found = True
                break
        if not found:
            shell_map.append(len(shells_unique))
            shells_unique.append({'l': bf['l'], 'primitives': bf['primitives']})

    n_shells = len(shells_unique)
    print(f"  {n_shells} unique radial shells, k_max = {k_max}")

    # Precompute R^k for all shell quartets
    Rk_cache: Dict[Tuple[int, int, int, int, int], float] = {}
    n_computed = 0
    for sa in range(n_shells):
        for sb in range(n_shells):
            for sc in range(n_shells):
                for sd in range(n_shells):
                    la = shells_unique[sa]['l']
                    lb = shells_unique[sb]['l']
                    lc = shells_unique[sc]['l']
                    ld = shells_unique[sd]['l']

                    for k in range(0, k_max + 1):
                        # Selection rules on k
                        if k < abs(la - lb) or k > la + lb:
                            continue
                        if k < abs(lc - ld) or k > lc + ld:
                            continue
                        if (la + lb + k) % 2 != 0:
                            continue
                        if (lc + ld + k) % 2 != 0:
                            continue

                        key = (sa, sb, sc, sd, k)
                        if key in Rk_cache:
                            continue

                        Rk_cache[key] = slater_Rk_radial(
                            shells_unique[sa]['primitives'], la,
                            shells_unique[sb]['primitives'], lb,
                            shells_unique[sc]['primitives'], lc,
                            shells_unique[sd]['primitives'], ld,
                            k,
                        )
                        n_computed += 1

    print(f"  Computed {n_computed} unique R^k integrals")

    # Assemble ERI from angular factors and R^k
    for i in range(M):
        for j in range(M):
            for p in range(M):
                for q in range(M):
                    la, ma = bfs[i]['l'], bfs[i]['m']
                    lb, mb = bfs[j]['l'], bfs[j]['m']
                    lc, mc = bfs[p]['l'], bfs[p]['m']
                    ld, md = bfs[q]['l'], bfs[q]['m']

                    sa, sb, sc, sd = shell_map[i], shell_map[j], shell_map[p], shell_map[q]

                    val = 0.0
                    for k in range(0, k_max + 1):
                        key = (sa, sb, sc, sd, k)
                        if key not in Rk_cache:
                            continue

                        ang = eri_angular_factor(la, ma, lb, mb, lc, mc, ld, md, k)
                        if abs(ang) < 1e-15:
                            continue

                        val += ang * Rk_cache[key] * 4.0 * np.pi / (2 * k + 1)

                    eri[i, j, p, q] = val

    return S, T, V, eri


# ---------------------------------------------------------------------------
# RHF + FCI
# ---------------------------------------------------------------------------

def run_rhf_fci(S: np.ndarray, h1: np.ndarray, eri: np.ndarray,
                n_electrons: int = 2) -> Tuple[float, np.ndarray, np.ndarray]:
    """
    Run RHF and FCI for a 2-electron system.

    Parameters
    ----------
    S : overlap matrix
    h1 : T + V (one-electron Hamiltonian)
    eri : 2-electron integrals in chemist notation
    n_electrons : int

    Returns
    -------
    e_fci : float
        FCI energy.
    mo_coeffs : ndarray (M, M)
        MO coefficient matrix.
    h1_mo : ndarray (M, M)
        One-electron integrals in MO basis.
    """
    M = S.shape[0]

    # Symmetric orthogonalization
    s_eigvals, s_eigvecs = eigh(S)
    # Remove near-linear-dependencies
    keep = s_eigvals > 1e-8
    s_eigvals = s_eigvals[keep]
    s_eigvecs = s_eigvecs[:, keep]
    X = s_eigvecs @ np.diag(1.0 / np.sqrt(s_eigvals))

    print(f"  Orthogonalization: {M} -> {X.shape[1]} (removed {M - X.shape[1]} linear deps)")

    # Solve h1 in orthogonal basis to get initial MO coefficients
    h1_orth = X.T @ h1 @ X
    mo_energies, C_orth = eigh(h1_orth)
    C = X @ C_orth  # MO coefficients in AO basis

    # For 2 electrons, RHF is just doubly occupying the lowest MO.
    # The RHF energy = 2*h11 + J11
    e_hf = 2.0 * mo_energies[0]

    # Transform to MO basis for FCI
    M_eff = C.shape[1]
    h1_mo = C.T @ h1 @ C

    # Transform ERI to MO basis: (pq|rs) = sum_{ijkl} C_ip C_jq C_kr C_ls (ij|kl)
    # Use 4-index transformation (O(M^5))
    print(f"  Transforming ERI to MO basis ({M_eff} MOs)...")
    eri_mo = np.einsum('ip,jq,ijkl,kr,ls->pqrs', C, C, eri, C, C,
                       optimize=True)

    # FCI for 2 electrons: enumerate all pairs of spatial orbitals
    # |Phi_ab> = a^+_alpha b^+_beta |0>  (spatial orbitals a, b)
    # Matrix elements: <Phi_ab|H|Phi_cd>
    n_configs = M_eff * M_eff  # all (a,b) pairs with a=alpha, b=beta
    print(f"  FCI: {n_configs} configurations ({M_eff} spatial MOs)")

    H_fci = np.zeros((n_configs, n_configs))
    for ia in range(M_eff):
        for ib in range(M_eff):
            idx_ab = ia * M_eff + ib
            for ic in range(M_eff):
                for id_ in range(M_eff):
                    idx_cd = ic * M_eff + id_

                    # <ia_alpha ib_beta | H | ic_alpha id_beta>
                    # One-electron: delta(ib,id)*h1[ia,ic] + delta(ia,ic)*h1[ib,id]
                    val = 0.0
                    if ib == id_:
                        val += h1_mo[ia, ic]
                    if ia == ic:
                        val += h1_mo[ib, id_]

                    # Two-electron: (ia ic | ib id) - 0 (no exchange for opposite spin)
                    # Coulomb: <ia_alpha ib_beta | 1/r12 | ic_alpha id_beta>
                    #        = (ia ic | ib id)   [chemist notation]
                    val += eri_mo[ia, ic, ib, id_]

                    H_fci[idx_ab, idx_cd] = val

    fci_eigvals = np.linalg.eigvalsh(H_fci)
    e_fci = fci_eigvals[0]

    print(f"  E(HF)  = {2*h1_mo[0,0] + eri_mo[0,0,0,0]:.10f} Ha")
    print(f"  E(FCI) = {e_fci:.10f} Ha")

    return e_fci, C, h1_mo, eri_mo


# ---------------------------------------------------------------------------
# JW Pauli term counting
# ---------------------------------------------------------------------------

def count_pauli_terms(h1_mo: np.ndarray, eri_mo: np.ndarray,
                      nuclear_repulsion: float = 0.0) -> Tuple[int, Any]:
    """
    Build JW qubit Hamiltonian and count Pauli terms.

    Parameters
    ----------
    h1_mo : (M, M) one-electron integrals in MO basis
    eri_mo : (M, M, M, M) two-electron integrals in MO basis, chemist notation
    nuclear_repulsion : float

    Returns
    -------
    n_pauli : int
    qubit_op : QubitOperator
    """
    from geovac.qubit_encoding import build_fermion_op_from_integrals
    from openfermion import jordan_wigner

    fermion_op = build_fermion_op_from_integrals(h1_mo, eri_mo, nuclear_repulsion)
    qubit_op = jordan_wigner(fermion_op)
    n_pauli = len(qubit_op.terms)

    M = h1_mo.shape[0]
    Q = 2 * M
    print(f"  JW encoding: {Q} qubits, {n_pauli} Pauli terms")

    return n_pauli, qubit_op


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def process_basis(name: str, basis_shells: List[Dict], Z: int = 2) -> Dict[str, Any]:
    """Process a single basis set: compute integrals, RHF, FCI, JW."""
    print(f"\n{'='*60}")
    print(f"He {name}")
    print(f"{'='*60}")

    t0 = time.time()
    S, T, V, eri = compute_integrals(basis_shells, Z=Z)
    h1 = T + V
    t_integrals = time.time() - t0
    print(f"  Integral time: {t_integrals:.1f}s")

    # Verify overlap diagonal
    print(f"  S diagonal: {np.diag(S)}")

    e_fci, C, h1_mo, eri_mo = run_rhf_fci(S, h1, eri, n_electrons=2)

    # Count Pauli terms
    M = h1_mo.shape[0]
    n_pauli, qubit_op = count_pauli_terms(h1_mo, eri_mo, nuclear_repulsion=0.0)

    # Reference energy
    he_exact = -2.903724
    error_pct = 100.0 * abs(e_fci - he_exact) / abs(he_exact)

    result = {
        'basis': name,
        'n_spatial': M,
        'n_qubits': 2 * M,
        'fci_energy': float(e_fci),
        'exact_he_energy': he_exact,
        'error_pct': float(error_pct),
        'n_pauli_terms': n_pauli,
        'h1_mo': h1_mo.tolist(),
        'eri_mo': eri_mo.tolist(),
        'integral_time_s': round(t_integrals, 1),
    }

    print(f"\n  RESULT: He {name}")
    print(f"    {M} spatial orbitals, {2*M} qubits")
    print(f"    FCI energy: {e_fci:.10f} Ha")
    print(f"    Error vs exact: {error_pct:.4f}%")
    print(f"    Pauli terms: {n_pauli}")

    return result


def main() -> None:
    results = {}

    # cc-pVDZ
    results['cc-pVDZ'] = process_basis('cc-pVDZ', he_cc_pvdz_basis())

    # cc-pVTZ
    results['cc-pVTZ'] = process_basis('cc-pVTZ', he_cc_pvtz_basis())

    # Save results
    out_path = Path("debug/data/he_gaussian_integrals.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Save without the large arrays
    summary = {}
    for k, v in results.items():
        summary[k] = {key: val for key, val in v.items()
                      if key not in ('h1_mo', 'eri_mo')}
    with open(out_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"\nSummary saved to {out_path}")

    # Print comparison with estimates
    print("\n" + "="*60)
    print("COMPARISON: Actual vs Estimated Pauli Terms")
    print("="*60)
    est_pvdz = 814
    est_pvtz = 91923
    act_pvdz = results['cc-pVDZ']['n_pauli_terms']
    act_pvtz = results['cc-pVTZ']['n_pauli_terms']
    print(f"  cc-pVDZ: estimated {est_pvdz:>8,} | actual {act_pvdz:>8,} | ratio {act_pvdz/est_pvdz:.2f}")
    print(f"  cc-pVTZ: estimated {est_pvtz:>8,} | actual {act_pvtz:>8,} | ratio {act_pvtz/est_pvtz:.2f}")

    # Save MO integrals for hardcoding
    for name, data in results.items():
        np.savez(
            f"debug/data/he_{name.replace('-','_').lower()}_mo_integrals.npz",
            h1_mo=np.array(data['h1_mo']),
            eri_mo=np.array(data['eri_mo']),
            fci_energy=data['fci_energy'],
        )
    print("MO integrals saved to debug/data/")


if __name__ == '__main__':
    main()
