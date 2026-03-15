"""
Numba-accelerated computational kernels for Hylleraas and prolate SCF.

All functions are @njit(cache=True) for nopython mode with disk caching.
Inputs are numpy arrays and scalars only — no Python objects.

The elliptic integral K(m) uses the Abramowitz & Stegun polynomial
approximation (17.3.34), accurate to |error| < 2e-8.

Author: GeoVac Development Team
Date: March 2026
"""

import numpy as np

try:
    from numba import njit
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False


# ======================================================================
# Guard: only define jitted functions if Numba is installed
# ======================================================================

if NUMBA_AVAILABLE:

    # ==================================================================
    # Elliptic integral K(m) — Abramowitz & Stegun 17.3.34
    # ==================================================================

    @njit(cache=True)
    def _ellipk_approx(m: float) -> float:
        """Complete elliptic integral K(m) via polynomial approximation.

        Abramowitz & Stegun 17.3.34, |error| < 2e-8.
        m is the parameter (not the modulus k; m = k^2).
        Valid for 0 <= m < 1.
        """
        if m < 0.0:
            m = 0.0
        if m >= 1.0:
            return 1e10  # K diverges at m=1

        m1 = 1.0 - m
        # Coefficients from A&S 17.3.34
        a0 = 1.3862943611198906
        a1 = 0.09666344259
        a2 = 0.03590092383
        a3 = 0.03742563713
        a4 = 0.01451196212

        b0 = 0.5
        b1 = 0.12498593597
        b2 = 0.06880248576
        b3 = 0.03328355346
        b4 = 0.00441787012

        ln_m1 = np.log(1.0 / max(m1, 1e-300))

        poly_a = a0 + m1 * (a1 + m1 * (a2 + m1 * (a3 + m1 * a4)))
        poly_b = b0 + m1 * (b1 + m1 * (b2 + m1 * (b3 + m1 * b4)))

        return poly_a + poly_b * ln_m1

    # ==================================================================
    # r12 computation (scalar)
    # ==================================================================

    @njit(cache=True)
    def _compute_r12_scalar(
        xi1: float, eta1: float, xi2: float, eta2: float,
        dphi: float, R: float,
    ) -> float:
        """Compute r12 including azimuthal angle difference (scalar)."""
        half_R = R / 2.0
        rho1_sq = max((xi1**2 - 1.0) * (1.0 - eta1**2), 0.0)
        rho2_sq = max((xi2**2 - 1.0) * (1.0 - eta2**2), 0.0)

        r12_sq = half_R**2 * (
            (xi1 * eta1 - xi2 * eta2)**2
            + rho1_sq + rho2_sq
            - 2.0 * np.cos(dphi) * np.sqrt(rho1_sq * rho2_sq)
        )
        return np.sqrt(max(r12_sq, 0.0))

    # ==================================================================
    # Basis function evaluation (scalar, unsymmetrized)
    # ==================================================================

    @njit(cache=True)
    def _eval_unsym(
        j: int, k: int, l: int, m: int, alpha: float,
        xi1: float, eta1: float, xi2: float, eta2: float,
    ) -> float:
        """Unsymmetrized basis: exp(-alpha*(xi1+xi2)) * xi1^j * xi2^k * eta1^l * eta2^m."""
        return (np.exp(-alpha * (xi1 + xi2))
                * xi1**j * xi2**k * eta1**l * eta2**m)

    @njit(cache=True)
    def _eval_sym(
        j: int, k: int, l: int, m: int, alpha: float, is_self: bool,
        xi1: float, eta1: float, xi2: float, eta2: float,
    ) -> float:
        """Symmetrized p=0 basis function."""
        direct = _eval_unsym(j, k, l, m, alpha, xi1, eta1, xi2, eta2)
        if is_self:
            return 2.0 * direct
        exchange = _eval_unsym(k, j, m, l, alpha, xi1, eta1, xi2, eta2)
        return direct + exchange

    @njit(cache=True)
    def _eval_sym_with_r12(
        j: int, k: int, l: int, m: int, p: int,
        alpha: float, is_self: bool,
        xi1: float, eta1: float, xi2: float, eta2: float,
        r12: float,
    ) -> float:
        """Symmetrized basis function with r12^p factor."""
        f = _eval_sym(j, k, l, m, alpha, is_self, xi1, eta1, xi2, eta2)
        if p > 0:
            f *= r12**p
        return f

    # ==================================================================
    # Basis function derivatives (scalar, for kinetic energy)
    # ==================================================================

    @njit(cache=True)
    def _eval_unsym_derivs(
        j: int, k: int, l: int, m: int, alpha: float,
        xi1: float, eta1: float, xi2: float, eta2: float,
    ) -> tuple:
        """Unsymmetrized basis and first derivatives (p=0).

        Returns (g, dg/dxi1, dg/deta1, dg/dxi2, dg/deta2).
        """
        exp_f = np.exp(-alpha * (xi1 + xi2))
        poly = xi1**j * xi2**k * eta1**l * eta2**m
        g = exp_f * poly

        # dg/dxi1
        if j > 0:
            dg_dxi1 = exp_f * xi2**k * eta1**l * eta2**m * (
                j * xi1**(j - 1) - alpha * xi1**j
            )
        else:
            dg_dxi1 = -alpha * g

        # dg/deta1
        if l > 0:
            dg_deta1 = exp_f * xi1**j * xi2**k * eta2**m * l * eta1**(l - 1)
        else:
            dg_deta1 = 0.0

        # dg/dxi2
        if k > 0:
            dg_dxi2 = exp_f * xi1**j * eta1**l * eta2**m * (
                k * xi2**(k - 1) - alpha * xi2**k
            )
        else:
            dg_dxi2 = -alpha * g

        # dg/deta2
        if m > 0:
            dg_deta2 = exp_f * xi1**j * xi2**k * eta1**l * m * eta2**(m - 1)
        else:
            dg_deta2 = 0.0

        return g, dg_dxi1, dg_deta1, dg_dxi2, dg_deta2

    @njit(cache=True)
    def _eval_sym_derivs(
        j: int, k: int, l: int, m: int, alpha: float, is_self: bool,
        xi1: float, eta1: float, xi2: float, eta2: float,
    ) -> tuple:
        """Symmetrized basis function and first derivatives (p=0).

        Returns (phi, dphi/dxi1, dphi/deta1, dphi/dxi2, dphi/deta2).
        """
        g_d, d1, d2, d3, d4 = _eval_unsym_derivs(
            j, k, l, m, alpha, xi1, eta1, xi2, eta2
        )
        if is_self:
            return 2.0 * g_d, 2.0 * d1, 2.0 * d2, 2.0 * d3, 2.0 * d4

        g_e, e1, e2, e3, e4 = _eval_unsym_derivs(
            k, j, m, l, alpha, xi1, eta1, xi2, eta2
        )
        return g_d + g_e, d1 + e1, d2 + e2, d3 + e3, d4 + e4

    # ==================================================================
    # Overlap matrix kernel (p=0): 4D integration
    # ==================================================================

    @njit(cache=True)
    def overlap_matrix_p0_kernel(
        bf_j: np.ndarray, bf_k: np.ndarray,
        bf_l: np.ndarray, bf_m: np.ndarray,
        bf_alpha: np.ndarray, bf_is_self: np.ndarray,
        xi: np.ndarray, eta: np.ndarray,
        w_xi: np.ndarray, w_eta: np.ndarray,
        R: float,
    ) -> np.ndarray:
        """Compute full overlap matrix for p=0 basis (4D quadrature)."""
        n_basis = len(bf_j)
        n_xi = len(xi)
        n_eta = len(eta)
        half_R = R / 2.0
        phi_factor = (2.0 * np.pi)**2
        S = np.zeros((n_basis, n_basis))

        for i in range(n_basis):
            ji, ki, li, mi = bf_j[i], bf_k[i], bf_l[i], bf_m[i]
            ai, si = bf_alpha[i], bf_is_self[i]
            for ij in range(i, n_basis):
                jj, kj, lj, mj = bf_j[ij], bf_k[ij], bf_l[ij], bf_m[ij]
                aj, sj = bf_alpha[ij], bf_is_self[ij]

                result = 0.0
                for a in range(n_xi):
                    xi1 = xi[a]
                    w_a = w_xi[a]
                    for c in range(n_xi):
                        xi2 = xi[c]
                        w_c = w_xi[c]
                        w_ac = w_a * w_c
                        for b in range(n_eta):
                            eta1 = eta[b]
                            J1 = half_R**3 * (xi1**2 - eta1**2)
                            w_b = w_eta[b]
                            for d in range(n_eta):
                                eta2 = eta[d]
                                J2 = half_R**3 * (xi2**2 - eta2**2)

                                fi = _eval_sym(
                                    ji, ki, li, mi, ai, si,
                                    xi1, eta1, xi2, eta2,
                                )
                                fj = _eval_sym(
                                    jj, kj, lj, mj, aj, sj,
                                    xi1, eta1, xi2, eta2,
                                )

                                result += (fi * fj * J1 * J2
                                           * w_ac * w_b * w_eta[d])

                S[i, ij] = result * phi_factor
                S[ij, i] = S[i, ij]

        return S

    # ==================================================================
    # Overlap matrix kernel (p>0): 5D integration
    # ==================================================================

    @njit(cache=True)
    def overlap_matrix_5d_kernel(
        bf_j: np.ndarray, bf_k: np.ndarray,
        bf_l: np.ndarray, bf_m: np.ndarray,
        bf_p: np.ndarray, bf_alpha: np.ndarray,
        bf_is_self: np.ndarray,
        xi: np.ndarray, eta: np.ndarray,
        w_xi: np.ndarray, w_eta: np.ndarray,
        dphi: np.ndarray, w_phi: np.ndarray,
        R: float,
    ) -> np.ndarray:
        """Compute full overlap matrix for general basis (5D quadrature)."""
        n_basis = len(bf_j)
        n_xi = len(xi)
        n_eta = len(eta)
        n_phi = len(dphi)
        half_R = R / 2.0
        S = np.zeros((n_basis, n_basis))

        for i in range(n_basis):
            ji, ki, li, mi, pi_ = bf_j[i], bf_k[i], bf_l[i], bf_m[i], bf_p[i]
            ai, si = bf_alpha[i], bf_is_self[i]
            for ij in range(i, n_basis):
                jj, kj, lj, mj, pj = (bf_j[ij], bf_k[ij], bf_l[ij],
                                       bf_m[ij], bf_p[ij])
                aj, sj = bf_alpha[ij], bf_is_self[ij]

                result = 0.0
                for a in range(n_xi):
                    xi1 = xi[a]
                    for c in range(n_xi):
                        xi2 = xi[c]
                        for b in range(n_eta):
                            eta1 = eta[b]
                            J1 = half_R**3 * (xi1**2 - eta1**2)
                            for dd in range(n_eta):
                                eta2 = eta[dd]
                                J2 = half_R**3 * (xi2**2 - eta2**2)

                                phi_int = 0.0
                                for ip in range(n_phi):
                                    r12 = _compute_r12_scalar(
                                        xi1, eta1, xi2, eta2, dphi[ip], R
                                    )
                                    r12 = max(r12, 1e-15)

                                    fi = _eval_sym_with_r12(
                                        ji, ki, li, mi, pi_, ai, si,
                                        xi1, eta1, xi2, eta2, r12,
                                    )
                                    fj = _eval_sym_with_r12(
                                        jj, kj, lj, mj, pj, aj, sj,
                                        xi1, eta1, xi2, eta2, r12,
                                    )
                                    phi_int += fi * fj * w_phi[ip]

                                result += (w_xi[a] * w_eta[b] * w_xi[c]
                                           * w_eta[dd] * J1 * J2
                                           * phi_int * 2.0 * np.pi)

                S[i, ij] = result
                S[ij, i] = result

        return S

    # ==================================================================
    # Hamiltonian matrix kernel (p=0): kinetic (IBP) + V_ne + V_ee
    # ==================================================================

    @njit(cache=True)
    def hamiltonian_p0_kernel(
        bf_j: np.ndarray, bf_k: np.ndarray,
        bf_l: np.ndarray, bf_m: np.ndarray,
        bf_alpha: np.ndarray, bf_is_self: np.ndarray,
        xi: np.ndarray, eta: np.ndarray,
        w_xi: np.ndarray, w_eta: np.ndarray,
        R: float, Z_A: float, Z_B: float,
    ) -> np.ndarray:
        """Compute full Hamiltonian matrix for p=0 basis.

        Kinetic energy via integration by parts, V_ne direct, V_ee via
        elliptic integral K.
        """
        n_basis = len(bf_j)
        n_xi = len(xi)
        n_eta = len(eta)
        half_R = R / 2.0
        phi_factor = (2.0 * np.pi)**2
        T_prefactor = R / 4.0
        H = np.zeros((n_basis, n_basis))

        for i_bf in range(n_basis):
            ji, ki, li, mi = bf_j[i_bf], bf_k[i_bf], bf_l[i_bf], bf_m[i_bf]
            ai, si = bf_alpha[i_bf], bf_is_self[i_bf]
            for j_bf in range(i_bf, n_basis):
                jj, kj, lj, mj = (bf_j[j_bf], bf_k[j_bf],
                                   bf_l[j_bf], bf_m[j_bf])
                aj, sj = bf_alpha[j_bf], bf_is_self[j_bf]

                T_ij = 0.0
                V_ij = 0.0

                for a in range(n_xi):
                    xi1 = xi[a]
                    xi1_sq_m1 = xi1**2 - 1.0

                    for c in range(n_xi):
                        xi2 = xi[c]
                        xi2_sq_m1 = xi2**2 - 1.0
                        w_ac = w_xi[a] * w_xi[c]

                        for b in range(n_eta):
                            eta1 = eta[b]
                            one_m_eta1_sq = 1.0 - eta1**2
                            J1 = half_R**3 * (xi1**2 - eta1**2)

                            # V_ne terms for electron 1 (scalar)
                            r1A = half_R * (xi1 + eta1)
                            r1B = half_R * abs(xi1 - eta1)

                            # V_ee cylindrical coords for electron 1
                            rho1 = half_R * np.sqrt(
                                max(xi1_sq_m1 * one_m_eta1_sq, 0.0)
                            )
                            z1 = half_R * xi1 * eta1

                            for d in range(n_eta):
                                eta2 = eta[d]
                                one_m_eta2_sq = 1.0 - eta2**2
                                J2 = half_R**3 * (xi2**2 - eta2**2)

                                # Basis derivatives
                                fi, dfi_xi1, dfi_eta1, dfi_xi2, dfi_eta2 = \
                                    _eval_sym_derivs(
                                        ji, ki, li, mi, ai, si,
                                        xi1, eta1, xi2, eta2,
                                    )
                                fj, dfj_xi1, dfj_eta1, dfj_xi2, dfj_eta2 = \
                                    _eval_sym_derivs(
                                        jj, kj, lj, mj, aj, sj,
                                        xi1, eta1, xi2, eta2,
                                    )

                                # Kinetic T1 + T2 (integration by parts)
                                T1 = (xi1_sq_m1 * dfi_xi1 * dfj_xi1
                                      + one_m_eta1_sq * dfi_eta1
                                      * dfj_eta1) * J2
                                T2 = (xi2_sq_m1 * dfi_xi2 * dfj_xi2
                                      + one_m_eta2_sq * dfi_eta2
                                      * dfj_eta2) * J1

                                w_total = w_ac * w_eta[b] * w_eta[d]
                                T_ij += (T1 + T2) * w_total

                                # V_ne
                                r2A = half_R * (xi2 + eta2)
                                r2B = half_R * abs(xi2 - eta2)

                                V_ne = (
                                    -Z_A / max(r1A, 1e-15)
                                    - Z_B / max(r1B, 1e-15)
                                    - Z_A / max(r2A, 1e-15)
                                    - Z_B / max(r2B, 1e-15)
                                )

                                # V_ee via elliptic K
                                rho2 = half_R * np.sqrt(
                                    max(xi2_sq_m1 * one_m_eta2_sq, 0.0)
                                )
                                z2 = half_R * xi2 * eta2

                                dz = z1 - z2
                                s = (rho1 + rho2)**2 + dz**2

                                if s > 1e-30:
                                    k2 = 4.0 * rho1 * rho2 / s
                                    k2 = min(max(k2, 0.0), 1.0 - 1e-15)
                                    K_val = _ellipk_approx(k2)
                                    V_ee = K_val / np.sqrt(s)
                                else:
                                    V_ee = 0.0

                                fifj = fi * fj
                                V_ne_contrib = fifj * V_ne * J1 * J2
                                V_ee_contrib = fifj * V_ee * J1 * J2

                                V_ij += w_total * (
                                    phi_factor * V_ne_contrib
                                    + 8.0 * np.pi * V_ee_contrib
                                )

                H[i_bf, j_bf] = T_prefactor * phi_factor * T_ij + V_ij
                H[j_bf, i_bf] = H[i_bf, j_bf]

        return H

    # ==================================================================
    # T + V_ne only kernel (no V_ee) — for Neumann V_ee replacement
    # ==================================================================

    @njit(cache=True)
    def hamiltonian_tvne_p0_kernel(
        bf_j: np.ndarray, bf_k: np.ndarray,
        bf_l: np.ndarray, bf_m: np.ndarray,
        bf_alpha: np.ndarray, bf_is_self: np.ndarray,
        xi: np.ndarray, eta: np.ndarray,
        w_xi: np.ndarray, w_eta: np.ndarray,
        R: float, Z_A: float, Z_B: float,
    ) -> np.ndarray:
        """Compute T + V_ne Hamiltonian (no V_ee) for p=0 basis.

        Same as hamiltonian_p0_kernel but without electron-electron
        repulsion. Used when V_ee is computed algebraically (Neumann).
        """
        n_basis = len(bf_j)
        n_xi = len(xi)
        n_eta = len(eta)
        half_R = R / 2.0
        phi_factor = (2.0 * np.pi)**2
        T_prefactor = R / 4.0
        H = np.zeros((n_basis, n_basis))

        for i_bf in range(n_basis):
            ji, ki, li, mi = bf_j[i_bf], bf_k[i_bf], bf_l[i_bf], bf_m[i_bf]
            ai, si = bf_alpha[i_bf], bf_is_self[i_bf]
            for j_bf in range(i_bf, n_basis):
                jj, kj, lj, mj = (bf_j[j_bf], bf_k[j_bf],
                                   bf_l[j_bf], bf_m[j_bf])
                aj, sj = bf_alpha[j_bf], bf_is_self[j_bf]

                T_ij = 0.0
                V_ij = 0.0

                for a in range(n_xi):
                    xi1 = xi[a]
                    xi1_sq_m1 = xi1**2 - 1.0

                    for c in range(n_xi):
                        xi2 = xi[c]
                        xi2_sq_m1 = xi2**2 - 1.0
                        w_ac = w_xi[a] * w_xi[c]

                        for b in range(n_eta):
                            eta1 = eta[b]
                            one_m_eta1_sq = 1.0 - eta1**2
                            J1 = half_R**3 * (xi1**2 - eta1**2)

                            r1A = half_R * (xi1 + eta1)
                            r1B = half_R * abs(xi1 - eta1)

                            for d in range(n_eta):
                                eta2 = eta[d]
                                one_m_eta2_sq = 1.0 - eta2**2
                                J2 = half_R**3 * (xi2**2 - eta2**2)

                                fi, dfi_xi1, dfi_eta1, dfi_xi2, dfi_eta2 = \
                                    _eval_sym_derivs(
                                        ji, ki, li, mi, ai, si,
                                        xi1, eta1, xi2, eta2,
                                    )
                                fj, dfj_xi1, dfj_eta1, dfj_xi2, dfj_eta2 = \
                                    _eval_sym_derivs(
                                        jj, kj, lj, mj, aj, sj,
                                        xi1, eta1, xi2, eta2,
                                    )

                                T1 = (xi1_sq_m1 * dfi_xi1 * dfj_xi1
                                      + one_m_eta1_sq * dfi_eta1
                                      * dfj_eta1) * J2
                                T2 = (xi2_sq_m1 * dfi_xi2 * dfj_xi2
                                      + one_m_eta2_sq * dfi_eta2
                                      * dfj_eta2) * J1

                                w_total = w_ac * w_eta[b] * w_eta[d]
                                T_ij += (T1 + T2) * w_total

                                r2A = half_R * (xi2 + eta2)
                                r2B = half_R * abs(xi2 - eta2)

                                V_ne = (
                                    -Z_A / max(r1A, 1e-15)
                                    - Z_B / max(r1B, 1e-15)
                                    - Z_A / max(r2A, 1e-15)
                                    - Z_B / max(r2B, 1e-15)
                                )

                                fifj = fi * fj
                                V_ne_contrib = fifj * V_ne * J1 * J2
                                V_ij += w_total * phi_factor * V_ne_contrib

                H[i_bf, j_bf] = T_prefactor * phi_factor * T_ij + V_ij
                H[j_bf, i_bf] = H[i_bf, j_bf]

        return H

    # ==================================================================
    # FD kinetic energy for p>0 (scalar)
    # ==================================================================

    @njit(cache=True)
    def _kinetic_fd_numba(
        j: int, k: int, l: int, m: int, p: int,
        alpha: float, is_self: bool,
        xi1: float, eta1: float, xi2: float, eta2: float,
        dphi: float, R: float, eps: float,
    ) -> float:
        """FD kinetic energy T|phi> at a single point for p>0 basis."""

        def _ev(x1: float, e1: float, x2: float, e2: float) -> float:
            r = _compute_r12_scalar(x1, e1, x2, e2, dphi, R)
            r = max(r, 1e-15)
            return _eval_sym_with_r12(j, k, l, m, p, alpha, is_self,
                                     x1, e1, x2, e2, r)

        f0 = _ev(xi1, eta1, xi2, eta2)
        T = 0.0

        # T1: xi1
        xi1_p = xi1 + eps
        xi1_m = max(xi1 - eps, 1.0 + 1e-8)
        h1 = xi1_p - xi1
        h2 = xi1 - xi1_m
        fp = _ev(xi1_p, eta1, xi2, eta2)
        fm = _ev(xi1_m, eta1, xi2, eta2)
        d2f = 2.0 * (fp * h2 + fm * h1 - f0 * (h1 + h2)) / (
            h1 * h2 * (h1 + h2))
        df = (fp - fm) / (h1 + h2)
        Lxi1 = (xi1**2 - 1.0) * d2f + 2.0 * xi1 * df

        # T1: eta1
        eta1_p = min(eta1 + eps, 1.0 - 1e-8)
        eta1_m = max(eta1 - eps, -1.0 + 1e-8)
        h1e = eta1_p - eta1
        h2e = eta1 - eta1_m
        fp = _ev(xi1, eta1_p, xi2, eta2)
        fm = _ev(xi1, eta1_m, xi2, eta2)
        d2f = 2.0 * (fp * h2e + fm * h1e - f0 * (h1e + h2e)) / (
            h1e * h2e * (h1e + h2e))
        df = (fp - fm) / (h1e + h2e)
        Leta1 = (1.0 - eta1**2) * d2f - 2.0 * eta1 * df

        w1 = xi1**2 - eta1**2
        T += -(2.0 / R**2) * (Lxi1 + Leta1) / max(w1, 1e-15)

        # T2: xi2
        xi2_p = xi2 + eps
        xi2_m = max(xi2 - eps, 1.0 + 1e-8)
        h1 = xi2_p - xi2
        h2 = xi2 - xi2_m
        fp = _ev(xi1, eta1, xi2_p, eta2)
        fm = _ev(xi1, eta1, xi2_m, eta2)
        d2f = 2.0 * (fp * h2 + fm * h1 - f0 * (h1 + h2)) / (
            h1 * h2 * (h1 + h2))
        df = (fp - fm) / (h1 + h2)
        Lxi2 = (xi2**2 - 1.0) * d2f + 2.0 * xi2 * df

        # T2: eta2
        eta2_p = min(eta2 + eps, 1.0 - 1e-8)
        eta2_m = max(eta2 - eps, -1.0 + 1e-8)
        h1e = eta2_p - eta2
        h2e = eta2 - eta2_m
        fp = _ev(xi1, eta1, xi2, eta2_p)
        fm = _ev(xi1, eta1, xi2, eta2_m)
        d2f = 2.0 * (fp * h2e + fm * h1e - f0 * (h1e + h2e)) / (
            h1e * h2e * (h1e + h2e))
        df = (fp - fm) / (h1e + h2e)
        Leta2 = (1.0 - eta2**2) * d2f - 2.0 * eta2 * df

        w2 = xi2**2 - eta2**2
        T += -(2.0 / R**2) * (Lxi2 + Leta2) / max(w2, 1e-15)

        return T

    # ==================================================================
    # Hamiltonian matrix kernel (p>0): 5D numerical integration
    # ==================================================================

    @njit(cache=True)
    def hamiltonian_general_kernel(
        bf_j: np.ndarray, bf_k: np.ndarray,
        bf_l: np.ndarray, bf_m: np.ndarray,
        bf_p: np.ndarray, bf_alpha: np.ndarray,
        bf_is_self: np.ndarray,
        xi: np.ndarray, eta: np.ndarray,
        w_xi: np.ndarray, w_eta: np.ndarray,
        dphi: np.ndarray, w_phi: np.ndarray,
        R: float, Z_A: float, Z_B: float,
    ) -> np.ndarray:
        """Compute full Hamiltonian for general basis (5D, FD kinetic)."""
        n_basis = len(bf_j)
        n_xi = len(xi)
        n_eta = len(eta)
        n_phi = len(dphi)
        half_R = R / 2.0
        eps_fd = 1e-4
        H = np.zeros((n_basis, n_basis))

        for i_bf in range(n_basis):
            ji, ki, li, mi, pi_ = (bf_j[i_bf], bf_k[i_bf], bf_l[i_bf],
                                    bf_m[i_bf], bf_p[i_bf])
            ai, si = bf_alpha[i_bf], bf_is_self[i_bf]

            for j_bf in range(i_bf, n_basis):
                jj, kj, lj, mj, pj = (bf_j[j_bf], bf_k[j_bf],
                                       bf_l[j_bf], bf_m[j_bf], bf_p[j_bf])
                aj, sj = bf_alpha[j_bf], bf_is_self[j_bf]

                val = 0.0
                for a in range(n_xi):
                    xi1 = xi[a]
                    for c in range(n_xi):
                        xi2 = xi[c]
                        for b in range(n_eta):
                            eta1 = eta[b]
                            J1 = half_R**3 * (xi1**2 - eta1**2)
                            for dd in range(n_eta):
                                eta2 = eta[dd]
                                J2 = half_R**3 * (xi2**2 - eta2**2)

                                phi_int = 0.0
                                for ip in range(n_phi):
                                    dp = dphi[ip]
                                    r12 = _compute_r12_scalar(
                                        xi1, eta1, xi2, eta2, dp, R
                                    )
                                    r12 = max(r12, 1e-15)

                                    fi = _eval_sym_with_r12(
                                        ji, ki, li, mi, pi_, ai, si,
                                        xi1, eta1, xi2, eta2, r12,
                                    )
                                    fj_val = _eval_sym_with_r12(
                                        jj, kj, lj, mj, pj, aj, sj,
                                        xi1, eta1, xi2, eta2, r12,
                                    )

                                    # V_ne
                                    r1A = half_R * (xi1 + eta1)
                                    r1B = half_R * abs(xi1 - eta1)
                                    r2A = half_R * (xi2 + eta2)
                                    r2B = half_R * abs(xi2 - eta2)
                                    v_ne = (
                                        -Z_A / max(r1A, 1e-15)
                                        - Z_B / max(r1B, 1e-15)
                                        - Z_A / max(r2A, 1e-15)
                                        - Z_B / max(r2B, 1e-15)
                                    )

                                    # V_ee
                                    v_ee = 1.0 / r12

                                    # Kinetic (FD)
                                    T_fj = _kinetic_fd_numba(
                                        jj, kj, lj, mj, pj, aj, sj,
                                        xi1, eta1, xi2, eta2,
                                        dp, R, eps_fd,
                                    )

                                    H_fj = T_fj + (v_ne + v_ee) * fj_val
                                    phi_int += fi * H_fj * w_phi[ip]

                                phi_int *= 2.0 * np.pi
                                val += (w_xi[a] * w_eta[b] * w_xi[c]
                                        * w_eta[dd] * J1 * J2 * phi_int)

                H[i_bf, j_bf] = val
                H[j_bf, i_bf] = val

        return H

    # ==================================================================
    # r₁₂² partial derivatives (analytical)
    # ==================================================================

    @njit(cache=True)
    def _compute_r12_and_derivs(
        xi1: float, eta1: float, xi2: float, eta2: float,
        dphi: float, R: float,
    ) -> tuple:
        """Compute r₁₂ and analytical partial derivatives of r₁₂².

        From r₁₂² = (R/2)² S where
          S = (ξ₁η₁-ξ₂η₂)² + ρ₁² + ρ₂² - 2ρ₁ρ₂cos(Δφ)
          ρᵢ² = (ξᵢ²-1)(1-ηᵢ²)

        Returns (r12, ∂(r₁₂²)/∂ξ₁, ∂(r₁₂²)/∂η₁,
                      ∂(r₁₂²)/∂ξ₂, ∂(r₁₂²)/∂η₂, ∂(r₁₂²)/∂Δφ).
        """
        half_R = R / 2.0
        half_R_sq = half_R * half_R
        cos_dp = np.cos(dphi)
        sin_dp = np.sin(dphi)

        diff = xi1 * eta1 - xi2 * eta2

        xi1_sq_m1 = xi1 * xi1 - 1.0
        ome1_sq = 1.0 - eta1 * eta1
        xi2_sq_m1 = xi2 * xi2 - 1.0
        ome2_sq = 1.0 - eta2 * eta2

        rho1_sq = max(xi1_sq_m1 * ome1_sq, 0.0)
        rho2_sq = max(xi2_sq_m1 * ome2_sq, 0.0)
        rho1 = np.sqrt(rho1_sq)
        rho2 = np.sqrt(rho2_sq)

        S = diff * diff + rho1_sq + rho2_sq - 2.0 * rho1 * rho2 * cos_dp
        r12 = half_R * np.sqrt(max(S, 0.0))

        # --- ∂S/∂ξ₁ ---
        # = 2η₁(ξ₁η₁-ξ₂η₂) + 2ξ₁(1-η₁²) - 2cos(Δφ)ρ₂·∂ρ₁/∂ξ₁
        # where ∂ρ₁/∂ξ₁ = ξ₁(1-η₁²)/ρ₁
        dS_dxi1 = 2.0 * eta1 * diff + 2.0 * xi1 * ome1_sq
        if rho1 > 1e-15:
            dS_dxi1 -= 2.0 * cos_dp * rho2 * xi1 * ome1_sq / rho1

        # --- ∂S/∂η₁ ---
        # = 2ξ₁(ξ₁η₁-ξ₂η₂) - 2η₁(ξ₁²-1) + 2cos(Δφ)ρ₂·η₁(ξ₁²-1)/ρ₁
        dS_deta1 = 2.0 * xi1 * diff - 2.0 * eta1 * xi1_sq_m1
        if rho1 > 1e-15:
            dS_deta1 += 2.0 * cos_dp * rho2 * eta1 * xi1_sq_m1 / rho1

        # --- ∂S/∂ξ₂ ---
        dS_dxi2 = -2.0 * eta2 * diff + 2.0 * xi2 * ome2_sq
        if rho2 > 1e-15:
            dS_dxi2 -= 2.0 * cos_dp * rho1 * xi2 * ome2_sq / rho2

        # --- ∂S/∂η₂ ---
        dS_deta2 = -2.0 * xi2 * diff - 2.0 * eta2 * xi2_sq_m1
        if rho2 > 1e-15:
            dS_deta2 += 2.0 * cos_dp * rho1 * eta2 * xi2_sq_m1 / rho2

        # --- ∂S/∂Δφ = 2ρ₁ρ₂sin(Δφ) ---
        dS_ddphi = 2.0 * rho1 * rho2 * sin_dp

        return (r12,
                half_R_sq * dS_dxi1, half_R_sq * dS_deta1,
                half_R_sq * dS_dxi2, half_R_sq * dS_deta2,
                half_R_sq * dS_ddphi)

    # ==================================================================
    # Full product-rule derivatives: φ = g(ξ,η) × r₁₂^p
    # ==================================================================

    @njit(cache=True)
    def _full_derivs_with_r12(
        j: int, k: int, l: int, m: int, p: int,
        alpha: float, is_self: bool,
        xi1: float, eta1: float, xi2: float, eta2: float,
        r12: float,
        dr12sq_dxi1: float, dr12sq_deta1: float,
        dr12sq_dxi2: float, dr12sq_deta2: float,
        dr12sq_ddphi: float,
    ) -> tuple:
        """Full derivatives of symmetrized basis function φ = g × r₁₂^p.

        Uses the product rule:
          ∂φ/∂x = (∂g/∂x) r₁₂^p + g × (p/2) r₁₂^{p-2} ∂(r₁₂²)/∂x

        Returns (φ, ∂φ/∂ξ₁, ∂φ/∂η₁, ∂φ/∂ξ₂, ∂φ/∂η₂, ∂φ/∂Δφ).
        """
        # Orbital part g and its derivatives
        g, dg_dxi1, dg_deta1, dg_dxi2, dg_deta2 = _eval_sym_derivs(
            j, k, l, m, alpha, is_self, xi1, eta1, xi2, eta2
        )

        if p == 0:
            return g, dg_dxi1, dg_deta1, dg_dxi2, dg_deta2, 0.0

        # r₁₂^p factor
        r12_p = r12**p

        # (p/2) × r₁₂^{p-2} for the coupling terms
        if r12 > 1e-15:
            r12_coupling = 0.5 * p * r12**(p - 2)
        else:
            r12_coupling = 0.0

        phi = g * r12_p

        # Product rule: ∂φ/∂x = (∂g/∂x)r₁₂^p + g·(p/2)r₁₂^{p-2}·∂(r₁₂²)/∂x
        g_coup = g * r12_coupling
        dphi_dxi1 = dg_dxi1 * r12_p + g_coup * dr12sq_dxi1
        dphi_deta1 = dg_deta1 * r12_p + g_coup * dr12sq_deta1
        dphi_dxi2 = dg_dxi2 * r12_p + g_coup * dr12sq_dxi2
        dphi_deta2 = dg_deta2 * r12_p + g_coup * dr12sq_deta2
        # g has no Δφ dependence, so ∂g/∂Δφ = 0
        dphi_ddphi = g_coup * dr12sq_ddphi

        return phi, dphi_dxi1, dphi_deta1, dphi_dxi2, dphi_deta2, dphi_ddphi

    # ==================================================================
    # Hamiltonian kernel (analytical IBP kinetic for all p)
    # ==================================================================

    @njit(cache=True)
    def hamiltonian_analytical_kernel(
        bf_j: np.ndarray, bf_k: np.ndarray,
        bf_l: np.ndarray, bf_m: np.ndarray,
        bf_p: np.ndarray, bf_alpha: np.ndarray,
        bf_is_self: np.ndarray,
        xi: np.ndarray, eta: np.ndarray,
        w_xi: np.ndarray, w_eta: np.ndarray,
        dphi: np.ndarray, w_phi: np.ndarray,
        R: float, Z_A: float, Z_B: float,
    ) -> np.ndarray:
        """Hamiltonian for general basis with analytical IBP kinetic energy.

        Replaces the FD kinetic energy in hamiltonian_general_kernel with
        exact analytical derivatives via the product rule for phi = g * r12^p.

        V_ee strategy:
        - p_i + p_j = 0: use elliptic integral K(k)/sqrt(s) for the
          azimuthal average of 1/r12 (avoids non-integrable singularity)
        - p_i + p_j > 0: use 1/r12 directly (r12^p regularizes it)

        IBP kinetic includes 3 terms per electron: xi, eta, and azimuthal
        (dphi) derivatives. The dphi term is required for p>0 because r12
        depends on dphi.
        """
        n_basis = len(bf_j)
        n_xi = len(xi)
        n_eta = len(eta)
        n_phi = len(dphi)
        half_R = R / 2.0
        T_prefactor = R / 4.0  # (2/R^2)(R/2)^3 = R/4
        phi_factor = (2.0 * np.pi)**2
        H = np.zeros((n_basis, n_basis))

        for i_bf in range(n_basis):
            ji, ki, li, mi, pi_ = (bf_j[i_bf], bf_k[i_bf], bf_l[i_bf],
                                    bf_m[i_bf], bf_p[i_bf])
            ai, si = bf_alpha[i_bf], bf_is_self[i_bf]

            for j_bf in range(i_bf, n_basis):
                jj, kj, lj, mj, pj = (bf_j[j_bf], bf_k[j_bf],
                                       bf_l[j_bf], bf_m[j_bf], bf_p[j_bf])
                aj, sj = bf_alpha[j_bf], bf_is_self[j_bf]

                q_total = pi_ + pj  # sum of r12 powers

                val = 0.0
                for a in range(n_xi):
                    xi1 = xi[a]
                    xi1_sq_m1 = xi1 * xi1 - 1.0

                    for c in range(n_xi):
                        xi2 = xi[c]
                        xi2_sq_m1 = xi2 * xi2 - 1.0
                        w_ac = w_xi[a] * w_xi[c]

                        for b in range(n_eta):
                            eta1 = eta[b]
                            ome1_sq = 1.0 - eta1 * eta1
                            J1 = half_R**3 * (xi1 * xi1 - eta1 * eta1)

                            r1A = half_R * (xi1 + eta1)
                            r1B = half_R * abs(xi1 - eta1)

                            # Cylindrical coords for V_ee elliptic K
                            rho1 = half_R * np.sqrt(
                                max(xi1_sq_m1 * ome1_sq, 0.0))
                            z1 = half_R * xi1 * eta1

                            for dd in range(n_eta):
                                eta2 = eta[dd]
                                ome2_sq = 1.0 - eta2 * eta2
                                J2 = half_R**3 * (xi2 * xi2 - eta2 * eta2)

                                r2A = half_R * (xi2 + eta2)
                                r2B = half_R * abs(xi2 - eta2)

                                v_ne = (-Z_A / max(r1A, 1e-15)
                                        - Z_B / max(r1B, 1e-15)
                                        - Z_A / max(r2A, 1e-15)
                                        - Z_B / max(r2B, 1e-15))

                                # phi denominators for IBP azimuthal
                                denom1 = xi1_sq_m1 * ome1_sq
                                denom2 = xi2_sq_m1 * ome2_sq
                                phi_fac1 = 0.0
                                phi_fac2 = 0.0
                                if denom1 > 1e-30:
                                    phi_fac1 = (xi1 * xi1 - eta1 * eta1) / denom1
                                if denom2 > 1e-30:
                                    phi_fac2 = (xi2 * xi2 - eta2 * eta2) / denom2

                                # --- 5D integration over dphi ---
                                T_phi_int = 0.0
                                Vne_phi_int = 0.0
                                Vee_phi_int = 0.0

                                for ip in range(n_phi):
                                    dp = dphi[ip]

                                    (r12, dr_dxi1, dr_deta1,
                                     dr_dxi2, dr_deta2,
                                     dr_ddphi) = _compute_r12_and_derivs(
                                        xi1, eta1, xi2, eta2, dp, R)
                                    r12 = max(r12, 1e-15)

                                    (phi_i, di_dxi1, di_deta1,
                                     di_dxi2, di_deta2,
                                     di_ddphi) = _full_derivs_with_r12(
                                        ji, ki, li, mi, pi_, ai, si,
                                        xi1, eta1, xi2, eta2, r12,
                                        dr_dxi1, dr_deta1,
                                        dr_dxi2, dr_deta2, dr_ddphi)

                                    (phi_j, dj_dxi1, dj_deta1,
                                     dj_dxi2, dj_deta2,
                                     dj_ddphi) = _full_derivs_with_r12(
                                        jj, kj, lj, mj, pj, aj, sj,
                                        xi1, eta1, xi2, eta2, r12,
                                        dr_dxi1, dr_deta1,
                                        dr_dxi2, dr_deta2, dr_ddphi)

                                    # IBP kinetic T1 + T2
                                    T1 = (xi1_sq_m1 * di_dxi1 * dj_dxi1
                                          + ome1_sq * di_deta1 * dj_deta1
                                          + phi_fac1 * di_ddphi * dj_ddphi)
                                    T1 *= J2

                                    T2 = (xi2_sq_m1 * di_dxi2 * dj_dxi2
                                          + ome2_sq * di_deta2 * dj_deta2
                                          + phi_fac2 * di_ddphi * dj_ddphi)
                                    T2 *= J1

                                    T_phi_int += (T1 + T2) * w_phi[ip]

                                    # V_ne (always 5D)
                                    fifj = phi_i * phi_j
                                    Vne_phi_int += fifj * v_ne * w_phi[ip]

                                    # V_ee: only via 5D when q > 0
                                    if q_total > 0:
                                        Vee_phi_int += fifj / r12 * w_phi[ip]

                                # Accumulate with proper phi factors
                                w_total = w_ac * w_eta[b] * w_eta[dd]

                                # T: prefactor * 2pi * integral
                                val += w_total * T_prefactor * 2.0 * np.pi * T_phi_int

                                # V_ne: 2pi * integral * J1*J2
                                val += w_total * 2.0 * np.pi * Vne_phi_int * J1 * J2

                                # V_ee
                                if q_total > 0:
                                    # Direct 1/r12 (smooth because r12^q regularizes)
                                    val += w_total * 2.0 * np.pi * Vee_phi_int * J1 * J2
                                else:
                                    # Elliptic K for p_i+p_j=0
                                    # phi_i, phi_j are dphi-independent,
                                    # so we need g_i * g_j (eval at any dphi)
                                    fi_0 = _eval_sym(
                                        ji, ki, li, mi, ai, si,
                                        xi1, eta1, xi2, eta2)
                                    fj_0 = _eval_sym(
                                        jj, kj, lj, mj, aj, sj,
                                        xi1, eta1, xi2, eta2)

                                    rho2 = half_R * np.sqrt(
                                        max(xi2_sq_m1 * ome2_sq, 0.0))
                                    z2 = half_R * xi2 * eta2
                                    dz = z1 - z2
                                    s = (rho1 + rho2)**2 + dz**2

                                    if s > 1e-30:
                                        k2 = 4.0 * rho1 * rho2 / s
                                        k2 = min(max(k2, 0.0), 1.0 - 1e-15)
                                        K_val = _ellipk_approx(k2)
                                        V_ee_K = K_val / np.sqrt(s)
                                    else:
                                        V_ee_K = 0.0

                                    # 8pi K(k)/sqrt(s) is the full
                                    # 2pi * int_0^{2pi} 1/r12 dphi result
                                    val += w_total * 8.0 * np.pi * fi_0 * fj_0 * V_ee_K * J1 * J2

                H[i_bf, j_bf] = val
                H[j_bf, i_bf] = val

        return H

    # ==================================================================
    # Coulomb potential kernel (for prolate_scf)
    # ==================================================================

    @njit(cache=True)
    def coulomb_potential_kernel(
        rho_w_flat: np.ndarray,
        rho_flat: np.ndarray,
        z_flat: np.ndarray,
        n_pts: int,
    ) -> np.ndarray:
        """Compute V_J[i] = sum_j rho_w[j] * 4*K(k2)/sqrt(s) for all i.

        Replaces batched numpy + scipy.ellipk with scalar Numba loop.
        """
        V_J = np.zeros(n_pts)

        for i in range(n_pts):
            rho1 = rho_flat[i]
            z1 = z_flat[i]
            vi = 0.0

            for j in range(n_pts):
                rho2 = rho_flat[j]
                z2 = z_flat[j]

                dz = z1 - z2
                s = (rho1 + rho2)**2 + dz**2

                if s > 1e-30:
                    k2 = 4.0 * rho1 * rho2 / s
                    k2 = min(max(k2, 0.0), 1.0 - 1e-15)
                    K_val = _ellipk_approx(k2)
                    vi += rho_w_flat[j] * 4.0 * K_val / np.sqrt(s)

            V_J[i] = vi

        return V_J

    # ==================================================================
    # V_ee integral kernel (for prolate_scf compute_vee_integral)
    # ==================================================================

    @njit(cache=True)
    def vee_integral_kernel(
        rho_ac_flat: np.ndarray,
        rho_bd_flat: np.ndarray,
        rho_flat: np.ndarray,
        z_flat: np.ndarray,
        n_pts: int,
    ) -> float:
        """Compute (ab|cd) two-electron integral via azimuthal averaging.

        result = sum_i sum_j rho_ac[i] * 8*pi*K(k2)/sqrt(s) * rho_bd[j]
        """
        result = 0.0

        for i in range(n_pts):
            rho1 = rho_flat[i]
            z1 = z_flat[i]
            w_ac_i = rho_ac_flat[i]

            if abs(w_ac_i) < 1e-30:
                continue

            for j in range(n_pts):
                rho2 = rho_flat[j]
                z2 = z_flat[j]

                dz = z1 - z2
                s = (rho1 + rho2)**2 + dz**2

                if s > 1e-30:
                    k2 = 4.0 * rho1 * rho2 / s
                    k2 = min(max(k2, 0.0), 1.0 - 1e-15)
                    K_val = _ellipk_approx(k2)
                    kernel = 8.0 * np.pi * K_val / np.sqrt(s)
                    result += w_ac_i * kernel * rho_bd_flat[j]

        return result

    # ==================================================================
    # Warmup: trigger JIT compilation on tiny inputs
    # ==================================================================

    def warmup_hylleraas_jit() -> None:
        """Compile all Hylleraas kernels on tiny inputs."""
        xi = np.array([1.5, 2.0])
        eta = np.array([-0.5, 0.5])
        w = np.array([0.5, 0.5])
        dphi = np.array([0.0, np.pi])
        w_phi = np.array([np.pi, np.pi])

        bf_j = np.array([0], dtype=np.int32)
        bf_k = np.array([0], dtype=np.int32)
        bf_l = np.array([0], dtype=np.int32)
        bf_m = np.array([0], dtype=np.int32)
        bf_p = np.array([0], dtype=np.int32)
        bf_alpha = np.array([1.0])
        bf_is_self = np.array([True])

        # p=0 kernels
        overlap_matrix_p0_kernel(
            bf_j, bf_k, bf_l, bf_m, bf_alpha, bf_is_self,
            xi, eta, w, w, 2.0,
        )
        hamiltonian_p0_kernel(
            bf_j, bf_k, bf_l, bf_m, bf_alpha, bf_is_self,
            xi, eta, w, w, 2.0, 1.0, 1.0,
        )

        # p>0 kernels
        bf_p_1 = np.array([1], dtype=np.int32)
        overlap_matrix_5d_kernel(
            bf_j, bf_k, bf_l, bf_m, bf_p_1, bf_alpha, bf_is_self,
            xi, eta, w, w, dphi, w_phi, 2.0,
        )
        hamiltonian_general_kernel(
            bf_j, bf_k, bf_l, bf_m, bf_p_1, bf_alpha, bf_is_self,
            xi, eta, w, w, dphi, w_phi, 2.0, 1.0, 1.0,
        )

        # Analytical IBP kernel
        hamiltonian_analytical_kernel(
            bf_j, bf_k, bf_l, bf_m, bf_p_1, bf_alpha, bf_is_self,
            xi, eta, w, w, dphi, w_phi, 2.0, 1.0, 1.0,
        )

        # Coulomb kernel
        rho_w = np.array([0.1, 0.2])
        rho = np.array([0.5, 1.0])
        z = np.array([0.0, 1.0])
        coulomb_potential_kernel(rho_w, rho, z, 2)
        vee_integral_kernel(rho_w, rho_w, rho, z, 2)
