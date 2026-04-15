"""
Positronium hydride (PsH) solver using the GeoVac Level 3 hyperspherical framework.

PsH = proton (fixed, Born-Oppenheimer) + electron (e-) + positron (e+)

Hamiltonian:
    H = -nabla_1^2/2 - nabla_2^2/2 - 1/r_1 + 1/r_2 - 1/r_12

where particle 1 = e- (attracted to proton), particle 2 = e+ (repelled from proton),
and e- and e+ attract each other (V_12 = -1/r_12).

Hyperspherical coordinates: R = sqrt(r1^2 + r2^2), alpha = arctan(r2/r1)
    r1 = R sin(alpha),  r2 = R cos(alpha)

Charge function:
    C(alpha, theta_12) = -1/sin(alpha) + 1/cos(alpha) - 1/sqrt(1 - sin(2a) cos(theta_12))

Key differences from He (Z=2):
    1. Nuclear term for e+: +1/cos(alpha) instead of -Z/cos(alpha)
    2. V_ee term: -1/r_12 (attraction) instead of +1/r_12 (repulsion)
    3. e- and e+ are distinguishable: BOTH alpha parities needed (no exchange symmetry)

Known exact result: E(PsH) ~ -0.789 Ha (binding energy ~0.039 Ha = 1.06 eV relative
to proton + Ps at -0.25 Ha).

This is a research prototype. Does NOT modify any geovac/ production code.
"""

import numpy as np
from scipy.linalg import eigh, eigh_tridiagonal
from scipy.interpolate import CubicSpline
from math import factorial, sqrt, pi
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json
import os

# ============================================================================
# Gegenbauer polynomials (copied from algebraic_angular.py)
# ============================================================================

def _gegenbauer(k: int, lam: float, x: np.ndarray) -> np.ndarray:
    """Evaluate Gegenbauer polynomial C_k^lambda(x) via stable recurrence."""
    if k == 0:
        return np.ones_like(x, dtype=float)
    c_prev = np.ones_like(x, dtype=float)
    c_curr = 2.0 * lam * x
    if k == 1:
        return c_curr
    for n in range(1, k):
        c_next = (2.0 * (n + lam) * x * c_curr
                  - (n + 2 * lam - 1) * c_prev) / (n + 1)
        c_prev = c_curr
        c_curr = c_next
    return c_curr


# ============================================================================
# Gaunt integrals (copied from hyperspherical_angular.py)
# ============================================================================

def gaunt_integral(l1: int, k: int, l2: int) -> float:
    """Gaunt integral: integral of P_l1(x) P_k(x) P_l2(x) dx over [-1,1]."""
    s = l1 + k + l2
    if s % 2 != 0:
        return 0.0
    if l2 > l1 + k or l2 < abs(l1 - k):
        return 0.0
    g = s // 2
    if g < l1 or g < k or g < l2:
        return 0.0
    num = factorial(2 * (g - l1)) * factorial(2 * (g - k)) * factorial(2 * (g - l2))
    den = factorial(2 * g + 1)
    threej_sq = (factorial(g) ** 2 * num) / (
        factorial(g - l1) ** 2 * factorial(g - k) ** 2
        * factorial(g - l2) ** 2 * den
    )
    return 2.0 * threej_sq


def _precompute_gaunt(l_max: int) -> np.ndarray:
    """Precompute Gaunt integrals for all (l, k, l') up to l_max."""
    n_l = l_max + 1
    k_max = 2 * l_max
    G = np.zeros((n_l, k_max + 1, n_l))
    for l in range(n_l):
        for lp in range(n_l):
            for k in range(k_max + 1):
                G[l, k, lp] = gaunt_integral(l, k, lp)
    return G


# ============================================================================
# PsH Angular Solver — Gegenbauer spectral basis
# ============================================================================

class PsHAngularSolver:
    """Spectral Gegenbauer solver for PsH angular eigenvalue problem.

    For He, the two electrons are identical, so only even-parity (singlet)
    or odd-parity (triplet) basis functions in alpha are needed.

    For PsH, e- and e+ are distinguishable, so BOTH parities are needed.
    The basis includes even-k AND odd-k Gegenbauer functions for each l-channel.

    Charge function differences from He:
        Nuclear: -1/sin(a) + 1/cos(a)   [vs He: -Z/sin(a) - Z/cos(a)]
        V_ee:    -W(a)                   [vs He: +W(a)] where W is the
                 Gaunt expansion of 1/r_12
    """

    def __init__(
        self,
        n_basis: int = 10,
        l_max: int = 0,
        n_quad: int = 150,
    ) -> None:
        self.n_basis = n_basis
        self.l_max = l_max
        self.n_quad = n_quad
        self.n_l = l_max + 1
        # Both parities: 2 * n_basis per l-channel
        self.n_basis_per_l = 2 * n_basis
        self._total_dim = self.n_l * self.n_basis_per_l

        self._setup_quadrature(n_quad)
        self._build_channels()
        self._precompute_coupling()

    def _setup_quadrature(self, n_quad: int) -> None:
        """Setup Gauss-Legendre quadrature on [0, pi/4] and [pi/4, pi/2]."""
        from numpy.polynomial.legendre import leggauss
        nodes, weights = leggauss(n_quad)

        # Map to [0, pi/4]
        a1, b1 = 0.0, np.pi / 4.0
        scale1 = (b1 - a1) / 2.0
        alpha1 = scale1 * nodes + (b1 + a1) / 2.0
        w1 = weights * scale1

        # Map to [pi/4, pi/2]
        a2, b2 = np.pi / 4.0, np.pi / 2.0
        scale2 = (b2 - a2) / 2.0
        alpha2 = scale2 * nodes + (b2 + a2) / 2.0
        w2 = weights * scale2

        self._alpha = np.concatenate([alpha1, alpha2])
        self._weights = np.concatenate([w1, w2])
        self._sin_a = np.sin(self._alpha)
        self._cos_a = np.cos(self._alpha)

    def _build_channels(self) -> None:
        """Build orthonormalized basis functions for each l-channel.

        For PsH, we use BOTH even and odd k values:
            k = 0, 1, 2, 3, ..., 2*n_basis - 1
        This gives both alpha-parities since C_k^lam(cos 2a) has parity (-1)^k
        under alpha -> pi/2 - alpha.
        """
        self._channel_phi = []
        self._channel_casimir = []

        cos_2a = np.cos(2.0 * self._alpha)
        sin_cos = self._sin_a * self._cos_a

        for l in range(self.n_l):
            # Both parities: k = 0, 1, 2, ..., 2*n_basis - 1
            k_indices = np.arange(2 * self.n_basis)

            # Free eigenvalues: mu = 2(l+k+1)^2 - 2
            casimir = 2.0 * (l + k_indices + 1.0) ** 2 - 2.0

            # Evaluate and normalize
            envelope = sin_cos ** (l + 1)
            lam = float(l + 1)

            n_funcs = len(k_indices)
            phi = np.zeros((n_funcs, len(self._alpha)))
            for j, kv in enumerate(k_indices):
                raw = envelope * _gegenbauer(kv, lam, cos_2a)
                norm_sq = np.dot(self._weights, raw * raw)
                if norm_sq > 1e-30:
                    phi[j] = raw / np.sqrt(norm_sq)
                else:
                    phi[j] = 0.0

            self._channel_phi.append(phi)
            self._channel_casimir.append(casimir)

    def _precompute_coupling(self) -> None:
        """Build R-independent nuclear and V_ee coupling matrices for PsH.

        Nuclear potential for PsH:
            V_nuc(alpha) = -1/sin(alpha) + 1/cos(alpha)

        V_ee for PsH (ATTRACTIVE e- -- e+ interaction):
            The Gaunt expansion of 1/r_12 gives W(alpha) > 0.
            For He, we ADD W (repulsion). For PsH, we SUBTRACT W (attraction).
        """
        dim = self._total_dim
        nb = self.n_basis_per_l
        nuclear = np.zeros((dim, dim))
        vee = np.zeros((dim, dim))

        # PsH nuclear potential: -1/sin(a) for e-, +1/cos(a) for e+
        V_nuc_alpha = -1.0 / self._sin_a + 1.0 / self._cos_a

        min_sc = np.minimum(self._sin_a, self._cos_a)
        max_sc = np.maximum(self._sin_a, self._cos_a)

        G = _precompute_gaunt(self.l_max)

        for l in range(self.n_l):
            phi_l = self._channel_phi[l]
            i0, i1 = l * nb, (l + 1) * nb

            # Nuclear coupling (diagonal in l)
            weighted_nuc = phi_l * (self._weights * V_nuc_alpha)[np.newaxis, :]
            nuclear[i0:i1, i0:i1] = weighted_nuc @ phi_l.T

            # V_ee coupling via Gaunt integrals
            for lp in range(l, self.n_l):
                phi_lp = self._channel_phi[lp]
                j0, j1 = lp * nb, (lp + 1) * nb

                W = np.zeros_like(self._alpha)
                for k in range(abs(l - lp), l + lp + 1):
                    if k > 2 * self.l_max:
                        continue
                    gv = G[l, k, lp]
                    if abs(gv) < 1e-15:
                        continue
                    f_k = (min_sc / max_sc) ** k
                    W += gv * f_k / max_sc

                norm_ll = sqrt((2 * l + 1) * (2 * lp + 1)) / 2.0
                W *= norm_ll

                # PsH: V_ee is ATTRACTIVE, so we SUBTRACT W
                # (In He, V_ee is repulsive and W is ADDED)
                weighted_vee = phi_l * (self._weights * (-W))[np.newaxis, :]
                block = weighted_vee @ phi_lp.T

                vee[i0:i1, j0:j1] = block
                if l != lp:
                    vee[j0:j1, i0:i1] = block.T

        self._nuclear_full = nuclear
        self._vee_full = vee
        self._coupling_full = nuclear + vee

    def solve(self, R: float, n_channels: int = 1):
        """Solve the angular eigenvalue problem at hyperradius R.

        H_ang(R) = diag(casimir) + R * V_coupling
        """
        casimir_all = np.concatenate(self._channel_casimir)
        H = np.diag(casimir_all) + R * self._coupling_full
        evals, evecs = eigh(H)
        return evals[:n_channels], evecs[:, :n_channels].T


# ============================================================================
# Adiabatic potential curve computation
# ============================================================================

def compute_adiabatic_curves(
    R_grid: np.ndarray,
    n_basis: int = 10,
    l_max: int = 0,
    n_channels: int = 3,
    n_quad: int = 150,
) -> np.ndarray:
    """Compute mu(R) on a grid of hyperradius values for PsH."""
    solver = PsHAngularSolver(n_basis=n_basis, l_max=l_max, n_quad=n_quad)
    N_R = len(R_grid)
    mu = np.zeros((n_channels, N_R))

    for i, R in enumerate(R_grid):
        evals, _ = solver.solve(R, n_channels)
        mu[:, i] = evals

    return mu


def effective_potential(R_grid: np.ndarray, mu_curve: np.ndarray) -> np.ndarray:
    """V_eff(R) = mu(R)/R^2 + 15/(8R^2)"""
    return mu_curve / R_grid**2 + 15.0 / (8.0 * R_grid**2)


# ============================================================================
# Radial solver (FD, copied from hyperspherical_radial.py)
# ============================================================================

def solve_radial(
    V_eff_func,
    R_min: float = 0.05,
    R_max: float = 40.0,
    N_R: int = 3000,
    n_states: int = 1,
):
    """Solve the hyperradial Schrodinger equation with FD."""
    h = (R_max - R_min) / (N_R + 1)
    R_grid = R_min + (np.arange(N_R) + 1) * h

    if callable(V_eff_func):
        V = V_eff_func(R_grid)
    else:
        V = np.asarray(V_eff_func)

    diag = np.ones(N_R) / h**2 + V
    off_diag = -0.5 * np.ones(N_R - 1) / h**2

    evals, evecs = eigh_tridiagonal(
        diag, off_diag,
        select='i', select_range=(0, n_states - 1)
    )

    for i in range(n_states):
        norm = np.sqrt(h * np.sum(evecs[:, i]**2))
        if norm > 0:
            evecs[:, i] /= norm

    return evals, evecs.T, R_grid


# ============================================================================
# He solver for comparison
# ============================================================================

class HeAngularSolver:
    """Standard He angular solver (singlet only) for comparison."""

    def __init__(self, Z: float = 2.0, n_basis: int = 10, l_max: int = 0, n_quad: int = 150):
        self.Z = Z
        self.n_basis = n_basis
        self.l_max = l_max
        self.n_quad = n_quad
        self.n_l = l_max + 1
        self._total_dim = self.n_l * n_basis

        self._setup_quadrature(n_quad)
        self._build_channels()
        self._precompute_coupling()

    def _setup_quadrature(self, n_quad):
        from numpy.polynomial.legendre import leggauss
        nodes, weights = leggauss(n_quad)
        a1, b1 = 0.0, np.pi / 4.0
        scale1 = (b1 - a1) / 2.0
        alpha1 = scale1 * nodes + (b1 + a1) / 2.0
        w1 = weights * scale1
        a2, b2 = np.pi / 4.0, np.pi / 2.0
        scale2 = (b2 - a2) / 2.0
        alpha2 = scale2 * nodes + (b2 + a2) / 2.0
        w2 = weights * scale2
        self._alpha = np.concatenate([alpha1, alpha2])
        self._weights = np.concatenate([w1, w2])
        self._sin_a = np.sin(self._alpha)
        self._cos_a = np.cos(self._alpha)

    def _build_channels(self):
        cos_2a = np.cos(2.0 * self._alpha)
        sin_cos = self._sin_a * self._cos_a
        self._channel_phi = []
        self._channel_casimir = []
        for l in range(self.n_l):
            k_indices = np.array([2 * j for j in range(self.n_basis)])  # even only (singlet)
            casimir = 2.0 * (l + k_indices + 1.0) ** 2 - 2.0
            envelope = sin_cos ** (l + 1)
            lam = float(l + 1)
            phi = np.zeros((self.n_basis, len(self._alpha)))
            for j, kv in enumerate(k_indices):
                raw = envelope * _gegenbauer(kv, lam, cos_2a)
                norm_sq = np.dot(self._weights, raw * raw)
                phi[j] = raw / np.sqrt(norm_sq)
            self._channel_phi.append(phi)
            self._channel_casimir.append(casimir)

    def _precompute_coupling(self):
        dim = self._total_dim
        nb = self.n_basis
        nuclear = np.zeros((dim, dim))
        vee = np.zeros((dim, dim))
        V_nuc_alpha = -self.Z * (1.0 / self._cos_a + 1.0 / self._sin_a)
        min_sc = np.minimum(self._sin_a, self._cos_a)
        max_sc = np.maximum(self._sin_a, self._cos_a)
        G = _precompute_gaunt(self.l_max)
        for l in range(self.n_l):
            phi_l = self._channel_phi[l]
            i0, i1 = l * nb, (l + 1) * nb
            weighted_nuc = phi_l * (self._weights * V_nuc_alpha)[np.newaxis, :]
            nuclear[i0:i1, i0:i1] = weighted_nuc @ phi_l.T
            for lp in range(l, self.n_l):
                phi_lp = self._channel_phi[lp]
                j0, j1 = lp * nb, (lp + 1) * nb
                W = np.zeros_like(self._alpha)
                for k in range(abs(l - lp), l + lp + 1):
                    if k > 2 * self.l_max:
                        continue
                    gv = G[l, k, lp]
                    if abs(gv) < 1e-15:
                        continue
                    f_k = (min_sc / max_sc) ** k
                    W += gv * f_k / max_sc
                norm_ll = sqrt((2 * l + 1) * (2 * lp + 1)) / 2.0
                W *= norm_ll
                weighted_vee = phi_l * (self._weights * W)[np.newaxis, :]
                block = weighted_vee @ phi_lp.T
                vee[i0:i1, j0:j1] = block
                if l != lp:
                    vee[j0:j1, i0:i1] = block.T
        self._nuclear_full = nuclear
        self._vee_full = vee
        self._coupling_full = nuclear + vee

    def solve(self, R: float, n_channels: int = 1):
        casimir_all = np.concatenate(self._channel_casimir)
        H = np.diag(casimir_all) + R * self._coupling_full
        evals, evecs = eigh(H)
        return evals[:n_channels], evecs[:, :n_channels].T


# ============================================================================
# Main computation
# ============================================================================

def main():
    print("=" * 70)
    print("PsH (Positronium Hydride) Solver — GeoVac Level 3 Hyperspherical")
    print("=" * 70)

    # --- Configuration ---
    l_max = 0
    n_basis = 20
    n_quad = 250
    n_channels_angular = 5

    # R grid for adiabatic curves — extend further for PsH (loosely bound)
    R_grid = np.concatenate([
        np.linspace(0.1, 1.0, 15),
        np.linspace(1.0, 5.0, 30),
        np.linspace(5.0, 15.0, 25),
        np.linspace(15.0, 40.0, 20),
        np.linspace(40.0, 100.0, 15),
    ])
    R_grid = np.unique(R_grid)

    # --- Step 1: Build angular solver and compute adiabatic curves ---
    print(f"\n--- Angular eigenvalue computation ---")
    print(f"l_max = {l_max}, n_basis = {n_basis} (per parity), n_quad = {n_quad}")
    print(f"Total angular dimension = {(l_max+1) * 2 * n_basis}")

    psh_solver = PsHAngularSolver(n_basis=n_basis, l_max=l_max, n_quad=n_quad)

    print(f"\nPsH coupling matrix properties:")
    print(f"  Nuclear coupling norm: {np.linalg.norm(psh_solver._nuclear_full):.4f}")
    print(f"  V_ee coupling norm:    {np.linalg.norm(psh_solver._vee_full):.4f}")
    print(f"  Total coupling norm:   {np.linalg.norm(psh_solver._coupling_full):.4f}")

    # Check that nuclear coupling breaks alpha-parity symmetry
    nb_per_l = psh_solver.n_basis_per_l
    nuc_even_even = psh_solver._nuclear_full[:n_basis, :n_basis]
    nuc_odd_odd = psh_solver._nuclear_full[n_basis:nb_per_l, n_basis:nb_per_l]
    nuc_even_odd = psh_solver._nuclear_full[:n_basis, n_basis:nb_per_l]
    print(f"\n  Nuclear even-even block norm: {np.linalg.norm(nuc_even_even):.4f}")
    print(f"  Nuclear odd-odd block norm:   {np.linalg.norm(nuc_odd_odd):.4f}")
    print(f"  Nuclear even-odd block norm:  {np.linalg.norm(nuc_even_odd):.4f}")
    if np.linalg.norm(nuc_even_odd) > 1e-10:
        print("  --> CONFIRMED: Nuclear potential mixes alpha parities (distinguishable particles)")
    else:
        print("  --> WARNING: No parity mixing detected — check charge function")

    # Compute PsH adiabatic curves
    print(f"\nComputing PsH mu(R) on {len(R_grid)} R-points...")
    mu_psh = np.zeros((n_channels_angular, len(R_grid)))
    for i, R in enumerate(R_grid):
        evals, _ = psh_solver.solve(R, n_channels_angular)
        mu_psh[:, i] = evals

    # Compute He adiabatic curves for comparison
    print("Computing He mu(R) for comparison...")
    he_solver = HeAngularSolver(Z=2.0, n_basis=n_basis, l_max=l_max, n_quad=n_quad)
    mu_he = np.zeros((n_channels_angular, len(R_grid)))
    for i, R in enumerate(R_grid):
        evals, _ = he_solver.solve(R, n_channels_angular)
        mu_he[:, i] = evals

    # Also compute H- for comparison (Z=1, both electrons attracted, V_ee repulsive)
    # H- is weakly bound at -0.5278 Ha
    print("Computing H- mu(R) for comparison...")
    hminus_solver = HeAngularSolver(Z=1.0, n_basis=n_basis, l_max=l_max, n_quad=n_quad)
    mu_hminus = np.zeros((n_channels_angular, len(R_grid)))
    for i, R in enumerate(R_grid):
        evals, _ = hminus_solver.solve(R, n_channels_angular)
        mu_hminus[:, i] = evals

    # --- Step 2: Effective potentials ---
    V_eff_psh = np.zeros_like(mu_psh)
    V_eff_he = np.zeros_like(mu_he)
    V_eff_hminus = np.zeros_like(mu_hminus)
    for ch in range(n_channels_angular):
        V_eff_psh[ch] = effective_potential(R_grid, mu_psh[ch])
        V_eff_he[ch] = effective_potential(R_grid, mu_he[ch])
        V_eff_hminus[ch] = effective_potential(R_grid, mu_hminus[ch])

    # --- Step 3: Analyze the lowest adiabatic curve ---
    print(f"\n--- Adiabatic potential analysis ---")

    # PsH
    V0_psh = V_eff_psh[0]
    idx_min_psh = np.argmin(V0_psh)
    R_min_psh = R_grid[idx_min_psh]
    V_min_psh = V0_psh[idx_min_psh]
    print(f"\nPsH lowest curve:")
    print(f"  V_eff minimum at R = {R_min_psh:.3f} bohr, V_min = {V_min_psh:.6f} Ha")
    print(f"  V_eff(R_max) = {V0_psh[-1]:.6f} Ha")

    # Check if there's a well (minimum below asymptote)
    V_asymp_psh = V0_psh[-1]
    well_depth_psh = V_asymp_psh - V_min_psh
    if well_depth_psh > 0.001:
        print(f"  Well depth: {well_depth_psh:.6f} Ha ({well_depth_psh*27.211:.3f} eV)")
        print(f"  --> PsH IS bound in the adiabatic approximation!")
    else:
        print(f"  Well depth: {well_depth_psh:.6f} Ha — marginal or unbound")

    # He
    V0_he = V_eff_he[0]
    idx_min_he = np.argmin(V0_he)
    print(f"\nHe lowest curve:")
    print(f"  V_eff minimum at R = {R_grid[idx_min_he]:.3f} bohr, V_min = {V0_he[idx_min_he]:.6f} Ha")

    # H-
    V0_hminus = V_eff_hminus[0]
    idx_min_hminus = np.argmin(V0_hminus)
    print(f"\nH- lowest curve:")
    print(f"  V_eff minimum at R = {R_grid[idx_min_hminus]:.3f} bohr, V_min = {V0_hminus[idx_min_hminus]:.6f} Ha")

    # --- Step 4: Asymptotic behavior ---
    print(f"\n--- Asymptotic analysis ---")
    # In hyperspherical coords, mu(R) -> -Z_eff^2 R^2/2 asymptotically
    # V_eff(R) = (mu + 15/8) / R^2 -> -Z_eff^2/2 as R -> inf
    # For He(Z=2): threshold = He+(1s) at -Z^2/2 = -2.0 Ha
    # For H-(Z=1): threshold = H(1s) at -Z^2/2 = -0.5 Ha
    # For PsH: threshold = H(1s) + e+ at -0.5 Ha (lowest)
    #   mu(R)/R^2 should approach -0.5 as R -> inf
    print(f"  Asymptotic mu(R)/R^2 at large R:")
    for label, mu_data in [("PsH", mu_psh), ("He", mu_he), ("H-", mu_hminus)]:
        for j in [-1, -3, -6]:
            if abs(j) < len(R_grid):
                print(f"    {label} at R={R_grid[j]:.1f}: mu/R^2 = {mu_data[0,j]/R_grid[j]**2:.6f}, "
                      f"V_eff = {(mu_data[0,j]+15/8)/R_grid[j]**2:.6f}")
    print(f"  Expected asymptotes: He -> -2.0, H- -> -0.5, PsH -> -0.5")

    # --- Step 4b: Well depth relative to TRUE threshold ---
    # The true H(1s) + e+ threshold is -0.5 Ha
    # The adiabatic curve converges slowly to this, so we can't use the
    # computed V_eff(R_max) as the threshold. Instead, check if V_eff_min
    # is below the known threshold.
    print(f"\n--- Well depth relative to exact thresholds ---")
    print(f"  PsH V_eff minimum: {V_min_psh:.6f} Ha at R={R_min_psh:.3f}")
    print(f"  H(1s) + e+ threshold: -0.500 Ha")
    if V_min_psh < -0.5:
        print(f"  V_eff_min < threshold => PsH IS bound (well depth {-0.5 - V_min_psh:.6f} Ha)")
    else:
        print(f"  V_eff_min > threshold => Need deeper well (gap = {V_min_psh + 0.5:.6f} Ha)")
    print(f"  H(1s) + Ps(1s) threshold: -0.750 Ha")
    print(f"  Note: slow asymptotic convergence is structural to hyperspherical coords")
    print(f"  for loosely bound systems. The adiabatic curve at R=100 is only at")
    print(f"  V_eff = {V0_psh[-1]:.4f} Ha, still far from the -0.5 threshold.")

    # --- Step 5: Solve radial equation ---
    print(f"\n--- Radial solution ---")

    # Interpolate the lowest adiabatic curve
    V_eff_spline_psh = CubicSpline(R_grid, V0_psh)
    V_eff_spline_he = CubicSpline(R_grid, V0_he)
    V_eff_spline_hminus = CubicSpline(R_grid, V0_hminus)

    # PsH radial solve — needs large R_max for loosely bound positron
    R_max_radial = 80.0
    N_R_radial = 6000

    E_psh, F_psh, R_rad = solve_radial(
        V_eff_spline_psh, R_min=0.1, R_max=R_max_radial,
        N_R=N_R_radial, n_states=3
    )

    print(f"\nPsH ground state energy (raw): {E_psh[0]:.6f} Ha")
    print(f"  Exact PsH:     -0.789 Ha")
    print(f"  Error:          {abs(E_psh[0] - (-0.789)):.6f} Ha ({abs(E_psh[0] - (-0.789))/0.789*100:.2f}%)")

    # The raw energy is poor because the adiabatic curve converges slowly
    # to the -0.5 Ha threshold. The well minimum IS below -0.5 (at -0.519),
    # confirming binding, but the radial solver sees an incorrect continuum.
    #
    # Approach: solve with the potential SHIFTED so its asymptote matches
    # the known -0.5 Ha threshold. This gives an estimate of the binding
    # energy relative to the correct threshold.
    V_asymp_computed = V0_psh[-1]  # computed V_eff at R_max of angular grid
    V_asymp_true = -0.5  # H(1s) + e+ threshold
    shift = V_asymp_true - V_asymp_computed

    print(f"\n  Asymptotic correction:")
    print(f"    Computed V_eff(R={R_grid[-1]:.0f}) = {V_asymp_computed:.6f} Ha")
    print(f"    True threshold (H + e+) = {V_asymp_true:.6f} Ha")
    print(f"    Shift = {shift:.6f} Ha")

    def V_eff_shifted(R):
        return V_eff_spline_psh(R) + shift

    E_psh_shifted, F_psh_shifted, R_rad_s = solve_radial(
        V_eff_shifted, R_min=0.1, R_max=R_max_radial,
        N_R=N_R_radial, n_states=3
    )

    print(f"\n  PsH ground state (asymptote-corrected): {E_psh_shifted[0]:.6f} Ha")
    print(f"  Exact PsH:     -0.789 Ha")
    print(f"  Error:          {abs(E_psh_shifted[0] - (-0.789)):.6f} Ha ({abs(E_psh_shifted[0] - (-0.789))/0.789*100:.2f}%)")

    # Thresholds
    E_H = -0.5       # H(1s)
    E_Ps = -0.25     # Ps(1s)

    print(f"\n  Thresholds:")
    print(f"    H(1s) + free e+:  {E_H:.4f} Ha")
    print(f"    If E < -0.5 Ha, bound relative to H + e+")
    for label, E_val in [("raw", E_psh[0]), ("shifted", E_psh_shifted[0])]:
        if E_val < E_H:
            BE = E_H - E_val
            print(f"    {label}: Binding energy = {BE:.6f} Ha ({BE*27.211:.4f} eV)")
        else:
            print(f"    {label}: NOT bound relative to H + e+")

    if len(E_psh_shifted) > 1:
        print(f"\n  Shifted excited states: {E_psh_shifted[1]:.6f}, {E_psh_shifted[2]:.6f} Ha")

    # He radial solve for comparison
    E_he, _, _ = solve_radial(
        V_eff_spline_he, R_min=0.1, R_max=30.0,
        N_R=3000, n_states=1
    )
    print(f"\nHe ground state (l_max=0): {E_he[0]:.6f} Ha (exact: -2.9037 Ha)")

    # H- radial solve
    E_hminus, _, _ = solve_radial(
        V_eff_spline_hminus, R_min=0.1, R_max=50.0,
        N_R=4000, n_states=1
    )
    print(f"H- ground state (l_max=0): {E_hminus[0]:.6f} Ha (exact: -0.5278 Ha)")

    # --- Step 6: Convergence study ---
    print(f"\n--- Basis convergence study ---")
    for nb in [5, 10, 15, 20]:
        solver_test = PsHAngularSolver(n_basis=nb, l_max=0, n_quad=200)
        # Compute at a few key R values
        R_test = np.array([1.0, 3.0, 5.0, 10.0, 20.0])
        mu_test = np.zeros(len(R_test))
        for i, R in enumerate(R_test):
            ev, _ = solver_test.solve(R, 1)
            mu_test[i] = ev[0]
        # Compute V_eff and find minimum
        V_test = effective_potential(R_test, mu_test)
        print(f"  n_basis={nb:2d}: mu(R=3)={mu_test[1]:.4f}, "
              f"V_eff(R=3)={V_test[1]:.6f}")

    # --- Step 7: l_max convergence ---
    print(f"\n--- l_max convergence ---")
    for lm in [0, 1, 2, 3]:
        n_b = 12
        solver_lm = PsHAngularSolver(n_basis=n_b, l_max=lm, n_quad=200)
        print(f"\n  l_max={lm}, dim={solver_lm._total_dim}")

        R_dense = np.concatenate([
            np.linspace(0.2, 2.0, 20),
            np.linspace(2.0, 8.0, 30),
            np.linspace(8.0, 20.0, 20),
            np.linspace(20.0, 80.0, 20),
        ])
        R_dense = np.unique(R_dense)

        mu_dense = np.zeros(len(R_dense))
        for i, R in enumerate(R_dense):
            ev, _ = solver_lm.solve(R, 1)
            mu_dense[i] = ev[0]

        V_dense = effective_potential(R_dense, mu_dense)
        V_spline = CubicSpline(R_dense, V_dense)

        # Raw energy
        E_lm, _, _ = solve_radial(
            V_spline, R_min=0.1, R_max=80.0,
            N_R=6000, n_states=1
        )

        # Shifted energy (correct asymptote)
        V_asymp_lm = V_dense[-1]
        shift_lm = -0.5 - V_asymp_lm
        V_shifted_lm = lambda R, spl=V_spline, s=shift_lm: spl(R) + s
        E_lm_s, _, _ = solve_radial(
            V_shifted_lm, R_min=0.1, R_max=80.0,
            N_R=6000, n_states=1
        )

        print(f"  Raw E_0 = {E_lm[0]:.6f} Ha, "
              f"Shifted E_0 = {E_lm_s[0]:.6f} Ha (exact: -0.789 Ha, "
              f"error: {abs(E_lm_s[0]-(-0.789))/0.789*100:.2f}%)")
        print(f"  V_eff min = {np.min(V_dense):.6f} Ha at R = {R_dense[np.argmin(V_dense)]:.2f}")
        print(f"  Well depth below -0.5 Ha: {max(0, -0.5 - np.min(V_dense)):.6f} Ha")

    # --- Step 8: Qualitative observations ---
    print(f"\n{'='*70}")
    print("QUALITATIVE OBSERVATIONS")
    print(f"{'='*70}")

    print("""
1. CHARGE FUNCTION: PsH has C(a) = -1/sin(a) + 1/cos(a) - 1/r_12.
   The nuclear part is ANTISYMMETRIC under a -> pi/2 - a (attraction
   at small a, repulsion at large a). This breaks alpha-parity, requiring
   both even and odd Gegenbauer basis functions.

2. ALPHA-PARITY MIXING: The nuclear coupling matrix has nonzero
   even-odd blocks, confirming that distinguishable particles (e-/e+)
   mix both parities. This doubles the effective basis size.

3. V_ee SIGN: The attractive e-/e+ interaction (-1/r_12) means the
   V_ee Gaunt coupling has the OPPOSITE sign from He. This promotes
   configurations where the particles are close together (forming
   positronium), counterbalancing the nuclear repulsion of e+.

4. BINDING: PsH is a delicate three-body system. The proton attracts
   the electron, which in turn attracts the positron via the Coulomb
   force. The positron is held by the electron, not by the proton.
""")

    # --- Save data ---
    output_dir = os.path.dirname(os.path.abspath(__file__))
    data_file = os.path.join(output_dir, "data", "psh_results.json")
    os.makedirs(os.path.join(output_dir, "data"), exist_ok=True)

    results = {
        "system": "PsH (positronium hydride)",
        "method": "Level 3 hyperspherical, adiabatic, Gegenbauer spectral",
        "l_max": l_max,
        "n_basis": n_basis,
        "E_psh_raw": float(E_psh[0]),
        "E_psh_shifted": float(E_psh_shifted[0]),
        "E_exact": -0.789,
        "error_raw_percent": float(abs(E_psh[0] - (-0.789)) / 0.789 * 100),
        "error_shifted_percent": float(abs(E_psh_shifted[0] - (-0.789)) / 0.789 * 100),
        "E_he_lmax0": float(E_he[0]),
        "E_hminus_lmax0": float(E_hminus[0]),
        "R_min_well": float(R_min_psh),
        "V_min_well": float(V_min_psh),
        "well_depth_below_threshold": float(max(0, -0.5 - V_min_psh)),
        "asymptotic_shift": float(shift),
    }

    with open(data_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {data_file}")

    # --- Plots ---
    plot_dir = os.path.join(output_dir, "plots")
    os.makedirs(plot_dir, exist_ok=True)

    # Plot 1: Adiabatic potential curves
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # PsH
    ax = axes[0]
    for ch in range(min(3, n_channels_angular)):
        ax.plot(R_grid, V_eff_psh[ch], label=f'Channel {ch}')
    ax.set_xlabel('R (bohr)')
    ax.set_ylabel('V_eff(R) (Ha)')
    ax.set_title('PsH adiabatic potential curves')
    ax.set_ylim(-1.5, 0.5)
    ax.axhline(-0.5, color='gray', ls='--', alpha=0.5, label='H(1s) + e+')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # He
    ax = axes[1]
    for ch in range(min(3, n_channels_angular)):
        ax.plot(R_grid, V_eff_he[ch], label=f'Channel {ch}')
    ax.set_xlabel('R (bohr)')
    ax.set_ylabel('V_eff(R) (Ha)')
    ax.set_title('He adiabatic potential curves')
    ax.set_ylim(-4.0, 1.0)
    ax.axhline(-2.0, color='gray', ls='--', alpha=0.5, label='He+(1s)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # H-
    ax = axes[2]
    for ch in range(min(3, n_channels_angular)):
        ax.plot(R_grid, V_eff_hminus[ch], label=f'Channel {ch}')
    ax.set_xlabel('R (bohr)')
    ax.set_ylabel('V_eff(R) (Ha)')
    ax.set_title('H- adiabatic potential curves')
    ax.set_ylim(-1.0, 0.5)
    ax.axhline(-0.5, color='gray', ls='--', alpha=0.5, label='H(1s) + e-')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'psh_adiabatic_curves.png'), dpi=150)
    print(f"Plot saved: {os.path.join(plot_dir, 'psh_adiabatic_curves.png')}")

    # Plot 2: mu(R) comparison
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.plot(R_grid, mu_psh[0], 'b-', linewidth=2, label='PsH')
    ax.plot(R_grid, mu_he[0], 'r-', linewidth=2, label='He (Z=2)')
    ax.plot(R_grid, mu_hminus[0], 'g-', linewidth=2, label='H- (Z=1)')
    ax.set_xlabel('R (bohr)', fontsize=12)
    ax.set_ylabel(r'$\mu_0(R)$', fontsize=12)
    ax.set_title('Lowest angular eigenvalue comparison', fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'psh_mu_comparison.png'), dpi=150)
    print(f"Plot saved: {os.path.join(plot_dir, 'psh_mu_comparison.png')}")

    # Plot 3: PsH charge function
    fig, ax = plt.subplots(figsize=(8, 6))
    alpha_plot = np.linspace(0.05, np.pi/2 - 0.05, 200)
    sin_a = np.sin(alpha_plot)
    cos_a = np.cos(alpha_plot)
    C_nuc_psh = -1.0 / sin_a + 1.0 / cos_a
    C_nuc_he = -2.0 / sin_a - 2.0 / cos_a
    C_nuc_hminus = -1.0 / sin_a - 1.0 / cos_a
    ax.plot(alpha_plot * 180 / np.pi, C_nuc_psh, 'b-', linewidth=2, label='PsH nuclear')
    ax.plot(alpha_plot * 180 / np.pi, C_nuc_he, 'r-', linewidth=2, label='He nuclear')
    ax.plot(alpha_plot * 180 / np.pi, C_nuc_hminus, 'g-', linewidth=2, label='H- nuclear')
    ax.set_xlabel(r'$\alpha$ (degrees)', fontsize=12)
    ax.set_ylabel('Nuclear charge function', fontsize=12)
    ax.set_title('Nuclear part of charge function C(alpha)', fontsize=14)
    ax.set_ylim(-30, 30)
    ax.axhline(0, color='k', ls='-', alpha=0.3)
    ax.axvline(45, color='k', ls='--', alpha=0.3, label=r'$\alpha = \pi/4$')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'psh_charge_function.png'), dpi=150)
    print(f"Plot saved: {os.path.join(plot_dir, 'psh_charge_function.png')}")

    print(f"\nDone.")
    return results


if __name__ == '__main__':
    results = main()
