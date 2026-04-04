"""
Sturmian CI solver — Coulomb Sturmian basis with variable orbital exponents.
=============================================================================

Track BU-1: Generalized Sturmian investigation for quantum simulation.

Coulomb Sturmians are hydrogenic orbitals with Z_eff = n*k where k is a
common scaling parameter. All Sturmians share the same exponential decay
rate exp(-k*r) and the same one-particle energy -k²/2, regardless of n.

Key properties:
  - Weighted orthogonality: ⟨a|1/r|b⟩ = (k/n_a) δ_{ab}
  - Nuclear attraction is diagonal in the Sturmian basis
  - Standard overlap ⟨a|b⟩ ≠ δ_{ab} for different n (same l)
  - Gaunt angular selection rules are preserved (angular part unchanged)

The solver Löwdin-orthogonalizes the Sturmian basis, transforms all
integrals, and runs standard FCI via Slater-Condon rules.

Author: GeoVac Development Team
Date: April 2026
"""

import numpy as np
from scipy.integrate import quad
from scipy.special import genlaguerre
from scipy.linalg import eigh, sqrtm, inv
from scipy.optimize import minimize_scalar
from itertools import combinations
from math import factorial as int_factorial
from typing import Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Hydrogenic radial wavefunctions
# ---------------------------------------------------------------------------

def hydrogenic_radial(r: np.ndarray, n: int, l: int, Z_eff: float) -> np.ndarray:
    """
    Normalized hydrogenic radial wavefunction R_{nl}(r; Z_eff).

    R_{nl} = N * (2*Z_eff*r/n)^l * exp(-Z_eff*r/n) * L_{n-l-1}^{2l+1}(2*Z_eff*r/n)

    Normalized so that ∫₀^∞ |R_{nl}|² r² dr = 1.
    """
    nr = n - l - 1
    rho = 2.0 * Z_eff * np.asarray(r, dtype=float) / n

    N_sq = (2.0 * Z_eff / n) ** 3 * int_factorial(nr) / (
        2.0 * n * int_factorial(n + l)
    )
    N = np.sqrt(N_sq)

    L_poly = genlaguerre(nr, 2 * l + 1)(rho)
    return N * rho ** l * np.exp(-rho / 2.0) * L_poly


# ---------------------------------------------------------------------------
# One-particle integral engine
# ---------------------------------------------------------------------------

def _radial_overlap(n1: int, l1: int, Z1: float,
                    n2: int, l2: int, Z2: float) -> float:
    """Overlap integral ∫ R_{n1,l1}(r;Z1) R_{n2,l2}(r;Z2) r² dr."""
    if l1 != l2:
        return 0.0

    r_max = max(n1 / Z1, n2 / Z2) * 60.0
    val, _ = quad(
        lambda r: (hydrogenic_radial(r, n1, l1, Z1)
                   * hydrogenic_radial(r, n2, l2, Z2) * r ** 2),
        0, r_max, limit=300
    )
    return val


def _radial_1_over_r(n1: int, l1: int, Z1: float,
                     n2: int, l2: int, Z2: float) -> float:
    """Matrix element ∫ R_{n1,l1}(r;Z1) (1/r) R_{n2,l2}(r;Z2) r² dr."""
    if l1 != l2:
        return 0.0

    r_max = max(n1 / Z1, n2 / Z2) * 60.0
    val, _ = quad(
        lambda r: (hydrogenic_radial(r, n1, l1, Z1)
                   * hydrogenic_radial(r, n2, l2, Z2) * r),
        0, r_max, limit=300
    )
    return val


# ---------------------------------------------------------------------------
# Two-electron Slater R^k integrals
# ---------------------------------------------------------------------------

def _slater_rk(n1: int, l1: int, Z1: float,
               n2: int, l2: int, Z2: float,
               n3: int, l3: int, Z3: float,
               n4: int, l4: int, Z4: float,
               k: int) -> float:
    """
    Slater R^k integral with variable Z_eff per orbital.

    R^k(1234) = ∫∫ R_1(r₁) R_3(r₁) (r_<^k / r_>^{k+1}) R_2(r₂) R_4(r₂) r₁² r₂² dr₁ dr₂

    Orbitals 1,3 share electron 1; orbitals 2,4 share electron 2.
    """
    # Use grid-based approach for robustness
    r_max1 = max(n1 / Z1, n3 / Z3) * 50.0
    r_max2 = max(n2 / Z2, n4 / Z4) * 50.0
    r_max = max(r_max1, r_max2)

    n_grid = 500
    r_grid = np.linspace(1e-12, r_max, n_grid)
    dr = r_grid[1] - r_grid[0]

    # Tabulate radial products
    f1 = hydrogenic_radial(r_grid, n1, l1, Z1) * hydrogenic_radial(r_grid, n3, l3, Z3) * r_grid ** 2
    f2 = hydrogenic_radial(r_grid, n2, l2, Z2) * hydrogenic_radial(r_grid, n4, l4, Z4) * r_grid ** 2

    # Build y_k(r1) = ∫ f2(r2) (r_<^k / r_>^{k+1}) r2² dr2
    # Split into r2 < r1 and r2 > r1
    yk = np.zeros(n_grid)
    for i in range(n_grid):
        r1 = r_grid[i]
        if r1 < 1e-30:
            continue
        # r2 < r1: r_<^k/r_>^{k+1} = r2^k / r1^{k+1}
        inner_lo = np.sum(f2[:i+1] * r_grid[:i+1] ** k) * dr / r1 ** (k + 1)
        # r2 > r1: r_<^k/r_>^{k+1} = r1^k / r2^{k+1}
        inner_hi = np.sum(f2[i:] / r_grid[i:] ** (k + 1)) * dr * r1 ** k
        yk[i] = inner_lo + inner_hi

    return float(np.sum(f1 * yk) * dr)


def _slater_rk_fast(n1: int, l1: int, Z1: float,
                    n2: int, l2: int, Z2: float,
                    n3: int, l3: int, Z3: float,
                    n4: int, l4: int, Z4: float,
                    k: int) -> float:
    """
    Fast Slater R^k integral using vectorized cumulative sums.

    Equivalent to _slater_rk but O(n_grid) instead of O(n_grid²),
    giving ~500x speedup for the generalized Sturmian solver where
    many more R^k integrals are needed (mixed Z_eff pairs).
    """
    r_max1 = max(n1 / Z1, n3 / Z3) * 50.0
    r_max2 = max(n2 / Z2, n4 / Z4) * 50.0
    r_max = max(r_max1, r_max2)

    n_grid = 500
    r_grid = np.linspace(1e-12, r_max, n_grid)
    dr = r_grid[1] - r_grid[0]

    f1 = (hydrogenic_radial(r_grid, n1, l1, Z1)
          * hydrogenic_radial(r_grid, n3, l3, Z3) * r_grid ** 2)
    f2 = (hydrogenic_radial(r_grid, n2, l2, Z2)
          * hydrogenic_radial(r_grid, n4, l4, Z4) * r_grid ** 2)

    # yk[i] = ∫₀^{r_i} f2(r2) r2^k dr2 / r_i^{k+1}
    #        + ∫_{r_i}^∞ f2(r2) / r2^{k+1} dr2 × r_i^k
    cumsum_lo = np.cumsum(f2 * r_grid ** k) * dr
    cumsum_hi = np.cumsum((f2 / r_grid ** (k + 1))[::-1])[::-1] * dr

    yk = cumsum_lo / r_grid ** (k + 1) + cumsum_hi * r_grid ** k
    yk[0] = 0.0  # r ≈ 0, skip (matches _slater_rk behavior)

    return float(np.sum(f1 * yk) * dr)


# ---------------------------------------------------------------------------
# Wigner 3j symbol and Gaunt coefficients
# ---------------------------------------------------------------------------

def _wigner3j(j1: int, j2: int, j3: int,
              m1: int, m2: int, m3: int) -> float:
    """Wigner 3j symbol for integer arguments (Racah formula)."""
    if m1 + m2 + m3 != 0:
        return 0.0
    if abs(j1 - j2) > j3 or j3 > j1 + j2:
        return 0.0
    if abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3:
        return 0.0

    def _tri(a: int, b: int, c: int) -> float:
        return (int_factorial(a + b - c) * int_factorial(a - b + c)
                * int_factorial(-a + b + c)
                / int_factorial(a + b + c + 1))

    pre = ((-1) ** (j1 - j2 - m3)
           * np.sqrt(_tri(j1, j2, j3)
                     * int_factorial(j1 + m1) * int_factorial(j1 - m1)
                     * int_factorial(j2 + m2) * int_factorial(j2 - m2)
                     * int_factorial(j3 + m3) * int_factorial(j3 - m3)))

    t_min = max(0, j2 - j3 - m1, j1 - j3 + m2)
    t_max = min(j1 + j2 - j3, j1 - m1, j2 + m2)

    s = 0.0
    for t in range(t_min, t_max + 1):
        s += ((-1) ** t
              / (int_factorial(t)
                 * int_factorial(j1 + j2 - j3 - t)
                 * int_factorial(j1 - m1 - t)
                 * int_factorial(j2 + m2 - t)
                 * int_factorial(j3 - j2 + m1 + t)
                 * int_factorial(j3 - j1 - m2 + t)))

    return pre * s


def _ck_coefficient(la: int, ma: int, lc: int, mc: int, k: int) -> float:
    """Gaunt angular coupling coefficient c^k(l,m,l',m')."""
    q = mc - ma
    pre = ((-1) ** ma * np.sqrt((2 * la + 1) * (2 * lc + 1)))
    w1 = _wigner3j(la, k, lc, 0, 0, 0)
    if abs(w1) < 1e-15:
        return 0.0
    w2 = _wigner3j(la, k, lc, -ma, q, mc)
    return pre * w1 * w2


# ---------------------------------------------------------------------------
# SturmianCI solver
# ---------------------------------------------------------------------------

class SturmianCI:
    """
    Variational CI solver using Coulomb Sturmian basis.

    Each orbital (n,l,m) uses Z_eff = n*k where k is a common scaling
    parameter. After Löwdin orthogonalization, standard FCI is performed.

    Parameters
    ----------
    Z : int
        Nuclear charge.
    n_electrons : int
        Number of electrons.
    max_n : int
        Maximum principal quantum number.
    """

    def __init__(self, Z: int, n_electrons: int, max_n: int) -> None:
        self.Z = Z
        self.n_electrons = n_electrons
        self.max_n = max_n

        # Build spatial states (same labeling as GeometricLattice)
        self.states: List[Tuple[int, int, int]] = []
        for n in range(1, max_n + 1):
            for ll in range(n):
                for m in range(-ll, ll + 1):
                    self.states.append((n, ll, m))
        self.n_spatial = len(self.states)
        self.n_sp = 2 * self.n_spatial

        # Build SD basis
        self.sd_basis = list(combinations(range(self.n_sp), n_electrons))
        self.n_sd = len(self.sd_basis)
        self._sd_index = {sd: i for i, sd in enumerate(self.sd_basis)}

    # ------------------------------------------------------------------
    # One-particle matrices
    # ------------------------------------------------------------------

    def _build_overlap(self, k: float) -> np.ndarray:
        """Spatial overlap matrix S[i,j] in the Sturmian basis."""
        n_sp = self.n_spatial
        S = np.eye(n_sp)
        for i in range(n_sp):
            ni, li, mi = self.states[i]
            Zi = ni * k
            for j in range(i + 1, n_sp):
                nj, lj, mj = self.states[j]
                if li != lj or mi != mj:
                    continue
                Zj = nj * k
                val = _radial_overlap(ni, li, Zi, nj, lj, Zj)
                S[i, j] = val
                S[j, i] = val
        return S

    def _build_h1_sturmian(self, k: float, S: np.ndarray) -> np.ndarray:
        """
        One-body Hamiltonian h1 = T + V_nuc in the Sturmian basis.

        Uses the analytical formula derived from the Sturmian equation:
          h1_{ab} = (k² - Z*k/n_a) δ_{ab} - k²/2 * S_{ab}

        Diagonal: k²/2 - Z*k/n_a
        Off-diagonal: -k²/2 * S_{ab}
        """
        n_sp = self.n_spatial
        Z = float(self.Z)

        # Start with -k²/2 * S
        h1 = -k ** 2 / 2.0 * S.copy()

        # Add diagonal correction: (k² - Z*k/n_a) on top of -k²/2
        for i in range(n_sp):
            ni = self.states[i][0]
            h1[i, i] += k ** 2 - Z * k / ni

        return h1

    # ------------------------------------------------------------------
    # Two-electron integrals
    # ------------------------------------------------------------------

    def _build_eri(self, k: float) -> Dict[Tuple[int, int, int, int], float]:
        """
        Build full ERI table ⟨ab|g|cd⟩ using Gaunt coefficients and
        Slater R^k integrals in the Sturmian basis.
        """
        states = self.states
        n_sp = self.n_spatial

        # Step 1: Pre-compute c^k angular coefficients
        ck_table: Dict[Tuple[int, int, int], float] = {}
        for a in range(n_sp):
            la, ma = states[a][1], states[a][2]
            for c in range(n_sp):
                lc, mc = states[c][1], states[c][2]
                k_max = la + lc
                for kk in range(0, k_max + 1):
                    if (la + lc + kk) % 2 != 0:
                        continue
                    val = _ck_coefficient(la, ma, lc, mc, kk)
                    if abs(val) > 1e-15:
                        ck_table[(a, c, kk)] = val

        # Group by (a,c)
        ac_k_map: Dict[Tuple[int, int], List[Tuple[int, float]]] = {}
        for (a, c, kk), val in ck_table.items():
            key = (a, c)
            if key not in ac_k_map:
                ac_k_map[key] = []
            ac_k_map[key].append((kk, val))

        # Step 2: Pre-compute needed R^k Slater integrals
        # Key: (na,la,nb,lb,nc,lc,nd,ld,k) -> value
        rk_cache: Dict[tuple, float] = {}

        def _get_rk(na: int, la: int, nb: int, lb: int,
                    nc: int, lc: int, nd: int, ld: int, kk: int) -> float:
            key = (na, la, nb, lb, nc, lc, nd, ld, kk)
            if key not in rk_cache:
                Za, Zb, Zc, Zd = na * k, nb * k, nc * k, nd * k
                rk_cache[key] = _slater_rk(na, la, Za, nb, lb, Zb,
                                           nc, lc, Zc, nd, ld, Zd, kk)
            return rk_cache[key]

        # Step 3: Assemble ERI table
        eri: Dict[Tuple[int, int, int, int], float] = {}
        for (a, c), ck_ac_list in ac_k_map.items():
            na, la, ma = states[a]
            nc, lc, mc = states[c]
            for (b, d), ck_bd_list in ac_k_map.items():
                nb, lb, mb = states[b]
                nd, ld, md = states[d]

                # m-selection rule
                if ma + mb != mc + md:
                    continue

                val = 0.0
                for k_ac, c_ac in ck_ac_list:
                    for k_bd, c_bd in ck_bd_list:
                        if k_ac != k_bd:
                            continue
                        rk_val = _get_rk(na, la, nb, lb, nc, lc, nd, ld, k_ac)
                        val += c_ac * c_bd * rk_val

                if abs(val) > 1e-15:
                    eri[(a, b, c, d)] = val

        return eri

    # ------------------------------------------------------------------
    # Löwdin orthogonalization and integral transformation
    # ------------------------------------------------------------------

    def _lowdin_transform(self, S: np.ndarray, h1: np.ndarray,
                          eri: Dict) -> Tuple[np.ndarray, np.ndarray]:
        """
        Löwdin-orthogonalize the basis and transform integrals.

        Returns (h1_ortho, eri_4d_ortho) in the orthonormalized basis.
        """
        n_sp = self.n_spatial

        # S^{-1/2} via eigendecomposition
        eigvals, eigvecs = eigh(S)
        # Ensure positive definite
        eigvals = np.maximum(eigvals, 1e-14)
        X = eigvecs @ np.diag(1.0 / np.sqrt(eigvals)) @ eigvecs.T

        # Transform h1
        h1_ortho = X @ h1 @ X

        # Transform ERIs: 4-index transformation
        eri_4d = np.zeros((n_sp, n_sp, n_sp, n_sp))
        for (a, b, c, d), val in eri.items():
            eri_4d[a, b, c, d] = val

        # 4-step transformation: a' = Σ_a X[a',a] ...
        tmp = np.einsum('pa,abcd->pbcd', X, eri_4d)
        tmp = np.einsum('qb,pbcd->pqcd', X, tmp)
        tmp = np.einsum('rc,pqcd->pqrd', X, tmp)
        eri_4d_ortho = np.einsum('sd,pqrd->pqrs', X, tmp)

        return h1_ortho, eri_4d_ortho

    # ------------------------------------------------------------------
    # FCI Hamiltonian assembly (Slater-Condon rules)
    # ------------------------------------------------------------------

    def _build_fci_hamiltonian(self, h1: np.ndarray,
                               eri_4d: np.ndarray) -> np.ndarray:
        """
        Build the FCI Hamiltonian matrix using Slater-Condon rules.

        Uses orthonormalized integrals (h1 and eri_4d).
        """
        n_sd = self.n_sd
        H = np.zeros((n_sd, n_sd))

        for I, sd_I in enumerate(self.sd_basis):
            # Diagonal
            val = 0.0
            for p in sd_I:
                sp_p = p >> 1
                val += h1[sp_p, sp_p]
            for idx_p in range(len(sd_I)):
                for idx_q in range(idx_p + 1, len(sd_I)):
                    p, q = sd_I[idx_p], sd_I[idx_q]
                    sp_p, sig_p = p >> 1, p & 1
                    sp_q, sig_q = q >> 1, q & 1
                    # Coulomb
                    val += eri_4d[sp_p, sp_q, sp_p, sp_q]
                    # Exchange (same spin only)
                    if sig_p == sig_q:
                        val -= eri_4d[sp_p, sp_q, sp_q, sp_p]
            H[I, I] = val

            # Off-diagonal: single and double excitations
            occ_set = set(sd_I)
            for J in range(I + 1, n_sd):
                sd_J = self.sd_basis[J]

                # Find differing orbitals
                diff_I = [x for x in sd_I if x not in sd_J]
                diff_J = [x for x in sd_J if x not in sd_I]

                n_diff = len(diff_I)

                if n_diff == 1:
                    # Single excitation: p -> r
                    p = diff_I[0]
                    r = diff_J[0]
                    sp_p, sig_p = p >> 1, p & 1
                    sp_r, sig_r = r >> 1, r & 1

                    if sig_p != sig_r:
                        continue

                    # Phase
                    phase = self._compute_phase(sd_I, sd_J, diff_I, diff_J)

                    # One-body
                    val = h1[sp_p, sp_r]

                    # Two-body: sum over common occupied
                    common = [x for x in sd_I if x in sd_J]
                    for q in common:
                        sp_q, sig_q = q >> 1, q & 1
                        val += eri_4d[sp_p, sp_q, sp_r, sp_q]
                        if sig_p == sig_q:
                            val -= eri_4d[sp_p, sp_q, sp_q, sp_r]

                    H[I, J] = phase * val
                    H[J, I] = phase * val

                elif n_diff == 2:
                    # Double excitation: p,q -> r,s
                    p, q = diff_I[0], diff_I[1]
                    r, s = diff_J[0], diff_J[1]
                    sp_p, sig_p = p >> 1, p & 1
                    sp_q, sig_q = q >> 1, q & 1
                    sp_r, sig_r = r >> 1, r & 1
                    sp_s, sig_s = s >> 1, s & 1

                    # Spin conservation
                    if sig_p + sig_q != sig_r + sig_s:
                        continue

                    phase = self._compute_phase(sd_I, sd_J, diff_I, diff_J)

                    # ⟨pq|rs⟩ - δ_spin ⟨pq|sr⟩
                    val = 0.0
                    if sig_p == sig_r and sig_q == sig_s:
                        val += eri_4d[sp_p, sp_q, sp_r, sp_s]
                    if sig_p == sig_s and sig_q == sig_r:
                        val -= eri_4d[sp_p, sp_q, sp_s, sp_r]

                    H[I, J] = phase * val
                    H[J, I] = phase * val

        return H

    @staticmethod
    def _compute_phase(sd_I: tuple, sd_J: tuple,
                       diff_I: list, diff_J: list) -> float:
        """Compute fermionic phase for excitation sd_I -> sd_J."""
        # Build intermediate determinant by removing diff_I orbitals
        # and inserting diff_J orbitals, counting transpositions
        sd_list = list(sd_I)
        n_swaps = 0

        for old, new in zip(sorted(diff_I), sorted(diff_J)):
            idx_old = sd_list.index(old)
            sd_list.pop(idx_old)
            n_swaps += idx_old

            # Find insertion point for new orbital
            ins_idx = 0
            for x in sd_list:
                if x < new:
                    ins_idx += 1
                else:
                    break
            sd_list.insert(ins_idx, new)
            n_swaps += ins_idx

        return (-1.0) ** n_swaps

    # ------------------------------------------------------------------
    # Solver
    # ------------------------------------------------------------------

    def solve(self, k: float, verbose: bool = False) -> Dict:
        """
        Solve CI at fixed k.

        Returns dict with energy, coefficients, orbital exponents, etc.
        """
        if verbose:
            print(f"  Sturmian CI: k={k:.4f}, Z={self.Z}, "
                  f"max_n={self.max_n}, n_spatial={self.n_spatial}, "
                  f"n_sd={self.n_sd}")

        # Build one-particle matrices
        S = self._build_overlap(k)
        h1 = self._build_h1_sturmian(k, S)

        # Check positive definiteness
        eigvals_S = np.linalg.eigvalsh(S)
        if eigvals_S[0] < 1e-8:
            if verbose:
                print(f"  WARNING: overlap near-singular, min eigenvalue = {eigvals_S[0]:.2e}")
            return {'energy': np.inf, 'k': k, 'converged': False}

        # Build ERIs
        eri = self._build_eri(k)

        # Löwdin orthogonalization
        h1_ortho, eri_4d_ortho = self._lowdin_transform(S, h1, eri)

        # Build and diagonalize FCI Hamiltonian
        H_fci = self._build_fci_hamiltonian(h1_ortho, eri_4d_ortho)
        eigvals, eigvecs = eigh(H_fci)

        E_gs = eigvals[0]
        c_gs = eigvecs[:, 0]

        # Orbital exponents for each state
        z_effs = {(n, l, m): n * k for n, l, m in self.states}

        # V_ee sparsity
        n_eri_total = self.n_spatial ** 4
        n_eri_nonzero = len(eri)
        vee_sparsity = n_eri_nonzero / n_eri_total if n_eri_total > 0 else 0

        if verbose:
            print(f"  E_gs = {E_gs:.6f} Ha")
            print(f"  V_ee nonzero fraction: {vee_sparsity:.4f} "
                  f"({n_eri_nonzero}/{n_eri_total})")

        return {
            'energy': E_gs,
            'k': k,
            'coefficients': c_gs,
            'eigenvalues': eigvals,
            'z_effs': z_effs,
            'overlap_matrix': S,
            'h1_sturmian': h1,
            'h1_ortho': h1_ortho,
            'eri_4d_ortho': eri_4d_ortho,
            'eri_count': n_eri_nonzero,
            'eri_total': n_eri_total,
            'vee_sparsity': vee_sparsity,
            'converged': True,
        }

    def optimize_k(self, k_range: Tuple[float, float] = (0.5, 5.0),
                   n_scan: int = 20,
                   verbose: bool = False) -> Dict:
        """
        Find optimal k that minimizes the ground state energy.

        Uses a coarse scan followed by Brent refinement.
        """
        if verbose:
            print(f"\nSturmian CI optimization: Z={self.Z}, "
                  f"n_e={self.n_electrons}, max_n={self.max_n}")
            print(f"  Scanning k in [{k_range[0]:.2f}, {k_range[1]:.2f}]")

        # Coarse scan
        k_values = np.linspace(k_range[0], k_range[1], n_scan)
        energies = []
        for kk in k_values:
            result = self.solve(kk)
            energies.append(result['energy'])
            if verbose:
                print(f"    k={kk:.3f}: E={result['energy']:.6f}")

        # Find bracket
        idx_min = int(np.argmin(energies))
        if idx_min == 0:
            k_lo, k_hi = k_values[0], k_values[1]
        elif idx_min == n_scan - 1:
            k_lo, k_hi = k_values[-2], k_values[-1]
        else:
            k_lo, k_hi = k_values[idx_min - 1], k_values[idx_min + 1]

        # Brent refinement
        def _objective(kk: float) -> float:
            r = self.solve(kk)
            return r['energy'] if r['converged'] else 1e10

        opt = minimize_scalar(_objective, bounds=(k_lo, k_hi), method='bounded')
        k_opt = opt.x

        # Final solve at optimal k
        result = self.solve(k_opt, verbose=verbose)
        result['k_scan'] = list(zip(k_values.tolist(), energies))

        if verbose:
            print(f"\n  Optimal k = {k_opt:.6f}")
            print(f"  E_opt = {result['energy']:.6f} Ha")

        return result


# ---------------------------------------------------------------------------
# Standard FCI baseline (for comparison)
# ---------------------------------------------------------------------------

class StandardFCI:
    """
    Standard FCI with hydrogenic basis (all Z_eff = Z_opt).

    This is the baseline: all orbitals use the same effective nuclear charge,
    optimized variationally.
    """

    def __init__(self, Z: int, n_electrons: int, max_n: int) -> None:
        self.Z = Z
        self.n_electrons = n_electrons
        self.max_n = max_n

        self.states: List[Tuple[int, int, int]] = []
        for n in range(1, max_n + 1):
            for ll in range(n):
                for m in range(-ll, ll + 1):
                    self.states.append((n, ll, m))
        self.n_spatial = len(self.states)
        self.n_sp = 2 * self.n_spatial

        self.sd_basis = list(combinations(range(self.n_sp), n_electrons))
        self.n_sd = len(self.sd_basis)
        self._sd_index = {sd: i for i, sd in enumerate(self.sd_basis)}

    def _build_h1(self, Z_eff: float) -> np.ndarray:
        """
        One-body Hamiltonian h1 = T + V_nuc in hydrogenic basis with Z_eff.

        Diagonal: ⟨nl|T+V_nuc|nl⟩ = Z_eff²/(2n²) - Z·Z_eff/n²
                = Z_eff(Z_eff - 2Z)/(2n²)

        Off-diagonal (same l, different n):
        ⟨n'l|T+V_nuc|nl⟩ = (Z_eff - Z) × ⟨n'l|1/r|nl⟩
        computed numerically.
        """
        Z = float(self.Z)
        n_sp = self.n_spatial
        h1 = np.zeros((n_sp, n_sp))

        for i, (ni, li, mi) in enumerate(self.states):
            # Diagonal: virial theorem gives T = Z_eff²/(2n²),
            # V_nuc = -Z·⟨1/r⟩ = -Z·Z_eff/n²
            h1[i, i] = Z_eff * (Z_eff - 2.0 * Z) / (2.0 * ni ** 2)

            # Off-diagonal: (Z_eff - Z) × ⟨n'l|1/r|nl⟩
            if abs(Z_eff - Z) > 1e-12:
                for j in range(i + 1, n_sp):
                    nj, lj, mj = self.states[j]
                    if li != lj or mi != mj:
                        continue
                    inv_r = _radial_1_over_r(ni, li, Z_eff, nj, lj, Z_eff)
                    val = (Z_eff - Z) * inv_r
                    h1[i, j] = val
                    h1[j, i] = val

        return h1

    def _build_eri(self, Z_eff: float) -> np.ndarray:
        """Build ERI table with uniform Z_eff for all orbitals."""
        states = self.states
        n_sp = self.n_spatial

        # Compute R^k integrals
        rk_cache: Dict[tuple, float] = {}
        unique_nl = sorted(set((n, l) for n, l, m in states))

        for n1, l1 in unique_nl:
            for n2, l2 in unique_nl:
                for n3, l3 in unique_nl:
                    for n4, l4 in unique_nl:
                        k_max = min(l1 + l3, l2 + l4)
                        for kk in range(0, k_max + 1):
                            if (l1 + l3 + kk) % 2 != 0:
                                continue
                            if (l2 + l4 + kk) % 2 != 0:
                                continue
                            key = (n1, l1, n2, l2, n3, l3, n4, l4, kk)
                            if key not in rk_cache:
                                rk_cache[key] = _slater_rk(
                                    n1, l1, Z_eff, n2, l2, Z_eff,
                                    n3, l3, Z_eff, n4, l4, Z_eff, kk
                                )

        # Assemble ERI
        ck_table: Dict[Tuple[int, int, int], float] = {}
        for a in range(n_sp):
            la, ma = states[a][1], states[a][2]
            for c in range(n_sp):
                lc, mc = states[c][1], states[c][2]
                for kk in range(0, la + lc + 1):
                    if (la + lc + kk) % 2 != 0:
                        continue
                    val = _ck_coefficient(la, ma, lc, mc, kk)
                    if abs(val) > 1e-15:
                        ck_table[(a, c, kk)] = val

        ac_k_map: Dict[Tuple[int, int], List[Tuple[int, float]]] = {}
        for (a, c, kk), val in ck_table.items():
            key = (a, c)
            if key not in ac_k_map:
                ac_k_map[key] = []
            ac_k_map[key].append((kk, val))

        eri_4d = np.zeros((n_sp, n_sp, n_sp, n_sp))
        for (a, c), ck_ac_list in ac_k_map.items():
            na, la, ma = states[a]
            nc, lc, mc = states[c]
            for (b, d), ck_bd_list in ac_k_map.items():
                nb, lb, mb = states[b]
                nd, ld, md = states[d]
                if ma + mb != mc + md:
                    continue
                val = 0.0
                for k_ac, c_ac in ck_ac_list:
                    for k_bd, c_bd in ck_bd_list:
                        if k_ac != k_bd:
                            continue
                        rk_key = (na, la, nb, lb, nc, lc, nd, ld, k_ac)
                        rk_val = rk_cache.get(rk_key, 0.0)
                        val += c_ac * c_bd * rk_val
                if abs(val) > 1e-15:
                    eri_4d[a, b, c, d] = val

        return eri_4d

    def solve(self, Z_eff: float, verbose: bool = False) -> Dict:
        """Solve standard FCI at fixed Z_eff."""
        h1 = self._build_h1(Z_eff)
        eri_4d = self._build_eri(Z_eff)

        # Build FCI Hamiltonian (reuse SturmianCI's method via standalone)
        solver = SturmianCI.__new__(SturmianCI)
        solver.sd_basis = self.sd_basis
        solver.n_sd = self.n_sd
        solver._sd_index = self._sd_index
        solver.states = self.states
        solver.n_spatial = self.n_spatial

        H_fci = solver._build_fci_hamiltonian(h1, eri_4d)
        eigvals, eigvecs = eigh(H_fci)

        n_eri_nonzero = int(np.count_nonzero(eri_4d))

        return {
            'energy': eigvals[0],
            'Z_eff': Z_eff,
            'eigenvalues': eigvals,
            'eri_count': n_eri_nonzero,
        }

    def optimize_zeff(self, z_range: Tuple[float, float] = (0.5, 3.0),
                      n_scan: int = 15,
                      verbose: bool = False) -> Dict:
        """Find optimal Z_eff."""
        z_values = np.linspace(z_range[0], z_range[1], n_scan)
        energies = []
        for zz in z_values:
            r = self.solve(zz)
            energies.append(r['energy'])

        idx_min = int(np.argmin(energies))
        if idx_min == 0:
            z_lo, z_hi = z_values[0], z_values[1]
        elif idx_min == n_scan - 1:
            z_lo, z_hi = z_values[-2], z_values[-1]
        else:
            z_lo, z_hi = z_values[idx_min - 1], z_values[idx_min + 1]

        opt = minimize_scalar(
            lambda z: self.solve(z)['energy'],
            bounds=(z_lo, z_hi), method='bounded'
        )

        result = self.solve(opt.x, verbose=verbose)
        if verbose:
            print(f"  Standard FCI: Z_eff_opt={opt.x:.6f}, E={result['energy']:.6f}")
        return result


# ---------------------------------------------------------------------------
# Generalized Sturmian CI (Track BU-1b)
# ---------------------------------------------------------------------------

class GeneralizedSturmianCI:
    """
    Generalized Sturmian CI with configuration-dependent orbital exponents.

    Each Slater determinant ν gets scaling β_ν from the isoenergetic condition:
    all configurations have total one-body energy E_trial. Orbitals within
    configuration ν use Z_eff = β_ν · Z.

    Self-consistency loop: E_trial → β_ν → (H, S) → diag → E_new → iterate.

    For 2 electrons, SD-level matrix elements use:
      S_μν = det(M)  where M_ij = ⟨a_i(Z_μ)|b_j(Z_ν)⟩
      H1_μν = Σ_ij h(a_i,b_j) · cofactor(M,i,j)
      V_ee_μν = g(a₁,a₂;b₁,b₂) - g(a₁,a₂;b₂,b₁)

    Parameters
    ----------
    Z : int
        Nuclear charge.
    n_electrons : int
        Number of electrons (currently supports 2).
    max_n : int
        Maximum principal quantum number.
    """

    def __init__(self, Z: int, n_electrons: int, max_n: int) -> None:
        if n_electrons != 2:
            raise NotImplementedError(
                "GeneralizedSturmianCI currently supports 2 electrons only"
            )
        self.Z = Z
        self.n_electrons = n_electrons
        self.max_n = max_n

        self.states: List[Tuple[int, int, int]] = []
        for n in range(1, max_n + 1):
            for ll in range(n):
                for m in range(-ll, ll + 1):
                    self.states.append((n, ll, m))
        self.n_spatial = len(self.states)
        self.n_sp = 2 * self.n_spatial

        self.sd_basis = list(combinations(range(self.n_sp), n_electrons))
        self.n_sd = len(self.sd_basis)

        # Precompute angular coefficients (β-independent)
        self._ac_k_map: Dict[Tuple[int, int], List[Tuple[int, float]]] = {}
        self._eri_angular_count = 0
        self._precompute_gaunt()

    def _precompute_gaunt(self) -> None:
        """Precompute Gaunt angular coupling coefficients."""
        states = self.states
        n_sp = self.n_spatial

        ck_table: Dict[Tuple[int, int, int], float] = {}
        for a in range(n_sp):
            la, ma = states[a][1], states[a][2]
            for c in range(n_sp):
                lc, mc = states[c][1], states[c][2]
                for kk in range(0, la + lc + 1):
                    if (la + lc + kk) % 2 != 0:
                        continue
                    val = _ck_coefficient(la, ma, lc, mc, kk)
                    if abs(val) > 1e-15:
                        ck_table[(a, c, kk)] = val

        for (a, c, kk), val in ck_table.items():
            key = (a, c)
            if key not in self._ac_k_map:
                self._ac_k_map[key] = []
            self._ac_k_map[key].append((kk, val))

        # Count angular ERI nonzeros (for sparsity comparison with BU-1)
        eri_nonzero: set = set()
        for (a, c), ck_ac in self._ac_k_map.items():
            for (b, d), ck_bd in self._ac_k_map.items():
                if states[a][2] + states[b][2] != states[c][2] + states[d][2]:
                    continue
                for k1, _ in ck_ac:
                    for k2, _ in ck_bd:
                        if k1 == k2:
                            eri_nonzero.add((a, b, c, d))
                            break
                    else:
                        continue
                    break
        self._eri_angular_count = len(eri_nonzero)

    def _compute_betas(self, E_trial: float) -> np.ndarray:
        """
        Compute β_ν for each SD from the isoenergetic condition.

        E_trial = -β²Z²/2 × Σ_i 1/n_i²  →  β = √(-2E/(Z²·Σ 1/n²))
        """
        betas = np.zeros(self.n_sd)
        Z = float(self.Z)
        for I, sd in enumerate(self.sd_basis):
            inv_n_sq_sum = 0.0
            for sp_idx in sd:
                n_val = self.states[sp_idx >> 1][0]
                inv_n_sq_sum += 1.0 / n_val ** 2
            betas[I] = np.sqrt(-2.0 * E_trial / (Z ** 2 * inv_n_sq_sum))
        return betas

    def _build_sd_matrices(self, betas: np.ndarray,
                           verbose: bool = False
                           ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Build SD-level overlap S and Hamiltonian H matrices.

        Uses determinantal formulas for 2-electron Slater determinants
        with configuration-dependent Z_eff = β_ν · Z.
        """
        n_sd = self.n_sd
        Z_nuc = float(self.Z)
        states = self.states
        S_mat = np.zeros((n_sd, n_sd))
        H_mat = np.zeros((n_sd, n_sd))

        # Integral caches (key: rounded numerical parameters)
        cache_ovlp: Dict[tuple, float] = {}
        cache_invr: Dict[tuple, float] = {}
        cache_rk: Dict[tuple, float] = {}

        def _ovlp(n1: int, l1: int, Z1: float,
                  n2: int, l2: int, Z2: float) -> float:
            if l1 != l2:
                return 0.0
            key = (n1, l1, round(Z1, 10), n2, l2, round(Z2, 10))
            if key not in cache_ovlp:
                if n1 == n2 and abs(Z1 - Z2) < 1e-10:
                    cache_ovlp[key] = 1.0
                else:
                    cache_ovlp[key] = _radial_overlap(n1, l1, Z1, n2, l2, Z2)
            return cache_ovlp[key]

        def _invr(n1: int, l1: int, Z1: float,
                  n2: int, l2: int, Z2: float) -> float:
            if l1 != l2:
                return 0.0
            key = (n1, l1, round(Z1, 10), n2, l2, round(Z2, 10))
            if key not in cache_invr:
                cache_invr[key] = _radial_1_over_r(n1, l1, Z1, n2, l2, Z2)
            return cache_invr[key]

        def _rk(n1: int, l1: int, Z1: float, n2: int, l2: int, Z2: float,
                n3: int, l3: int, Z3: float, n4: int, l4: int, Z4: float,
                kk: int) -> float:
            key = (n1, l1, round(Z1, 10), n2, l2, round(Z2, 10),
                   n3, l3, round(Z3, 10), n4, l4, round(Z4, 10), kk)
            if key not in cache_rk:
                cache_rk[key] = _slater_rk_fast(
                    n1, l1, Z1, n2, l2, Z2,
                    n3, l3, Z3, n4, l4, Z4, kk
                )
            return cache_rk[key]

        def _sp_overlap(a: int, Z_a: float, b: int, Z_b: float) -> float:
            """Spin-orbital overlap ⟨a(Z_a)|b(Z_b)⟩."""
            if (a & 1) != (b & 1):
                return 0.0
            sa, sb = a >> 1, b >> 1
            na, la, ma = states[sa]
            nb, lb, mb = states[sb]
            if ma != mb:
                return 0.0
            return _ovlp(na, la, Z_a, nb, lb, Z_b)

        def _sp_h1(a: int, Z_a: float, b: int, Z_b: float) -> float:
            """⟨a(Z_a)|T + V_nuc|b(Z_b)⟩ using ket eigenvalue identity."""
            if (a & 1) != (b & 1):
                return 0.0
            sa, sb = a >> 1, b >> 1
            na, la, ma = states[sa]
            nb, lb, mb = states[sb]
            if la != lb or ma != mb:
                return 0.0
            # From (T - Z_b/r)|ket⟩ = -Z_b²/(2n_b²)|ket⟩:
            # ⟨bra|T - Z_nuc/r|ket⟩ = -Z_b²/(2n_b²)·S + (Z_b - Z_nuc)·⟨1/r⟩
            s_val = _ovlp(na, la, Z_a, nb, lb, Z_b)
            invr_val = _invr(na, la, Z_a, nb, lb, Z_b)
            return -Z_b ** 2 / (2.0 * nb ** 2) * s_val + (Z_b - Z_nuc) * invr_val

        def _spatial_coulomb(p1: int, p2: int, q1: int, q2: int,
                             Z_bra: float, Z_ket: float) -> float:
            """⟨p1(Z_bra) p2(Z_bra)|1/r₁₂|q1(Z_ket) q2(Z_ket)⟩."""
            if states[p1][2] + states[p2][2] != states[q1][2] + states[q2][2]:
                return 0.0
            ck_13 = self._ac_k_map.get((p1, q1), [])
            ck_24 = self._ac_k_map.get((p2, q2), [])
            if not ck_13 or not ck_24:
                return 0.0
            np1, lp1 = states[p1][0], states[p1][1]
            np2, lp2 = states[p2][0], states[p2][1]
            nq1, lq1 = states[q1][0], states[q1][1]
            nq2, lq2 = states[q2][0], states[q2][1]
            val = 0.0
            for k1, c1 in ck_13:
                for k2, c2 in ck_24:
                    if k1 != k2:
                        continue
                    rk_val = _rk(np1, lp1, Z_bra, np2, lp2, Z_bra,
                                 nq1, lq1, Z_ket, nq2, lq2, Z_ket, k1)
                    val += c1 * c2 * rk_val
            return val

        # ---- Assemble SD-level matrices ----
        for I in range(n_sd):
            a1, a2 = self.sd_basis[I]
            Z_I = betas[I] * Z_nuc

            for J in range(I, n_sd):
                b1, b2 = self.sd_basis[J]
                Z_J = betas[J] * Z_nuc

                # 2×2 orbital overlap matrix M
                m11 = _sp_overlap(a1, Z_I, b1, Z_J)
                m12 = _sp_overlap(a1, Z_I, b2, Z_J)
                m21 = _sp_overlap(a2, Z_I, b1, Z_J)
                m22 = _sp_overlap(a2, Z_I, b2, Z_J)

                # Overlap = det(M)
                s_val = m11 * m22 - m12 * m21

                # One-body: Σ_ij h_ij · cofactor(M,i,j)
                h11 = _sp_h1(a1, Z_I, b1, Z_J)
                h12 = _sp_h1(a1, Z_I, b2, Z_J)
                h21 = _sp_h1(a2, Z_I, b1, Z_J)
                h22 = _sp_h1(a2, Z_I, b2, Z_J)
                h1_val = (h11 * m22 - h12 * m21
                          - h21 * m12 + h22 * m11)

                # Two-body: g_direct - g_exchange
                vee_val = 0.0
                sp_a1, sp_a2 = a1 >> 1, a2 >> 1
                sp_b1, sp_b2 = b1 >> 1, b2 >> 1

                # Direct: ⟨a₁a₂|V₁₂|b₁b₂⟩
                if (a1 & 1) == (b1 & 1) and (a2 & 1) == (b2 & 1):
                    vee_val += _spatial_coulomb(
                        sp_a1, sp_a2, sp_b1, sp_b2, Z_I, Z_J)

                # Exchange: -⟨a₁a₂|V₁₂|b₂b₁⟩
                if (a1 & 1) == (b2 & 1) and (a2 & 1) == (b1 & 1):
                    vee_val -= _spatial_coulomb(
                        sp_a1, sp_a2, sp_b2, sp_b1, Z_I, Z_J)

                h_val = h1_val + vee_val

                S_mat[I, J] = s_val
                S_mat[J, I] = s_val
                H_mat[I, J] = h_val
                H_mat[J, I] = h_val

        # Symmetrize for numerical stability
        S_mat = (S_mat + S_mat.T) / 2.0
        H_mat = (H_mat + H_mat.T) / 2.0

        if verbose:
            print(f"    Caches: {len(cache_ovlp)} overlaps, "
                  f"{len(cache_invr)} 1/r, {len(cache_rk)} R^k")

        return S_mat, H_mat

    def solve(self, E_trial: float, max_iter: int = 20,
              tol: float = 1e-6, damping: float = 0.5,
              verbose: bool = False) -> Dict:
        """
        Self-consistent generalized Sturmian CI.

        Parameters
        ----------
        E_trial : float
            Initial energy guess (must be < 0).
        max_iter : int
            Maximum self-consistency iterations.
        tol : float
            Convergence tolerance on energy (Ha).
        damping : float
            E_next = damping * E_computed + (1 - damping) * E_old.
        verbose : bool
            Print convergence info.
        """
        if E_trial >= 0:
            raise ValueError("E_trial must be negative")

        if verbose:
            print(f"\n  Generalized Sturmian CI: Z={self.Z}, max_n={self.max_n}, "
                  f"n_spatial={self.n_spatial}, n_sd={self.n_sd}")
            print(f"  E_trial = {E_trial:.6f}")

        E = E_trial
        convergence_log: List[Dict] = []
        converged = False

        for it in range(max_iter):
            betas = self._compute_betas(E)

            if verbose:
                unique_b = sorted(set(np.round(betas, 6)))
                print(f"\n  iter {it + 1}: E = {E:.6f}, "
                      f"{len(unique_b)} unique beta")

            S, H = self._build_sd_matrices(betas, verbose=verbose)

            # Check overlap conditioning
            eig_S = np.linalg.eigvalsh(S)
            n_neg = int(np.sum(eig_S < 0))
            min_eig = float(eig_S[0])

            if verbose:
                print(f"    S: min_eig = {min_eig:.2e}, "
                      f"cond = {eig_S[-1] / max(abs(min_eig), 1e-30):.1f}")

            if min_eig < 1e-12:
                # Löwdin with threshold for robustness
                eig_v, eig_U = np.linalg.eigh(S)
                keep = eig_v > 1e-8
                n_keep = int(np.sum(keep))
                if verbose:
                    print(f"    Löwdin truncation: {n_keep}/{len(eig_v)} kept")
                eig_v = eig_v[keep]
                eig_U = eig_U[:, keep]
                X = eig_U @ np.diag(1.0 / np.sqrt(eig_v))
                H_orth = X.T @ H @ X
                eigvals_h, eigvecs_h = eigh(H_orth)
                E_new = float(eigvals_h[0])
                c_orth = eigvecs_h[:, 0]
                # Map back to original basis
                c_full = X @ c_orth
            else:
                try:
                    eigvals_h, eigvecs_h = eigh(H, S)
                    E_new = float(eigvals_h[0])
                    c_full = eigvecs_h[:, 0]
                except np.linalg.LinAlgError:
                    if verbose:
                        print("    ERROR: eigenvalue solve failed")
                    return {'energy': np.inf, 'converged': False,
                            'convergence': convergence_log}

            dE = abs(E_new - E)
            convergence_log.append({
                'iteration': it + 1, 'E_trial': E,
                'E_computed': E_new, 'dE': dE
            })

            if verbose:
                print(f"    E_computed = {E_new:.6f}, dE = {dE:.2e}")

            if dE < tol:
                converged = True
                E = E_new
                if verbose:
                    print(f"  CONVERGED at iteration {it + 1}")
                break

            E = damping * E_new + (1.0 - damping) * E

        # Final solve at converged energy
        betas_final = self._compute_betas(E)
        S_final, H_final = self._build_sd_matrices(betas_final)

        eig_S_final = np.linalg.eigvalsh(S_final)
        if eig_S_final[0] < 1e-12:
            eig_v, eig_U = np.linalg.eigh(S_final)
            keep = eig_v > 1e-8
            eig_v, eig_U = eig_v[keep], eig_U[:, keep]
            X = eig_U @ np.diag(1.0 / np.sqrt(eig_v))
            H_orth = X.T @ H_final @ X
            eigvals_final, eigvecs_final = eigh(H_orth)
            c_gs = X @ eigvecs_final[:, 0]
        else:
            eigvals_final, eigvecs_final = eigh(H_final, S_final)
            c_gs = eigvecs_final[:, 0]

        E_gs = float(eigvals_final[0])

        # β analysis: group by (n₁, n₂) configuration type
        beta_by_config: Dict[Tuple[int, ...], float] = {}
        for I, sd in enumerate(self.sd_basis):
            spatial_ns = tuple(sorted([self.states[sp >> 1][0] for sp in sd]))
            if spatial_ns not in beta_by_config:
                beta_by_config[spatial_ns] = betas_final[I]

        # Detailed β info for top-weight SDs
        beta_info: List[Dict] = []
        Z_nuc = float(self.Z)
        for I, sd in enumerate(self.sd_basis):
            orbs = []
            for sp_idx in sd:
                n, l, m = self.states[sp_idx >> 1]
                spin = 'a' if sp_idx % 2 == 0 else 'b'
                orbs.append(f"{n}{'spdf'[l]}({spin})")
            beta_info.append({
                'sd_index': I, 'orbitals': orbs,
                'beta': betas_final[I], 'Z_eff': betas_final[I] * Z_nuc,
                'weight': c_gs[I] ** 2,
            })

        return {
            'energy': E_gs,
            'betas': betas_final,
            'beta_info': beta_info,
            'beta_by_config': beta_by_config,
            'coefficients': c_gs,
            'eigenvalues': eigvals_final,
            'convergence': convergence_log,
            'overlap_matrix': S_final,
            'converged': converged,
            'n_iterations': len(convergence_log),
            'eri_count': self._eri_angular_count,
        }


# ---------------------------------------------------------------------------
# Convenience functions
# ---------------------------------------------------------------------------

def compare_he_generalized(max_n: int = 2, verbose: bool = True) -> Dict:
    """
    Compare all three He methods: Standard FCI, Coulomb Sturmian, Generalized.

    Returns comparison dict with energies, errors, beta analysis, etc.
    """
    E_exact = -2.903724  # Ha (NIST)
    Z, n_e = 2, 2

    if verbose:
        print(f"\n{'=' * 65}")
        print(f"He 3-Way Comparison (max_n={max_n})")
        print(f"{'=' * 65}")
        print(f"Exact: {E_exact:.6f} Ha")

    # Standard FCI
    if verbose:
        print(f"\n--- Standard FCI ---")
    std = StandardFCI(Z, n_e, max_n)
    std_r = std.optimize_zeff(z_range=(1.0, 3.0), n_scan=15, verbose=verbose)
    E_std = std_r['energy']
    err_std = 100 * abs(E_std - E_exact) / abs(E_exact)

    # Coulomb Sturmian (BU-1)
    if verbose:
        print(f"\n--- Coulomb Sturmian ---")
    sturm = SturmianCI(Z, n_e, max_n)
    sturm_r = sturm.optimize_k(k_range=(0.5, 4.0), n_scan=20, verbose=verbose)
    E_sturm = sturm_r['energy']
    err_sturm = 100 * abs(E_sturm - E_exact) / abs(E_exact)

    # Generalized Sturmian (BU-1b)
    if verbose:
        print(f"\n--- Generalized Sturmian ---")
    gen = GeneralizedSturmianCI(Z, n_e, max_n)
    gen_r = gen.solve(E_sturm, max_iter=20, tol=1e-6, damping=0.5,
                      verbose=verbose)
    E_gen = gen_r['energy']
    err_gen = 100 * abs(E_gen - E_exact) / abs(E_exact)

    if verbose:
        print(f"\n{'=' * 65}")
        print(f"SUMMARY (max_n={max_n}):")
        print(f"  Standard FCI:       E = {E_std:.6f}, err = {err_std:.4f}%")
        print(f"  Coulomb Sturmian:   E = {E_sturm:.6f}, err = {err_sturm:.4f}%")
        print(f"  Generalized Sturm:  E = {E_gen:.6f}, err = {err_gen:.4f}%")
        print(f"  Converged: {gen_r['converged']} ({gen_r['n_iterations']} iter)")
        print(f"\n  beta by config type:")
        for cfg, beta in sorted(gen_r['beta_by_config'].items()):
            z_eff = beta * Z
            print(f"    n={cfg}: beta={beta:.4f}, Z_eff={z_eff:.4f}")
        print(f"\n  Angular ERI count: {gen_r['eri_count']} "
              f"(BU-1: {sturm_r['eri_count']})")
        print(f"{'=' * 65}")

    return {
        'max_n': max_n,
        'E_exact': E_exact,
        'E_standard': E_std, 'err_standard': err_std,
        'E_sturmian': E_sturm, 'err_sturmian': err_sturm,
        'E_generalized': E_gen, 'err_generalized': err_gen,
        'converged': gen_r['converged'],
        'n_iterations': gen_r['n_iterations'],
        'beta_by_config': gen_r['beta_by_config'],
        'eri_sturmian': sturm_r['eri_count'],
        'eri_generalized': gen_r['eri_count'],
        'convergence': gen_r['convergence'],
        'beta_info': gen_r['beta_info'],
    }


def compare_he(max_n: int = 2, verbose: bool = True) -> Dict:
    """
    Run Sturmian CI and standard FCI for He at given max_n.

    Returns comparison dict with energies, errors, sparsity, etc.
    """
    E_exact = -2.903724  # Ha (NIST)
    Z = 2
    n_e = 2

    if verbose:
        print(f"\n{'='*60}")
        print(f"He Comparison: Sturmian vs Standard FCI (max_n={max_n})")
        print(f"{'='*60}")
        print(f"Exact energy: {E_exact:.6f} Ha")

    # Standard FCI with Z_eff optimization
    if verbose:
        print(f"\n--- Standard FCI (Z_eff optimized) ---")
    std = StandardFCI(Z, n_e, max_n)
    std_result = std.optimize_zeff(z_range=(1.0, 3.0), n_scan=15, verbose=verbose)
    E_std = std_result['energy']
    err_std = 100 * abs(E_std - E_exact) / abs(E_exact)

    if verbose:
        print(f"  Energy: {E_std:.6f} Ha")
        print(f"  Error:  {err_std:.4f}%")

    # Sturmian CI with k optimization
    if verbose:
        print(f"\n--- Sturmian CI (k optimized) ---")
    sturm = SturmianCI(Z, n_e, max_n)
    sturm_result = sturm.optimize_k(k_range=(0.5, 4.0), n_scan=20, verbose=verbose)
    E_sturm = sturm_result['energy']
    err_sturm = 100 * abs(E_sturm - E_exact) / abs(E_exact)

    k_opt = sturm_result['k']
    if verbose:
        print(f"  Energy: {E_sturm:.6f} Ha")
        print(f"  Error:  {err_sturm:.4f}%")
        print(f"  k_opt:  {k_opt:.6f}")
        print(f"\n  Orbital exponents (Z_eff = n*k):")
        for (n, l, m), zeff in sturm_result['z_effs'].items():
            if m == 0:  # Print one per (n,l) shell
                print(f"    n={n}, l={l}: Z_eff = {zeff:.4f}")

    # Sparsity comparison
    if verbose:
        print(f"\n--- Sparsity ---")
        print(f"  Sturmian ERI nonzero: {sturm_result['eri_count']}")
        print(f"  Standard ERI nonzero: {std_result['eri_count']}")
        print(f"  Sturmian V_ee fraction: {sturm_result['vee_sparsity']:.4f}")

    # Overlap matrix diagnostics
    S = sturm_result['overlap_matrix']
    off_diag_max = np.max(np.abs(S - np.eye(len(S))))
    if verbose:
        print(f"\n  Overlap max off-diagonal: {off_diag_max:.6f}")
        print(f"  Overlap condition number: {np.linalg.cond(S):.2f}")

    # Summary
    if verbose:
        print(f"\n{'='*60}")
        print(f"SUMMARY (max_n={max_n}):")
        print(f"  Standard FCI: E = {E_std:.6f}, error = {err_std:.4f}%")
        print(f"  Sturmian CI:  E = {E_sturm:.6f}, error = {err_sturm:.4f}%")
        winner = "Sturmian" if E_sturm < E_std else "Standard"
        print(f"  Winner: {winner} (dE = {abs(E_sturm - E_std)*1000:.3f} mHa)")
        print(f"{'='*60}")

    return {
        'max_n': max_n,
        'E_exact': E_exact,
        'E_standard': E_std,
        'E_sturmian': E_sturm,
        'err_standard': err_std,
        'err_sturmian': err_sturm,
        'k_opt': k_opt,
        'Z_eff_opt': std_result['Z_eff'],
        'eri_nonzero_sturmian': sturm_result['eri_count'],
        'eri_nonzero_standard': std_result['eri_count'],
        'vee_sparsity_sturmian': sturm_result['vee_sparsity'],
        'overlap_cond': float(np.linalg.cond(S)),
    }
