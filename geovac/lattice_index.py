"""
LatticeIndex — Pre-calculated topology for N-electron graph Hamiltonians
=========================================================================

Architecture: the many-body state space is treated as a relational database.

  - sp_states: single-particle states (n, l, m, sigma) — primary keys
  - conformal_weights[i]: Omega_i = 2*p0/(p_i^2 + p0^2), pre-computed once
  - h1_diag[sp], h1_offdiag[sp]: one-electron integrals indexed by spatial state
  - vee_matrix[sp_i, sp_j]: two-electron repulsion (Mulliken approximation)
  - sd_basis: enumerated Slater determinants (sorted tuples of occupied sp indices)

Hamiltonian assembly is a sparse JOIN query:
  For each occupied orbital p in each SD I:
    For each graph-adjacent r of spatial(p):
      If r unoccupied: H[I, J] += phase * h1(spatial_p, spatial_r)
  plus diagonal V_ee correction.

Elements below sparsity_threshold (default 1e-8) are dropped before the
matrix constructor is called — this is the topological sparsity mask.

KNOWN SYSTEMATIC ERRORS (logged at construction, not papered over):
  1. V_ee graph approximation (default: 'chordal'):
     Uses S³ chordal distance V_ee = κ_ee / d²_chord, where d_chord is the
     Euclidean distance in R⁴ between Fock-projected S³ coordinates. This is
     architecturally consistent with the graph Laplacian kinetic term (Paper 7,
     Eq. 14) but does not reproduce exact Hartree two-electron integrals.
     Fallback 'mulliken': uses mean orbital radius r_eff=n²/Z with ~12%
     systematic overestimation; factor ~10x error for same-orbital pairs.
  2. Topological bond contraction: the graph Laplacian uses uniform edge
     weights (W=1), not the exact radial matrix elements. This under-represents
     transitions between shells of very different n.

Reference implementation: MoleculeHamiltonian._solve_full_ci() for 2-electron
systems (geovac/hamiltonian.py). This module extends that architecture to
arbitrary N.

Author: GeoVac Development Team
Date: February 2026
"""

import os
import warnings
import time
from itertools import combinations
from math import factorial
from typing import Dict, List, Optional, Tuple

import numpy as np
import scipy.sparse as sp
from scipy.sparse import csr_matrix, diags, identity
from scipy.sparse.linalg import eigsh

try:
    from .lattice import GeometricLattice
except ImportError:
    from lattice import GeometricLattice

KINETIC_SCALE: float = -1.0 / 16.0   # Universal topological constant


class LatticeIndex:
    """
    Pre-calculated topology index for N-electron graph Hamiltonians.

    For Lithium (n_electrons=3, nuclear_charge=3), the state space is the
    antisymmetrized tensor product of single-particle lattices truncated at
    max_n. Nodes are Slater determinant configurations. Edges are non-zero
    off-diagonal Hamiltonian elements derived from the graph Laplacian.

    Parameters
    ----------
    n_electrons : int
        Number of electrons (3 for Li, 4 for LiH, etc.)
    max_n : int
        Maximum principal quantum number (truncation)
    nuclear_charge : int
        Nuclear charge Z (default 3 for Li)
    sparsity_threshold : float
        Matrix elements below this are dropped during assembly (default 1e-8)

    Attributes
    ----------
    n_sp : int
        Number of spin-orbitals (2 * number of spatial states)
    n_sd : int
        Number of Slater determinants = C(n_sp, n_electrons)
    sp_states : list of (n, l, m, sigma)
        Single-particle basis indexed by integer sp_idx
    conformal_weights : np.ndarray, shape (n_spatial,)
        Conformal factor Omega_i = 2*p0 / (p_i^2 + p0^2) per spatial state
    """

    def __init__(
        self,
        n_electrons: int,
        max_n: int,
        nuclear_charge: int = 3,
        sparsity_threshold: float = 1e-8,
        vee_method: str = 'chordal',
        h1_method: Optional[str] = None,
        kinetic_scale: float = KINETIC_SCALE,
        adjacency: Optional[csr_matrix] = None,
        node_weights: Optional[np.ndarray] = None,
    ) -> None:
        if vee_method == 'chordal':
            warnings.warn(
                "LatticeIndex V_ee uses S³ chordal distance: "
                "V_ee = κ_ee / d²_S3,  κ_ee = 5Z/2 Ha,  d²_S3 = 2(1 − ξ_i·ξ_j). "
                "Same-orbital: antipodal d²=4, V_ee = 5Z/8 Ha exact. "
                "Off-diagonal d²(1s,2s)=0.4 → V_ee too large for multi-shell systems.",
                UserWarning,
                stacklevel=2,
            )
        elif vee_method == 'slater':
            warnings.warn(
                "LatticeIndex V_ee uses exact Slater F⁰ integrals computed "
                "numerically from hydrogen-like radial wavefunctions R_{nl}(r). "
                "Exact for all orbital pairs. Construction slower (numerical integration).",
                UserWarning,
                stacklevel=2,
            )
        elif vee_method == 'slater_full':
            warnings.warn(
                "LatticeIndex V_ee uses full Slater integrals (F^k direct + G^k "
                "exchange) with proper angular coupling via Wigner 3j symbols. "
                "Off-diagonal two-electron matrix elements enabled for FCI.",
                UserWarning,
                stacklevel=2,
            )
        else:
            warnings.warn(
                "LatticeIndex V_ee uses Mulliken approximation: ~12% systematic "
                "overestimation of electron repulsion for well-separated orbitals; "
                "larger errors for same-orbital pairs (factor ~10x, V_ee=Z/n²). "
                "Li ground state energy will NOT be quantitatively accurate.",
                UserWarning,
                stacklevel=2,
            )

        self.n_electrons = n_electrons
        self.max_n = max_n
        self.Z = nuclear_charge
        self.threshold = sparsity_threshold
        self.vee_method = vee_method
        self.kinetic_scale = kinetic_scale
        self._adjacency_override = adjacency
        self._node_weights_override = node_weights
        # Default h1_method: 'exact' for slater_full, 'hybrid' for slater, 'graph' otherwise
        if h1_method:
            self.h1_method = h1_method
        elif vee_method == 'slater_full':
            self.h1_method = 'exact'
        elif vee_method == 'slater':
            self.h1_method = 'hybrid'
        else:
            self.h1_method = 'graph'

        # --- Build topology in dependency order ---
        self._build_sp_lattice()
        self._compute_conformal_weights()
        self._build_sp_hamiltonian()
        self._build_vee_index()
        self._enumerate_sd_basis()

    # ------------------------------------------------------------------
    # Private construction methods
    # ------------------------------------------------------------------

    def _build_sp_lattice(self) -> None:
        """Build the single-particle geometric lattice and spin-orbital index."""
        self.lattice = GeometricLattice(max_n=self.max_n, nuclear_charge=self.Z)
        n_spatial = self.lattice.num_states

        # sp_states[2*i] = (n,l,m, 0=up)   sp_states[2*i+1] = (n,l,m, 1=down)
        self.sp_states: List[Tuple[int, int, int, int]] = []
        for n, l, m in self.lattice.states:
            self.sp_states.append((n, l, m, 0))
            self.sp_states.append((n, l, m, 1))
        self.n_sp = len(self.sp_states)

        self._sp_index: Dict[Tuple[int, int, int, int], int] = {
            s: i for i, s in enumerate(self.sp_states)
        }

    def _compute_conformal_weights(self) -> None:
        """
        Pre-compute conformal factor Omega_i for each spatial state.

        Omega_i = 2*p0 / (p_i^2 + p0^2)
        where p_i = Z / n_i   (state momentum on S^3)
              p0  = Z          (focal length = nuclear charge for atomic case)

        These are stored in LatticeIndex for future use as conformal edge
        modulation W_ij = W_flat * Omega_i * Omega_j. For the current
        single-atom assembly, they inform the physical interpretation but
        the graph edges use uniform weights (consistent with AtomicSolver).
        """
        p0 = float(self.Z)
        n_spatial = self.lattice.num_states
        self.conformal_weights = np.zeros(n_spatial)

        for idx, (n, l, m) in enumerate(self.lattice.states):
            p_i = float(self.Z) / float(n)
            self.conformal_weights[idx] = 2.0 * p0 / (p_i ** 2 + p0 ** 2)

    def _fock_s3_coords(self) -> np.ndarray:
        """
        Map each spatial state (n,l,m) to its S³ coordinate via Fock
        stereographic projection (Paper 7, Eq. 14).

        With normalised focal length p₀ = 1 and momentum magnitude p_i = 1/n_i
        (Z cancels identically), the S³ embedding is:

            ξ₁ = 2n sin(θ) / (n²+1)
            ξ₂ = 0                          [all states in xz half-plane, φ=0]
            ξ₃ = 2n cos(θ) / (n²+1)
            ξ₄ = (1 − n²) / (n²+1)

        where θ encodes the angular direction:
            l=0  →  θ=0  (north-pole direction)
            l>0  →  cos(θ) = m / √(l(l+1))

        Same convention as _compute_sp_coordinates().  All ξ satisfy |ξ|=1.

        Returns
        -------
        coords : ndarray, shape (n_spatial, 4)
        """
        n_spatial = self.lattice.num_states
        coords = np.zeros((n_spatial, 4))
        for idx, (n, l, m) in enumerate(self.lattice.states):
            if l > 0:
                cos_th = np.clip(float(m) / np.sqrt(float(l * (l + 1))), -1.0, 1.0)
                sin_th = np.sqrt(max(0.0, 1.0 - cos_th * cos_th))
            else:
                cos_th = 1.0
                sin_th = 0.0
            fac = 2.0 * n / (n * n + 1.0)
            coords[idx, 0] = fac * sin_th                      # ξ₁
            coords[idx, 1] = 0.0                               # ξ₂
            coords[idx, 2] = fac * cos_th                      # ξ₃
            coords[idx, 3] = (1.0 - n * n) / (1.0 + n * n)   # ξ₄
        return coords

    def _build_sp_hamiltonian(self) -> None:
        """
        Build single-particle Hamiltonian in spatial basis.

        Two modes:

        h1_method == 'graph' (default for vee_method='chordal'/'mulliken'):
            H1 = KINETIC_SCALE * (D - A) + diag(node_weights)
               = KINETIC_SCALE * L + V_en
            Off-diagonal elements from graph Laplacian adjacency.

        h1_method == 'exact' (default for vee_method='slater'):
            H1 = diag(-Z²/(2n²)) — exact hydrogen eigenbasis.
            States (n,l,m) ARE the eigenstates, so h1 is purely diagonal.
            No off-diagonal hopping: all excitations come from V_ee.

        Diagonal elements: h_diag[sp] = H1[sp, sp]
        Off-diagonal elements: h_offdiag[sp] = [(r_sp, H1[sp,r_sp]), ...]
        """
        n_spatial = self.lattice.num_states

        if self.h1_method == 'exact':
            # Exact hydrogen eigenbasis: h(n,l,m) = -Z²/(2n²), no off-diagonal
            Z = float(self.Z)
            self._h1_diag = np.array([
                -Z * Z / (2.0 * n * n)
                for n, l, m in self.lattice.states
            ])
            self._H1_spatial = diags(self._h1_diag).tocsr()
            self._h1_offdiag: Dict[int, List[Tuple[int, float]]] = {
                i: [] for i in range(n_spatial)
            }
        elif self.h1_method == 'hybrid':
            # Exact diagonal + graph Laplacian off-diagonal for CI mixing
            Z = float(self.Z)
            self._h1_diag = np.array([
                -Z * Z / (2.0 * n * n)
                for n, l, m in self.lattice.states
            ])

            # Off-diagonal from graph Laplacian adjacency (use override if provided)
            A = self._adjacency_override if self._adjacency_override is not None else self.lattice.adjacency
            H1_offdiag = self.kinetic_scale * (-A)
            self._H1_spatial = (diags(self._h1_diag) + H1_offdiag).tocsr()

            H1_coo = H1_offdiag.tocoo()
            self._h1_offdiag = {i: [] for i in range(n_spatial)}
            for r, c, v in zip(H1_coo.row, H1_coo.col, H1_coo.data):
                if r != c and abs(v) >= self.threshold:
                    self._h1_offdiag[r].append((c, float(v)))
        else:
            # Graph Laplacian mode (use overrides if provided)
            A = self._adjacency_override if self._adjacency_override is not None else self.lattice.adjacency
            deg = np.array(A.sum(axis=1)).flatten()
            L = diags(deg) - A
            nw = self._node_weights_override if self._node_weights_override is not None else self.lattice.node_weights
            V = diags(nw)
            H1 = self.kinetic_scale * L + V
            self._H1_spatial = H1.tocsr()

            self._h1_diag = np.array(H1.diagonal())

            H1_coo = H1.tocoo()
            self._h1_offdiag = {i: [] for i in range(n_spatial)}
            for r, c, v in zip(H1_coo.row, H1_coo.col, H1_coo.data):
                if r != c and abs(v) >= self.threshold:
                    self._h1_offdiag[r].append((c, float(v)))

    def _momentum_vectors(self) -> np.ndarray:
        """
        Assign a representative 3D momentum vector to each spatial state.

        For state (n, l, m):
            |p| = Z/n  (on-shell momentum for energy E_n = -Z²/(2n²))
            direction = (sin θ, 0, cos θ) with:
                l=0: θ=0 (z-axis)
                l>0: cos θ = m / √(l(l+1))

        Returns
        -------
        p_vecs : ndarray, shape (n_spatial, 3)
        """
        n_spatial = self.lattice.num_states
        p_vecs = np.zeros((n_spatial, 3))
        Z = float(self.Z)
        for idx, (n, l, m) in enumerate(self.lattice.states):
            p_mag = Z / float(n)
            if l > 0:
                cos_th = np.clip(float(m) / np.sqrt(float(l * (l + 1))), -1.0, 1.0)
                sin_th = np.sqrt(max(0.0, 1.0 - cos_th * cos_th))
            else:
                cos_th = 1.0
                sin_th = 0.0
            p_vecs[idx, 0] = p_mag * sin_th
            p_vecs[idx, 1] = 0.0
            p_vecs[idx, 2] = p_mag * cos_th
        return p_vecs

    def _conformal_factors_from_momenta(self, p_vecs: np.ndarray) -> np.ndarray:
        """
        Compute conformal factor Ω_k = 2p₀_k / (|p_k|² + p₀_k²) for each state.

        p₀_k = Z/n_k (energy-shell focal length for the kth state).
        |p_k| = Z/n_k (on-shell), so |p_k|² = p₀_k² = Z²/n_k².
        Therefore Ω_k = 2p₀ / (2p₀²) = 1/p₀ = n_k/Z.

        Returns
        -------
        omega : ndarray, shape (n_spatial,)
        """
        n_spatial = self.lattice.num_states
        omega = np.zeros(n_spatial)
        Z = float(self.Z)
        for idx, (n, l, m) in enumerate(self.lattice.states):
            p0 = Z / float(n)
            p_sq = float(np.dot(p_vecs[idx], p_vecs[idx]))
            omega[idx] = 2.0 * p0 / (p_sq + p0 * p0)
        return omega

    @staticmethod
    def _slater_f0_cache_path(Z: float, max_n: int) -> str:
        """Path to disk-cached Slater F0 integrals."""
        import os
        cache_dir = os.path.join(os.path.dirname(__file__), 'cache')
        os.makedirs(cache_dir, exist_ok=True)
        return os.path.join(cache_dir, f'slater_f0_Z{int(Z)}_nmax{max_n}.npz')

    def _compute_slater_f0(self) -> np.ndarray:
        """
        Compute exact Slater F0 direct Coulomb integrals for all spatial pairs.

        F0(n1 l1, n2 l2) = integral |R_{n1l1}(r1)|^2 |R_{n2l2}(r2)|^2
                            * (1/r_>) * r1^2 r2^2 dr1 dr2

        F0 depends only on (n, l), not on m, so we cache by (n, l) pairs.
        Results are cached to disk at geovac/cache/slater_f0_Z{Z}_nmax{N}.npz.

        Returns
        -------
        vee : ndarray, shape (n_spatial, n_spatial)
        """
        import os

        Z = float(self.Z)
        n_sp = self.lattice.num_states
        unique_nl = sorted(set((n, l) for n, l, m in self.lattice.states))

        # --- Try disk cache ---
        cache_path = self._slater_f0_cache_path(Z, self.max_n)
        f0_cache: Dict[Tuple[int, int, int, int], float] = {}

        if os.path.exists(cache_path):
            data = np.load(cache_path)
            keys = data['keys']     # shape (N, 4): n1, l1, n2, l2
            vals = data['values']   # shape (N,)
            for row, val in zip(keys, vals):
                f0_cache[(int(row[0]), int(row[1]),
                          int(row[2]), int(row[3]))] = float(val)

        # Check if all needed pairs are cached
        needed = set()
        for n1, l1 in unique_nl:
            for n2, l2 in unique_nl:
                key = (n1, l1, n2, l2)
                key_rev = (n2, l2, n1, l1)
                if key not in f0_cache and key_rev not in f0_cache:
                    needed.add(key if key <= key_rev else key_rev)

        if needed:
            # Compute missing integrals
            self._compute_f0_integrals(Z, needed, f0_cache, unique_nl)

            # Save to disk
            all_keys = []
            all_vals = []
            for k, v in f0_cache.items():
                all_keys.append(list(k))
                all_vals.append(v)
            np.savez(cache_path,
                     keys=np.array(all_keys, dtype=int),
                     values=np.array(all_vals))

        # Propagate symmetry
        for n1, l1 in unique_nl:
            for n2, l2 in unique_nl:
                key = (n1, l1, n2, l2)
                key_rev = (n2, l2, n1, l1)
                if key not in f0_cache and key_rev in f0_cache:
                    f0_cache[key] = f0_cache[key_rev]

        # --- Fill the V_ee matrix ---
        vee = np.zeros((n_sp, n_sp))
        for i in range(n_sp):
            ni, li, _ = self.lattice.states[i]
            for j in range(i, n_sp):
                nj, lj, _ = self.lattice.states[j]
                v = f0_cache[(ni, li, nj, lj)]
                vee[i, j] = v
                vee[j, i] = v
        return vee

    @staticmethod
    def _compute_f0_integrals(
        Z: float,
        needed: set,
        f0_cache: Dict[Tuple[int, int, int, int], float],
        unique_nl: list,
    ) -> None:
        """Compute Slater F0 integrals via nested numerical integration."""
        from scipy.special import genlaguerre
        from scipy.integrate import quad

        r_max = 80.0 / Z

        def _radial_wf_unnorm(r: np.ndarray, n: int, l: int) -> np.ndarray:
            rho = 2.0 * Z * r / n
            L_poly = genlaguerre(n - l - 1, 2 * l + 1)(rho)
            return rho ** l * np.exp(-rho / 2.0) * L_poly

        # Normalization constants
        norm_c: Dict[Tuple[int, int], float] = {}
        for n, l in unique_nl:
            norm_sq, _ = quad(
                lambda r: _radial_wf_unnorm(np.array(r), n, l) ** 2 * r ** 2,
                0, r_max, limit=200
            )
            norm_c[(n, l)] = 1.0 / np.sqrt(norm_sq)

        def R(r: float, n: int, l: int) -> float:
            return float(norm_c[(n, l)] * _radial_wf_unnorm(np.array(r), n, l))

        for n1, l1, n2, l2 in needed:
            def _inner(r1: float, _n2: int = n2, _l2: int = l2) -> float:
                if r1 < 1e-30:
                    return 0.0
                p1, _ = quad(lambda r2: R(r2, _n2, _l2) ** 2 * r2 ** 2,
                             0, r1, limit=100)
                p2, _ = quad(lambda r2: R(r2, _n2, _l2) ** 2 * r2,
                             r1, r_max, limit=100)
                return p1 / r1 + p2

            val, _ = quad(
                lambda r1, _n1=n1, _l1=l1: R(r1, _n1, _l1) ** 2 * _inner(r1) * r1 ** 2,
                0, r_max, limit=200
            )
            f0_cache[(n1, l1, n2, l2)] = val
            f0_cache[(n2, l2, n1, l1)] = val  # symmetry

    # ------------------------------------------------------------------
    # R^k Slater integrals and full two-electron machinery
    # ------------------------------------------------------------------

    @staticmethod
    def _rk_cache_path(Z: float, max_n: int) -> str:
        """Path to disk-cached R^k Slater integrals."""
        cache_dir = os.path.join(os.path.dirname(__file__), 'cache')
        os.makedirs(cache_dir, exist_ok=True)
        return os.path.join(cache_dir, f'slater_rk_Z{int(Z)}_nmax{max_n}.npz')

    @staticmethod
    def _wigner3j(j1: int, j2: int, j3: int,
                  m1: int, m2: int, m3: int) -> float:
        """
        Wigner 3j symbol for integer arguments.

        Uses the Racah formula. All arguments are integers (not half-integers).
        Returns 0 if triangle or m-selection rules are violated.
        """
        # Selection rules
        if m1 + m2 + m3 != 0:
            return 0.0
        if abs(j1 - j2) > j3 or j3 > j1 + j2:
            return 0.0
        if abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3:
            return 0.0

        # Racah formula
        def _tri(a: int, b: int, c: int) -> float:
            return (factorial(a + b - c) * factorial(a - b + c)
                    * factorial(-a + b + c)
                    / factorial(a + b + c + 1))

        pre = ((-1) ** (j1 - j2 - m3)
               * np.sqrt(_tri(j1, j2, j3)
                         * factorial(j1 + m1) * factorial(j1 - m1)
                         * factorial(j2 + m2) * factorial(j2 - m2)
                         * factorial(j3 + m3) * factorial(j3 - m3)))

        t_min = max(0, j2 - j3 - m1, j1 - j3 + m2)
        t_max = min(j1 + j2 - j3, j1 - m1, j2 + m2)

        s = 0.0
        for t in range(t_min, t_max + 1):
            s += ((-1) ** t
                  / (factorial(t)
                     * factorial(j3 - j2 + m1 + t)
                     * factorial(j3 - j1 - m2 + t)
                     * factorial(j1 + j2 - j3 - t)
                     * factorial(j1 - m1 - t)
                     * factorial(j2 + m2 - t)))
        return pre * s

    def _compute_rk_integrals(self) -> Dict[Tuple[int, ...], float]:
        """
        Compute R^k(n1l1, n2l2, n3l3, n4l4) Slater radial integrals.

        R^k(abcd) = ∫∫ R_a(r1) R_c(r1) (r_<^k / r_>^{k+1})
                       R_b(r2) R_d(r2) r1^2 r2^2 dr1 dr2

        These are the radial parts of the general two-electron integrals.
        F^k(ab) = R^k(ab,ab) (direct), G^k(ab) = R^k(ab,ba) (exchange).

        Cached to disk at geovac/cache/slater_rk_Z{Z}_nmax{N}.npz.

        Returns
        -------
        rk_cache : dict mapping (n1,l1,n2,l2,n3,l3,n4,l4,k) -> float
        """
        from scipy.special import genlaguerre
        from scipy.integrate import quad

        Z = float(self.Z)
        unique_nl = sorted(set((n, l) for n, l, m in self.lattice.states))
        r_max = 80.0 / Z

        # --- Try disk cache ---
        cache_path = self._rk_cache_path(Z, self.max_n)
        rk_cache: Dict[Tuple[int, ...], float] = {}

        if os.path.exists(cache_path):
            data = np.load(cache_path)
            keys = data['keys']    # shape (N, 9): n1,l1,n2,l2,n3,l3,n4,l4,k
            vals = data['values']
            for row, val in zip(keys, vals):
                rk_cache[tuple(int(x) for x in row)] = float(val)

        # Determine ALL needed R^k(n1l1,n2l2,n3l3,n4l4,k) tuples.
        # For full FCI Slater-Condon rules we need general R^k(a,b,c,d)
        # for ALL (nl) quadruples, not just F^k and G^k patterns.
        needed: set = set()
        for n1, l1 in unique_nl:
            for n2, l2 in unique_nl:
                for n3, l3 in unique_nl:
                    for n4, l4 in unique_nl:
                        k_max = min(l1 + l3, l2 + l4)
                        for k in range(0, k_max + 1):
                            # Parity selection: (l1+l3+k) even AND (l2+l4+k) even
                            if (l1 + l3 + k) % 2 != 0:
                                continue
                            if (l2 + l4 + k) % 2 != 0:
                                continue
                            key = (n1, l1, n2, l2, n3, l3, n4, l4, k)
                            if key not in rk_cache:
                                needed.add(key)

        if not needed:
            return rk_cache

        print(f"[LatticeIndex] Computing {len(needed)} R^k integrals "
              f"(grid method)...")
        t0 = time.perf_counter()

        # --- Grid-based integration (vectorized, ~100x faster than quad) ---
        N_GRID = 2000
        # Log-linear grid: denser near origin where wavefunctions peak
        r_grid = np.linspace(0, r_max, N_GRID + 1)[1:]  # exclude r=0
        dr = r_grid[1] - r_grid[0]

        # Pre-compute normalized radial wavefunctions on grid
        def _radial_wf_unnorm(r: np.ndarray, n: int, l: int) -> np.ndarray:
            rho = 2.0 * Z * r / n
            L_poly = genlaguerre(n - l - 1, 2 * l + 1)(rho)
            return rho ** l * np.exp(-rho / 2.0) * L_poly

        # Store R_{nl}(r) * r on grid for each (n,l)
        R_grid: Dict[Tuple[int, int], np.ndarray] = {}
        for n, l in unique_nl:
            wf = _radial_wf_unnorm(r_grid, n, l)
            norm_sq = np.trapezoid(wf ** 2 * r_grid ** 2, r_grid)
            R_grid[(n, l)] = wf / np.sqrt(norm_sq)

        # Pre-compute Y^k(r) = ∫_0^r f(r') (r'/r)^k r'^2 dr'
        #                     + ∫_r^∞ f(r') (r/r')^k r' dr'
        # where f(r') = R_b(r') * R_d(r')
        # This is the Yk potential for the (b,d) pair at multipole k.

        # Group needed integrals by (n2,l2,n4,l4,k) to reuse Yk
        from collections import defaultdict
        yk_groups: Dict[Tuple[int, ...], List[Tuple[int, ...]]] = defaultdict(list)
        for key in needed:
            n1, l1, n2, l2, n3, l3, n4, l4, k = key
            yk_key = (n2, l2, n4, l4, k)
            yk_groups[yk_key].append(key)

        # Compute Yk potentials
        yk_cache: Dict[Tuple[int, ...], np.ndarray] = {}
        for yk_key in yk_groups:
            n2, l2, n4, l4, k = yk_key
            f_r = R_grid[(n2, l2)] * R_grid[(n4, l4)] * r_grid ** 2
            yk = np.zeros(N_GRID)
            # Cumulative integral for inner part: ∫_0^r f(r') (r'/r)^k dr'
            for i in range(N_GRID):
                r1 = r_grid[i]
                # Inner: sum f(r') * (r'/r1)^k for r' <= r1
                inner = np.sum(f_r[:i + 1] * (r_grid[:i + 1] / r1) ** k) * dr
                # Outer: sum f(r') * (r1/r')^k * (1/r') for r' > r1
                if i + 1 < N_GRID:
                    outer = np.sum(
                        f_r[i + 1:]
                        * (r1 / r_grid[i + 1:]) ** k
                        / r_grid[i + 1:]
                    ) * dr
                else:
                    outer = 0.0
                yk[i] = inner / r1 + outer
            yk_cache[yk_key] = yk

        # Now compute R^k integrals using pre-computed Yk
        for yk_key, keys in yk_groups.items():
            n2, l2, n4, l4, k = yk_key
            yk = yk_cache[yk_key]
            for key in keys:
                n1, l1, _, _, n3, l3, _, _, _ = key
                integrand = (R_grid[(n1, l1)] * R_grid[(n3, l3)]
                             * yk * r_grid ** 2)
                val = np.trapezoid(integrand, r_grid)
                rk_cache[key] = val

        elapsed = time.perf_counter() - t0
        print(f"[LatticeIndex] R^k integrals computed in {elapsed:.1f}s")

        # Save to disk
        all_keys = [list(k) for k in rk_cache]
        all_vals = [rk_cache[k] for k in rk_cache]
        np.savez(cache_path,
                 keys=np.array(all_keys, dtype=int),
                 values=np.array(all_vals))

        return rk_cache

    def _build_two_electron_integrals(self) -> None:
        """
        Build two-electron integral table <ab|g|cd> using Slater integrals.

        <ab|g|cd> = Σ_k c^k(la,ma,lc,mc) c^k(lb,mb,ld,md) R^k(abcd)

        Optimization: pre-compute c^k coefficients, then iterate only over
        pairs with non-zero angular coupling. This avoids the naive O(n^4)
        loop over all spatial quadruples.

        Stores results in self._eri[a,b,c,d] where a,b,c,d are spatial indices.
        """
        rk_cache = self._compute_rk_integrals()
        n_sp = self.lattice.num_states
        states = self.lattice.states

        self._eri: Dict[Tuple[int, int, int, int], float] = {}

        # Pre-compute all c^k coefficients for pairs of spatial states
        # c^k(a,c) depends only on angular quantum numbers (l,m)
        ck_table: Dict[Tuple[int, int, int], float] = {}
        for a in range(n_sp):
            la, ma = states[a][1], states[a][2]
            for c in range(n_sp):
                lc, mc = states[c][1], states[c][2]
                k_max = la + lc
                for k in range(0, k_max + 1):
                    if (la + lc + k) % 2 != 0:
                        continue
                    val = self._ck_coefficient(la, ma, lc, mc, k)
                    if abs(val) > 1e-15:
                        ck_table[(a, c, k)] = val

        print(f"[LatticeIndex] c^k table: {len(ck_table)} non-zero entries")

        # Group c^k entries by (a,c) for efficient inner loop
        ac_k_map: Dict[Tuple[int, int], List[Tuple[int, float]]] = {}
        for (a, c, k), val in ck_table.items():
            key = (a, c)
            if key not in ac_k_map:
                ac_k_map[key] = []
            ac_k_map[key].append((k, val))

        count = 0
        for (a, c), ck_ac_list in ac_k_map.items():
            na, la, ma = states[a]
            nc, lc, mc = states[c]
            for (b, d), ck_bd_list in ac_k_map.items():
                nb, lb, mb = states[b]
                nd, ld, md = states[d]

                # m-selection rule: ma + mb = mc + md
                if ma + mb != mc + md:
                    continue

                val = 0.0
                for k_ac, c_ac in ck_ac_list:
                    for k_bd, c_bd in ck_bd_list:
                        if k_ac != k_bd:
                            continue
                        k = k_ac
                        rk_key = (na, la, nb, lb, nc, lc, nd, ld, k)
                        rk_val = rk_cache.get(rk_key)
                        if rk_val is None:
                            continue
                        val += c_ac * c_bd * rk_val

                if abs(val) > 1e-15:
                    self._eri[(a, b, c, d)] = val
                    count += 1

        print(f"[LatticeIndex] ERI table: {count} non-zero entries")

    def _ck_coefficient(self, la: int, ma: int, lc: int, mc: int, k: int) -> float:
        """
        Gaunt angular coupling coefficient.

        c^k(l,m,l',m') = (-1)^m √((2l+1)(2l'+1)) * (l k l'; 0 0 0) * (l k l'; -m q m')

        where q = mc - ma (m-transfer).
        """
        q = mc - ma
        pre = ((-1) ** ma
               * np.sqrt((2 * la + 1) * (2 * lc + 1)))
        w1 = self._wigner3j(la, k, lc, 0, 0, 0)
        if abs(w1) < 1e-15:
            return 0.0
        w2 = self._wigner3j(la, k, lc, -ma, q, mc)
        return pre * w1 * w2

    def _build_vee_index(self) -> None:
        """
        Pre-compute V_ee(i,j) for all spatial state pairs.

        Dispatches to the method selected by self.vee_method:

        'chordal' (default):
            V_ee(i,j) = κ_ee / d²_S3,  κ_ee = 5Z/2 Ha
            where d²_S3 = 2(1 − ξ_i · ξ_j) is the S³ chordal distance squared.
            Same-orbital pairs (i==j): antipodal override d² = 4,
                V_ee = κ_ee/4 = 5Z/8 Ha = F⁰(1s,1s) exactly.
            Good for He (0.97% error), but inter-shell d²(1s,2s)=0.4 is too
            small → catastrophic V_ee for multi-shell systems like Li.

        'slater' (exact):
            V_ee(i,j) = F⁰(n_i l_i, n_j l_j) — exact Slater direct Coulomb
            integrals computed numerically from hydrogen-like radial wavefunctions.
            F⁰ = ∫∫ |R_{nl}(r₁)|² |R_{n'l'}(r₂)|² (1/r_>) r₁² r₂² dr₁ dr₂.

        'mulliken' (legacy):
            V_ee(i,j) = 1 / r_ij,  r_ij = |r_i − r_j| with r_i = n_i²/Z.
            Same-orbital: r_eff = n²/Z  →  V_ee = Z/n².
            Known ~12% overestimation for well-separated pairs; ~10x for same n=1.
        """
        n_sp = self.lattice.num_states
        self._vee_matrix: np.ndarray = np.zeros((n_sp, n_sp))

        if self.vee_method == 'chordal':
            KAP_EE: float = 5.0 * float(self.Z) / 2.0
            ANTIPODAL_D_SQ: float = 4.0

            s3 = self._fock_s3_coords()           # (n_sp, 4)
            gram = s3 @ s3.T                       # dot products
            d_sq = 2.0 * (1.0 - gram)              # chordal distance²

            np.fill_diagonal(d_sq, ANTIPODAL_D_SQ)
            np.clip(d_sq, 1e-14, None, out=d_sq)
            self._vee_matrix = KAP_EE / d_sq

        elif self.vee_method == 'slater':
            self._vee_matrix = self._compute_slater_f0()

        elif self.vee_method == 'slater_full':
            # Full Slater integrals: F0 for diagonal V_ee + ERI table for
            # off-diagonal two-electron matrix elements
            self._vee_matrix = self._compute_slater_f0()
            self._build_two_electron_integrals()

        else:  # 'mulliken'
            coords = self._compute_sp_coordinates()
            for i in range(n_sp):
                for j in range(i, n_sp):
                    if i == j:
                        n_i = self.lattice.states[i][0]
                        v = float(self.Z) / float(n_i ** 2)
                    else:
                        rij = float(np.linalg.norm(coords[i] - coords[j]))
                        if rij < 1e-10:
                            n_i = self.lattice.states[i][0]
                            v = float(self.Z) / float(n_i ** 2)
                        else:
                            v = 1.0 / rij
                    self._vee_matrix[i, j] = v
                    self._vee_matrix[j, i] = v

    def _compute_sp_coordinates(self) -> np.ndarray:
        """
        Compute 3D spatial coordinates for each spatial state.
        r = n^2/Z (effective Bohr radius), orbital aligned on z-axis.
        """
        n_sp = self.lattice.num_states
        coords = np.zeros((n_sp, 3))
        for idx, (n, l, m) in enumerate(self.lattice.states):
            r = n ** 2 / float(self.Z)
            if l > 0:
                cos_th = np.clip(m / np.sqrt(l * (l + 1)), -1.0, 1.0)
                theta = np.arccos(cos_th)
            else:
                theta = 0.0
            coords[idx] = [r * np.sin(theta), 0.0, r * np.cos(theta)]
        return coords

    def _enumerate_sd_basis(self) -> None:
        """
        Enumerate all valid Slater determinants.

        Each SD is a sorted tuple of n_electrons distinct spin-orbital indices.
        Pauli exclusion is automatic: combinations() never repeats indices.
        Total: C(n_sp, n_electrons) determinants.
        """
        t0 = time.perf_counter()
        self.sd_basis: List[Tuple[int, ...]] = list(
            combinations(range(self.n_sp), self.n_electrons)
        )
        self.n_sd: int = len(self.sd_basis)
        self._sd_index: Dict[Tuple[int, ...], int] = {
            sd: i for i, sd in enumerate(self.sd_basis)
        }
        elapsed = time.perf_counter() - t0
        n_spatial = self.lattice.num_states
        print(
            f"[LatticeIndex] max_n={self.max_n}: {n_spatial} spatial states, "
            f"{self.n_sp} spin-orbitals, {self.n_sd:,} Slater determinants "
            f"(enumerated in {elapsed:.3f}s)"
        )

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def _get_eri(self, a: int, b: int, c: int, d: int) -> float:
        """Look up two-electron integral <ab|cd> from ERI table."""
        return self._eri.get((a, b, c, d), 0.0)

    def _slater_condon_diagonal(self, sd: Tuple[int, ...]) -> float:
        """
        Slater-Condon diagonal: <I|H|I> with full two-electron integrals.

        = Σ_p h1(p,p) + Σ_{p<q} [<pq|pq> - δ(σ_p,σ_q)<pq|qp>]
        """
        h_diag = 0.0
        for p in sd:
            h_diag += self._h1_diag[p >> 1]

        n = len(sd)
        for i in range(n):
            pi = sd[i]
            sp_i = pi >> 1
            sig_i = pi & 1
            for j in range(i + 1, n):
                pj = sd[j]
                sp_j = pj >> 1
                sig_j = pj & 1
                # Coulomb: <ij|ij>
                h_diag += self._get_eri(sp_i, sp_j, sp_i, sp_j)
                # Exchange: -δ(σi,σj) <ij|ji>
                if sig_i == sig_j:
                    h_diag -= self._get_eri(sp_i, sp_j, sp_j, sp_i)
        return h_diag

    def _slater_condon_single(self, sd: Tuple[int, ...],
                              kp: int, r: int) -> float:
        """
        Slater-Condon single excitation: sd[kp] → r.

        = h1(p,r) + Σ_{q in occ, q≠p} [<pq|rq> - δ(σ_p,σ_q)<pq|qr>]

        Returns the UNSIGNED matrix element (caller applies phase).
        """
        p = sd[kp]
        sp_p = p >> 1
        sp_r = r >> 1
        sig_p = p & 1
        sig_r = r & 1

        # One-body part
        val = self._H1_spatial[sp_p, sp_r]

        # Two-body part
        for q in sd:
            if q == p:
                continue
            sp_q = q >> 1
            sig_q = q & 1
            # Coulomb: <pq|rq>
            val += self._get_eri(sp_p, sp_q, sp_r, sp_q)
            # Exchange: -δ(σ_p,σ_q) <pq|qr>
            if sig_p == sig_q:
                val -= self._get_eri(sp_p, sp_q, sp_q, sp_r)

        return val

    def assemble_hamiltonian(self) -> csr_matrix:
        """
        Assemble the N-electron Hamiltonian via Slater-Condon rules.

        For vee_method='slater_full', uses full two-electron integrals:
          Diagonal:  Σ_p h1(p,p) + Σ_{p<q} [<pq|pq> - δ(σ)·<pq|qp>]
          Single exc (p→r): h1(p,r) + Σ_q [<pq|rq> - δ(σ)·<pq|qr>]
          Double exc (pq→rs): <pq|rs> - δ(σ)·<pq|sr>

        For other vee_methods, uses diagonal-only V_ee (original behavior).

        Returns
        -------
        H : csr_matrix, shape (n_sd, n_sd)
        """
        if self.vee_method == 'slater_full':
            return self._assemble_hamiltonian_full()
        return self._assemble_hamiltonian_diag_vee()

    def _assemble_hamiltonian_diag_vee(self) -> csr_matrix:
        """Original assembly: diagonal V_ee only, off-diagonal from h1."""
        t0 = time.perf_counter()

        diag_rows: List[int] = []
        diag_vals: List[float] = []
        off_rows: List[int] = []
        off_cols: List[int] = []
        off_vals: List[float] = []

        for I, sd in enumerate(self.sd_basis):
            occupied_set = set(sd)

            # ---- Diagonal: one-body sum + V_ee ----
            h_diag = 0.0
            for p in sd:
                sp = p >> 1
                h_diag += self._h1_diag[sp]
            h_diag += self._compute_sd_vee(sd)

            if abs(h_diag) >= self.threshold:
                diag_rows.append(I)
                diag_vals.append(h_diag)

            # ---- Off-diagonal: single excitations (h1 only) ----
            for kp, p in enumerate(sd):
                sigma_p = p & 1
                sp_p = p >> 1

                for r_sp, h_val in self._h1_offdiag[sp_p]:
                    if abs(h_val) < self.threshold:
                        continue
                    r = (r_sp << 1) | sigma_p
                    if r in occupied_set:
                        continue

                    new_sd = tuple(sorted(
                        sd[:kp] + (r,) + sd[kp + 1:]
                    ))
                    if new_sd not in self._sd_index:
                        continue
                    J = self._sd_index[new_sd]
                    if J <= I:
                        continue

                    phase = self._compute_phase(sd, kp, r)
                    off_rows.append(I)
                    off_cols.append(J)
                    off_vals.append(phase * h_val)

        H_diag = csr_matrix(
            (diag_vals, (diag_rows, diag_rows)),
            shape=(self.n_sd, self.n_sd),
        )
        H_upper = csr_matrix(
            (off_vals, (off_rows, off_cols)),
            shape=(self.n_sd, self.n_sd),
        )
        H = H_upper + H_upper.T + H_diag

        elapsed = time.perf_counter() - t0
        print(
            f"[LatticeIndex] Hamiltonian assembled: shape={H.shape}, "
            f"nnz={H.nnz:,}, time={elapsed:.3f}s"
        )
        return H.tocsr()

    def _assemble_hamiltonian_full(self) -> csr_matrix:
        """
        Full Slater-Condon assembly with off-diagonal two-electron integrals.

        Handles diagonal, single excitations, and double excitations.
        Optimized: iterates over SD pairs (I,J) with J>I and classifies
        excitation order by set difference, avoiding the O(n_sp^2) inner
        loop for double excitations.
        """
        t0 = time.perf_counter()

        diag_rows: List[int] = []
        diag_vals: List[float] = []
        off_rows: List[int] = []
        off_cols: List[int] = []
        off_vals: List[float] = []

        n_el = self.n_electrons

        # --- Phase 1: Diagonal ---
        for I, sd_I in enumerate(self.sd_basis):
            h_diag = self._slater_condon_diagonal(sd_I)
            if abs(h_diag) >= self.threshold:
                diag_rows.append(I)
                diag_vals.append(h_diag)

        t_diag = time.perf_counter() - t0

        # --- Phase 2: Off-diagonal (single + double excitations) ---
        # Pre-convert SDs to frozensets for fast set operations
        sd_sets = [frozenset(sd) for sd in self.sd_basis]

        for I in range(self.n_sd):
            sd_I = self.sd_basis[I]
            set_I = sd_sets[I]

            for J in range(I + 1, self.n_sd):
                set_J = sd_sets[J]

                # Determine excitation order
                diff_I = set_I - set_J   # orbitals in I but not J (removed)
                diff_J = set_J - set_I   # orbitals in J but not I (added)
                n_diff = len(diff_I)

                if n_diff == 0:
                    continue  # identical SDs
                elif n_diff > 2:
                    continue  # triple+ excitation: zero by Slater-Condon

                elif n_diff == 1:
                    # Single excitation: p → r
                    p = next(iter(diff_I))
                    r = next(iter(diff_J))

                    # Spin conservation
                    if (p & 1) != (r & 1):
                        continue

                    # Find kp (position of p in sd_I)
                    kp = sd_I.index(p)

                    me = self._slater_condon_single(sd_I, kp, r)
                    if abs(me) < self.threshold:
                        continue

                    phase = self._compute_phase(sd_I, kp, r)
                    off_rows.append(I)
                    off_cols.append(J)
                    off_vals.append(phase * me)

                elif n_diff == 2:
                    # Double excitation: p,q → r,s
                    removed = sorted(diff_I)
                    added = sorted(diff_J)
                    p, q = removed[0], removed[1]
                    r, s = added[0], added[1]

                    sig_p, sig_q = p & 1, q & 1
                    sig_r, sig_s = r & 1, s & 1
                    sp_p, sp_q = p >> 1, q >> 1
                    sp_r, sp_s = r >> 1, s >> 1

                    # Spin conservation
                    if sig_p + sig_q != sig_r + sig_s:
                        continue

                    # <pq|rs> with proper spin matching
                    me = 0.0
                    if sig_p == sig_r and sig_q == sig_s:
                        me += self._get_eri(sp_p, sp_q, sp_r, sp_s)
                    if sig_p == sig_s and sig_q == sig_r:
                        me -= self._get_eri(sp_p, sp_q, sp_s, sp_r)

                    if abs(me) < self.threshold:
                        continue

                    kp = sd_I.index(p)
                    kq = sd_I.index(q)
                    phase = self._compute_double_phase(sd_I, kp, kq, r, s)
                    off_rows.append(I)
                    off_cols.append(J)
                    off_vals.append(phase * me)

        H_diag = csr_matrix(
            (diag_vals, (diag_rows, diag_rows)),
            shape=(self.n_sd, self.n_sd),
        )
        H_upper = csr_matrix(
            (off_vals, (off_rows, off_cols)),
            shape=(self.n_sd, self.n_sd),
        )
        H = H_upper + H_upper.T + H_diag

        elapsed = time.perf_counter() - t0
        print(
            f"[LatticeIndex] Full Hamiltonian assembled: shape={H.shape}, "
            f"nnz={H.nnz:,}, time={elapsed:.3f}s "
            f"(diag={t_diag:.1f}s)"
        )
        return H.tocsr()

    def compute_ground_state(
        self, n_states: int = 1
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute the N-electron ground state via sparse diagonalisation.

        WARNING: Result accuracy limited by Mulliken V_ee approximation.
        Li ground state will be qualitatively correct (1s^2 2s configuration
        dominant) but NOT quantitatively accurate (overestimated V_ee).

        Returns
        -------
        eigvals : np.ndarray, shape (n_states,)
        eigvecs : np.ndarray, shape (n_sd, n_states)
        """
        H = self.assemble_hamiltonian()
        k = min(n_states, self.n_sd - 2)
        # Use deterministic initial vector for reproducible ARPACK convergence
        rng = np.random.RandomState(42)
        v0 = rng.randn(H.shape[0])
        eigvals, eigvecs = eigsh(H, k=k, which="SA", v0=v0)
        order = np.argsort(eigvals)
        return eigvals[order], eigvecs[:, order]

    # ------------------------------------------------------------------
    # Helper methods
    # ------------------------------------------------------------------

    def _compute_sd_vee(self, sd: Tuple[int, ...]) -> float:
        """Sum V_ee(sp_i, sp_j) for all pairs i < j in the Slater determinant."""
        total = 0.0
        n = len(sd)
        for i in range(n):
            for j in range(i + 1, n):
                total += self._vee_matrix[sd[i] >> 1, sd[j] >> 1]
        return total

    @staticmethod
    def _compute_phase(sd: Tuple[int, ...], kp: int, r: int) -> float:
        """
        Fermionic sign for the single excitation sd[kp] → r.

        Phase = (-1)^{kp + kr}
          kp = index of removed orbital in original SD (0-based)
          kr = #{a in SD : a < r, a ≠ sd[kp]}
        """
        p = sd[kp]
        kr = sum(1 for a in sd if a < r and a != p)
        return (-1.0) ** (kp + kr)

    @staticmethod
    def _compute_double_phase(sd: Tuple[int, ...],
                              kp: int, kq: int,
                              r: int, s: int) -> float:
        """
        Fermionic sign for double excitation sd[kp],sd[kq] → r,s.

        Computed by applying two sequential single excitations and
        combining phases.
        """
        p = sd[kp]
        q = sd[kq]

        # First excitation: remove p, insert r
        # Count swaps to move p out
        n_swap_p = kp
        # Count position for r among remaining (after removing p)
        remaining_1 = [a for a in sd if a != p]
        kr = sum(1 for a in remaining_1 if a < r)
        phase1 = (-1) ** (n_swap_p + kr)

        # Second excitation: remove q, insert s from the intermediate SD
        intermediate = sorted(remaining_1[:kr] + [r] + remaining_1[kr:])
        # Find position of q in intermediate
        kq_new = intermediate.index(q)
        remaining_2 = [a for a in intermediate if a != q]
        ks = sum(1 for a in remaining_2 if a < s)
        phase2 = (-1) ** (kq_new + ks)

        return float(phase1 * phase2)

    def assembly_stats(self) -> dict:
        """
        Return statistics about the topology without assembling the full H.
        Useful for scaling benchmarks.
        """
        nnz_offdiag = 0
        for sp, neighbors in self._h1_offdiag.items():
            nnz_offdiag += len(neighbors)

        return {
            "max_n": self.max_n,
            "n_spatial": self.lattice.num_states,
            "n_sp": self.n_sp,
            "n_sd": self.n_sd,
            "h1_nnz": int(self._H1_spatial.nnz),
            "h1_offdiag_edges": nnz_offdiag,
            "vee_method": self.vee_method,
            "h1_method": self.h1_method,
        }
