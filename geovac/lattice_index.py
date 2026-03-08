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
from scipy.integrate import quad

try:
    from .lattice import GeometricLattice
except ImportError:
    from lattice import GeometricLattice

KINETIC_SCALE: float = -1.0 / 16.0   # Universal topological constant


# ---------------------------------------------------------------------------
# S³ density-overlap V_ee (Paper 7, Section VI)
# ---------------------------------------------------------------------------

def _phi_s_orbital(n: int, t: float) -> float:
    """
    Fock-projected density for s-orbital (l=0) with principal quantum number n.

    Phi_n(t) = rho_tilde_{ns,ns}(2Z*t) evaluated with Z=1 (Z scaling handled
    separately by the master formula prefactor 4Z/pi).

    Parameters
    ----------
    n : int
        Principal quantum number (n >= 1)
    t : float
        Dimensionless momentum transfer t = q/(2Z)

    Returns
    -------
    float
        Projected density value at t
    """
    if n == 1:
        # Phi_1s(t) = 1/(1+t^2)^2
        return 1.0 / (1.0 + t * t) ** 2
    elif n == 2:
        # Phi_2s(t) = (1-12t^2+32t^4)/(4t^2+1)^4
        t2 = t * t
        return (1.0 - 12.0 * t2 + 32.0 * t2 * t2) / (4.0 * t2 + 1.0) ** 4
    else:
        # General n: compute from hydrogen density FT with p0=2 (Z=1 units)
        # rho_ns(q) for hydrogen, then substitute q=2t
        # Use numerical evaluation via position-space radial wavefunctions
        return _phi_s_orbital_general(n, t)


def _phi_s_orbital_general(n: int, t: float) -> float:
    """
    General s-orbital projected density via numerical FT of |R_{n0}(r)|^2.

    rho_ns(q) = integral_0^inf |R_{n0}(r)|^2 * sin(qr)/(qr) * r^2 * 4pi dr
    Then Phi_n(t) = rho_ns(2t) with Z=1.

    For n=1,2 the closed-form expressions in _phi_s_orbital are preferred.
    """
    from scipy.special import assoc_laguerre
    q = 2.0 * t  # momentum in Z=1 units, p0=2
    # R_{n0}(r) = 2*(1/n)^{3/2} * L_{n-1}^1(2r/n) * exp(-r/n) for Z=1
    # |R_{n0}|^2 r^2 is the radial probability density
    def integrand(r: float) -> float:
        x = 2.0 * r / n
        lag = assoc_laguerre(x, n - 1, 1)
        R = 2.0 * (1.0 / n) ** 1.5 * lag * np.exp(-r / n)
        rho_r = R * R * r * r
        if q < 1e-14:
            return rho_r  # sin(qr)/(qr) -> 1
        qr = q * r
        return rho_r * np.sin(qr) / qr
    val, _ = quad(integrand, 0, np.inf, limit=300)
    return val * 4.0 * np.pi


def _form_factor_nl(n: int, l: int, Z: float, q: float) -> float:
    """
    Spherically averaged charge form factor for hydrogenic (n,l) orbital.

    rho_tilde(q) = integral_0^inf |R_{nl,Z}(r)|^2 * sin(qr)/(qr) * r^2 dr

    For s-orbitals (l=0) with n<=2, uses closed-form expressions.
    Otherwise, uses numerical integration of hydrogenic radial wavefunctions.

    The form factor satisfies rho_tilde(0) = 1 (normalization) and
    rho_tilde(q) -> 0 as q -> inf.

    Parameters
    ----------
    n : int
        Principal quantum number
    l : int
        Angular momentum quantum number
    Z : float
        Nuclear charge
    q : float
        Momentum transfer in atomic units (1/bohr)

    Returns
    -------
    float
        Form factor value (dimensionless)
    """
    # Closed-form for n=1,2 s-orbitals (verified against known Slater F0)
    if l == 0 and n <= 2:
        return _phi_s_orbital(n, q / (2.0 * Z))

    # General (n, l, Z) via numerical integration of |R_{nl,Z}(r)|²
    from scipy.special import assoc_laguerre
    from math import factorial

    # Hydrogenic radial wavefunction R_{nl,Z}(r):
    # R_{nl}(r) = N * (2Zr/n)^l * L_{n-l-1}^{2l+1}(2Zr/n) * exp(-Zr/n)
    # N = (2Z/n)^{3/2} * sqrt((n-l-1)! / (2n * ((n+l)!)^3)) * (n+l)!
    nr = n - l - 1  # radial quantum number
    N_sq = (2.0 * Z / n) ** 3 * factorial(nr) / (2.0 * n * (factorial(n + l)) ** 3)
    N_sq *= factorial(n + l) ** 2  # from the associated Laguerre normalization
    # Full: N² = (2Z/n)³ × (n-l-1)! / (2n × (n+l)!) — standard Griffiths convention

    def integrand(r: float) -> float:
        if r < 1e-30:
            return 0.0
        x = 2.0 * Z * r / n
        lag = assoc_laguerre(x, nr, 2 * l + 1)
        # |R_{nl}|² = N² × x^{2l} × L² × exp(-x)
        R_sq = N_sq * x ** (2 * l) * lag * lag * np.exp(-x)
        rho_r = R_sq * r * r  # |R_{nl}|² r²
        if q < 1e-14:
            return rho_r  # sin(qr)/(qr) -> 1
        qr = q * r
        return rho_r * np.sin(qr) / qr

    val, _ = quad(integrand, 0, np.inf, limit=300)
    return val


def compute_vee_s3_overlap(
    n_a: int,
    l_a: int,
    n_b: int,
    l_b: int,
    Z: float,
    p0: float | None = None,
) -> float:
    """
    Compute the direct Coulomb (Slater F0) integral via S3 density-overlap.

    V_ee is a NODE property on S3: each orbital contributes a density profile
    Phi(t) and F0 = (4Z/pi) * integral_0^inf Phi_a(t) * Phi_b(t) dt.

    This is NOT an edge property -- do not use chordal distance d2_chord for
    l>0 orbitals (overestimates 1s-2s by ~30x).

    Currently implemented for s-orbital pairs (l_a = l_b = 0).
    For l>0, the full 3D angular convolution with spherical harmonic multipoles
    is required (deferred to future work -- see Paper 7 Section VI.E).

    Parameters
    ----------
    n_a, l_a : int
        Quantum numbers for orbital a
    n_b, l_b : int
        Quantum numbers for orbital b
    Z : float
        Nuclear charge
    p0 : float or None
        Fock projection scale (default: 2Z for s-orbitals)

    Returns
    -------
    float
        Slater F0 direct Coulomb integral in Hartree atomic units

    Raises
    ------
    NotImplementedError
        If l_a > 0 or l_b > 0 (angular extension pending, see Paper 7 Sec VI.E)
    """
    if l_a > 0 or l_b > 0:
        raise NotImplementedError(
            f"S3 density-overlap V_ee not implemented for l>0 orbitals "
            f"(got l_a={l_a}, l_b={l_b}). The full 3D angular convolution "
            f"with spherical harmonic multipoles is required. "
            f"See Paper 7 Section VI.E. Do NOT fall back to chordal d2 -- "
            f"it overestimates inter-shell integrals by ~30x."
        )

    def integrand(t: float) -> float:
        return _phi_s_orbital(n_a, t) * _phi_s_orbital(n_b, t)

    integral, _ = quad(integrand, 0, np.inf, limit=200)
    return (4.0 * Z / np.pi) * integral


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
        fci_method: str = 'auto',
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
        elif vee_method == 's3_overlap':
            warnings.warn(
                "LatticeIndex V_ee uses S3 density-overlap formula: "
                "F0(a,b) = (4Z/pi) int Phi_a(t) Phi_b(t) dt (Paper 7 Sec VI). "
                "V_ee is a NODE property (density overlap), NOT an edge property. "
                "s-orbital pairs only; l>0 raises NotImplementedError.",
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
        # Default h1_method: 'exact' for slater_full/s3_overlap, 'hybrid' for slater, 'graph' otherwise
        if h1_method:
            self.h1_method = h1_method
        elif vee_method in ('slater_full', 's3_overlap'):
            self.h1_method = 'exact'
        elif vee_method == 'slater':
            self.h1_method = 'hybrid'
        else:
            self.h1_method = 'graph'

        self.fci_method = fci_method

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

            ξ₁ = fac × sin(θ) × cos(φ)
            ξ₂ = fac × sin(θ) × sin(φ)
            ξ₃ = fac × cos(θ)
            ξ₄ = (1 − n²) / (n²+1)

        where  fac = 2n / (n²+1)  and the momentum direction (θ, φ) is:

            θ  encodes the z-projection:
                l=0  →  θ=0
                l>0  →  cos(θ) = m / √(l(l+1))

            φ  encodes the orbital angular momentum magnitude:
                l=0  →  φ=0  (irrelevant: sin(θ)=0)
                l>0  →  φ = 2π × l / n

        The azimuthal angle φ breaks the degeneracy between states with the
        same (n, m=0) but different l, which previously all mapped to
        θ=π/2, φ=0 → identical S³ coordinates.  Distributing l-subshells
        around the azimuthal circle ensures each state gets a unique S³ point.

        All ξ satisfy |ξ|=1.

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
                phi = 2.0 * np.pi * l / n
            else:
                cos_th = 1.0
                sin_th = 0.0
                phi = 0.0
            fac = 2.0 * n / (n * n + 1.0)
            coords[idx, 0] = fac * sin_th * np.cos(phi)        # ξ₁
            coords[idx, 1] = fac * sin_th * np.sin(phi)        # ξ₂
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

    def _compute_s3_overlap(self) -> np.ndarray:
        """
        Compute V_ee via S3 density-overlap formula (Paper 7, Section VI).

        F0(a,b) = (4Z/pi) * integral_0^inf Phi_a(t) * Phi_b(t) dt

        where t = q/(2Z) is dimensionless momentum transfer and Phi_a(t) is the
        Fock-projected orbital density on S3.

        V_ee is a NODE property (density overlap on S3), NOT an edge property
        (pairwise chordal distance). The chordal ansatz kappa/d2 overestimates
        inter-shell integrals (e.g. 1s-2s) by ~30x.

        Currently supports s-orbital pairs only (l=0). For l>0 pairs, falls
        back to the position-space Slater F0 integral (_compute_f0_integrals).

        Returns
        -------
        vee : ndarray, shape (n_spatial, n_spatial)
        """
        Z = float(self.Z)
        n_sp = self.lattice.num_states
        vee = np.zeros((n_sp, n_sp))
        unique_nl = sorted(set((n, l) for n, l, m in self.lattice.states))

        # Pre-compute F0 for each unique (n1,l1)-(n2,l2) pair
        f0_cache: Dict[Tuple[int, int, int, int], float] = {}
        for n1, l1 in unique_nl:
            for n2, l2 in unique_nl:
                key = (n1, l1, n2, l2)
                key_rev = (n2, l2, n1, l1)
                if key in f0_cache or key_rev in f0_cache:
                    continue
                if l1 == 0 and l2 == 0:
                    # Use S3 density-overlap formula
                    val = compute_vee_s3_overlap(n1, l1, n2, l2, Z)
                else:
                    # l>0: fall back to position-space Slater F0
                    computed = self._compute_f0_integrals(Z, {key})
                    val = computed[key]
                f0_cache[key] = val
                f0_cache[key_rev] = val

        # Fill matrix
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

        's3_overlap' (Paper 7, Sec VI):
            F⁰(a,b) = (4Z/π) ∫₀^∞ Φ_a(t)·Φ_b(t) dt, where Φ_a is the
            Fock-projected orbital density on S³. NODE property, not edge.
            s-orbital pairs via closed-form densities; l>0 falls back to
            position-space Slater F⁰.

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

        elif self.vee_method == 's3_overlap':
            self._vee_matrix = self._compute_s3_overlap()

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

        Uses the method selected by ``fci_method``:
        - ``'auto'``: direct CI for N_SD >= 5000, matrix otherwise
        - ``'direct'``: excitation-driven assembly (O(N_SD x N_connected))
        - ``'matrix'``: explicit pairwise assembly (O(N^2_SD), reference)

        Returns
        -------
        eigvals : np.ndarray, shape (n_states,)
        eigvecs : np.ndarray, shape (n_sd, n_states)
        """
        method = self.fci_method
        if method == 'auto':
            method = 'direct' if self.n_sd >= 5000 else 'matrix'

        if method == 'direct' and self.vee_method == 'slater_full':
            from .direct_ci import DirectCISolver
            solver = DirectCISolver(self)
            return solver.solve(n_states=n_states)

        # Fall back to matrix method (original implementation)
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


# ===========================================================================
# MolecularLatticeIndex — Heteronuclear diatomic FCI
# ===========================================================================


class MolecularLatticeIndex:
    """
    FCI solver for diatomic molecules using two atomic hydrogenic lattices.

    Constructs the molecular Hamiltonian as:
        H_mol = H1_mol + V_ee(same-atom + cross-atom s-orbital) + V_NN

    The one-electron Hamiltonian combines:
    - Exact atomic eigenvalues -Z²/(2n²) on each atom (diagonal)
    - Cross-nuclear attraction via electrostatic potential (s-orbitals only)
    - Inter-atomic kinetic hopping via graph bridges (off-diagonal)

    V_ee combines exact same-atom Slater integrals with cross-atom direct
    Coulomb J_AB via Fourier convolution of S³ momentum-space densities.
    Only s-orbital cross-atom pairs are included (l>0 requires angular-
    dependent form factors for variational rigour).

    Parameters
    ----------
    Z_A, Z_B : int
        Nuclear charges of atoms A and B
    nmax_A, nmax_B : int
        Basis size for each atom (can differ)
    R : float
        Internuclear distance in Bohr
    n_electrons : int
        Total electron count (= Z_A + Z_B for neutral molecule)
    n_bridges : int
        Number of inter-atomic bridge connections (default 20)
    vee_method : str
        Two-electron integral method ('slater_full' recommended)
    fci_method : str
        FCI assembly method ('auto', 'direct', or 'matrix')
    kinetic_scale : float
        Universal kinetic scale factor (default -1/16)
    """

    # Reuse LatticeIndex static methods for fermionic phase computation
    _compute_phase = staticmethod(LatticeIndex._compute_phase)
    _compute_double_phase = staticmethod(LatticeIndex._compute_double_phase)

    def __init__(
        self,
        Z_A: int,
        Z_B: int,
        nmax_A: int,
        nmax_B: int,
        R: float,
        n_electrons: int,
        n_bridges: int = 20,
        vee_method: str = 'slater_full',
        fci_method: str = 'auto',
        kinetic_scale: float = KINETIC_SCALE,
        sparsity_threshold: float = 1e-8,
    ) -> None:
        self.Z_A = Z_A
        self.Z_B = Z_B
        self.nmax_A = nmax_A
        self.nmax_B = nmax_B
        self.R = R
        self.n_electrons = n_electrons
        self.n_bridges = n_bridges
        self.vee_method = vee_method
        self.fci_method = fci_method
        self.kinetic_scale = kinetic_scale
        self.threshold = sparsity_threshold

        # Ghost atom support: Z=0 means basis functions are present
        # but carry no nuclear attraction (Boys-Bernardi counterpoise).
        self._ghost_A = (Z_A == 0)
        self._ghost_B = (Z_B == 0)

        # --- Build atomic lattice indices for ERI computation ---
        # For ghost atoms (Z=0), we still need the lattice topology (states,
        # adjacency) but set Z=1 internally so GeometricLattice doesn't fail.
        # The h1_diag will be overwritten to 0 in _build_molecular_h1.
        Z_A_eff = Z_A if Z_A > 0 else 1
        Z_B_eff = Z_B if Z_B > 0 else 1
        label = "ghost" if (self._ghost_A or self._ghost_B) else "LiH"
        print(f"[MolecularLatticeIndex] {label}: Z_A={Z_A}, Z_B={Z_B}, "
              f"nmax_A={nmax_A}, nmax_B={nmax_B}, R={R:.2f} bohr, "
              f"Ne={n_electrons}")

        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self._li_A = LatticeIndex(
                n_electrons=1, max_n=nmax_A, nuclear_charge=Z_A_eff,
                vee_method=vee_method, h1_method='exact',
                fci_method='matrix', kinetic_scale=kinetic_scale,
            )
            self._li_B = LatticeIndex(
                n_electrons=1, max_n=nmax_B, nuclear_charge=Z_B_eff,
                vee_method=vee_method, h1_method='exact',
                fci_method='matrix', kinetic_scale=kinetic_scale,
            )

        self._n_spatial_A = self._li_A.lattice.num_states
        self._n_spatial_B = self._li_B.lattice.num_states
        self._n_spatial = self._n_spatial_A + self._n_spatial_B
        self.n_sp = 2 * self._n_spatial

        # --- Build molecular data structures ---
        self._build_combined_states()
        self._build_combined_adjacency()
        self._build_molecular_h1()
        self._build_molecular_vee()
        self._enumerate_sd_basis()

        # Nuclear repulsion (exact; zero for ghost atoms)
        self.V_NN = float(Z_A * Z_B) / R if (Z_A > 0 and Z_B > 0) else 0.0
        print(f"[MolecularLatticeIndex] V_NN = {self.V_NN:.6f} Ha")

    # ------------------------------------------------------------------
    # Construction methods
    # ------------------------------------------------------------------

    def _build_combined_states(self) -> None:
        """Build combined spin-orbital state list from both atoms."""
        self.sp_states: List[Tuple[int, int, int, int]] = []
        self._spatial_atom: List[int] = []  # which atom each spatial state belongs to

        for n, l, m in self._li_A.lattice.states:
            self.sp_states.append((n, l, m, 0))
            self.sp_states.append((n, l, m, 1))
            self._spatial_atom.append(0)

        for n, l, m in self._li_B.lattice.states:
            self.sp_states.append((n, l, m, 0))
            self.sp_states.append((n, l, m, 1))
            self._spatial_atom.append(1)

    def _build_combined_adjacency(self) -> None:
        """
        Build combined adjacency matrix: block-diagonal intra-atom edges
        plus inter-atomic bridges with STO overlap and conformal weighting.
        """
        from scipy.sparse import block_diag as sp_block_diag, coo_matrix

        A_A = self._li_A.lattice.adjacency
        A_B = self._li_B.lattice.adjacency
        nA = self._n_spatial_A
        N = self._n_spatial

        # Block-diagonal intra-atom adjacency
        block = sp_block_diag([A_A, A_B], format='coo')
        rows = list(block.row)
        cols = list(block.col)
        data = list(block.data.astype(np.float64))

        # Inter-atomic bridges (replicates MoleculeHamiltonian bridge logic)
        lat_A = self._li_A.lattice
        lat_B = self._li_B.lattice

        boundary_A = lat_A._get_boundary_states_prioritized()
        boundary_B = lat_B._get_boundary_states_prioritized()
        n_actual = 0  # reset; only set if bridges are actually created

        # Skip bridges if either atom is a ghost (Z=0): no physical hopping
        n_candidate = min(len(boundary_A), len(boundary_B), self.n_bridges)
        if n_candidate > 0 and self.R > 1e-10 and not self._ghost_A and not self._ghost_B:
            R_AB = self.R
            S_R = (1.0 + R_AB + R_AB**2 / 3.0) * np.exp(-R_AB)

            V_nn = float(self.Z_A * self.Z_B) / R_AB
            p0_A = np.sqrt(float(self.Z_A)**2 + V_nn)
            p0_B = np.sqrt(float(self.Z_B)**2 + V_nn)

            local_A = np.array(boundary_A[:n_candidate], dtype=np.intp)
            local_B = np.array(boundary_B[:n_candidate], dtype=np.intp)

            n_vals_A = np.array([lat_A.states[k][0] for k in local_A],
                                dtype=np.float64)
            n_vals_B = np.array([lat_B.states[k][0] for k in local_B],
                                dtype=np.float64)

            p_A = float(self.Z_A) / n_vals_A
            p_B = float(self.Z_B) / n_vals_B
            omega_A = 2.0 * p0_A / (p_A**2 + p0_A**2)
            omega_B = 2.0 * p0_B / (p_B**2 + p0_B**2)

            weights = S_R * omega_A * omega_B

            mask = np.abs(weights) >= 1e-8
            local_A = local_A[mask]
            local_B_offset = local_B[mask] + nA
            weights = weights[mask]
            n_actual = int(mask.sum())

            # Symmetric bridges
            rows.extend(local_A.tolist())
            rows.extend(local_B_offset.tolist())
            cols.extend(local_B_offset.tolist())
            cols.extend(local_A.tolist())
            data.extend(weights.tolist())
            data.extend(weights.tolist())

        self._adjacency_combined = coo_matrix(
            (data, (rows, cols)), shape=(N, N)
        ).tocsr()
        self._n_bridges_actual = n_actual
        print(f"[MolecularLatticeIndex] Bridges: {n_actual} "
              f"(requested {self.n_bridges})")

    def _build_molecular_h1(self) -> None:
        """
        Build one-electron Hamiltonian for the molecule.

        Diagonal: exact atomic eigenvalues -Z²/(2n²) + cross-nuclear attraction
        (electrostatic potential, s-orbitals only for variational balance).
        Off-diagonal: kinetic hopping from combined adjacency (with bridges).
        """
        n = self._n_spatial
        nA = self._n_spatial_A

        # Diagonal: exact atomic eigenvalues (zero for ghost atoms)
        h1_diag = np.zeros(n)
        if not self._ghost_A:
            for i, (ni, li, mi) in enumerate(self._li_A.lattice.states):
                h1_diag[i] = -float(self.Z_A)**2 / (2.0 * ni**2)
        if not self._ghost_B:
            for j, (nj, lj, mj) in enumerate(self._li_B.lattice.states):
                h1_diag[nA + j] = -float(self.Z_B)**2 / (2.0 * nj**2)

        # Cross-nuclear attraction via electrostatic potential.
        # Only s-orbitals (l=0) get cross-nuclear attraction to match the
        # scope of cross-atom V_ee (also s-orbital only).  Including
        # cross-nuclear attraction for p-orbitals WITHOUT compensating
        # cross-atom V_ee would violate the variational bound.
        # Skip entirely if either atom is a ghost (Z=0 → no attraction).
        if not self._ghost_A and not self._ghost_B:
            for i, (ni, li, mi) in enumerate(self._li_A.lattice.states):
                if li == 0:
                    h1_diag[i] += self._fourier_cross_attraction(
                        ni, li, self.Z_A, self.Z_B, self.R)
            for j, (nj, lj, mj) in enumerate(self._li_B.lattice.states):
                if lj == 0:
                    h1_diag[nA + j] += self._fourier_cross_attraction(
                        nj, lj, self.Z_B, self.Z_A, self.R)

        self._h1_diag = h1_diag

        # Off-diagonal: kinetic hopping from combined graph Laplacian
        H1_offdiag = self.kinetic_scale * (-self._adjacency_combined)
        self._H1_spatial = (diags(h1_diag) + H1_offdiag).tocsr()

        # Build offdiag dict for matrix assembly path
        H1_coo = H1_offdiag.tocoo()
        self._h1_offdiag: Dict[int, List[Tuple[int, float]]] = {
            i: [] for i in range(n)
        }
        for r, c, v in zip(H1_coo.row, H1_coo.col, H1_coo.data):
            if r != c and abs(v) >= self.threshold:
                self._h1_offdiag[r].append((c, float(v)))

        # Report H1 statistics
        n_cross_A = sum(1 for i in range(nA) if h1_diag[i] !=
                        -float(self.Z_A)**2 / (2.0 * self._li_A.lattice.states[i][0]**2))
        print(f"[MolecularLatticeIndex] H1: {n} spatial states, "
              f"off-diag nnz={H1_offdiag.nnz}")

    @staticmethod
    def _mulliken_cross_attraction(
        n: int, l: int, Z_self: int, Z_other: int, R: float
    ) -> float:
        """
        Cross-nuclear attraction via Mulliken approximation with variational cap.

        Replicates MoleculeHamiltonian._build_cross_nuclear_attraction logic:
            V_cross = max(-Z_other/R * S_eff, -Z_other * Z_self / n²)
        where the cap prevents variational collapse for diffuse orbitals.
        """
        if R < 1e-10:
            return 0.0
        R_eff = R * float(Z_self) / (n**2)
        S_1s = np.exp(-R_eff) * (1.0 + R_eff + R_eff**2 / 3.0)
        ang = 1.0 / (2 * l + 1)
        S_eff = min(S_1s * ang, 1.0)
        v_mulliken = (-float(Z_other) / R) * S_eff
        v_limit = -float(Z_other) * float(Z_self) / (n**2)
        return max(v_mulliken, v_limit)

    @staticmethod
    def _fourier_cross_attraction(
        n: int, l: int, Z_self: int, Z_other: int, R: float
    ) -> float:
        """
        Cross-nuclear attraction via direct electrostatic potential.

        For the spherically averaged orbital density, the potential of
        nucleus B at distance R from the electron's own nucleus is:

            V = -Z_other × [(1/R) ∫₀^R |R_{nl}(r)|² r² dr
                           + ∫_R^∞ |R_{nl}(r)|² r dr]

        Exact for the spherically averaged density. Correct limits:
        - V → -Z_other/R as R → ∞ (since ∫|R|²r²dr = 1)
        - Smooth, finite at all R > 0

        Uses direct r-space integration (fast, single quadrature).

        Parameters
        ----------
        n, l : int
            Quantum numbers of the orbital on atom A
        Z_self : int
            Nuclear charge of atom A (where the electron sits)
        Z_other : int
            Nuclear charge of atom B (the attracting nucleus)
        R : float
            Internuclear distance in Bohr

        Returns
        -------
        float
            Cross-nuclear attraction energy in Hartree (negative)
        """
        if R < 1e-10:
            return 0.0

        from scipy.integrate import quad
        from scipy.special import assoc_laguerre
        from math import factorial

        Z = float(Z_self)

        # Build |R_{nl,Z}(r)|² using closed-form for common cases
        if n == 1 and l == 0:
            def R_sq(r: float) -> float:
                return 4.0 * Z**3 * np.exp(-2.0 * Z * r)
        elif n == 2 and l == 0:
            def R_sq(r: float) -> float:
                x = Z * r
                return (Z**3 / 8.0) * (2.0 - x)**2 * np.exp(-x)
        elif n == 2 and l == 1:
            def R_sq(r: float) -> float:
                x = Z * r
                return (Z**3 / 24.0) * x**2 * np.exp(-x)
        else:
            nr = n - l - 1
            N_sq = (2.0 * Z / n)**3 * factorial(nr) / (2.0 * n * factorial(n + l))

            def R_sq(r: float) -> float:
                x = 2.0 * Z * r / n
                lag = assoc_laguerre(x, nr, 2 * l + 1)
                return N_sq * x**(2 * l) * lag * lag * np.exp(-x)

        # Electrostatic potential: Φ(R) = (1/R)∫₀^R ρr²dr + ∫_R^∞ ρr dr
        inner, _ = quad(lambda r: R_sq(r) * r * r, 0, R, limit=100)
        outer, _ = quad(lambda r: R_sq(r) * r, R, np.inf, limit=100)

        return -float(Z_other) * (inner / R + outer)

    def _build_molecular_vee(self) -> None:
        """
        Build combined V_ee matrix and ERI table.

        Same-atom Slater integrals are exact (cached). Cross-atom two-electron
        integrals use Fourier convolution of momentum-space densities with
        sin(qR)/(qR) structure factor (Paper 7, Section VI):

            J_AB(a,b,R) = (2/π) ∫₀^∞ ρ̃_A(q) ρ̃_B(q) sin(qR)/(qR) dq

        where ρ̃(q) is the spherically averaged charge form factor.
        Note: no q² factor — the q² volume element cancels 1/q² from Coulomb.

        Only s-orbital (l=0) cross-atom pairs are included. The l>0 pairs
        require angular-dependent form factors; spherical averaging
        underestimates repulsion and violates the variational bound.

        Only direct Coulomb J is included; exchange K is neglected
        (overestimates repulsion slightly, preserving variational bound).

        Known limitation: cross-nuclear attraction is also restricted to
        s-orbitals for balance (see _build_molecular_h1). The resulting
        Hamiltonian has a ~0.2% variational violation at equilibrium, likely
        due to incomplete screening of s-orbital cross-nuclear attraction
        by cross-atom V_ee in the FCI ground state.
        """
        n = self._n_spatial
        nA = self._n_spatial_A

        self._vee_matrix = np.zeros((n, n))
        # Only include same-atom V_ee for real (non-ghost) atoms
        if not self._ghost_A:
            self._vee_matrix[:nA, :nA] = self._li_A._vee_matrix
        if not self._ghost_B:
            self._vee_matrix[nA:, nA:] = self._li_B._vee_matrix

        # ERI table: reindex atom B entries with offset
        self._eri: Dict[Tuple[int, int, int, int], float] = {}

        if not self._ghost_A and hasattr(self._li_A, '_eri'):
            for key, val in self._li_A._eri.items():
                self._eri[key] = val

        if not self._ghost_B and hasattr(self._li_B, '_eri'):
            for (a, b, c, d), val in self._li_B._eri.items():
                self._eri[(a + nA, b + nA, c + nA, d + nA)] = val

        # Cross-atom V_ee via Fourier convolution (skip for ghost atoms)
        n_cross = 0
        if not self._ghost_A and not self._ghost_B:
            n_cross = self._build_cross_atom_vee()

        n_eri = len(self._eri)
        print(f"[MolecularLatticeIndex] V_ee: {n_eri} ERI entries "
              f"({n_cross} cross-atom)")

    def _build_cross_atom_vee(self) -> int:
        """
        Compute cross-atom direct Coulomb integrals J_AB(a,b,R) via
        Fourier convolution of momentum-space form factors.

        J_AB(a,b,R) = (2/π) ∫₀^∞ ρ̃_A(q) ρ̃_B(q) sin(qR)/(qR) dq

        Uses adaptive quadrature (scipy.integrate.quad) with closed-form
        form factors for s-orbitals. Only s-orbital pairs (l=0) are included;
        l>0 cross-atom V_ee requires angular-dependent form factors that
        preserve variational rigour (spherical averaging violates the
        variational bound — see Paper 7 Section VI.E).

        Only direct Coulomb J is computed; exchange K is neglected
        (overestimates repulsion slightly, preserves variational bound).

        Returns
        -------
        int
            Number of cross-atom ERI entries added.
        """
        nA = self._n_spatial_A
        states_A = self._li_A.lattice.states
        states_B = self._li_B.lattice.states
        R = self.R
        count = 0

        # Cache J by (n_a, n_b) since l=0 for all pairs, m doesn't matter
        j_cache: Dict[Tuple[int, int], float] = {}

        for i_a, (n_a, l_a, m_a) in enumerate(states_A):
            if l_a > 0:
                continue  # s-orbital pairs only
            for i_b, (n_b, l_b, m_b) in enumerate(states_B):
                if l_b > 0:
                    continue  # s-orbital pairs only
                j_b = i_b + nA

                cache_key = (n_a, n_b)
                if cache_key in j_cache:
                    j_ab = j_cache[cache_key]
                else:
                    def integrand(q: float, _na=n_a, _nb=n_b) -> float:
                        rho_a = _phi_s_orbital(_na, q / (2.0 * self.Z_A))
                        rho_b = _phi_s_orbital(_nb, q / (2.0 * self.Z_B))
                        qR = q * R
                        if qR < 1e-10:
                            sinc = 1.0
                        else:
                            sinc = np.sin(qR) / qR
                        return rho_a * rho_b * sinc

                    result, _ = quad(integrand, 0, np.inf, limit=300)
                    j_ab = (2.0 / np.pi) * result
                    j_cache[cache_key] = j_ab

                if abs(j_ab) < 1e-15:
                    continue

                self._vee_matrix[i_a, j_b] = j_ab
                self._vee_matrix[j_b, i_a] = j_ab
                self._eri[(i_a, j_b, i_a, j_b)] = j_ab
                self._eri[(j_b, i_a, j_b, i_a)] = j_ab
                count += 2

        return count

    def _enumerate_sd_basis(self) -> None:
        """Enumerate all N-electron Slater determinants over combined basis."""
        t0 = time.perf_counter()
        self.sd_basis: List[Tuple[int, ...]] = list(
            combinations(range(self.n_sp), self.n_electrons)
        )
        self.n_sd: int = len(self.sd_basis)
        self._sd_index: Dict[Tuple[int, ...], int] = {
            sd: i for i, sd in enumerate(self.sd_basis)
        }
        elapsed = time.perf_counter() - t0
        print(
            f"[MolecularLatticeIndex] {self._n_spatial} spatial, "
            f"{self.n_sp} spin-orbitals, {self.n_sd:,} SDs "
            f"(enumerated in {elapsed:.3f}s)"
        )

    # ------------------------------------------------------------------
    # ERI access (duck-type compatible with LatticeIndex)
    # ------------------------------------------------------------------

    def _get_eri(self, a: int, b: int, c: int, d: int) -> float:
        """Look up two-electron integral <ab|cd> from ERI table."""
        return self._eri.get((a, b, c, d), 0.0)

    # ------------------------------------------------------------------
    # Hamiltonian assembly (matrix method, for small systems)
    # ------------------------------------------------------------------

    def assemble_hamiltonian(self) -> csr_matrix:
        """
        Assemble the molecular FCI Hamiltonian via Slater-Condon rules.

        Uses full two-electron integrals (same-atom approximation).
        """
        if self.vee_method == 'slater_full':
            return self._assemble_hamiltonian_full()
        return self._assemble_hamiltonian_diag_vee()

    def _assemble_hamiltonian_full(self) -> csr_matrix:
        """Full Slater-Condon assembly with off-diagonal two-electron integrals."""
        t0 = time.perf_counter()

        diag_rows: List[int] = []
        diag_vals: List[float] = []
        off_rows: List[int] = []
        off_cols: List[int] = []
        off_vals: List[float] = []

        # Diagonal
        for I, sd_I in enumerate(self.sd_basis):
            h_diag = 0.0
            for p in sd_I:
                h_diag += self._h1_diag[p >> 1]
            n = len(sd_I)
            for i in range(n):
                pi = sd_I[i]
                sp_i = pi >> 1
                sig_i = pi & 1
                for j in range(i + 1, n):
                    pj = sd_I[j]
                    sp_j = pj >> 1
                    h_diag += self._get_eri(sp_i, sp_j, sp_i, sp_j)
                    if (pj & 1) == sig_i:
                        h_diag -= self._get_eri(sp_i, sp_j, sp_j, sp_i)
            if abs(h_diag) >= self.threshold:
                diag_rows.append(I)
                diag_vals.append(h_diag)

        t_diag = time.perf_counter() - t0

        # Off-diagonal (single + double excitations)
        sd_sets = [frozenset(sd) for sd in self.sd_basis]
        for I in range(self.n_sd):
            sd_I = self.sd_basis[I]
            set_I = sd_sets[I]
            for J in range(I + 1, self.n_sd):
                set_J = sd_sets[J]
                diff_I = set_I - set_J
                diff_J = set_J - set_I
                n_diff = len(diff_I)
                if n_diff == 0 or n_diff > 2:
                    continue
                if n_diff == 1:
                    p = next(iter(diff_I))
                    r = next(iter(diff_J))
                    if (p & 1) != (r & 1):
                        continue
                    kp = sd_I.index(p)
                    sp_p = p >> 1
                    sp_r = r >> 1
                    me = self._H1_spatial[sp_p, sp_r]
                    for q in sd_I:
                        if q == p:
                            continue
                        sp_q = q >> 1
                        me += self._get_eri(sp_p, sp_q, sp_r, sp_q)
                        if (p & 1) == (q & 1):
                            me -= self._get_eri(sp_p, sp_q, sp_q, sp_r)
                    if abs(me) < self.threshold:
                        continue
                    phase = self._compute_phase(sd_I, kp, r)
                    off_rows.append(I)
                    off_cols.append(J)
                    off_vals.append(phase * me)
                else:
                    removed = sorted(diff_I)
                    added = sorted(diff_J)
                    p, q = removed
                    r, s = added
                    if (p & 1) + (q & 1) != (r & 1) + (s & 1):
                        continue
                    sp_p, sp_q = p >> 1, q >> 1
                    sp_r, sp_s = r >> 1, s >> 1
                    me = 0.0
                    if (p & 1) == (r & 1) and (q & 1) == (s & 1):
                        me += self._get_eri(sp_p, sp_q, sp_r, sp_s)
                    if (p & 1) == (s & 1) and (q & 1) == (r & 1):
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
            (diag_vals, (diag_rows, diag_rows)), shape=(self.n_sd, self.n_sd))
        H_upper = csr_matrix(
            (off_vals, (off_rows, off_cols)), shape=(self.n_sd, self.n_sd))
        H = H_upper + H_upper.T + H_diag

        elapsed = time.perf_counter() - t0
        print(f"[MolecularLatticeIndex] H assembled: shape={H.shape}, "
              f"nnz={H.nnz:,}, time={elapsed:.3f}s")
        return H.tocsr()

    def _assemble_hamiltonian_diag_vee(self) -> csr_matrix:
        """Diagonal-only V_ee assembly (for non-slater_full methods)."""
        t0 = time.perf_counter()
        diag_rows: List[int] = []
        diag_vals: List[float] = []
        off_rows: List[int] = []
        off_cols: List[int] = []
        off_vals: List[float] = []

        for I, sd in enumerate(self.sd_basis):
            occupied_set = set(sd)
            h_diag = 0.0
            for p in sd:
                h_diag += self._h1_diag[p >> 1]
            # Diagonal V_ee
            n = len(sd)
            for i in range(n):
                for j in range(i + 1, n):
                    h_diag += self._vee_matrix[sd[i] >> 1, sd[j] >> 1]
            if abs(h_diag) >= self.threshold:
                diag_rows.append(I)
                diag_vals.append(h_diag)
            for kp, p in enumerate(sd):
                sigma_p = p & 1
                sp_p = p >> 1
                for r_sp, h_val in self._h1_offdiag[sp_p]:
                    if abs(h_val) < self.threshold:
                        continue
                    r = (r_sp << 1) | sigma_p
                    if r in occupied_set:
                        continue
                    new_sd = tuple(sorted(sd[:kp] + (r,) + sd[kp + 1:]))
                    J = self._sd_index.get(new_sd)
                    if J is None or J <= I:
                        continue
                    phase = self._compute_phase(sd, kp, r)
                    off_rows.append(I)
                    off_cols.append(J)
                    off_vals.append(phase * h_val)

        H_diag = csr_matrix(
            (diag_vals, (diag_rows, diag_rows)), shape=(self.n_sd, self.n_sd))
        H_upper = csr_matrix(
            (off_vals, (off_rows, off_cols)), shape=(self.n_sd, self.n_sd))
        H = H_upper + H_upper.T + H_diag
        elapsed = time.perf_counter() - t0
        print(f"[MolecularLatticeIndex] H assembled: shape={H.shape}, "
              f"nnz={H.nnz:,}, time={elapsed:.3f}s")
        return H.tocsr()

    # ------------------------------------------------------------------
    # Solver
    # ------------------------------------------------------------------

    def compute_ground_state(
        self, n_states: int = 1
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute the molecular ground state via FCI.

        Adds nuclear repulsion V_NN = Z_A * Z_B / R to eigenvalues.

        Returns
        -------
        eigvals : np.ndarray, shape (n_states,)
            Total energies (electronic + nuclear repulsion)
        eigvecs : np.ndarray, shape (n_sd, n_states)
        """
        method = self.fci_method
        if method == 'auto':
            method = 'direct' if self.n_sd >= 5000 else 'matrix'

        if method == 'direct' and self.vee_method == 'slater_full':
            from .direct_ci import DirectCISolver
            solver = DirectCISolver(self)
            eigvals, eigvecs = solver.solve(n_states=n_states)
        else:
            H = self.assemble_hamiltonian()
            k = min(n_states, self.n_sd - 2)
            rng = np.random.RandomState(42)
            v0 = rng.randn(H.shape[0])
            eigvals, eigvecs = eigsh(H, k=k, which="SA", v0=v0)
            order = np.argsort(eigvals)
            eigvals, eigvecs = eigvals[order], eigvecs[:, order]

        # Add nuclear repulsion
        eigvals = eigvals + self.V_NN
        print(f"[MolecularLatticeIndex] E_elec = {eigvals[0] - self.V_NN:.6f} Ha, "
              f"V_NN = {self.V_NN:.6f} Ha, E_total = {eigvals[0]:.6f} Ha")
        return eigvals, eigvecs


# ---------------------------------------------------------------------------
# Boys-Bernardi Counterpoise Correction for BSSE
# ---------------------------------------------------------------------------

def compute_bsse_correction(
    Z_A: int,
    Z_B: int,
    nmax_A: int,
    nmax_B: int,
    R: float,
    n_electrons_A: int,
    n_electrons_B: int,
    vee_method: str = 'slater_full',
    fci_method: str = 'auto',
) -> dict:
    """
    Boys-Bernardi counterpoise correction for Basis Set Superposition Error.

    Computes four energies:
      E_A_own   = atom A with its own basis (nmax_A)
      E_B_own   = atom B with its own basis (nmax_B)
      E_A_ghost = atom A in full molecular basis (nmax_A + ghost nmax_B)
      E_B_ghost = atom B in full molecular basis (ghost nmax_A + nmax_B)

    BSSE = (E_A_ghost - E_A_own) + (E_B_ghost - E_B_own)
         = artificial lowering due to basis borrowing (negative number)

    Parameters
    ----------
    Z_A, Z_B : int
        Nuclear charges
    nmax_A, nmax_B : int
        Basis sizes for each atom
    R : float
        Internuclear distance in Bohr (defines ghost orbital placement)
    n_electrons_A, n_electrons_B : int
        Electron counts for each atom
    vee_method : str
        Two-electron integral method (default 'slater_full')
    fci_method : str
        FCI assembly method (default 'auto')

    Returns
    -------
    dict
        Keys: E_A_own, E_B_own, E_A_ghost, E_B_ghost, BSSE, BSSE_A, BSSE_B
    """
    import warnings

    # Own-basis energies
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        li_A = LatticeIndex(
            n_electrons=n_electrons_A, max_n=nmax_A, nuclear_charge=Z_A,
            vee_method=vee_method, h1_method='exact', fci_method=fci_method,
        )
        E_A_own = li_A.compute_ground_state(n_states=1)[0][0]

        li_B = LatticeIndex(
            n_electrons=n_electrons_B, max_n=nmax_B, nuclear_charge=Z_B,
            vee_method=vee_method, h1_method='exact', fci_method=fci_method,
        )
        E_B_own = li_B.compute_ground_state(n_states=1)[0][0]

    # Ghost-basis energies: atom A with ghost B orbitals
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mol_A_ghost = MolecularLatticeIndex(
            Z_A=Z_A, Z_B=0, nmax_A=nmax_A, nmax_B=nmax_B,
            R=R, n_electrons=n_electrons_A,
            vee_method=vee_method, fci_method=fci_method,
        )
        E_A_ghost = mol_A_ghost.compute_ground_state(n_states=1)[0][0]

    # Ghost-basis energies: atom B with ghost A orbitals
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mol_B_ghost = MolecularLatticeIndex(
            Z_A=0, Z_B=Z_B, nmax_A=nmax_A, nmax_B=nmax_B,
            R=R, n_electrons=n_electrons_B,
            vee_method=vee_method, fci_method=fci_method,
        )
        E_B_ghost = mol_B_ghost.compute_ground_state(n_states=1)[0][0]

    BSSE_A = E_A_ghost - E_A_own
    BSSE_B = E_B_ghost - E_B_own
    BSSE = BSSE_A + BSSE_B

    return {
        'E_A_own': E_A_own,
        'E_B_own': E_B_own,
        'E_A_ghost': E_A_ghost,
        'E_B_ghost': E_B_ghost,
        'BSSE': BSSE,
        'BSSE_A': BSSE_A,
        'BSSE_B': BSSE_B,
    }
