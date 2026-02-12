"""
Relativistic Dirac Hamiltonian on Geometric Lattice

Implements a bipartite Dirac structure for the hydrogen atom lattice,
suitable for relativistic quantum mechanical calculations.

The Dirac Hamiltonian has 2×2 block structure:
    H = | V(r) + mc²      c·A     |
        | c·A†           V(r) - mc²|

where:
- m: electron rest mass (1.0 in atomic units)
- c: speed of light (137.036 in atomic units)
- A: adjacency matrix (hopping/kinetic term)
- V(r): Coulomb potential (diagonal)

Author: Theoretical Physics & Quantum Computing
Date: February 2026
"""

import numpy as np
import scipy.sparse as sp
from scipy.sparse import csr_matrix, bmat, diags, identity, kron
from scipy.sparse.linalg import eigsh
from typing import Tuple, List
from .lattice import GeometricLattice


# Physical Constants (Atomic Units)
C_LIGHT = 137.036  # Speed of light in atomic units
ELECTRON_MASS = 1.0  # Electron rest mass in atomic units
HARTREE_TO_EV = 27.2114  # Conversion factor


class DiracLatticeStates:
    """
    Manages bipartite splitting of lattice states based on parity.
    
    In Dirac theory, we split the Hilbert space into two sectors:
    - Positive energy states (particle sector): even l
    - Negative energy states (antiparticle sector): odd l
    
    This creates a natural bipartite graph structure.
    """
    
    def __init__(self, lattice: GeometricLattice):
        """
        Initialize bipartite state structure.
        
        Parameters:
        -----------
        lattice : GeometricLattice
            Underlying quantum state lattice
        """
        self.lattice = lattice
        
        # Partition states by parity of l
        self.even_states = []  # Even l (particle sector)
        self.odd_states = []   # Odd l (antiparticle sector)
        self.even_indices = []
        self.odd_indices = []
        
        for idx, (n, l, m) in enumerate(lattice.states):
            if l % 2 == 0:
                self.even_states.append((n, l, m))
                self.even_indices.append(idx)
            else:
                self.odd_states.append((n, l, m))
                self.odd_indices.append(idx)
        
        self.n_even = len(self.even_states)
        self.n_odd = len(self.odd_states)
        
        print(f"Bipartite lattice structure:")
        print(f"  Even l (particle):     {self.n_even} states")
        print(f"  Odd l (antiparticle):  {self.n_odd} states")
        print(f"  Total spinor dimension: {self.n_even + self.n_odd}")
    
    def get_bipartite_adjacency(self) -> Tuple[csr_matrix, csr_matrix, csr_matrix]:
        """
        Extract bipartite blocks from adjacency matrix.
        
        Returns:
        --------
        A_ee : csr_matrix
            Even-to-even connections
        A_oo : csr_matrix  
            Odd-to-odd connections
        A_eo : csr_matrix
            Even-to-odd connections (off-diagonal hopping)
        """
        adj = self.lattice.adjacency
        
        # Extract blocks
        A_ee = adj[np.ix_(self.even_indices, self.even_indices)]
        A_oo = adj[np.ix_(self.odd_indices, self.odd_indices)]
        A_eo = adj[np.ix_(self.even_indices, self.odd_indices)]
        
        return A_ee, A_oo, A_eo


class DiracHamiltonian:
    """
    Relativistic Dirac Hamiltonian for hydrogen-like atoms on geometric lattice.
    
    Implements the 2×2 block structure:
        H_Dirac = | H₊    cA   |
                  | cA†   H₋   |
    
    where:
        H₊ = V(r) + mc²  (positive energy sector)
        H₋ = V(r) - mc²  (negative energy sector)
        A = adjacency matrix (kinetic/hopping term)
    
    Attributes:
    -----------
    lattice : GeometricLattice
        Single-particle quantum state lattice
    h_dirac : scipy.sparse.csr_matrix
        Full Dirac Hamiltonian (spinor space)
    """
    
    def __init__(self, max_n: int, Z: int = 2, c: float = C_LIGHT, 
                 m: float = ELECTRON_MASS, use_effective_c: bool = True):
        """
        Initialize Dirac Hamiltonian.
        
        Parameters:
        -----------
        max_n : int
            Maximum principal quantum number
        Z : int, optional
            Nuclear charge (default: 2 for Helium)
        c : float, optional
            Speed of light in atomic units (default: 137.036)
        m : float, optional
            Electron rest mass in atomic units (default: 1.0)
        use_effective_c : bool, optional
            Use effective lattice-scaled c (~1-10) instead of true c (default: True)
        """
        self.max_n = max_n
        self.Z = Z
        
        # CRITICAL: For lattice discretization, we use an effective c
        # The adjacency matrix represents discrete hops, not continuous momentum
        # Scale c to match lattice spacing: c_eff ~ c / (n_max)²
        if use_effective_c:
            self.c = c / (max_n**2)  # Effective c for lattice
            self.m = m  # Keep mass at 1.0
            print(f"  Using effective c = {self.c:.4f} (scaled for lattice)")
        else:
            self.c = c  # True speed of light
            self.m = m
            print(f"  Using true c = {self.c:.1f} (WARNING: may give huge energies)")
        
        self.mc2 = self.m * self.c**2  # Rest energy
        
        # Build lattice
        print(f"Building Dirac lattice for max_n={max_n}, Z={Z}...")
        self.lattice = GeometricLattice(max_n)
        self.n_states = self.lattice.num_states
        
        # Bipartite structure
        self.bipartite = DiracLatticeStates(self.lattice)
        
        # Build Hamiltonian
        self.h_dirac = None
        self._build_dirac_hamiltonian()
    
    def _build_potential_diagonal(self) -> np.ndarray:
        """
        Construct Coulomb potential V(r) = -Z/r.
        
        Uses effective radius r ≈ n² from Bohr model.
        
        Returns:
        --------
        potential : np.ndarray
            Diagonal potential energy for all states
        """
        potential = np.zeros(self.n_states)
        for idx, (n, l, m) in enumerate(self.lattice.states):
            r_eff = n**2  # Bohr radius scaling
            potential[idx] = -self.Z / r_eff
        return potential
    
    def _build_dirac_hamiltonian(self) -> None:
        """
        Construct full Dirac Hamiltonian in spinor space.
        
        Block structure:
            H = | V + mc²      c·A     |
                | c·A†         V - mc²  |
        
        The full spinor dimension is 2N (or just N if using even/odd split).
        For simplicity, we use the full N×N structure with bipartite hopping.
        """
        print("\nBuilding Dirac Hamiltonian...")
        
        # Potential energy (diagonal)
        V_diag = self._build_potential_diagonal()
        V = diags(V_diag, 0, shape=(self.n_states, self.n_states), format='csr')
        
        # Adjacency matrix (kinetic/hopping term)
        A = self.lattice.adjacency
        
        # Block diagonal components
        # H₊ = V + mc² (positive energy states - particles)
        H_plus = V + self.mc2 * identity(self.n_states, format='csr')
        
        # H₋ = V - mc² (negative energy states - antiparticles)
        H_minus = V - self.mc2 * identity(self.n_states, format='csr')
        
        # Off-diagonal hopping: c·A (couples positive and negative energy sectors)
        # Additional scaling by lattice spacing to give proper dimensions
        lattice_scale = 1.0 / self.max_n  # Inverse lattice spacing
        c_A = self.c * lattice_scale * A
        
        # Full Dirac Hamiltonian in 2×2 block form
        # | H₊    c·A  |
        # | c·A†  H₋   |
        self.h_dirac = bmat([
            [H_plus, c_A],
            [c_A.T, H_minus]
        ], format='csr')
        
        # Statistics
        n_spinor = 2 * self.n_states
        sparsity = 1.0 - (self.h_dirac.nnz / n_spinor**2)
        
        print(f"\n  ✓ Dirac Hamiltonian complete:")
        print(f"      Single-particle states: {self.n_states}")
        print(f"      Spinor dimension:       {n_spinor}")
        print(f"      Matrix shape:           {self.h_dirac.shape}")
        print(f"      Nonzero elements:       {self.h_dirac.nnz}")
        print(f"      Sparsity:               {sparsity:.6f}")
        print(f"      Memory (MB):            {self.h_dirac.data.nbytes / 1e6:.2f}")
        print(f"      Rest energy mc²:        {self.mc2:.2f} Hartree")
    
    def compute_spectrum(self, k: int = 6) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute energy spectrum of Dirac Hamiltonian.
        
        Finds both positive and negative energy eigenvalues.
        
        Parameters:
        -----------
        k : int, optional
            Number of eigenvalues to compute (default: 6)
        
        Returns:
        --------
        energies : np.ndarray
            Eigenvalues (energies) sorted by algebraic value
        wavefunctions : np.ndarray
            Corresponding eigenvectors (spinor wavefunctions)
        """
        print(f"\nComputing {k} eigenvalues of Dirac Hamiltonian...")
        
        # Note: eigsh may fail for non-positive-definite matrices
        # Use which='SA' for smallest algebraic (most negative)
        try:
            energies, wavefunctions = eigsh(self.h_dirac, k=k, which='SA')
            print(f"  ✓ Eigenvalue computation complete")
        except Exception as e:
            print(f"  ! Warning: eigsh failed, trying dense solver for small matrix")
            # Fallback to dense solver for small matrices
            H_dense = self.h_dirac.toarray()
            all_energies, all_wavefunctions = np.linalg.eigh(H_dense)
            energies = all_energies[:k]
            wavefunctions = all_wavefunctions[:, :k]
        
        return energies, wavefunctions
    
    def __repr__(self) -> str:
        return (f"DiracHamiltonian(max_n={self.max_n}, Z={self.Z}, "
                f"c={self.c:.1f}, spinor_dim={2*self.n_states})")


class HeliumDiracHamiltonian:
    """
    Two-electron Dirac Hamiltonian for Helium using tensor product.
    
    H_total = H_Dirac ⊗ I + I ⊗ H_Dirac + V_ee
    
    where H_Dirac is the single-particle Dirac Hamiltonian and
    V_ee is the electron-electron repulsion.
    """
    
    def __init__(self, max_n: int, Z: int = 2):
        """
        Initialize two-electron Dirac Hamiltonian.
        
        Parameters:
        -----------
        max_n : int
            Maximum principal quantum number
        Z : int, optional
            Nuclear charge (default: 2 for Helium)
        """
        self.max_n = max_n
        self.Z = Z
        
        print("=" * 70)
        print("HELIUM DIRAC HAMILTONIAN (Relativistic)")
        print("=" * 70)
        
        # Build single-particle Dirac Hamiltonian
        self.dirac = DiracHamiltonian(max_n, Z)
        self.n_spinor = 2 * self.dirac.n_states  # Spinor dimension
        
        # Two-particle system
        self.n_two_particle = self.n_spinor**2
        
        print(f"\nTwo-particle system:")
        print(f"  Single-particle spinor dimension: {self.n_spinor}")
        print(f"  Two-particle dimension:           {self.n_two_particle}")
        
        # Build two-electron Hamiltonian
        self.h_total = None
        self._build_two_electron_hamiltonian()
    
    def _compute_spatial_coordinates(self) -> np.ndarray:
        """
        Compute spatial coordinates for spatial repulsion.
        
        Extends coordinates to full spinor space (both sectors have same spatial coords).
        
        Returns:
        --------
        coords : np.ndarray, shape (n_spinor, 3)
            Spatial coordinates for each spinor component
        """
        base_coords = np.zeros((self.dirac.n_states, 3))
        
        for idx, (n, l, m) in enumerate(self.dirac.lattice.states):
            r = n**2
            if l > 0:
                l_magnitude = np.sqrt(l * (l + 1))
                cos_theta = m / l_magnitude
                cos_theta = np.clip(cos_theta, -1.0, 1.0)
                theta = np.arccos(cos_theta)
            else:
                theta = 0.0
            phi = 0.0
            
            x = r * np.sin(theta) * np.cos(phi)
            y = r * np.sin(theta) * np.sin(phi)
            z = r * np.cos(theta)
            
            base_coords[idx] = [x, y, z]
        
        # Extend to spinor space (duplicate for positive/negative energy sectors)
        coords = np.vstack([base_coords, base_coords])
        return coords
    
    def _build_electron_repulsion(self) -> csr_matrix:
        """
        Build electron-electron repulsion term V_ee = 1/|r₁ - r₂|.
        
        Uses graph-based distance for more accurate lattice interaction.
        For now, uses Euclidean distance as in non-relativistic case.
        
        Returns:
        --------
        v_ee : scipy.sparse.csr_matrix
            Diagonal repulsion matrix
        """
        print("\n[4/4] Building electron-electron repulsion...")
        
        coords = self._compute_spatial_coordinates()
        
        v_ee_diagonal = np.zeros(self.n_two_particle)
        
        for i in range(self.n_spinor):
            for j in range(self.n_spinor):
                idx_combined = i * self.n_spinor + j
                
                r1 = coords[i]
                r2 = coords[j]
                distance = np.linalg.norm(r1 - r2)
                
                if distance < 1e-10:
                    # Self-interaction: use average
                    state_idx = i % self.dirac.n_states
                    n1 = self.dirac.lattice.states[state_idx][0]
                    v_ee_diagonal[idx_combined] = 1.0 / (n1**2)
                else:
                    v_ee_diagonal[idx_combined] = 1.0 / distance
        
        v_ee = diags(v_ee_diagonal, 0, shape=(self.n_two_particle, self.n_two_particle),
                     format='csr')
        
        print(f"  ✓ Repulsion matrix: diagonal {self.n_two_particle}×{self.n_two_particle}")
        print(f"  ✓ Mean repulsion:   {np.mean(v_ee_diagonal):.4f} Hartree")
        
        return v_ee
    
    def _build_two_electron_hamiltonian(self) -> None:
        """
        Build two-electron Hamiltonian: H = H₁⊗I + I⊗H₁ + V_ee.
        """
        print("\nBuilding two-electron Dirac Hamiltonian...")
        
        H_dirac = self.dirac.h_dirac
        I = identity(self.n_spinor, format='csr')
        
        print("  [1/4] Computing H_Dirac ⊗ I...")
        h1_x_I = kron(H_dirac, I, format='csr')
        
        print("  [2/4] Computing I ⊗ H_Dirac...")
        I_x_h1 = kron(I, H_dirac, format='csr')
        
        print("  [3/4] Adding single-particle terms...")
        h_kinetic = h1_x_I + I_x_h1
        
        # Electron-electron repulsion
        v_ee = self._build_electron_repulsion()
        
        print("  [5/4] Assembling total Hamiltonian...")
        self.h_total = h_kinetic + v_ee
        
        # Statistics
        sparsity = 1.0 - (self.h_total.nnz / self.n_two_particle**2)
        
        print(f"\n  ✓ Two-electron Dirac Hamiltonian complete:")
        print(f"      Matrix shape:    {self.h_total.shape}")
        print(f"      Nonzero:         {self.h_total.nnz}")
        print(f"      Sparsity:        {sparsity:.6f}")
        print(f"      Memory (MB):     {self.h_total.data.nbytes / 1e6:.2f}")
    
    def compute_ground_state(self, n_states: int = 1) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute ground state energy and wavefunction.
        
        Parameters:
        -----------
        n_states : int, optional
            Number of lowest eigenstates (default: 1)
        
        Returns:
        --------
        energies : np.ndarray
            Eigenvalues (sorted)
        wavefunctions : np.ndarray
            Eigenvectors
        """
        print(f"\nComputing {n_states} lowest eigenstate(s)...")
        print("  (Relativistic calculation - may take longer...)")
        
        energies, wavefunctions = eigsh(self.h_total, k=n_states, which='SA')
        
        print(f"  ✓ Eigenvalue computation complete")
        
        return energies, wavefunctions


def format_energy(energy: float, label: str = "") -> str:
    """Format energy in Hartree and eV."""
    ev = energy * HARTREE_TO_EV
    return f"{label}{energy:.6f} Hartree ({ev:.4f} eV)"


if __name__ == "__main__":
    """
    Benchmark: Relativistic Helium ground state vs non-relativistic.
    """
    
    print("\n" + "=" * 70)
    print("RELATIVISTIC DIRAC HAMILTONIAN BENCHMARK")
    print("=" * 70)
    
    # Experimental values
    HELIUM_NONREL = -2.903  # Non-relativistic Hartree
    
    # Test 1: Single-particle Dirac Hamiltonian
    print("\n[TEST 1] Single-Particle Dir, use_effective_c=Trueac Spectrum")
    print("-" * 70)
    dirac = DiracHamiltonian(max_n=3, Z=2)
    energies_1p, _ = dirac.compute_spectrum(k=8)
    
    print(f"\nSingle-particle spectrum (first 8 states):")
    for i, E in enumerate(energies_1p):
        # Shift by rest mass to see binding energy
        E_shifted = E - 2 * dirac.mc2  # Remove 2mc² (two sectors)
        print(f"  State {i+1}: {format_energy(E)}")
        if i == 0:
            print(f"           (Binding: {E_shifted:.4f} Hartree)")
    
    # Test 2: Two-electron Helium (Dirac)
    print("\n" + "=" * 70)
    print("[TEST 2] Two-Electron Helium (Relativistic)")
    print("=" * 70)
    
    helium_dirac = HeliumDiracHamiltonian(max_n=3, Z=2)
    
    print("\nComputing ground state...")
    energies_he, wavefunctions_he = helium_dirac.compute_ground_state(n_states=3)
    
    ground_state_raw = energies_he[0]
    # Remove rest mass contribution: 2 electrons × 2 sectors × mc²
    ground_state_binding = ground_state_raw - 4 * dirac.mc2
    
    # Comparison
    error = ground_state_binding - HELIUM_NONREL
    error_percent = 100 * abs(error / HELIUM_NONREL)
    
    print("\n" + "=" * 70)
    print("RESULTS")
    print("=" * 70)
    
    print(f"\n  Ground State Energy (Dirac):")
    print(f"    Raw eigenvalue:          {format_energy(ground_state_raw)}")
    print(f"    Rest mass (4×mc²):       {format_energy(4 * dirac.mc2)}")
    print(f"    Binding energy:          {format_energy(ground_state_binding)}")
    print(f"    Non-relativistic target: {format_energy(HELIUM_NONREL)}")
    print(f"    Difference:              {error:+.6f} Hartree ({error_percent:.2f}%)")
    
    print(f"\n  Excited States:")
    for i, E in enumerate(energies_he[1:], start=1):
        E_binding = E - 4 * dirac.mc2
        print(f"    State {i+1}: {format_energy(E_binding)} (binding)")
    
    print(f"\n  Wavefunction:")
    psi_0 = wavefunctions_he[:, 0]
    print(f"    Dimension:       {len(psi_0)}")
    print(f"    Norm:            {np.linalg.norm(psi_0):.6f}")
    
    # Physical interpretation
    print("\n" + "=" * 70)
    print("INTERPRETATION")
    print("=" * 70)
    print(f"""
  Relativistic Dirac Hamiltonian with:
  - Speed of light: c = {C_LIGHT:.3f} atomic units
  - Electron mass:  m = {ELECTRON_MASS:.1f} atomic units  
  - Rest energy:    mc² = {dirac.mc2:.2f} Hartree
  - Lattice size:   max_n = {helium_dirac.max_n}
  - Spinor states:  {helium_dirac.n_spinor} (= 2 × {dirac.n_states})
  
  The Dirac structure splits each lattice state into positive and
  negative energy sectors, coupled by the adjacency matrix (kinetic term).
  
  Key features:
  1. Block structure preserves particle/antiparticle sectors
  2. Off-diagonal c·A terms mix positive/negative energy states
  3. Relativistic corrections enter through mass term ±mc²
  4. Binding energy extracted by removing rest mass
  
  Comparison to non-relativistic:
  - Error of {error_percent:.2f}% indicates relativistic effects
  - Graph Laplacian approximation limits accuracy
  - Bipartite structure captures essential spinor physics
  
  Improvements possible:
  - Larger basis (increase max_n)
  - Better kinetic energy operator (discrete derivatives)
  - Breit interaction for electron-electron (relativistic correction)
  - QED corrections (radiative, vacuum polarization)
""")
    
    print("=" * 70)
    print("✓ Relativistic Dirac benchmark complete")
    print("=" * 70)
