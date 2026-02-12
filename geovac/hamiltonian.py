"""
Helium Atom Hamiltonian on Geometric Lattice

Implements a 2-electron Hamiltonian using tensor product space over
the discrete quantum state lattice (n, l, m).

H_total = H₁ ⊗ I + I ⊗ H₁ + V_ee

where H₁ is the single-particle Hamiltonian and V_ee is electron-electron repulsion.

Author: Computational Quantum Physics
Date: February 2026
"""

import numpy as np
import scipy.sparse as sp
from scipy.sparse import csr_matrix, diags, identity, kron
from scipy.sparse.linalg import eigsh
from typing import Tuple
from .lattice import GeometricLattice


# Physical Constants (Atomic Units)
HARTREE_TO_EV = 27.2114  # Conversion factor
HELIUM_Z = 2  # Nuclear charge for Helium


class HeliumHamiltonian:
    """
    Two-electron Hamiltonian for Helium atom on geometric lattice.
    
    Uses tensor product space: |ψ⟩ = |n₁,l₁,m₁⟩ ⊗ |n₂,l₂,m₂⟩
    
    Attributes:
    -----------
    lattice : GeometricLattice
        Single-particle lattice
    h1 : scipy.sparse.csr_matrix
        Single-particle Hamiltonian
    h2 : scipy.sparse.csr_matrix
        Two-particle Hamiltonian (full system)
    """
    
    def __init__(self, max_n: int, Z: int = HELIUM_Z, kinetic_scale: float = 1.0):
        """
        Initialize Helium Hamiltonian.
        
        Parameters:
        -----------
        max_n : int
            Maximum principal quantum number for lattice
        Z : int, optional
            Nuclear charge (default: 2 for Helium)
        kinetic_scale : float, optional
            Scaling factor for kinetic energy (default: 1.0)
            Used to calibrate graph Laplacian to physical kinetic energy
        """
        self.max_n = max_n
        self.Z = Z
        self.kinetic_scale = kinetic_scale
        
        # Build single-particle lattice
        print(f"Building lattice for max_n={max_n}...")
        self.lattice = GeometricLattice(max_n)
        self.n_states = self.lattice.num_states
        
        print(f"  → {self.n_states} single-particle states")
        print(f"  → {self.n_states**2} two-particle states")
        
        # Build Hamiltonians
        self.h1 = None
        self.h2 = None
        
        self._build_single_particle_hamiltonian()
        self._build_two_particle_hamiltonian()
    
    def _build_single_particle_hamiltonian(self) -> None:
        """
        Construct single-particle Hamiltonian: H₁ = T + V
        
        Kinetic Energy (T):
        - Use graph Laplacian from adjacency matrix
        - T = -½ ∇² ≈ -½ L where L is graph Laplacian
        - L = D - A (degree matrix - adjacency matrix)
        
        Potential Energy (V):
        - Coulomb potential: V = -Z/r
        - Approximate r ≈ n² (Bohr radius scaling)
        - V_ii = -Z/n²
        """
        print("\nBuilding single-particle Hamiltonian...")
        
        # === Kinetic Energy: Graph Laplacian ===
        # Degree matrix (diagonal): D_ii = sum of row i in adjacency
        adjacency = self.lattice.adjacency
        degree = np.array(adjacency.sum(axis=1)).flatten()
        D = diags(degree, 0, shape=(self.n_states, self.n_states), format='csr')
        
        # Graph Laplacian: L = D - A
        laplacian = D - adjacency
        
        # Kinetic energy operator: T = -½ L × kinetic_scale
        # (Factor of ½ is the reduced mass in atomic units)
        # kinetic_scale calibrates graph Laplacian to physical kinetic energy
        T = -0.5 * self.kinetic_scale * laplacian
        
        # === Potential Energy: Coulomb Attraction ===
        # V = -Z/r with r ≈ n²
        potential = np.zeros(self.n_states)
        for idx, (n, l, m) in enumerate(self.lattice.states):
            r_eff = n**2  # Effective radius
            potential[idx] = -self.Z / r_eff
        
        V = diags(potential, 0, shape=(self.n_states, self.n_states), format='csr')
        
        # Total single-particle Hamiltonian
        self.h1 = T + V
        
        print(f"  ✓ Kinetic energy (graph Laplacian): {self.n_states}×{self.n_states}")
        print(f"  ✓ Potential energy (Coulomb): diagonal")
        print(f"  ✓ H₁ matrix: {self.h1.shape}, {self.h1.nnz} nonzero")
    
    def _compute_spatial_coordinates(self) -> np.ndarray:
        """
        Compute 3D spatial coordinates for each quantum state.
        
        Uses spherical coordinates:
        - r = n² (Bohr radius scaling)
        - θ = arccos(m/√(l(l+1))) if l>0, else 0
        - φ = 0 (mean field approximation)
        
        Returns:
        --------
        coords : np.ndarray, shape (n_states, 3)
            Cartesian coordinates (x, y, z)
        """
        coords = np.zeros((self.n_states, 3))
        
        for idx, (n, l, m) in enumerate(self.lattice.states):
            r = n**2  # Radial distance
            
            # Polar angle from magnetic quantum number
            if l > 0:
                # θ = arccos(m/√(l(l+1)))
                l_magnitude = np.sqrt(l * (l + 1))
                cos_theta = m / l_magnitude
                # Clamp to [-1, 1] to avoid numerical issues
                cos_theta = np.clip(cos_theta, -1.0, 1.0)
                theta = np.arccos(cos_theta)
            else:
                theta = 0.0
            
            # Azimuthal angle (mean field approximation)
            phi = 0.0
            
            # Convert to Cartesian
            x = r * np.sin(theta) * np.cos(phi)
            y = r * np.sin(theta) * np.sin(phi)
            z = r * np.cos(theta)
            
            coords[idx] = [x, y, z]
        
        return coords
    
    def _build_electron_repulsion(self) -> csr_matrix:
        """
        Construct electron-electron repulsion term: V_ee = 1/|r₁ - r₂|
        
        In tensor product space, this is a diagonal operator acting on
        combined states |n₁,l₁,m₁⟩ ⊗ |n₂,l₂,m₂⟩.
        
        The interaction energy is estimated as:
        V_ee(i,j) = 1/|r_i - r_j|
        
        where r_i and r_j are the spatial coordinates of states i and j.
        
        Returns:
        --------
        v_ee : scipy.sparse.csr_matrix
            Diagonal matrix of electron-electron repulsion
        """
        print("\nBuilding electron-electron repulsion...")
        
        # Get spatial coordinates for all states
        coords = self._compute_spatial_coordinates()
        
        # Build diagonal interaction matrix
        n_two_particle = self.n_states**2
        v_ee_diagonal = np.zeros(n_two_particle)
        
        # Iterate over all pairs of single-particle states
        for i in range(self.n_states):
            for j in range(self.n_states):
                # Combined state index in tensor product space
                idx_combined = i * self.n_states + j
                
                # Euclidean distance between the two electrons
                r1 = coords[i]
                r2 = coords[j]
                distance = np.linalg.norm(r1 - r2)
                
                # Avoid division by zero (same state)
                if distance < 1e-10:
                    # Self-interaction: use average self-repulsion
                    # Approximate as 1/n² for state n
                    n1 = self.lattice.states[i][0]
                    v_ee_diagonal[idx_combined] = 1.0 / (n1**2)
                else:
                    # Coulomb repulsion: 1/r₁₂
                    v_ee_diagonal[idx_combined] = 1.0 / distance
        
        v_ee = diags(v_ee_diagonal, 0, shape=(n_two_particle, n_two_particle), 
                     format='csr')
        
        print(f"  ✓ Electron-electron repulsion: {n_two_particle}×{n_two_particle} diagonal")
        print(f"  ✓ Mean repulsion energy: {np.mean(v_ee_diagonal):.4f} Hartree")
        
        return v_ee
    
    def _build_two_particle_hamiltonian(self) -> None:
        """
        Construct two-particle Hamiltonian using tensor products.
        
        H₂ = H₁ ⊗ I + I ⊗ H₁ + V_ee
        
        where:
        - H₁ ⊗ I: First electron kinetic + potential
        - I ⊗ H₁: Second electron kinetic + potential
        - V_ee: Electron-electron repulsion
        
        Uses scipy.sparse.kron for efficient sparse tensor products.
        """
        print("\nBuilding two-particle Hamiltonian...")
        
        # Identity matrix for tensor products
        I = identity(self.n_states, format='csr')
        
        # First electron Hamiltonian: H₁ ⊗ I
        print("  → Computing H₁ ⊗ I...")
        h1_x_I = kron(self.h1, I, format='csr')
        
        # Second electron Hamiltonian: I ⊗ H₁
        print("  → Computing I ⊗ H₁...")
        I_x_h1 = kron(I, self.h1, format='csr')
        
        # Electron-electron repulsion
        print("  → Computing V_ee...")
        v_ee = self._build_electron_repulsion()
        
        # Total Hamiltonian
        print("  → Assembling total Hamiltonian...")
        self.h2 = h1_x_I + I_x_h1 + v_ee
        
        # Statistics
        n_total = self.n_states**2
        sparsity = 1.0 - (self.h2.nnz / (n_total**2))
        
        print(f"\n  ✓ Two-particle Hamiltonian complete:")
        print(f"      Shape:      {self.h2.shape}")
        print(f"      Nonzero:    {self.h2.nnz}")
        print(f"      Sparsity:   {sparsity:.6f}")
        print(f"      Memory:     {self.h2.data.nbytes / 1e6:.2f} MB")
    
    def compute_ground_state(self, n_states: int = 1) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute ground state energy and wavefunction.
        
        Uses sparse eigenvalue solver (Lanczos algorithm) to find
        the lowest energy eigenstate.
        
        Parameters:
        -----------
        n_states : int, optional
            Number of lowest eigenstates to compute (default: 1)
        
        Returns:
        --------
        energies : np.ndarray
            Eigenvalues (energies) in ascending order
        wavefunctions : np.ndarray
            Corresponding eigenvectors (columns)
        """
        print(f"\nComputing {n_states} lowest eigenstate(s)...")
        print("  (This may take a moment for large matrices...)")
        
        # Use eigsh for symmetric/Hermitian matrices
        # which='SA' finds smallest algebraic (most negative)
        energies, wavefunctions = eigsh(self.h2, k=n_states, which='SA')
        
        print(f"  ✓ Eigenvalue computation complete")
        
        return energies, wavefunctions
    
    def __repr__(self) -> str:
        return (f"HeliumHamiltonian(max_n={self.max_n}, Z={self.Z}, "
                f"single_states={self.n_states}, "
                f"two_particle_dim={self.n_states**2})")


def format_energy(energy: float, label: str = "") -> str:
    """Format energy value in both Hartree and eV."""
    ev = energy * HARTREE_TO_EV
    return f"{label}{energy:.6f} Hartree ({ev:.4f} eV)"


if __name__ == "__main__":
    """
    Benchmark: Compute Helium ground state energy and compare to experiment.
    """
    print("=" * 70)
    print("HELIUM ATOM HAMILTONIAN - GEOMETRIC LATTICE BENCHMARK")
    print("=" * 70)
    
    # Experimental reference value
    HELIUM_EXPERIMENT = -2.903  # Hartree
    
    # Build Hamiltonian for max_n=3
    print("\n[1/3] Building Hamiltonian...")
    print("-" * 70)
    hamiltonian = HeliumHamiltonian(max_n=3, Z=2)
    
    # Matrix properties
    print("\n[2/3] Matrix Properties:")
    print("-" * 70)
    n_single = hamiltonian.n_states
    n_two = n_single**2
    print(f"  Single-particle states:  {n_single}")
    print(f"  Two-particle states:     {n_two}")
    print(f"  Matrix shape:            {hamiltonian.h2.shape}")
    print(f"  Matrix nonzero:          {hamiltonian.h2.nnz}")
    print(f"  Matrix sparsity:         {1 - hamiltonian.h2.nnz/n_two**2:.4f}")
    
    # Compute ground state
    print("\n[3/3] Computing Ground State Energy...")
    print("-" * 70)
    energies, wavefunctions = hamiltonian.compute_ground_state(n_states=3)
    
    # Results
    ground_state_energy = energies[0]
    error = ground_state_energy - HELIUM_EXPERIMENT
    error_percent = 100 * abs(error / HELIUM_EXPERIMENT)
    
    print("\n" + "=" * 70)
    print("RESULTS")
    print("=" * 70)
    print(f"\n  Ground State Energy:")
    print(f"    Calculated:      {format_energy(ground_state_energy)}")
    print(f"    Experimental:    {format_energy(HELIUM_EXPERIMENT)}")
    print(f"    Error:           {error:+.6f} Hartree ({error_percent:.2f}%)")
    
    print(f"\n  Excited States:")
    for i, E in enumerate(energies[1:], start=1):
        print(f"    State {i+1}:         {format_energy(E)}")
    
    print(f"\n  Wavefunction:")
    psi_0 = wavefunctions[:, 0]
    print(f"    Dimension:       {len(psi_0)}")
    print(f"    Norm:            {np.linalg.norm(psi_0):.6f}")
    print(f"    Max amplitude:   {np.max(np.abs(psi_0)):.6f}")
    
    # Physical interpretation
    print("\n" + "=" * 70)
    print("INTERPRETATION")
    print("=" * 70)
    print(f"""
  This calculation uses a geometric lattice approximation where:
  - Lattice size: max_n = {hamiltonian.max_n}
  - Graph Laplacian approximates kinetic energy operator
  - Bohr radius scaling: r ≈ n²
  - Spherical coordinate mapping for electron-electron distance
  
  The error of {error_percent:.2f}% is expected due to:
  1. Small basis set (only {n_single} single-particle states)
  2. Graph Laplacian approximation for kinetic energy
  3. Mean field approximation for azimuthal angle
  4. Lack of correlation corrections
  
  To improve accuracy:
  - Increase max_n (more basis states)
  - Add correlation corrections
  - Use better kinetic energy operator
  - Include electron correlation explicitly
""")
    
    print("=" * 70)
    print("✓ Benchmark complete")
    print("=" * 70)
