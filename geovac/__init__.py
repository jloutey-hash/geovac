"""
GeoVac: Topological Quantum Chemistry Solver
=============================================

Solves the wave function using Spectral Graph Theory (Sparse Graph Laplacians)
rather than continuous integration. Universal kinetic scale K_vac = -1/16.

**Three Laws of Isoelectronic Scaling (v0.4.1):**
- Law 1 (Conformal): Kinetic energy scales as Z²
- Law 2 (Coulomb): Potential energy scales as Z
- Law 3 (Torsion): Lattice torsion gamma = mu * (Z - Z_ref), mu = 1/4

**Performance:**
- O(N) complexity scaling with >97% matrix sparsity
- Isoelectronic accuracy: Li+ 0.03%, Be2+ 0.15%
- H2+ topological control: <0.1% error
- Ultra-fast: Single atoms in <10ms, molecules in <50ms

**Key Classes:**
- AtomicSolver: Unified atomic solver with isoelectronic scaling
- GeometricLattice: Sparse graph lattice with (n,l,m) quantum state nodes
- MoleculeHamiltonian: Molecular bonding via sparse bridges
- DiracHamiltonian: Relativistic spinor-based solver (experimental)

Quick Start:
-----------
>>> from geovac import AtomicSolver, UNIVERSAL_KINETIC_SCALE
>>> solver = AtomicSolver(max_n=10, Z=1)
>>> E, psi = solver.compute_ground_state()
>>> print(f"Hydrogen: {E[0]:.6f} Ha")

>>> # Isoelectronic scaling (Li+, Be2+, etc.)
>>> solver = AtomicSolver(max_n=10, Z=3)
>>> solver.apply_isoelectronic_scaling()
>>> E, psi = solver.compute_ground_state()
"""

__version__ = '0.5.0'
__author__ = 'J. Loutey'
__license__ = 'MIT'

# Core imports - expose main classes at top level
from .lattice import GeometricLattice
from .hamiltonian import HeliumHamiltonian, HeliumPackingSolver, MoleculeHamiltonian
from .dirac_hamiltonian import DiracHamiltonian, DiracLatticeStates
from .atomic_solver import AtomicSolver, solve_hydrogen, solve_atom

# Holographic/AdS-CFT modules live in ADSCFT/ package:
#   from ADSCFT import MuonicHydrogenSolver, compute_holographic_entropy, etc.

# Define public API
__all__ = [
    'GeometricLattice',
    'HeliumHamiltonian',
    'MoleculeHamiltonian',
    'HeliumPackingSolver',
    'DiracHamiltonian',
    'DiracLatticeStates',
    'AtomicSolver',
    'solve_hydrogen',
    'solve_atom',
    '__version__',
]

# Package metadata
UNIVERSAL_KINETIC_SCALE = -1/16  # Universal topological constant (-0.0625)
CALIBRATED_KINETIC_SCALE = -0.10298808  # Helium-specific (for backward compatibility)
EXPERIMENTAL_HELIUM_ENERGY = -2.90338583  # Hartree (NIST reference)

# Physical constants validated by H/He+/H2+ convergence analysis
HYDROGEN_GROUND_STATE = -0.5  # Hartree (exact)
H2_PLUS_USES_UNIVERSAL_SCALE = True  # 0% error confirms topology is correct
H2_CORRELATION_ERROR = 0.17  # 17% missing from mean-field (expected)

# Convenience function for quick calculations
def solve_helium(max_n=3, kinetic_scale=CALIBRATED_KINETIC_SCALE):
    """
    Convenience function to solve Helium ground state with calibrated parameters.
    
    Parameters
    ----------
    max_n : int, optional
        Maximum principal quantum number for lattice (default: 3)
    kinetic_scale : float, optional
        Calibration factor for graph Laplacian kinetic energy
        (default: -0.103, matches experimental ground state)
    
    Returns
    -------
    energy : float
        Ground state energy in Hartree atomic units
    wavefunction : ndarray
        Ground state wavefunction (normalized)
    
    Examples
    --------
    >>> energy, psi = solve_helium(max_n=3)
    >>> print(f"E₀ = {energy:.6f} Ha")
    E₀ = -2.903000 Ha
    """
    h = HeliumHamiltonian(max_n=max_n, Z=2, kinetic_scale=kinetic_scale)
    return h.compute_ground_state()
