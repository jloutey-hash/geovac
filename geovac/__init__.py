"""
GeoVac: Computational Quantum Chemistry via Spectral Graph Theory
=================================================================

Solves electronic structure using sparse graph Laplacians over hydrogenic
quantum numbers. Universal kinetic scale K_vac = -1/16.

**Key Classes:**
- AtomicSolver: Single-electron atoms with Z^2 scaling
- LatticeIndex: N-electron FCI via relational lattice index (He, Li validated)
- MoleculeHamiltonian: Molecular bonding via sparse bridges (H2 Full CI)
- GeometricLattice: Core sparse graph lattice with (n,l,m) quantum state nodes
- TimePropagator: Crank-Nicolson unitary time evolution

**FCI Accuracy (slater_full):**
- He (2e): 0.35% at max_n=5 (hybrid h1, monotonic convergence)
- Li (3e): 1.10% at max_n=4 (exact h1, monotonic convergence)
- H (1e): 0.57% at max_n=30

Quick Start:
-----------
>>> from geovac import LatticeIndex
>>> idx = LatticeIndex(n_electrons=2, max_n=3, nuclear_charge=2,
...                    vee_method='slater_full', h1_method='hybrid')
>>> H = idx.assemble_hamiltonian()
"""

__version__ = '1.2.0'
__author__ = 'J. Loutey'
__license__ = 'MIT'

# Core imports - expose main classes at top level
from .lattice import GeometricLattice
from .hamiltonian import MoleculeHamiltonian
from .dirac_hamiltonian import DiracHamiltonian, DiracLatticeStates
from .atomic_solver import AtomicSolver, solve_hydrogen, solve_atom
from .dynamics import TimePropagator
from .lattice_index import LatticeIndex, MolecularLatticeIndex, compute_vee_s3_overlap, compute_bsse_correction, compute_cross_atom_J, compute_cross_atom_K, compute_overlap_element, compute_atomic_p0, compute_exact_cross_nuclear
from .direct_ci import DirectCISolver
from .frozen_core import FrozenCoreLatticeIndex
from .locked_shell import LockedShellMolecule
from .qubit_encoding import JordanWignerEncoder, PauliAnalysis
from .aimd import VelocityVerlet, LangevinThermostat, run_lih_aimd, run_li_nve
from .benchmark import run_unitarity_test, run_scaling_benchmark, run_li_energy_audit

# Holographic/AdS-CFT modules live in ADSCFT/ package:
#   from ADSCFT import MuonicHydrogenSolver, compute_holographic_entropy, etc.

# Define public API
__all__ = [
    'GeometricLattice',
    'MoleculeHamiltonian',
    'DiracHamiltonian',
    'DiracLatticeStates',
    'AtomicSolver',
    'solve_hydrogen',
    'solve_atom',
    'TimePropagator',
    'LatticeIndex',
    'MolecularLatticeIndex',
    'JordanWignerEncoder',
    'PauliAnalysis',
    'FrozenCoreLatticeIndex',
    'LockedShellMolecule',
    'compute_bsse_correction',
    'compute_cross_atom_K',
    'VelocityVerlet',
    'LangevinThermostat',
    'run_lih_aimd',
    'run_li_nve',
    'run_unitarity_test',
    'run_scaling_benchmark',
    'run_li_energy_audit',
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

