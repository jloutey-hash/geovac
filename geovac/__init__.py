"""
GeoVac: spectral-graph discretization for quantum chemistry
============================================================

Sparse Hamiltonians built on the (n, l, m) label set of the S^3
hyperspherical harmonics (Fock's projection; Papers 0, 7), with natural
geometries for molecules (prolate spheroidal, hyperspherical, composed;
Papers 11, 12, 13, 15, 17) and qubit-Hamiltonian export (Paper 14).

The atomic default is the physical diagonal -Z^2/(2 n^2).  The legacy graph
model kappa*(D - A) matches the ground state only in the band-edge limit
(kappa = -1/16, an Observation); its excited states are not hydrogen's, and the old
bridged-graph (LCAO concatenation) molecular picture is superseded
(FCI-M guardrail: R-independent kinetic term, no equilibrium).

Quick Start (Atoms):
-------------------
>>> from geovac import AtomicSolver
>>> solver = AtomicSolver(max_n=10, Z=1)

Helium: use the Level-3 hyperspherical solver
``geovac.hyperspherical_radial.solve_helium`` (zero fitted parameters).

Archived (2026-09-26): the legacy HeliumHamiltonian / HeliumPackingSolver /
MoleculeHamiltonian classes, the top-level ``solve_helium`` wrapper and
``CALIBRATED_KINETIC_SCALE`` (-0.10298808, FITTED to the experimental helium
energy) now live in ``geovac._archive.superseded.hamiltonian``. The old
top-level names still resolve, with a DeprecationWarning.

Available Classes:
-----------------
- GeometricLattice: Sparse graph lattice with (n,l,m) quantum state nodes
- AtomicSolver: single-electron atomic solver
- DiracHamiltonian: Relativistic spinor-based solver (experimental)

Status: measured results live in the papers and docs/claims_register.md.
"""

__version__ = '0.4.0'
__author__ = 'J. Loutey'
__license__ = 'MIT'

# Core imports - expose main classes at top level
from .lattice import GeometricLattice
from .dirac_hamiltonian import DiracHamiltonian, DiracLatticeStates
from .atomic_solver import AtomicSolver, solve_hydrogen, solve_atom

# Archived modules (moved to old_research_archive/specialized_solvers/):
# - muonic_hydrogen
# - holographic_analysis
# - fundamental_constants
# - hyperfine_contact
# - proton_radius_improved

# Define public API
__all__ = [
    'GeometricLattice',
    'DiracHamiltonian',
    'DiracLatticeStates',
    'AtomicSolver',
    'solve_hydrogen',
    'solve_atom',
    '__version__',
]

# Package metadata
UNIVERSAL_KINETIC_SCALE = -1/16  # matched energy scale kappa (an Observation, Papers 0/7; -0.0625)
EXPERIMENTAL_HELIUM_ENERGY = -2.90338583  # Hartree (NIST reference)
HYDROGEN_GROUND_STATE = -0.5  # Hartree (exact)

# ---------------------------------------------------------------------------
# Deprecated legacy names (archived 2026-09-26 to
# geovac._archive.superseded.hamiltonian). They carried a kinetic scale
# -0.10298808 FITTED to the experimental helium ground-state energy -- an
# empirical parameter, against the no-fitted-parameters rule. Resolved lazily
# (PEP 562) so that plain ``import geovac`` emits no warning; any USE warns.
# ---------------------------------------------------------------------------
import warnings as _warnings

_LEGACY_FITTED_MSG = (
    "geovac.{name} is archived (superseded): the legacy helium/molecule "
    "classes use a kinetic scale -0.10298808 FITTED to the experimental "
    "helium energy (an empirical parameter). Import it from "
    "geovac._archive.superseded.hamiltonian if you need the old behaviour; "
    "for helium use geovac.hyperspherical_radial.solve_helium (Level-3 "
    "hyperspherical solver, zero fitted parameters)."
)

_LEGACY_NAMES = frozenset({
    'HeliumHamiltonian',
    'HeliumPackingSolver',
    'MoleculeHamiltonian',
    'CALIBRATED_KINETIC_SCALE',
})


def __getattr__(name):
    if name in _LEGACY_NAMES:
        _warnings.warn(_LEGACY_FITTED_MSG.format(name=name),
                       DeprecationWarning, stacklevel=2)
        from ._archive.superseded import hamiltonian as _legacy
        return getattr(_legacy, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def solve_helium(max_n=3, kinetic_scale=None):
    """DEPRECATED legacy wrapper -- uses a kinetic scale FITTED to experiment.

    Delegates, bit-identically, to the archived
    ``geovac._archive.superseded.hamiltonian.solve_helium`` (which builds
    ``HeliumHamiltonian(max_n, Z=2, kinetic_scale=-0.10298808)`` and returns
    ``compute_ground_state()``), after a DeprecationWarning. The default
    -0.10298808 was fitted to the experimental helium ground-state energy.

    For helium use ``geovac.hyperspherical_radial.solve_helium`` (the Level-3
    hyperspherical solver; different signature and return type).
    """
    _warnings.warn(_LEGACY_FITTED_MSG.format(name='solve_helium'),
                   DeprecationWarning, stacklevel=2)
    from ._archive.superseded import hamiltonian as _legacy
    if kinetic_scale is None:
        kinetic_scale = _legacy.CALIBRATED_KINETIC_SCALE
    return _legacy.solve_helium(max_n=max_n, kinetic_scale=kinetic_scale)
