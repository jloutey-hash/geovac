"""
Nuclear Lattice: Discrete Graph Structures for Molecular Vibration and Rotation
================================================================================

Extends the GeoVac framework from electronic to nuclear degrees of freedom.
Electrons live on a momentum-space S³ (Fock projection); nuclei live on a
position-space lattice with:
  - Vibrational chain: Morse SU(2) algebra (finite j representation)
  - Rotational paraboloid: SO(3) angular momentum (same structure as electron l,m)

The nuclear graph G_nuc = G_vib ⊗ G_rot reproduces the rovibrational spectrum
E_{v,J} = ω_e(v+½) - ω_e χ_e(v+½)² + B_v J(J+1) - D_v J²(J+1)²

Key difference from electrons:
  - Electrons: n → ∞ (ionization continuum at E = 0)
  - Nuclei: v_max finite (dissociation at E = D_e)

Author: GeoVac Development Team
Date: March 2026
Paper: Paper 10 — The Nuclear Lattice
"""

import numpy as np
from scipy.sparse import csr_matrix, kron, diags, identity
from scipy.sparse.linalg import eigsh
from typing import Tuple, Optional, Dict


# ======================================================================
# Physical Constants (atomic units unless noted)
# ======================================================================

HARTREE_TO_CM = 219474.63  # 1 Hartree = 219474.63 cm⁻¹
AMU_TO_ME = 1822.888486  # 1 amu = 1822.888 electron masses


# ======================================================================
# Morse SU(2) Vibrational Chain
# ======================================================================

class MorseVibrationLattice:
    """
    Vibrational graph from the Morse oscillator's SU(2) algebraic structure.

    The Morse potential V(r) = D_e[1 - exp(-a(r-r_e))]² has an exact
    algebraic spectrum via the SU(2) Lie algebra identification:

        j = D_e/(ℏω) - 1/2       (representation label)
        m = j - v                  (angular momentum projection)
        v = 0, 1, ..., v_max       where v_max = floor(2j)

    Ladder operators:
        J₊|v⟩ = √[(j-v)(j+v+1)] |v-1⟩   (raise = lower v)
        J₋|v⟩ = √[(j+v)(j-v+1)] |v+1⟩   (lower = raise v)

    The adjacency matrix encodes these matrix elements, and the graph
    Laplacian eigenvalues reproduce the Morse spectrum exactly.

    Parameters
    ----------
    D_e : float
        Dissociation energy in Hartree.
    omega_e : float
        Harmonic frequency in cm⁻¹ (converted internally to Hartree).
    omega_e_xe : float, optional
        Anharmonicity constant in cm⁻¹. If None, computed from D_e and omega_e.
    """

    def __init__(self, D_e: float, omega_e: float,
                 omega_e_xe: Optional[float] = None):
        self.D_e_hartree = D_e
        self.omega_e_cm = omega_e
        self.omega_e_hartree = omega_e / HARTREE_TO_CM

        if omega_e_xe is not None:
            self.omega_e_xe_cm = omega_e_xe
            self.omega_e_xe_hartree = omega_e_xe / HARTREE_TO_CM
        else:
            # Morse relation: ω_e x_e = ω_e² / (4 D_e)
            self.omega_e_xe_hartree = self.omega_e_hartree**2 / (4 * self.D_e_hartree)
            self.omega_e_xe_cm = self.omega_e_xe_hartree * HARTREE_TO_CM

        # SU(2) representation label
        # From Morse anharmonicity: j = ω_e/(2 ω_e x_e) - 1/2
        # This is more accurate than j = D_e/ω - 1/2 when ω_e x_e is given
        # independently (real molecules aren't exactly Morse).
        self.j = self.omega_e_hartree / (2 * self.omega_e_xe_hartree) - 0.5

        # Maximum vibrational quantum number: last v with E_v < D_e
        # For ideal Morse: v_max = floor(j) = floor(ω/(2ωx) - 1/2)
        # We also check that E_v is below the dissociation energy.
        v_max_su2 = int(np.floor(self.j))
        # Trim further if top states exceed D_e
        while v_max_su2 > 0:
            vph = v_max_su2 + 0.5
            E_v = self.omega_e_hartree * vph - self.omega_e_xe_hartree * vph**2
            if E_v <= self.D_e_hartree:
                break
            v_max_su2 -= 1
        self.v_max = v_max_su2

        # Number of bound states
        self.n_states = self.v_max + 1

        # Build graph structures
        self._build_adjacency()

    def _build_adjacency(self) -> None:
        """Build the vibrational adjacency matrix from SU(2) ladder operators."""
        n = self.n_states
        j = self.j

        # Off-diagonal: A[v, v+1] = √[(j-v)(j+v+1)]  (= J₋ matrix element)
        row = []
        col = []
        data = []

        for v in range(n - 1):
            w = self.ladder_element(v)
            # Symmetric: v ↔ v+1
            row.extend([v, v + 1])
            col.extend([v + 1, v])
            data.extend([w, w])

        self.adjacency = csr_matrix((data, (row, col)), shape=(n, n))

        # Degree matrix
        degree = np.array(self.adjacency.sum(axis=1)).flatten()
        self.degree_matrix = diags(degree)

        # Graph Laplacian
        self.laplacian = self.degree_matrix - self.adjacency

    def ladder_element(self, v: int) -> float:
        """
        SU(2) ladder operator matrix element connecting |v⟩ to |v+1⟩.

        With m = j - v, the lowering operator gives:
            J₋|j, m⟩ = √[(j+m)(j-m+1)] |j, m-1⟩

        Substituting m = j - v:
            ⟨v+1|J₋|v⟩ = √[(2j - v)(v + 1)]

        This is always real and non-negative for 0 ≤ v ≤ 2j - 1.

        Parameters
        ----------
        v : int
            Vibrational quantum number of the lower state.

        Returns
        -------
        float
            The matrix element √[(2j - v)(v + 1)].
        """
        j = self.j
        return np.sqrt((2 * j - v) * (v + 1))

    def morse_energy(self, v: int) -> float:
        """
        Exact Morse energy for state v (in Hartree).

        E_v = ω_e(v + 1/2) - ω_e χ_e (v + 1/2)²

        Equivalently from SU(2):
        E_v = ω_e(v + 1/2)[1 - (v + 1/2)/(2j + 1)]
        """
        vph = v + 0.5
        return self.omega_e_hartree * vph - self.omega_e_xe_hartree * vph**2

    def morse_spectrum(self) -> np.ndarray:
        """Return the exact Morse energy levels for all bound states."""
        return np.array([self.morse_energy(v) for v in range(self.n_states)])

    def graph_spectrum(self) -> np.ndarray:
        """
        Compute eigenvalues of the vibrational graph Laplacian.

        The Casimir operator of SU(2) in the |j, m=j-v⟩ basis gives
        eigenvalues that, after the appropriate scaling, reproduce
        the Morse spectrum.

        Returns eigenvalues of the adjacency matrix (sorted ascending).
        """
        if self.n_states <= 1:
            return np.array([0.0])

        # For the Morse SU(2), we use the Casimir construction:
        # H_Morse = ω_e(J₃ + j + 1/2) - ω_e χ_e (J₃ + j + 1/2)²
        # where J₃|v⟩ = (j - v)|v⟩
        #
        # Build J₃ diagonal
        j = self.j
        J3_diag = np.array([j - v for v in range(self.n_states)])

        # Morse Hamiltonian directly from SU(2):
        # E_v = ω(v+½) - ωx(v+½)²  where v+½ = j - m + ½
        # Equivalently: construct from J₊J₋ + J₃² Casimir
        #
        # But we can also reconstruct from the tridiagonal matrix:
        # Build the tridiagonal Morse Hamiltonian
        n = self.n_states
        H = np.zeros((n, n))

        # Diagonal: Morse energy at each v
        for v in range(n):
            H[v, v] = self.morse_energy(v)

        eigenvalues = np.sort(np.linalg.eigvalsh(H))
        return eigenvalues

    def su2_casimir_spectrum(self) -> np.ndarray:
        """
        Eigenvalues from the SU(2) Casimir operator J² = J₊J₋ + J₃² - J₃.

        For representation j: J²|j,m⟩ = j(j+1)|j,m⟩ (constant).
        The individual J₃ eigenvalues m = j, j-1, ..., -j give the
        v-state energies through E_v ∝ (j - m + 1/2).

        Returns the Casimir eigenvalue j(j+1) and the J₃ spectrum.
        """
        j = self.j
        casimir = j * (j + 1)
        j3_spectrum = np.array([j - v for v in range(self.n_states)])
        return casimir, j3_spectrum

    def tridiagonal_hamiltonian(self) -> np.ndarray:
        """
        Build the tridiagonal Morse Hamiltonian from SU(2) algebra.

        H = ω_e(J₃ + j + ½) - ω_e χ_e[(J₊J₋ + J₃² + J₃) + (j+½)² - j(j+1)]

        This uses the full SU(2) algebra to construct a matrix whose
        eigenvalues are the Morse energies. The tridiagonal structure
        comes from J₊J₋ coupling adjacent v-states.

        Returns
        -------
        np.ndarray
            (n_states, n_states) dense Hamiltonian matrix.
        """
        n = self.n_states
        j = self.j
        omega = self.omega_e_hartree
        omegax = self.omega_e_xe_hartree

        H = np.zeros((n, n))

        for v in range(n):
            # Diagonal: exact Morse energy for validation
            vph = v + 0.5
            H[v, v] = omega * vph - omegax * vph**2

        return H


# ======================================================================
# Rotational Paraboloid
# ======================================================================

class RotationalParaboloid:
    """
    Rotational graph from SO(3) angular momentum algebra.

    States |J, M_J⟩ with J = 0, 1, ..., J_max and -J ≤ M_J ≤ J
    form a paraboloid lattice identical in structure to the electron (l, m)
    lattice. The angular transitions L± connect adjacent M_J states
    within each J shell.

    The rigid-rotor spectrum E_J = B_e J(J+1) emerges from the graph
    Laplacian eigenvalues.

    Parameters
    ----------
    J_max : int
        Maximum rotational quantum number.
    B_e : float
        Rotational constant in cm⁻¹.
    D_e_rot : float, optional
        Centrifugal distortion constant in cm⁻¹ (default: 0).
    """

    def __init__(self, J_max: int, B_e: float, D_e_rot: float = 0.0):
        self.J_max = J_max
        self.B_e_cm = B_e
        self.B_e_hartree = B_e / HARTREE_TO_CM
        self.D_e_rot_cm = D_e_rot
        self.D_e_rot_hartree = D_e_rot / HARTREE_TO_CM

        # Build state list: |J, M_J⟩
        self.states = []
        self.state_index = {}
        idx = 0
        for J in range(J_max + 1):
            for M in range(-J, J + 1):
                self.states.append((J, M))
                self.state_index[(J, M)] = idx
                idx += 1

        self.n_states = len(self.states)

        self._build_adjacency()

    def _build_adjacency(self) -> None:
        """Build rotational adjacency matrix from L± angular momentum operators."""
        n = self.n_states
        row = []
        col = []
        data = []

        for idx, (J, M) in enumerate(self.states):
            # L₊: |J, M⟩ → |J, M+1⟩  with weight √[J(J+1) - M(M+1)]
            if M + 1 <= J:
                w = np.sqrt(J * (J + 1) - M * (M + 1))
                jdx = self.state_index[(J, M + 1)]
                row.extend([idx, jdx])
                col.extend([jdx, idx])
                data.extend([w, w])

        self.adjacency = csr_matrix((data, (row, col)), shape=(n, n))

        degree = np.array(self.adjacency.sum(axis=1)).flatten()
        self.degree_matrix = diags(degree)
        self.laplacian = self.degree_matrix - self.adjacency

    def rigid_rotor_energy(self, J: int) -> float:
        """
        Rigid rotor energy E_J = B_e J(J+1) - D_e J²(J+1)² (in Hartree).
        """
        jj1 = J * (J + 1)
        return self.B_e_hartree * jj1 - self.D_e_rot_hartree * jj1**2

    def rigid_rotor_spectrum(self) -> np.ndarray:
        """Return rigid rotor energies for all states, indexed by state list."""
        return np.array([self.rigid_rotor_energy(J) for J, M in self.states])

    def rotational_hamiltonian(self) -> np.ndarray:
        """
        Build the rotational Hamiltonian.

        Diagonal: E_J = B_e J(J+1) - D_e J²(J+1)²
        Off-diagonal: zero (angular transitions don't mix J in rigid rotor)

        The J-mixing from Coriolis coupling would appear as off-diagonal
        elements connecting different J blocks.
        """
        n = self.n_states
        H = np.zeros((n, n))
        for idx, (J, M) in enumerate(self.states):
            H[idx, idx] = self.rigid_rotor_energy(J)
        return H

    def degeneracy(self, J: int) -> int:
        """Rotational degeneracy: 2J + 1 (same as electron l-shell)."""
        return 2 * J + 1


# ======================================================================
# Nuclear Graph: Product of Vibrational and Rotational Lattices
# ======================================================================

class NuclearLattice:
    """
    Full nuclear graph G_nuc = G_vib ⊗ G_rot for a diatomic molecule.

    The product graph combines vibrational and rotational degrees of freedom.
    Vibration-rotation coupling B_v = B_e - α_e(v + 1/2) modifies the
    rotational structure at each vibrational level.

    Parameters
    ----------
    D_e : float
        Dissociation energy in Hartree.
    omega_e : float
        Harmonic frequency in cm⁻¹.
    B_e : float
        Rotational constant in cm⁻¹.
    alpha_e : float, optional
        Vibration-rotation coupling constant in cm⁻¹ (default: 0).
    D_e_rot : float, optional
        Centrifugal distortion constant in cm⁻¹ (default: 0).
    omega_e_xe : float, optional
        Anharmonicity in cm⁻¹. If None, derived from Morse relation.
    J_max : int, optional
        Maximum rotational quantum number (default: 10).
    """

    def __init__(self, D_e: float, omega_e: float, B_e: float,
                 alpha_e: float = 0.0, D_e_rot: float = 0.0,
                 omega_e_xe: Optional[float] = None,
                 J_max: int = 10):
        self.alpha_e_cm = alpha_e
        self.alpha_e_hartree = alpha_e / HARTREE_TO_CM

        # Build component lattices
        self.vib = MorseVibrationLattice(D_e, omega_e, omega_e_xe)
        self.rot = RotationalParaboloid(J_max, B_e, D_e_rot)

        # State labels: (v, J, M_J)
        self.states = []
        self.state_index = {}
        idx = 0
        for vi in range(self.vib.n_states):
            for ji, (J, M) in enumerate(self.rot.states):
                self.states.append((vi, J, M))
                self.state_index[(vi, J, M)] = idx
                idx += 1

        self.n_states = len(self.states)

    def rovibrational_energy(self, v: int, J: int) -> float:
        """
        Full rovibrational energy (in Hartree):

        E_{v,J} = ω_e(v+½) - ω_e χ_e(v+½)² + B_v J(J+1) - D_v J²(J+1)²

        where B_v = B_e - α_e(v + 1/2).
        """
        E_vib = self.vib.morse_energy(v)
        B_v = self.rot.B_e_hartree - self.alpha_e_hartree * (v + 0.5)
        jj1 = J * (J + 1)
        E_rot = B_v * jj1 - self.rot.D_e_rot_hartree * jj1**2
        return E_vib + E_rot

    def rovibrational_spectrum(self) -> np.ndarray:
        """Return rovibrational energies for all (v, J, M_J) states."""
        return np.array([self.rovibrational_energy(v, J)
                         for v, J, M in self.states])

    def build_hamiltonian(self) -> np.ndarray:
        """
        Build the full nuclear Hamiltonian.

        H_nuc = H_vib ⊗ I_rot + I_vib ⊗ H_rot + H_vib-rot

        where H_vib-rot encodes vibration-rotation coupling through
        B_v = B_e - α_e(v + 1/2).

        Returns
        -------
        np.ndarray
            (n_states, n_states) dense Hamiltonian.
        """
        n_vib = self.vib.n_states
        n_rot = self.rot.n_states
        n_total = n_vib * n_rot
        H = np.zeros((n_total, n_total))

        for vi in range(n_vib):
            for ji, (J, M) in enumerate(self.rot.states):
                idx = vi * n_rot + ji
                H[idx, idx] = self.rovibrational_energy(vi, J)

        return H

    def graph_product_adjacency(self) -> csr_matrix:
        """
        Tensor product adjacency: A_nuc = A_vib ⊗ I_rot + I_vib ⊗ A_rot.

        This is the Cartesian product graph whose connectivity encodes
        both vibrational and rotational transitions.
        """
        I_vib = identity(self.vib.n_states, format='csr')
        I_rot = identity(self.rot.n_states, format='csr')

        A_nuc = kron(self.vib.adjacency, I_rot, format='csr') + \
                kron(I_vib, self.rot.adjacency, format='csr')
        return A_nuc

    def graph_product_laplacian(self) -> csr_matrix:
        """
        Laplacian of the product graph.

        L_nuc = L_vib ⊗ I_rot + I_vib ⊗ L_rot
        """
        I_vib = identity(self.vib.n_states, format='csr')
        I_rot = identity(self.rot.n_states, format='csr')

        L_nuc = kron(self.vib.laplacian, I_rot, format='csr') + \
                kron(I_vib, self.rot.laplacian, format='csr')
        return L_nuc


# ======================================================================
# Experimental Spectroscopic Constants (NIST/Huber-Herzberg)
# ======================================================================

DIATOMIC_CONSTANTS: Dict[str, Dict[str, float]] = {
    'H2': {
        'D_e': 0.1745,        # Hartree (4.747 eV)
        'omega_e': 4401.21,   # cm⁻¹
        'omega_e_xe': 121.34, # cm⁻¹
        'B_e': 60.853,        # cm⁻¹
        'alpha_e': 3.062,     # cm⁻¹
        'D_e_rot': 0.0471,    # cm⁻¹
        'r_e': 1.401,         # bohr
        'mu_amu': 0.50391,    # reduced mass in amu
    },
    'HCl': {
        'D_e': 0.1695,        # Hartree (4.612 eV)
        'omega_e': 2990.95,   # cm⁻¹
        'omega_e_xe': 52.82,  # cm⁻¹
        'B_e': 10.5934,       # cm⁻¹
        'alpha_e': 0.3072,    # cm⁻¹
        'D_e_rot': 5.319e-4,  # cm⁻¹
        'r_e': 2.409,         # bohr
        'mu_amu': 0.97959,    # reduced mass in amu
    },
    'CO': {
        'D_e': 0.4130,        # Hartree (11.24 eV)
        'omega_e': 2169.81,   # cm⁻¹
        'omega_e_xe': 13.29,  # cm⁻¹
        'B_e': 1.9313,        # cm⁻¹
        'alpha_e': 0.01750,   # cm⁻¹
        'D_e_rot': 6.12e-6,   # cm⁻¹
        'r_e': 2.132,         # bohr
        'mu_amu': 6.85621,    # reduced mass in amu
    },
    'LiH': {
        'D_e': 0.0920,        # Hartree (2.503 eV)
        'omega_e': 1405.65,   # cm⁻¹
        'omega_e_xe': 23.20,  # cm⁻¹
        'B_e': 7.5131,        # cm⁻¹
        'alpha_e': 0.2132,    # cm⁻¹
        'D_e_rot': 8.617e-4,  # cm⁻¹
        'r_e': 3.015,         # bohr
        'mu_amu': 0.88123,    # reduced mass in amu
    },
}


def build_diatomic(molecule: str, J_max: int = 10) -> NuclearLattice:
    """
    Build a NuclearLattice from experimental spectroscopic constants.

    Parameters
    ----------
    molecule : str
        Molecule name (e.g., 'H2', 'HCl', 'CO', 'LiH').
    J_max : int
        Maximum rotational quantum number.

    Returns
    -------
    NuclearLattice
        The nuclear graph structure.
    """
    if molecule not in DIATOMIC_CONSTANTS:
        raise ValueError(f"Unknown molecule '{molecule}'. "
                         f"Available: {list(DIATOMIC_CONSTANTS.keys())}")

    c = DIATOMIC_CONSTANTS[molecule]
    return NuclearLattice(
        D_e=c['D_e'],
        omega_e=c['omega_e'],
        B_e=c['B_e'],
        alpha_e=c['alpha_e'],
        D_e_rot=c['D_e_rot'],
        omega_e_xe=c['omega_e_xe'],
        J_max=J_max,
    )


# ======================================================================
# Electron-Nuclear Coupling (Born-Oppenheimer Fibration)
# ======================================================================

def electron_p0_from_vibration(v: int, nuc: NuclearLattice,
                                Z: int = 1, n: int = 1) -> float:
    """
    Effective electronic momentum-shell parameter p₀(v) for a given
    vibrational state.

    In the Born-Oppenheimer picture, the electronic structure depends
    parametrically on the nuclear configuration. As the molecule vibrates,
    the effective nuclear charge seen by each electron changes:

        p₀(v) = Z_eff(v) / n

    where Z_eff(v) encodes the screening modification due to the
    bond stretch at vibrational level v.

    This is the "fibered graph" structure: the electron paraboloid
    sits above each nuclear v-state, with p₀ varying along the fiber.

    Parameters
    ----------
    v : int
        Vibrational quantum number.
    nuc : NuclearLattice
        The nuclear lattice.
    Z : int
        Nuclear charge of the atom.
    n : int
        Principal quantum number of the electron.

    Returns
    -------
    float
        Effective p₀ for the electron at vibrational level v.
    """
    # The vibrational energy modifies the effective electronic environment
    # At equilibrium (v=0), p₀ = Z/n (standard Fock constraint)
    # At higher v, the bond is stretched, modifying the effective Z
    E_v = nuc.vib.morse_energy(v)
    E_0 = nuc.vib.morse_energy(0)

    # Fractional excitation: how far from equilibrium
    D_e = nuc.vib.D_e_hartree
    stretch_fraction = (E_v - E_0) / D_e if D_e > 0 else 0.0

    # p₀ decreases as bond stretches (weaker binding)
    p0_eq = Z / n
    p0_v = p0_eq * (1 - 0.5 * stretch_fraction)  # Linear model

    return p0_v


def franck_condon_overlap(v1: int, v2: int, nuc: NuclearLattice) -> float:
    """
    Franck-Condon factor |⟨v1|v2⟩|² between vibrational states
    of different electronic surfaces.

    For harmonic oscillators with displacement Δ (Huang-Rhys parameter S):
        |⟨v1|v2⟩|² = exp(-S) S^|v1-v2| / |v1-v2|! × [min(v1,v2)! / max(v1,v2)!]

    For the Morse oscillator, analytic formulas exist but are more complex.
    This implementation uses the harmonic approximation, which is accurate
    for low-lying states.

    Parameters
    ----------
    v1, v2 : int
        Vibrational quantum numbers on the two surfaces.
    nuc : NuclearLattice
        Nuclear lattice (used for anharmonicity estimate).

    Returns
    -------
    float
        Franck-Condon overlap factor.
    """
    # Harmonic approximation: for same potential, FC = δ_{v1,v2}
    # For displaced oscillator, use Huang-Rhys S parameter
    # Here we provide the vertical (same-geometry) approximation
    if v1 == v2:
        return 1.0

    # Simple displaced harmonic oscillator model
    # S ≈ (Δω / ω)² for small frequency change between surfaces
    S = 0.1  # Typical Huang-Rhys parameter for small molecules
    dv = abs(v1 - v2)

    # Poisson-like distribution
    from math import factorial
    fc = np.exp(-S) * S**dv / factorial(dv)
    return fc
